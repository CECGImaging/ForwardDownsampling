/*************************************************************************************************
 * barycentricInterp                                                                             *
 * MATLAB extension to calculate an operator matrix for linear interpolation                     *
 * of values from one triangle mesh to another using barycentric coordinates                     *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
 * Authors: Steffen Schuler                                                                      *
 *          Institute of Biomedical Engineering                                                  *
 *          Karlsruhe Institute of Technology                                                    *
 *          www.ibt.kit.edu                                                                      *
 * Date:    October 2018                                                                         *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
 * Before first usage:                                                                           *
 * Compile the MEX using the script mexCompile_barycentricInterp.m                               *
 *                                                                                               *
 * Usage inside MATLAB:                                                                          *
 * barycentricInterp(struct sourceMesh, struct targetMesh, double normalLength, [bool verbose])  *
 *                                                                                               *
 * The function creates a matrix 'interpMat' in the workspace. This matrix can be used to        *
 * linearly interpolate values from all points in sourceMesh to all points in targetMesh.        *
 *                                                                                               *
 * Additionally, two debug files 'barycentricInterp_intersections.vtp' and                       *
 * 'barycentricInterp_target.vtp' are created in the current folder.                             *
 *                                                                                               *
 * sourceMesh and targetMesh must have the following fields (as created by readVTK):             *
 *     points:    [numPoints x 3 double]                                                         *
 *     cells:     [numCells x 3 int32]                                                           *
 *                                                                                               *
 * normalLength defines the maximum search distance in positive and negative normal direction.   *
 *                                                                                               *
 * verbose is false on default                                                                   *
 *************************************************************************************************/

// vtk
#include <vtkSmartPointer.h>
#include <vtkErrorCode.h>
#include <vtkPolyData.h>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkPointLocator.h>
#include <vtkOBBTree.h>
#include <vtkIdList.h>
#include <vtkTriangle.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>

// mex
#include "mex.h"
#include "matrix.h"

// std lib
#include <string>
#include <sstream>
#include <map>

// math_pack
#include "Vector3.h"

using namespace math_pack;

bool verbose_;
void printVerbose(const char* msg)
{
    if(verbose_)
    {
        mexPrintf(msg);
        mexEvalString("drawnow;");
    }
}

void GetPointNormal(vtkPolyData* mesh, vtkIdType pointId, double normalLength, Vector3<double>& normal, double &concavity)
{
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> neighborPointIds = vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(pointId, cellIds);
    
    normal = Vector3<double>(0, 0, 0);
    for(vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
    {
        vtkIdType cellId = cellIds->GetId(i);
        mesh->GetCellPoints(cellId, pointIds);
        if(pointIds->GetNumberOfIds() != 3)
            throw std::runtime_error("getPointNormal(): Cell " + std::to_string(cellId) + " is not a triangle.");
        double v[9];
        for(int i = 0; i < 3; i++)
        {
            vtkIdType neighborPointId = pointIds->GetId(i);
            double* vv = v+(3*i);
            mesh->GetPoint(neighborPointId, vv);
            if(neighborPointId != pointId)
                neighborPointIds->InsertUniqueId(neighborPointId);
        }
        Vector3<double> n;
        vtkTriangle::ComputeNormal(v, v+3, v+6, n.GetArray());
        normal += n;
    }
    normal.Normalize();
    normal *= normalLength;
    
    Vector3<double> p(0, 0, 0);
    vtkIdType numNeighborPoints = neighborPointIds->GetNumberOfIds();
    for(int i = 0; i < numNeighborPoints; i++)
    {
        Vector3<double> q;
        mesh->GetPoint(neighborPointIds->GetId(i), q.GetArray());
        p += q;
    }
    Vector3<double> u;
    mesh->GetPoint(pointId, u.GetArray());
    p = p/numNeighborPoints - u;
    concavity = p * normal/normal.Norm();
}

void GetIntersection(vtkPolyData* mesh, vtkIdList* pointIds, Vector3<double> p, Vector3<double> n, Vector3<double>& x, double* barycentricCoords)
{
    Vector3<double> a, b, c, triangleNormal;
    mesh->GetPoint(pointIds->GetId(0), a.GetArray());
    mesh->GetPoint(pointIds->GetId(1), b.GetArray());
    mesh->GetPoint(pointIds->GetId(2), c.GetArray());
    vtkTriangle::ComputeNormal(a.GetArray(), b.GetArray(), c.GetArray(), triangleNormal.GetArray());
    
    double d = (a-p)*triangleNormal / (n*triangleNormal);
    x = p + d*n;
    
    Vector3<double> v0 = b-a, v1 = c-a, v2 = x-a;
    double d00 = v0*v0;
    double d01 = v0*v1;
    double d11 = v1*v1;
    double d20 = v2*v0;
    double d21 = v2*v1;
    double denom = d00*d11 - d01*d01;
    barycentricCoords[1] = (d11*d20 - d01*d21)/denom;
    barycentricCoords[2] = (d00*d21 - d01*d20)/denom;
    
    double tol = 1e-10;
    if(barycentricCoords[1] < -tol || barycentricCoords[1] > 1+tol || barycentricCoords[2] < -tol || barycentricCoords[2] > 1+tol)
        throw std::runtime_error("Barycentric coordinates out of range [0,1].");
    
    barycentricCoords[1] = (barycentricCoords[1] < 1) ? barycentricCoords[1] : 1;
    barycentricCoords[1] = (barycentricCoords[1] > 0) ? barycentricCoords[1] : 0;
    barycentricCoords[2] = (barycentricCoords[2] < 1) ? barycentricCoords[2] : 1;
    barycentricCoords[2] = (barycentricCoords[2] > 0) ? barycentricCoords[2] : 0;
    
    barycentricCoords[0] = 1.0 - barycentricCoords[1] - barycentricCoords[2];
}

struct linearCombination
{
    vtkIdType pointIds[3] = {0};
    double coefficients[3] = {0};
};

/* MATLAB entry function
 * nlhs/nrhs contain the number of left/right-hand-side arguments to this function
 * plhs/prhs are arrays of pointers to the arguments in MATLAB data format */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nlhs > 0)
        mexErrMsgTxt("Too many output arguments. There is no return value. Instead, a variable 'interpMat' is created in the workspace.");
    
    if(nrhs < 3)
        mexErrMsgTxt("Not enough input arguments. Syntax: barycentricInterp(struct sourceMesh, struct targetMesh, double normalLength, [bool verbose])");
    if(nrhs > 4)
        mexErrMsgTxt("Too many input arguments. Syntax: barycentricInterp(struct sourceMesh, struct targetMesh, double normalLength, [bool verbose])");
    
    double normalLength = mxGetScalar(prhs[2]); // Defines the maximum search distance in positive and negative normal direction
    
    if(nrhs == 4)
        verbose_ = mxGetScalar(prhs[3]);
    else
        verbose_ = false;
    
    ///// Read source mesh /////
    
    printVerbose("Reading source mesh...\n");
    
    const mxArray* sourceMesh = prhs[0];
    if(!mxIsStruct(sourceMesh))
        mexErrMsgTxt("First argument must be a struct. Syntax: barycentricInterp(struct sourceMesh, struct targetMesh, double normalLength, [bool verbose])");
    
    mxArray* points = mxGetField(sourceMesh, 0, "points");
    if(points == NULL)
        mexErrMsgTxt("Field 'points' of sourceMesh could not be read.");
    if(!mxIsDouble(points))
        mexErrMsgTxt("Field 'points' of sourceMesh must be of class double. You could try points = double(points)");

    mxArray* cells = mxGetField(sourceMesh, 0, "cells");
    if(cells == NULL)
        mexErrMsgTxt("Field 'cells' of sourceMesh could not be read.");
    if(!mxIsClass(cells, "int32"))
        mexErrMsgTxt("Field 'cells' of sourceMesh must be of class int32. You could try cells = int32(cells)");
    
    // Process cells
    
    printVerbose("Processing cells...\n");
    
    mwSize numCells = mxGetM(cells);
    if(mxGetN(cells) != 3)
        mexErrMsgTxt("sourceMesh.cells must have exactly 3 columns. Make shure the mesh is a pure triangle mesh.");
    vtkIdType pointIds[3];
    
    vtkSmartPointer<vtkPolyData> sourcePolyData = vtkSmartPointer<vtkPolyData>::New();
    sourcePolyData->Allocate(numCells);
    for(unsigned long i = 0; i < numCells; i++)
    {
        for(unsigned long j = 0; j < 3; j++)
        {
            pointIds[j] = ((int32_t*)mxGetData(cells))[numCells*j+i]-1;
        }
        sourcePolyData->InsertNextCell(VTK_TRIANGLE, 3, pointIds);
    }
    
    // Process points
    
    printVerbose("Processing points...\n");
    
    vtkSmartPointer<vtkPoints> sourcePoints = vtkSmartPointer<vtkPoints>::New();
    unsigned long numSourcePoints = mxGetM(points);
    sourcePoints->SetNumberOfPoints(numSourcePoints);
    for(unsigned long i = 0; i < numSourcePoints; i++)
    {
        double pointCoords[3];
        for(int j = 0; j < 3; j++)
        {
            pointCoords[j] = mxGetPr(points)[numSourcePoints*j+i];
        }
        sourcePoints->SetPoint(i, pointCoords);
    }
    sourcePoints->Modified();
    sourcePolyData->SetPoints(sourcePoints);

    ///// Read target mesh /////
    
    printVerbose("Reading target mesh...\n");
    
    const mxArray* targetMesh = prhs[1];
    if(!mxIsStruct(targetMesh))
        mexErrMsgTxt("First argument must be a struct. Syntax: barycentricInterp(struct targetMesh, struct targetMesh, double normalLength, [bool verbose])");
    
    points = mxGetField(targetMesh, 0, "points");
    if(points == NULL)
        mexErrMsgTxt("Field 'points' of targetMesh could not be read.");
    if(!mxIsDouble(points))
        mexErrMsgTxt("Field 'points' of targetMesh must be of class double. You could try points = double(points)");
    
    cells = mxGetField(targetMesh, 0, "cells");
    if(cells == NULL)
        mexErrMsgTxt("Field 'cells' of targetMesh could not be read.");
    if(!mxIsClass(cells, "int32"))
        mexErrMsgTxt("Field 'cells' of targetMesh must be of class int32. You could try cells = int32(cells)");
    
    // Process cells
    
    printVerbose("Processing cells...\n");
    
    numCells = mxGetM(cells);
    if(mxGetN(cells) != 3)
        mexErrMsgTxt("targetMesh.cells must have exactly 3 columns. Make shure the mesh is a pure triangle mesh.");
    
    vtkSmartPointer<vtkPolyData> targetPolyData = vtkSmartPointer<vtkPolyData>::New();
    targetPolyData->Allocate(numCells);
    for(unsigned long i = 0; i < numCells; i++)
    {
        for(unsigned long j = 0; j < 3; j++)
        {
            pointIds[j] = ((int32_t*)mxGetData(cells))[numCells*j+i]-1;
        }
        targetPolyData->InsertNextCell(VTK_TRIANGLE, 3, pointIds);
    }
    
    // Process points
    
    printVerbose("Processing points...\n");
    
    vtkSmartPointer<vtkPoints> targetPoints = vtkSmartPointer<vtkPoints>::New();
    unsigned long numTargetPoints = mxGetM(points);
    targetPoints->SetNumberOfPoints(numTargetPoints);
    for(unsigned long i = 0; i < numTargetPoints; i++)
    {
        double pointCoords[3];
        for(int j = 0; j < 3; j++)
        {
            pointCoords[j] = mxGetPr(points)[numTargetPoints*j+i];
        }
        targetPoints->SetPoint(i, pointCoords);
    }
    targetPoints->Modified();
    targetPolyData->SetPoints(targetPoints);
    
    ///// For each point in targetPolyData: Find triangle in sourcePolyData that is intersected by the line in normal direction passing through this point /////
    
    vtkSmartPointer<vtkDoubleArray> normalsArray = vtkSmartPointer<vtkDoubleArray>::New();
    normalsArray->SetNumberOfComponents(3);
    normalsArray->SetNumberOfTuples(targetPolyData->GetNumberOfPoints());
    normalsArray->SetName("normals");
    targetPolyData->GetPointData()->AddArray(normalsArray);
    
    vtkSmartPointer<vtkDoubleArray> concavityArray = vtkSmartPointer<vtkDoubleArray>::New();
    concavityArray->SetNumberOfComponents(1);
    concavityArray->SetNumberOfTuples(targetPolyData->GetNumberOfPoints());
    concavityArray->SetName("concavity");
    targetPolyData->GetPointData()->AddArray(concavityArray);
    
    vtkSmartPointer<vtkIntArray> notFoundArray = vtkSmartPointer<vtkIntArray>::New();
    notFoundArray->SetNumberOfComponents(1);
    notFoundArray->SetNumberOfTuples(targetPolyData->GetNumberOfPoints());
    notFoundArray->SetName("notFound");
    targetPolyData->GetPointData()->AddArray(notFoundArray);
    
    vtkSmartPointer<vtkIntArray> foundIdArray = vtkSmartPointer<vtkIntArray>::New();
    foundIdArray->SetNumberOfComponents(1);
    foundIdArray->SetNumberOfTuples(targetPolyData->GetNumberOfPoints());
    foundIdArray->SetName("foundId");
    targetPolyData->GetPointData()->AddArray(foundIdArray);
    
    vtkSmartPointer<vtkOBBTree> obbTree = vtkSmartPointer<vtkOBBTree>::New();
    obbTree->SetDataSet(sourcePolyData);
    obbTree->BuildLocator();
    //obbTree->SetTolerance(tolerance);
    
    vtkSmartPointer<vtkIdList> sourcePointIds = vtkSmartPointer<vtkIdList>::New();
    std::map<vtkIdType, linearCombination> sourceToTargetMapping;
    vtkIdType numNotFound = 0;
    
    vtkSmartPointer<vtkPolyData> intersectionMesh = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> intersectionPoints = vtkSmartPointer<vtkPoints>::New();
    intersectionMesh->SetPoints(intersectionPoints);
    
    std::stringstream ss;
    
    for(vtkIdType id = 0; id < numTargetPoints; id++)
    {
        Vector3<double> n, p1, p2, p3;
        targetPolyData->GetPoint(id, p1.GetArray());
        double concavity;
        GetPointNormal(targetPolyData, id, normalLength, n, concavity);
        normalsArray->SetTuple(id, n.GetArray());
        concavityArray->SetTuple1(id, concavity);
        p2 = p1 + n;
        p3 = p1 - n;
        double t, x[3], pCoords[3];
        int subId;
        vtkIdType sourceCellId;
        int found = obbTree->IntersectWithLine(p1.GetArray(), p2.GetArray(), 0.0, t, x, pCoords, subId, sourceCellId);
        if(!found)
            found = obbTree->IntersectWithLine(p1.GetArray(), p3.GetArray(), 0.0, t, x, pCoords, subId, sourceCellId);
        if(!found)
        {
            notFoundArray->SetTuple1(id, 1);
            foundIdArray->SetTuple1(id, -1);
            ss << "target point " << id << " --> NO CORRESPONDING SOURCE CELL FOUND." << std::endl;
            
            linearCombination lc;
            sourceToTargetMapping.insert(std::pair<vtkIdType, linearCombination>(id, lc));
            
            numNotFound++;
        }
        else // found
        {
            notFoundArray->SetTuple1(id, 0);
            foundIdArray->SetTuple1(id, sourceCellId);
            ss << "target point " << id << " --> source cell " << sourceCellId;
            
            sourcePolyData->GetCellPoints(sourceCellId, sourcePointIds);
            Vector3<double> x;
            linearCombination lc;
            GetIntersection(sourcePolyData, sourcePointIds, p1, n, x, lc.coefficients);
            intersectionPoints->InsertNextPoint(x.GetArray());
            
            ss << "          " << id << "  =  " << std::fixed << std::setprecision(2);
            for(vtkIdType i = 0; i < 3; i++)
            {
                lc.pointIds[i] = sourcePointIds->GetId(i);
                ss << lc.coefficients[i] << " * " << lc.pointIds[i];
                if(i < 2)
                    ss << "  +  ";
            }
            ss << std::endl;
            sourceToTargetMapping.insert(std::pair<vtkIdType, linearCombination>(id, lc));
        }
        //printVerbose(ss.str().c_str());
        ss.str("");
    }
    printVerbose(("Number of target nodes, for which a corresponding source triangle could be found:     " + std::to_string(sourceToTargetMapping.size()) + "\nNumber of target nodes, for which a corresponding source triangle could not be found: " + std::to_string(numNotFound) + "\n").c_str());
    
    ///// Write files /////
    
    printVerbose("Writing debug files...\n");
    
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("barycentricInterp_target.vtp");
    writer->SetInputData(targetPolyData);
    writer->Write();
    
    writer->SetFileName("barycentricInterp_intersections.vtp");
    writer->SetInputData(intersectionMesh);
    writer->Write();
    
    //////////
    
    std::stringstream sstream;
    sstream << "interpMat = sparse(" + std::to_string(numTargetPoints) + "," + std::to_string(numSourcePoints) + ");\n";
    for(auto sourceToTarget : sourceToTargetMapping)
    {
        std::string iStr = std::to_string(sourceToTarget.first+1);
        linearCombination lc = sourceToTarget.second;
        sstream << "interpMat(" + iStr + "," + std::to_string(lc.pointIds[0]+1) + ")=" + std::to_string(lc.coefficients[0]) + "; ";
        sstream << "interpMat(" + iStr + "," + std::to_string(lc.pointIds[1]+1) + ")=" + std::to_string(lc.coefficients[1]) + "; ";
        sstream << "interpMat(" + iStr + "," + std::to_string(lc.pointIds[2]+1) + ")=" + std::to_string(lc.coefficients[2]) + ";\n";
    }
    mexEvalString(sstream.str().c_str());
}
