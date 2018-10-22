/*******************************************************************************
 * readVTK.cxx                                                                 *
 * MATLAB extension to read VTK files (legacy: vtk or XML: vtp, vtu)           *
 * Also supported are: ply, stl, obj, (off not yet)                            *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
 * Author: Steffen Schuler                                                     *
 *         Institute of Biomedical Engineering                                 *
 *         Karlsruhe Institute of Technology                                   *
 *         www.ibt.kit.edu                                                     *
 * Date:   June 2017                                                           *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
 * Before first usage:                                                         *
 * Compile the MEX using the script mexCompile_readVTK.m                       *
 *                                                                             *
 * Usage inside MATLAB:                                                        *
 * outStruct = readVTK(char filename, [bool verbose]);                         *
 *                                                                             *
 * verbose is false on default                                                 *
 *******************************************************************************/

// general
#include <vtkSmartPointer.h>

// legacy format reader
#include <vtkDataSetReader.h>

// XML format readers
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>

// readers for other common polygon mesh formats
#include <vtkPLYReader.h>
#include <vtkSTLReader.h>
#include <vtkOBJReader.h>
//#include "additionalReaders/vtkOFFReader.h"

// data sets
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

// data set attributes
#include <vtkPointData.h>
#include <vtkCellData.h>

// mex
#include "mex.h"
#include "matrix.h"

// This function is adapted from vtkMatlabMexAdapter
mxClassID GetMatlabDataType(vtkDataArray* inArray)
{
    int dataType = inArray->GetDataType();
    switch(dataType)
    {
        case VTK_BIT:
            return(mxUINT8_CLASS);
        case VTK_CHAR:
            return(mxINT8_CLASS);
        case VTK_SIGNED_CHAR:
            return(mxCHAR_CLASS);
        case VTK_UNSIGNED_CHAR:
            return(mxUINT8_CLASS);
        case VTK_SHORT:
            return(mxINT16_CLASS);
        case VTK_UNSIGNED_SHORT:
            return(mxUINT16_CLASS);
        case VTK_INT:
            return(mxINT32_CLASS);
        case VTK_ID_TYPE:
            return(mxINT32_CLASS);
        case VTK_UNSIGNED_INT:
            return(mxUINT32_CLASS);
        case VTK_LONG:
            return(mxINT64_CLASS);
        case VTK_UNSIGNED_LONG:
            return(mxUINT64_CLASS);
        case VTK_LONG_LONG:
            return(mxINT64_CLASS);
        case VTK_UNSIGNED_LONG_LONG:
            return(mxUINT64_CLASS);
        case VTK_FLOAT:
            return(mxSINGLE_CLASS);
        case VTK_DOUBLE:
            return(mxDOUBLE_CLASS);
        default:
            mexWarnMsgTxt((std::string("Data array type ") + std::to_string(dataType) + " is not supported. Using double.").c_str());
            return(mxDOUBLE_CLASS);
    }
}

// This function is adapted from vtkMatlabMexAdapter
mxArray* vtkDataArrayToMxArray(vtkDataArray* inArray)
{
    mxClassID classID = GetMatlabDataType(inArray);

    mwSize numTuples = inArray->GetNumberOfTuples();
    mwSize numComponents = inArray->GetNumberOfComponents();
    
    mxArray* outArray = mxCreateNumericMatrix(numTuples, numComponents, classID, mxREAL);
    mwSize numBytes = mxGetElementSize(outArray);
    
    if(numBytes != inArray->GetElementComponentSize())
        mexErrMsgTxt("Data size mismatch between Matlab and VTK.");
    
    unsigned char* source;
    unsigned char* dest = (unsigned char*)mxGetData(outArray);
    for(int i = 0; i < numTuples; i++)
    {
        for(int j = 0; j < numComponents; j++)
        {
            source = (unsigned char*)inArray->GetVoidPointer(i*numComponents + j);
            for(int k = 0; k < numBytes; k++)
                dest[(numTuples*numBytes)*j + numBytes*i + k] = source[k];
        }
    }
    return(outArray);
}

int getNumberOfCellPoints(uint8_t cellType)
{
    switch(cellType)
    {
        case VTK_VERTEX:
            return(1);
        case VTK_LINE:
            return(2);
        case VTK_TRIANGLE:
            return(3);
        case VTK_PIXEL:
            return(4);
        case VTK_QUAD:
            return(4);
        case VTK_TETRA:
            return(4);
        case VTK_VOXEL:
            return(8);
        case VTK_HEXAHEDRON:
            return(8);
        case VTK_WEDGE:
            return(6);
        case VTK_PYRAMID:
            return(5);
        case VTK_PENTAGONAL_PRISM:
            return(10);
        case VTK_HEXAGONAL_PRISM:
            return(12);
        case VTK_QUADRATIC_EDGE:
            return(3);
        case VTK_QUADRATIC_TRIANGLE:
            return(6);
        case VTK_QUADRATIC_QUAD:
            return(8);
        case VTK_QUADRATIC_TETRA:
            return(10);
        case VTK_QUADRATIC_HEXAHEDRON:
            return(20);
        case VTK_QUADRATIC_WEDGE:
            return(15);
        case VTK_QUADRATIC_PYRAMID:
            return(13);
        default:
            mexErrMsgTxt(("Cell type " + std::string(vtkCellTypes::GetClassNameFromTypeId(cellType)) + " is not supported.").c_str());
    }
}

bool verbose_;
void printVerbose(const char* msg)
{
    if(verbose_)
    {
        mexPrintf(msg);
        mexEvalString("drawnow;");
    }
}

/* MATLAB entry function
 * nlhs/nrhs contain the number of left/right-hand-side arguments to this function
 * plhs/prhs are arrays of pointers to the arguments in MATLAB data format */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs < 1)
        mexErrMsgTxt("Not enough input arguments. Syntax: outStruct = readVTK(char filename, [bool verbose])");
    if(nrhs > 2)
        mexErrMsgTxt("Too many input arguments. Syntax: outStruct = readVTK(char filename, [bool verbose])");
    if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments. Syntax: outStruct = readVTK(char filename, [bool verbose])");
    
    if(nrhs == 2)
        verbose_ = mxGetScalar(prhs[1]);
    else
        verbose_ = false;
    
    ///// Read file /////
    
    printVerbose("Reading file...\n");
    
    std::string path(mxArrayToString(prhs[0]));
    long pos = path.find_last_of(".");
    std::string extension = path.substr(pos+1, path.size()-pos);
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    
    vtkSmartPointer<vtkPointSet> data;
    
    if(extension == "vtk")
    {
        vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
        reader->SetFileName(path.c_str());
        if(!reader->OpenVTKFile())
            mexErrMsgTxt("File cannot be opened by vtkDataSetReader. Does it exist?");
        reader->Update();
        
        if(reader->IsFilePolyData())
            data = reader->GetPolyDataOutput();
        else if(reader->IsFileUnstructuredGrid())
            data = reader->GetUnstructuredGridOutput();
        else
            mexErrMsgTxt("VTK internal dataformat is not supported.");
    }
    else if(extension == "vtp")
    {
        vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        if(!reader->CanReadFile(path.c_str()))
            mexErrMsgTxt("File cannot be read by vtkXMLPolyDataReader. Does it exist?");
        reader->SetFileName(path.c_str());
        reader->Update();
        data = reader->GetOutput();
    }
    else if(extension == "vtu")
    {
        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        if(!reader->CanReadFile(path.c_str()))
            mexErrMsgTxt("File cannot be read by vtkXMLUnstructuredGridReader. Does it exist?");
        reader->SetFileName(path.c_str());
        reader->Update();
        data = reader->GetOutput();
    }
    else if(extension == "ply")
    {
        vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
        if(!reader->CanReadFile(path.c_str()))
            mexErrMsgTxt("File cannot be read by vtkPLYReader. Does it exist?");
        reader->SetFileName(path.c_str());
        reader->Update();
        data = reader->GetOutput();
    }
    else if(extension == "stl")
    {
        vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
        reader->SetFileName(path.c_str());
        reader->Update();
        data = reader->GetOutput();
    }
    else if(extension == "obj")
    {
        vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
        reader->SetFileName(path.c_str());
        reader->Update();
        data = reader->GetOutput();
    }
//    else if(extension == "off")
//    {
//        vtkSmartPointer<vtkOFFReader> reader = vtkSmartPointer<vtkOFFReader>::New();
//        reader->SetFileName(path.c_str());
//        reader->Update();
//        data = reader->GetOutput();
//    }
    else
        mexErrMsgTxt((std::string("Unknown extension \"") + extension + "\".").c_str());
    
    ///// Process points and cells /////
    
    mwSize numPoints = data->GetNumberOfPoints();   
    mwSize numCells = data->GetNumberOfCells();
    
    vtkSmartPointer<vtkCellTypes> typesOfCells = vtkSmartPointer<vtkCellTypes>::New();
    data->GetCellTypes(typesOfCells);
    mwSize maxPointsPerCell = 0;
    for(int i = 0; i < typesOfCells->GetNumberOfTypes(); i++)
    {
        unsigned char cellType = typesOfCells->GetCellType(i);
        mwSize n = getNumberOfCellPoints(cellType);
        maxPointsPerCell = std::max(n, maxPointsPerCell);
    }
    
    mxArray* points    = mxCreateDoubleMatrix(numPoints, 3, mxREAL);
    mxArray* cells     = mxCreateNumericMatrix(numCells, maxPointsPerCell, mxINT32_CLASS, mxREAL);
    mxArray* cellTypes = mxCreateNumericMatrix(numCells, 1, mxUINT8_CLASS, mxREAL);
    
    printVerbose("Processing points...\n");
    
    for(vtkIdType i = 0; i < numPoints; i++)
    {
        for(vtkIdType j = 0; j < 3; j++)
            mxGetPr(points)[numPoints*j+i] = data->GetPoint(i)[j];
    }

    printVerbose("Processing cells...\n");
    
    for(vtkIdType i = 0; i < numCells; i++)
    {
        vtkCell* cell = data->GetCell(i);
        
        vtkIdType numPoints = cell->GetNumberOfPoints();
        for(int j = 0; j < numPoints; j++)
            ((int32_t*)mxGetData(cells))[numCells*j+i] = (int32_t)cell->GetPointId(j)+1;
        
        ((uint8_t*)mxGetData(cellTypes))[i] = cell->GetCellType();
    }
    
    ///// Process point data /////
    
    printVerbose("Processing point data...\n");
    
    vtkPointData* pointData = data->GetPointData();
    int numPointDataArrays = pointData->GetNumberOfArrays();
    
    const char** pointDataArrayNames = new const char*[numPointDataArrays];
    for(int i = 0; i < numPointDataArrays; i++)
        pointDataArrayNames[i] = pointData->GetArrayName(i);
    
    mxArray* pointDataStruct = mxCreateStructMatrix(1, 1, numPointDataArrays, pointDataArrayNames);
    for(int i = 0; i < numPointDataArrays; i++)
    {
        mxArray* mxPointDataArray = vtkDataArrayToMxArray(pointData->GetArray(i));
        mxSetFieldByNumber(pointDataStruct, 0, i, mxPointDataArray);
    }
    
    ///// Process cell data /////
    
    printVerbose("Processing cell data...\n");
    
    vtkCellData* cellData = data->GetCellData();
    int numCellDataArrays = cellData->GetNumberOfArrays();
    
    const char** cellDataArrayNames = new const char*[numCellDataArrays];
    for(int i = 0; i < numCellDataArrays; i++)
        cellDataArrayNames[i] = cellData->GetArrayName(i);
    
    mxArray* cellDataStruct = mxCreateStructMatrix(1, 1, numCellDataArrays, cellDataArrayNames);
    for(int i = 0; i < numCellDataArrays; i++)
    {
        mxArray* mxCellDataArray = vtkDataArrayToMxArray(cellData->GetArray(i));
        mxSetFieldByNumber(cellDataStruct, 0, i, mxCellDataArray);
    }
    
    ///// Process field data /////
    
    printVerbose("Processing field data...\n");
    
    vtkFieldData* fieldData = data->GetFieldData();
    int numFieldDataArrays = fieldData->GetNumberOfArrays();

    const char** fieldDataArrayNames = new const char*[numFieldDataArrays];
    for(int i = 0; i < numFieldDataArrays; i++)
        fieldDataArrayNames[i] = fieldData->GetArrayName(i);
    
    mxArray* fieldDataStruct = mxCreateStructMatrix(1, 1, numFieldDataArrays, fieldDataArrayNames);
    for(int i = 0; i < numFieldDataArrays; i++)
    {
        mxArray* mxFieldDataArray = vtkDataArrayToMxArray(fieldData->GetArray(i));
        mxSetFieldByNumber(fieldDataStruct, 0, i, mxFieldDataArray);
    }
    
    ///// Assemble output struct /////
    
    printVerbose("Assembling output struct...\n");
    
    const char* fieldnames[6] = {"points", "cells", "cellTypes", "pointData", "cellData", "fieldData"};
    mxArray* outStruct = mxCreateStructMatrix(1, 1, 6, fieldnames);
    mxSetField(outStruct, 0, "points", points);
    mxSetField(outStruct, 0, "cells", cells);
    mxSetField(outStruct, 0, "cellTypes", cellTypes);
    mxSetField(outStruct, 0, "pointData", pointDataStruct);
    mxSetField(outStruct, 0, "cellData", cellDataStruct);
    mxSetField(outStruct, 0, "fieldData", fieldDataStruct);
    
    plhs[0] = outStruct;
}
