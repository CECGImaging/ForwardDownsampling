/*******************************************************************************
 * writeVTK.cxx                                                                *
 * MATLAB extension to write VTK files (legacy: vtk or XML: vtp, vtu)          *
 * Also supported are: ply, stl, obj, off (only geometry, no data arrays)      *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
 * Authors: Nils Daub, Steffen Schuler                                         *
 *          Institute of Biomedical Engineering                                *
 *          Karlsruhe Institute of Technology                                  *
 *          www.ibt.kit.edu                                                    *
 * Date:    November 2017                                                      *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
 * Before first usage:                                                         *
 * Compile the MEX using the script mexCompile_writeVTK.m                      *
 *                                                                             *
 * Usage inside MATLAB:                                                        *
 * writeVTK(struct inStruct, char filename, [bool verbose], [char dataMode]);  *
 *                                                                             *
 * inStruct must have the following fields:                                    *
 *     points:    [numPoints x 3 double]                                       *
 *     cells:     [numCells x maxNumPointsPerCell int32]                       *
 *     cellTypes: [numCells x 1 uint8]                                         *
 *     pointData: [1 x 1 struct of point data arrays] (optional)               *
 *     cellData:  [1 x 1 struct of cell data arrays] (optional)                *
 *     fieldData: [1 x 1 struct of field data arrays] (optional)               *
 *                                                                             *
 * verbose is false on default                                                 *
 *                                                                             *
 * dataMode can be chosen from {'ascii', 'binary'} for vtk, ply, stl           *
 *                         and {'ascii', 'binary', 'appended'} for vtp, vtu    *
 *                         and {'ascii'} for obj, off                          *
 *******************************************************************************/

// general
#include <vtkSmartPointer.h>
#include <vtkErrorCode.h>

// legacy format writer
#include <vtkDataSetWriter.h>

// XML format writers
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

// writers for other common polygon mesh formats
#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>
#include "additionalWriters/vtkOBJWriter.h"
#include "additionalWriters/vtkOFFWriter.h"

// data sets
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

// data set attributes
#include <vtkPointData.h>
#include <vtkCellData.h>

// data array types
#include <vtkCharArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkTypeInt8Array.h>
#include <vtkTypeUInt8Array.h>
#include <vtkTypeInt16Array.h>
#include <vtkTypeUInt16Array.h>
#include <vtkTypeInt32Array.h>
#include <vtkTypeUInt32Array.h>
#include <vtkTypeInt64Array.h>
#include <vtkTypeUInt64Array.h>

// mex
#include "mex.h"
#include "matrix.h"

// This function is adapted from vtkMatlabMexAdapter
vtkDataArray* GetVtkDataType(mxClassID cid)
{
    switch(cid)
    {
        case mxCHAR_CLASS:
            return(vtkCharArray::New());
        case mxLOGICAL_CLASS:
            return(vtkUnsignedShortArray::New());
        case mxDOUBLE_CLASS:
            return(vtkDoubleArray::New());
        case mxSINGLE_CLASS:
            return(vtkFloatArray::New());
        case mxINT8_CLASS:
            return(vtkTypeInt8Array::New());
        case mxUINT8_CLASS:
            return(vtkTypeUInt8Array::New());
        case mxINT16_CLASS:
            return(vtkTypeInt16Array::New());
        case mxUINT16_CLASS:
            return(vtkTypeUInt16Array::New());
        case mxINT32_CLASS:
            return(vtkTypeInt32Array::New());
        case mxUINT32_CLASS:
            return(vtkTypeUInt32Array::New());
        case mxINT64_CLASS:
            return(vtkTypeInt64Array::New());
        case mxUINT64_CLASS:
            return(vtkTypeUInt64Array::New());
        default:
            mexWarnMsgTxt((std::string("Data array type ") + std::to_string(cid) + " is not supported. Using double.").c_str());
            return(vtkDoubleArray::New());
    }
}

// This function is adapted from vtkMatlabMexAdapter
vtkDataArray* mxArrayToVtkDataArray(const mxArray* inArray)
{
    if(inArray == NULL)
    {
        mexErrMsgTxt("NULL input to mxArrayTovtkDataArray().");
        return(NULL);
    }
    if(mxGetNumberOfDimensions(inArray) > 2)
    {
        mexErrMsgTxt("Input to mxArrayTovtkDataArray() has more than two dimensions, cannot convert to vtkDataArray.");
        return(NULL);
    }
    if(mxIsCell(inArray))
    {
        mexErrMsgTxt("Input to mxArrayTovtkDataArray() is a Cell Array, cannot convert to vtkDataArray.");
        return(NULL);
    }
    if(mxIsSparse(inArray))
    {
        mexErrMsgTxt("Input to mxArrayTovtkDataArray() is a Sparse Array, cannot convert to vtkDataArray.");
        return(NULL);
    }
    
    size_t numRows = mxGetM(inArray);
    size_t numColumns = mxGetN(inArray);
    
    vtkDataArray* outArray = GetVtkDataType(mxGetClassID(inArray));
    
    size_t numBytes = mxGetElementSize(inArray);
    if(numBytes != outArray->GetElementComponentSize())
    {
        outArray->Delete();
        mexErrMsgTxt("Data size mismatch between Matlab and VTK.");
        return(NULL);
    }
    
    outArray->SetNumberOfTuples(numRows);
    outArray->SetNumberOfComponents((int)numColumns);
    
    double* tuple = (double*) mxMalloc(sizeof(double)*numColumns);
    unsigned char* source = (unsigned char*)mxGetData(inArray);
    unsigned char* dest;
    
    for(size_t i = 0; i < numRows; i++)
    {
        outArray->InsertTuple(i,tuple);
        for(size_t j = 0; j < numColumns; j++)
        {
            dest = (unsigned char*) outArray->GetVoidPointer(i*numColumns + j);
            for(size_t k = 0; k < numBytes; k++)
                dest[k] = source[j*(numRows*numBytes) + i*numBytes + k];
        }
    }
    mxFree(tuple);
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

int isPolyDataCell(uint8_t cellType)
{
    switch(cellType)
    {
        case VTK_VERTEX:
            return(true);
        case VTK_POLY_VERTEX:
            return(true);
        case VTK_LINE:
            return(true);
        case VTK_POLY_LINE:
            return(true);
        case VTK_TRIANGLE:
            return(true);
        case VTK_QUAD:
            return(true);
        case VTK_POLYGON:
            return(true);
        case VTK_TRIANGLE_STRIP:
            return(true);
        default:
            return(false);
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
    if(nrhs < 2)
        mexErrMsgTxt("Not enough input arguments. Syntax: writeVTK(struct inStruct, char filename, [bool verbose])");
    if(nrhs > 4)
        mexErrMsgTxt("Too many input arguments. Syntax: writeVTK(struct inStruct, char filename, [bool verbose])");
    if(nlhs > 0)
        mexErrMsgTxt("Too many output arguments. There is no return value. Syntax: writeVTK(struct inStruct, char filename, [bool verbose])");
    
    std::string dataModeString;
    if(nrhs >= 3)
    {
        if(nrhs == 4)
        {
            char* dataModeChar = mxArrayToString(prhs[3]);
            if(dataModeChar == NULL)
                mexErrMsgTxt("Data mode argument could not be read. Must be one of {'ascii', 'binary', 'appended'}.");
            dataModeString = dataModeChar;
        }
        verbose_ = mxGetScalar(prhs[2]);
    }
    else
        verbose_ = false;
    
    ///// Read input struct /////
    
    printVerbose("Reading input struct...\n");
    
    const mxArray* inStruct = prhs[0];
    if(!mxIsStruct(inStruct))
        mexErrMsgTxt("First argument must be a struct. Syntax: writeVTK(struct inStruct, char filename, [bool verbose])");
    
    mxArray* points = mxGetField(inStruct, 0, "points");
    if(points == NULL)
        mexErrMsgTxt("Field 'points' of input struct could not be read.");
    if(!mxIsDouble(points))
        mexErrMsgTxt("Field 'points' of input struct must be of class double. You could try points = double(points)");

    mxArray* cells = mxGetField(inStruct, 0, "cells");
    if(cells == NULL)
        mexErrMsgTxt("Field 'cells' of input struct could not be read.");
    if(!mxIsClass(cells, "int32"))
        mexErrMsgTxt("Field 'cells' of input struct must be of class int32. You could try cells = int32(cells)");
    
    mxArray* cellTypes = mxGetField(inStruct, 0, "cellTypes");
    if(cellTypes == NULL)
        mexErrMsgTxt("Field 'cellTypes' of input struct could not be read.");
    if(!mxIsClass(cellTypes, "uint8"))
        mexErrMsgTxt("Field 'cellTypes' of input struct must be of class uint8. You could try cellTypes = uint8(cellTypes)");
    
    mxArray* pointDataStruct = mxGetField(inStruct, 0, "pointData");
    mxArray* cellDataStruct = mxGetField(inStruct, 0, "cellData");
    mxArray* fieldDataStruct = mxGetField(inStruct, 0, "fieldData");
    
    ///// Check file format compatibility /////
    
    printVerbose("Checking file format compatibility...\n");
    
    std::string path(mxArrayToString(prhs[1]));
    long pos = path.find_last_of(".");
    std::string extension = path.substr(pos+1, path.size()-pos);
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    
    if(extension != "vtk" && extension != "vtp" && extension != "vtu" && extension != "ply" && extension != "stl" && extension != "obj" && extension != "off")
        mexErrMsgTxt(("File extension " + extension + " is not supported (must be vtk, vtp, vtu, ply, stl, obj or off).").c_str());
    
    bool requirePolyData = false;
    if(extension != "vtu" && extension != "vtk")
        requirePolyData = true;
    
    bool isPolyData = false;
    if(extension == "vtk" || requirePolyData)
    {
        // Check, whether there are only cells compatible with vtkPolyData
        unsigned long nonPolyDataCells = 0;
        for(size_t i = 0; i < mxGetNumberOfElements(cellTypes); i++)
        {
            uint8_t cellType = ((uint8_t*)mxGetData(cellTypes))[i];
            if(!isPolyDataCell(cellType))
                nonPolyDataCells++;
        }
        if(nonPolyDataCells == 0)
            isPolyData = true;
        else if(requirePolyData)
            mexErrMsgTxt(("vtkPolyData only supports linear surface cells (vertices, lines or triangles). " + std::to_string(nonPolyDataCells) + " incompatible cells found.").c_str());
    }
    
    ///// Process cells /////
    
    printVerbose("Processing cells...\n");
    
    mwSize numCells = mxGetM(cells);
    if(mxGetM(cellTypes)*mxGetN(cellTypes) != numCells)
        mexErrMsgTxt("The number of elements in 'cellTypes' must be equal to the number of rows in 'cells'.");
        
    mwSize maxNumCellPoints = mxGetN(cells);
    vtkIdType* pointIds = new vtkIdType[maxNumCellPoints];
    
    vtkPointSet* pointSet;
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    
    if(isPolyData)
    {
        polyData->Allocate(numCells);
        for(unsigned long i = 0; i < numCells; i++)
        {
            uint8_t cellType = ((uint8_t*)mxGetData(cellTypes))[i];
            int numCellPoints = getNumberOfCellPoints(cellType);
            for(unsigned long j = 0; j < numCellPoints; j++)
            {
                pointIds[j] = ((int32_t*)mxGetData(cells))[numCells*j+i]-1;
            }
            polyData->InsertNextCell(cellType, numCellPoints, pointIds);
        }
        pointSet = polyData;
    }
    else
    {
        unstructuredGrid->Allocate(numCells);
        for(unsigned long i = 0; i < numCells; i++)
        {
            uint8_t cellType = ((uint8_t*)mxGetData(cellTypes))[i];
            int numCellPoints = getNumberOfCellPoints(cellType);
            for(unsigned long j = 0; j < numCellPoints; j++)
            {
                pointIds[j] = ((int32_t*)mxGetData(cells))[numCells*j+i]-1;
            }
            unstructuredGrid->InsertNextCell(cellType, numCellPoints, pointIds);
        }
        pointSet = unstructuredGrid;
    }
    
    ///// Process points /////
    
    printVerbose("Processing points...\n");
    
    vtkSmartPointer<vtkPoints> vtkpoints = vtkSmartPointer<vtkPoints>::New();
    unsigned long numPoints = mxGetM(points);
    vtkpoints->SetNumberOfPoints(numPoints);
    for(unsigned long i = 0; i < numPoints; i++)
    {
        double pointCoords[3];
        for(int j = 0; j < 3; j++)
        {
            pointCoords[j] = mxGetPr(points)[numPoints*j+i];
        }
        vtkpoints->SetPoint(i, pointCoords);
    }
    vtkpoints->Modified();
    pointSet->SetPoints(vtkpoints);
    
    ///// Process data arrays /////
    
    if(extension == "vtk" || extension == "vtp" || extension == "vtu")
    {
        ///// Process point data /////
        
        if(pointDataStruct != NULL)
        {
            printVerbose("Processing point data...\n");
            
            int numPointDataArrays = mxGetNumberOfFields(pointDataStruct);
            for(int i = 0; i < numPointDataArrays; i++)
            {
                vtkDataArray* pointDataArray = mxArrayToVtkDataArray(mxGetFieldByNumber(pointDataStruct, 0, i));
                pointDataArray->SetName(mxGetFieldNameByNumber(pointDataStruct, i));
                pointSet->GetPointData()->AddArray(pointDataArray);
            }
        }
        
        ///// Process cell data /////
        
        if(cellDataStruct != NULL)
        {
            printVerbose("Processing cell data...\n");
            
            int numCellDataArrays = mxGetNumberOfFields(cellDataStruct);
            for(int i = 0; i < numCellDataArrays; i++)
            {
                vtkDataArray* cellDataArray = mxArrayToVtkDataArray(mxGetFieldByNumber(cellDataStruct, 0, i));
                cellDataArray->SetName(mxGetFieldNameByNumber(cellDataStruct, i));
                pointSet->GetCellData()->AddArray(cellDataArray);
            }
        }
        
        ///// Process field data /////
        
        if(fieldDataStruct != NULL)
        {
            printVerbose("Processing field data...\n");
            
            int numFieldDataArrays = mxGetNumberOfFields(fieldDataStruct);
            for(int i = 0; i < numFieldDataArrays; i++)
            {
                vtkDataArray* fieldDataArray = mxArrayToVtkDataArray(mxGetFieldByNumber(fieldDataStruct, 0, i));
                fieldDataArray->SetName(mxGetFieldNameByNumber(fieldDataStruct, i));
                pointSet->GetFieldData()->AddArray(fieldDataArray);
            }
        }
    }
    
    ///// Write file /////
    
    printVerbose("Writing file...\n");
    
    int dataMode = 0;
    if(!dataModeString.empty())
    {
        if(dataModeString == "ascii")
            dataMode = 1;
        else if(dataModeString == "binary")
            dataMode = 2;
        else if(dataModeString == "appended")
        {
            if(extension == "vtp" || extension == "vtu")
                dataMode = 3;
            else
                mexWarnMsgTxt("Data mode 'appended' only supported for vtp and vtu. Using VTK default.");
        }
        else
            mexWarnMsgTxt(("Unknown data mode '" + std::string(dataModeString) + "'. Using VTK default.").c_str());
    }
    
    if(extension == "vtk")
    {
        vtkSmartPointer<vtkDataSetWriter> writer = vtkSmartPointer<vtkDataSetWriter>::New();
        writer->SetInputData(pointSet);
        writer->SetFileName(path.c_str());
        if(dataMode)
            writer->SetFileType(dataMode);
        writer->Write();
    }
    else if(extension == "vtp")
    {
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetInputData(pointSet);
        writer->SetFileName(path.c_str());
        if(dataMode)
            writer->SetDataMode(dataMode-1);
        writer->Write();
    }
    else if(extension == "vtu")
    {
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetInputData(pointSet);
        writer->SetFileName(path.c_str());
        if(dataMode)
            writer->SetDataMode(dataMode-1);
        writer->Write();
    }
    else if(extension == "ply")
    {
        vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
        writer->SetInputData(pointSet);
        writer->SetFileName(path.c_str());
        if(dataMode)
            writer->SetFileType(dataMode);
        writer->Write();
    }
    else if(extension == "stl")
    {
        vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
        writer->SetInputData(pointSet);
        writer->SetFileName(path.c_str());
        if(dataMode)
            writer->SetFileType(dataMode);
        writer->Write();
    }
    else if(extension == "obj")
    {
        vtkSmartPointer<vtkOBJWriter> writer = vtkSmartPointer<vtkOBJWriter>::New();
        writer->SetInputData(pointSet);
        writer->SetFileName(path.c_str());
        if(dataMode > 1)
            mexWarnMsgTxt("Data mode is always 'ascii' for obj.");
        writer->Write();
    }
    else if(extension == "off")
    {
        vtkSmartPointer<vtkOFFWriter> writer = vtkSmartPointer<vtkOFFWriter>::New();
        writer->SetInputData(pointSet);
        writer->SetFileName(path.c_str());
        if(dataMode > 1)
            mexWarnMsgTxt("Data mode is always 'ascii' for off.");
        writer->Write();
    }
}
