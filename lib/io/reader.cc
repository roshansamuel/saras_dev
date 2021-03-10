/********************************************************************************************************************************************
 * Saras
 * 
 * Copyright (C) 2019, Mahendra K. Verma
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     1. Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *     2. Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *     3. Neither the name of the copyright holder nor the
 *        names of its contributors may be used to endorse or promote products
 *        derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ********************************************************************************************************************************************
 */
/*! \file reader.cc
 *
 *  \brief Definitions for functions of class reader
 *  \sa reader.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "reader.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the reader class
 *
 *          The constructor initializes the variables and parameters for parallel file reading through HDF5
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   wField is a vector of fields to be read into
 ********************************************************************************************************************************************
 */
reader::reader(const grid &mesh, std::vector<field> &rFields): mesh(mesh), rFields(rFields) {
    /** Initialize the common global and local limits for file writing */
    initLimits();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to initialize the global and local limits for setting file views
 *
 *          All the necessary limits of the local arrays with respect to the global array for setting the
 *          dataspace views for HDF5 are appropriately set here.
 ********************************************************************************************************************************************
 */
void reader::initLimits() {
    hid_t sDSpace;
    hid_t tDSpace;

    herr_t status;

    blitz::TinyVector<int, 3> gloSize, sdStart, locSize;

#ifdef PLANAR
    hsize_t dimsf[2];           /* dataset dimensions */
    hsize_t offset[2];          /* offset of hyperslab */
#else
    hsize_t dimsf[3];           /* dataset dimensions */
    hsize_t offset[3];          /* offset of hyperslab */
#endif

    for (unsigned int i=0; i < rFields.size(); i++) {
        gloSize = mesh.globalSize;
        locSize = mesh.coreSize;

#ifdef PLANAR
        gloSize(1) = 1;
        locSize(1) = 1;
#endif

        // Since only the last rank along X and Y directions include the extra point (shared across processors), subArrayStarts are same for all ranks
        sdStart = mesh.subarrayStarts;

        // Create a dataspace representing the full limits of the array - this is the target dataspace
#ifdef PLANAR
        dimsf[0] = locSize(0);
        dimsf[1] = locSize(2);
        tDSpace = H5Screate_simple(2, dimsf, NULL);
#else
        dimsf[0] = locSize(0);
        dimsf[1] = locSize(1);
        dimsf[2] = locSize(2);
        tDSpace = H5Screate_simple(3, dimsf, NULL);
#endif

        // Modify the view of the *target* dataspace by using a hyperslab - *this view will be used to read from memory*

        // ************* NOTE: Check if the hyperslab view must be changed to take care of the common point at the the MPI decomposed sub-domain boundaries ****************//
#ifdef PLANAR
        dimsf[0] = locSize(0);
        dimsf[1] = locSize(2);
        offset[0] = 0;
        offset[1] = 0;
#else
        dimsf[0] = locSize(0);
        dimsf[1] = locSize(1);
        dimsf[2] = locSize(2);
        offset[0] = 0;
        offset[1] = 0;
        offset[2] = 0;
#endif

        status = H5Sselect_hyperslab(tDSpace, H5S_SELECT_SET, offset, NULL, dimsf, NULL);
        if (status) {
            if (mesh.rankData.rank == 0) {
                std::cout << "Error in creating hyperslab while writing data. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }

        // Create a dataspace representing the full limits of the global array - i.e. the dataspace for input file
#ifdef PLANAR
        dimsf[0] = gloSize(0);
        dimsf[1] = gloSize(2);
        sDSpace = H5Screate_simple(2, dimsf, NULL);
#else
        dimsf[0] = gloSize(0);
        dimsf[1] = gloSize(1);
        dimsf[2] = gloSize(2);
        sDSpace = H5Screate_simple(3, dimsf, NULL);
#endif

        // Modify the view of the *source* dataspace by using a hyperslab according to its position in the global file dataspace
#ifdef PLANAR
        dimsf[0] = locSize(0);
        dimsf[1] = locSize(2);
        offset[0] = sdStart[0];
        offset[1] = sdStart[2];
#else
        dimsf[0] = locSize(0);
        dimsf[1] = locSize(1);
        dimsf[2] = locSize(2);
        offset[0] = sdStart[0];
        offset[1] = sdStart[1];
        offset[2] = sdStart[2];
#endif

        status = H5Sselect_hyperslab(sDSpace, H5S_SELECT_SET, offset, NULL, dimsf, NULL);
        if (status) {
            if (mesh.rankData.rank == 0) {
                std::cout << "Error in creating hyperslab while writing data. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }

        localSize.push_back(locSize);
        sourceDSpace.push_back(sDSpace);
        targetDSpace.push_back(tDSpace);
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to read files in HDF5 format in parallel
 *
 *          It opens a file in the output folder and all the processors read in parallel from the file
 *
 ********************************************************************************************************************************************
 */
real reader::readData() {
    hid_t plist_id;

    hid_t fileHandle;

    hid_t dataSet;

    herr_t status;

    real time;

    // Create a property list for collectively opening a file by all processors
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    // First create a file handle with the path to the input file
    H5E_BEGIN_TRY {
        fileHandle = H5Fopen("output/restartFile.h5", H5F_ACC_RDONLY, plist_id);
    } H5E_END_TRY;

    // Abort if file doesn't exist
    if (fileHandle < 0) {
        if (mesh.rankData.rank == 0) {
            std::cout << "ERROR: Restart flag is true, but could not open restart file. Aborting" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }

    // Close the property list for later reuse
    H5Pclose(plist_id);

    // Check the restart file for consistency with input parameters
    restartCheck(fileHandle);

    // Read the scalar value containing the time from the restart file
    hid_t timeDSpace = H5Screate(H5S_SCALAR);
    dataSet = H5Dopen2(fileHandle, "Time", H5P_DEFAULT);
    status = H5Dread(dataSet, H5T_NATIVE_REAL, timeDSpace, timeDSpace, H5P_DEFAULT, &time);

    // Close dataset for future use and dataspace for clearing resources
    H5Dclose(dataSet);
    H5Sclose(timeDSpace);

    // Create a property list to use collective data read
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    for (unsigned int i=0; i < rFields.size(); i++) {
#ifdef PLANAR
        fieldData.resize(blitz::TinyVector<int, 2>(localSize[i](0), localSize[i](2)));
#else
        fieldData.resize(localSize[i]);
#endif

        // Create the dataset *for the array in memory*, linking it to the file handle.
        // Correspondingly, it will use the *core* dataspace, as only the core has to be written excluding the pads
        dataSet = H5Dopen2(fileHandle, rFields[i].fieldName.c_str(), H5P_DEFAULT);

        // Write the dataset. Most important thing to note is that the 3rd and 4th arguments represent the *source* and *destination* dataspaces.
        // The source here is the sourceDSpace pointing to the file. Note that its view has been adjusted using hyperslab.
        // The destination is the targetDSpace. Though the targetDSpace is smaller than the sourceDSpace,
        // only the appropriate hyperslab within the sourceDSpace is transferred to the destination.

        // Note that the targetDSpace and sourceDSpace have switched positions
        // This is another point where the reader differs from the writer
        status = H5Dread(dataSet, H5T_NATIVE_REAL, targetDSpace[i], sourceDSpace[i], plist_id, fieldData.dataFirst());
        if (status) {
            if (mesh.rankData.rank == 0) {
                std::cout << "Error in reading input from HDF file. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }

        //Read data
        copyData(rFields[i]);

        H5Dclose(dataSet);
    }

    // CLOSE/RELEASE RESOURCES
    H5Pclose(plist_id);
    H5Fclose(fileHandle);

    return time;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to copy data from blitz array without pads into solver variables
 *
 *          In order to simplify the file views while reading from disk to memory,
 *          the variables are copied into a local blitz array without the pads.
 *          These variables are then written into the solver variables supplied to the reader class.
 *
 ********************************************************************************************************************************************
 */
void reader::copyData(field &outField) {
#ifdef PLANAR
    for (int i=0; i < fieldData.shape()[0]; i++) {
        for (int k=0; k < fieldData.shape()[1]; k++) {
            outField.F(i, 0, k) = fieldData(i, k);
        }
    }
#else
    for (int i=0; i < fieldData.shape()[0]; i++) {
        for (int j=0; j < fieldData.shape()[1]; j++) {
            for (int k=0; k < fieldData.shape()[2]; k++) {
                outField.F(i, j, k) = fieldData(i, j, k);
            }
        }
    }
#endif
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to check compatibility of restart file with input parameters
 *
 *          Firstly, the dimensions of arrays in restart file must match with those in the parameters.
 *          Secondly, the array limits of the data are verified with the sizes specified in the input parameters.
 *
 ********************************************************************************************************************************************
 */
void reader::restartCheck(hid_t fHandle) {
    // Use the pressure data to get size of dataset
    hid_t pData = H5Dopen2(fHandle, "P", H5P_DEFAULT);
    hid_t pSpace = H5Dget_space(pData);
    const int ndims = H5Sget_simple_extent_ndims(pSpace);
#ifdef PLANAR
    if (ndims != 2) {
        if (mesh.rankData.rank == 0) {
            std::cout << "ERROR: Dimensionality of restart file conflicts with solver parameters. Aborting" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }
#else
    if (ndims != 3) {
        if (mesh.rankData.rank == 0) {
            std::cout << "ERROR: Dimensionality of restart file conflicts with solver parameters. Aborting" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }
#endif
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(pSpace, dims, NULL);

    // Abort if size of dataset doesn't match input parameters
    if (int(dims[0]) != mesh.globalSize(0)) {
        if (mesh.rankData.rank == 0) {
            std::cout << "ERROR: Array limits in restart file conflicts with solver parameters. Aborting" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }
#ifdef PLANAR
    if (int(dims[1]) != mesh.globalSize(2)) {
        if (mesh.rankData.rank == 0) {
            std::cout << "ERROR: Array limits in restart file conflicts with solver parameters. Aborting" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }
#else
    if (int(dims[1]) != mesh.globalSize(1)) {
        if (mesh.rankData.rank == 0) {
            std::cout << "ERROR: Array limits in restart file conflicts with solver parameters. Aborting" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }

    if (int(dims[2]) != mesh.globalSize(2)) {
        if (mesh.rankData.rank == 0) {
            std::cout << "ERROR: Array limits in restart file conflicts with solver parameters. Aborting" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }
#endif

    // Close dataset and dataspace
    H5Dclose(pData);
    H5Sclose(pSpace);
}

reader::~reader() { }
