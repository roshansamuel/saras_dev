#include "writer.h"
#include <h5si.h>
/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the writer class
 *
 *          The constructor initializes the variables and parameters for parallel file writing through HDF5
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   wField is a vector of sfields to be written
 ********************************************************************************************************************************************
 */
writer::writer(const grid &mesh, std::vector<field> &wFields): mesh(mesh), wFields(wFields) {
    /** Initialize the common global and local limits for file writing */
    initLimits();

    /** Create output directory if it doesn't exist */
    outputCheck();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to initialize the global and local limits for setting file views
 *
 *          All the necessary limits of the local arrays with respect to the global array for setting the
 *          dataspace views for HDF5 are appropriately set here.
 ********************************************************************************************************************************************
 */
void writer::initLimits() {
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

    for (unsigned int i=0; i < wFields.size(); i++) {
        gloSize = mesh.globalSize;
        if (not wFields[i].xStag) {
            gloSize(0) -= 1;
        }

#ifndef PLANAR
        if (not wFields[i].yStag) {
            gloSize(1) -= 1;
        }
#else
        gloSize(1) = 1;
#endif

        if (not wFields[i].zStag) {
            gloSize(2) -= 1;
        }

        locSize = mesh.collocCoreSize;
        if (wFields[i].xStag) {
            // All subdomains exclude the last point (which is shared across 2 processors), except those at the last rank along X direction, thereby capturing the full domain
            locSize(0) = mesh.staggrCoreSize(0) - 1;
            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                locSize(0) += 1;
            }
        }

#ifndef PLANAR
        if (wFields[i].yStag) {
            // As with X direction, all subdomains exclude the last point (which is shared across 2 processors), except those at the last rank along Y direction
            locSize(1) = mesh.staggrCoreSize(1) - 1;
            if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
                locSize(1) += 1;
            }
        }
#else
        locSize(1) = 1;
#endif

        if (wFields[i].zStag) {
            locSize(2) = mesh.staggrCoreSize(2);
        }

        // Since only the last rank along X and Y directions include the extra point (shared across processors), subArrayStarts are same for all ranks
        sdStart = mesh.subarrayStarts;

        // Create a dataspace representing the full limits of the array - this is the source dataspace
#ifdef PLANAR
        dimsf[0] = locSize(0);
        dimsf[1] = locSize(2);
        sDSpace = H5Screate_simple(2, dimsf, NULL);
#else
        dimsf[0] = locSize(0);
        dimsf[1] = locSize(1);
        dimsf[2] = locSize(2);
        sDSpace = H5Screate_simple(3, dimsf, NULL);
#endif

        // Modify the view of the *source* dataspace by using a hyperslab - *this view will be used to read from memory*

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

        status = H5Sselect_hyperslab(sDSpace, H5S_SELECT_SET, offset, NULL, dimsf, NULL);
        if (status) {
            if (mesh.rankData.rank == 0) {
                std::cout << "Error in creating hyperslab while writing data. Aborting" << std::endl;
            }
            exit(0);
        }

        // Create a dataspace representing the full limits of the global array - i.e. the dataspace for output file
#ifdef PLANAR
        dimsf[0] = gloSize(0);
        dimsf[1] = gloSize(2);
        tDSpace = H5Screate_simple(2, dimsf, NULL);
#else
        dimsf[0] = gloSize(0);
        dimsf[1] = gloSize(1);
        dimsf[2] = gloSize(2);
        tDSpace = H5Screate_simple(3, dimsf, NULL);
#endif

        // Modify the view of the *target* dataspace by using a hyperslab according to its position in the global file dataspace
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

        status = H5Sselect_hyperslab(tDSpace, H5S_SELECT_SET, offset, NULL, dimsf, NULL);
        if (status) {
            if (mesh.rankData.rank == 0) {
                std::cout << "Error in creating hyperslab while writing data. Aborting" << std::endl;
            }
            exit(0);
        }

        localSize.push_back(locSize);
        sourceDSpace.push_back(sDSpace);
        targetDSpace.push_back(tDSpace);
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create output folder if it does not exist
 *
 *          The function uses methods available in the sys/stat.h header to create an output folder if it does not exist
 *          Without this, the code used to give segmentation fault while attempting to write solution files
 ********************************************************************************************************************************************
 */
void writer::outputCheck() {
    struct stat info;
    int createStatus;

    if (mesh.rankData.rank == 0) {
        // Check if output directory exists
        if (stat("output", &info) != 0) {
            createStatus = mkdir("output", S_IRWXU | S_IRWXG);

            // Raise error if the filesystem is read-only or something
            if (createStatus) {
                std::cout << "Error in while attempting to create output directory. Aborting" << std::endl;
                exit(0);
            }
        }
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to write output in HDF5 format in parallel
 *
 *          It opens a file in the output folder and all the processors write in parallel into the file
 *
 * \param   time is a double precision value containing the time to be used for naming the file
 ********************************************************************************************************************************************
 */

void writer::writeData(double time) {
    hid_t plist_id;

    hid_t fileHandle;

    hid_t dataSet;

    herr_t status;

    std::ostringstream constFile;

    char* fileName;

    // Create a property list for collectively opening a file by all processors
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    // Generate the filename corresponding to the solution file
    fileName = new char[100];
    constFile.str(std::string());
    constFile << "output/Soln_" << std::fixed << std::setfill('0') << std::setw(9) << std::setprecision(4) << time << ".h5";
    strcpy(fileName, constFile.str().c_str());

    // First create a file handle with the path to the output file
    fileHandle = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    // Close the property list for later reuse
    H5Pclose(plist_id);

    // Create a property list to use collective data write
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    for (unsigned int i=0; i < wFields.size(); i++) {
#ifdef PLANAR
        fieldData.resize(blitz::TinyVector<int, 2>(localSize[i](0), localSize[i](2)));
#else
        fieldData.resize(localSize[i]);
#endif

        //Write data
        copyData(wFields[i]);

        // Create the dataset *for the file*, linking it to the file handle.
        // Correspondingly, it will use the *core* dataspace, as only the core has to be written excluding the pads
        dataSet = H5Dcreate2(fileHandle, wFields[i].fieldName.c_str(), H5T_NATIVE_DOUBLE, targetDSpace[i], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Write the dataset. Most important thing to note is that the 3rd and 4th arguments represent the *source* and *destination* dataspaces.
        // The source here is the sourceDSpace pointing to the memory buffer. Note that its view has been adjusted using hyperslab.
        // The destination is the targetDSpace. Though the targetDSpace is smaller than the sourceDSpace,
        // only the appropriate hyperslab within the sourceDSpace is transferred to the destination.

        status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, sourceDSpace[i], targetDSpace[i], plist_id, fieldData.dataFirst());
        if (status) {
            if (mesh.rankData.rank == 0) {
                std::cout << "Error in writing output to HDF file. Aborting" << std::endl;
            }
            exit(0);
        }

        H5Dclose(dataSet);
    }

    // CLOSE/RELEASE RESOURCES
    H5Pclose(plist_id);
    H5Fclose(fileHandle);

    delete fileName;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to copy data to a blitz array without pads
 *
 *          In order to simplify the file views while writing to disk from memory,
 *          the variables are copied into a local blitz array without the pads.
 *
 ********************************************************************************************************************************************
 */
void writer::copyData(field &outField) {
#ifdef PLANAR
    for (int i=0; i < fieldData.shape()[0]; i++) {
        for (int k=0; k < fieldData.shape()[1]; k++) {
            fieldData(i, k) = outField.F(i, 0, k);
        }
    }
#else
    for (int i=0; i < fieldData.shape()[0]; i++) {
        for (int j=0; j < fieldData.shape()[1]; j++) {
            for (int k=0; k < fieldData.shape()[2]; k++) {
                fieldData(i, j, k) = outField.F(i, j, k);
            }
        }
    }
#endif
}

writer::~writer() { }
