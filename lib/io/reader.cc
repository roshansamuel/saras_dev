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
        if (not rFields[i].xStag) {
            gloSize(0) -= 1;
        }

#ifndef PLANAR
        if (not rFields[i].yStag) {
            gloSize(1) -= 1;
        }
#else
        gloSize(1) = 1;
#endif

        if (not rFields[i].zStag) {
            gloSize(2) -= 1;
        }

        locSize = mesh.collocCoreSize;
        if (rFields[i].xStag) {
            // Though the last point (which is shared across 2 processors) was excluded for some subdomains in the writer class, while reading, these overlapping points have to be considered
            // This is one point where reader differs from writer
            locSize(0) = mesh.staggrCoreSize(0);
        }

#ifndef PLANAR
        if (rFields[i].yStag) {
            // As with X direction, all subdomains include both the boundary points (which may be shared across 2 processors)
            locSize(1) = mesh.staggrCoreSize(1);
        }
#else
        locSize(1) = 1;
#endif

        if (rFields[i].zStag) {
            locSize(2) = mesh.staggrCoreSize(2);
        }

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
            exit(0);
        }

        localSize.push_back(locSize);
        sourceDSpace.push_back(sDSpace);
        targetDSpace.push_back(tDSpace);
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to scan files in the output folder and identify the last written file
 *
 *          This function uses dirent.h for POSIX systems to list all files in the directory.
 *          This part of the code may restrict cross-platform capabilities
 *
 ********************************************************************************************************************************************
 */
double reader::getLastFile() {
    DIR *dir;
    struct stat info;
    struct dirent *dirp;
    std::vector<double> timeArray;
    std::vector<double>::iterator pos;

    // Check if output directory exists
    if (stat("output", &info) != 0) {
        if (mesh.rankData.rank == 0) {
            std::cout << "Error output directory does not exist for solver in restart mode. Aborting" << std::endl;
        }
        exit(0);
    }

    dir = opendir("output/");

    while ((dirp = readdir(dir))!=NULL) {
        std::string fName(dirp->d_name);

        if (fName.find("Soln_") != std::string::npos) {
            double timeVal = std::atof(fName.substr(5, 9).c_str());
            timeArray.push_back(timeVal);
        }
    }

    if (timeArray.empty()) {
        if (mesh.rankData.rank == 0) {
            std::cout << "Output directory has no solution files (Soln_XXXX.XXXX.h5 files) to read in restart mode. Aborting" << std::endl;
        }
        exit(0);
    }

    pos = std::max_element(timeArray.begin(), timeArray.end());

    return *pos;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to read files in HDF5 format in parallel
 *
 *          It opens a file in the output folder and all the processors read in parallel from the file
 *
 ********************************************************************************************************************************************
 */
double reader::readData() {
    hid_t plist_id;

    hid_t fileHandle;

    hid_t dataSet;

    herr_t status;

    std::ostringstream constFile;

    char* fileName;

    double time;

    time = getLastFile();

    // Create a property list for collectively opening a file by all processors
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    // Generate the filename corresponding to the solution file
    fileName = new char[100];
    constFile.str(std::string());
    constFile << "output/Soln_" << std::fixed << std::setfill('0') << std::setw(9) << std::setprecision(4) << time << ".h5";
    strcpy(fileName, constFile.str().c_str());

    // First create a file handle with the path to the input file
    fileHandle = H5Fopen(fileName, H5F_ACC_RDONLY, plist_id);

    // Close the property list for later reuse
    H5Pclose(plist_id);

    // Create a property list to use collective data write
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
        status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, targetDSpace[i], sourceDSpace[i], plist_id, fieldData.dataFirst());
        if (status) {
            if (mesh.rankData.rank == 0) {
                std::cout << "Error in reading input from HDF file. Aborting" << std::endl;
            }
            exit(0);
        }

        //Read data
        copyData(rFields[i]);

        H5Dclose(dataSet);
    }

    // CLOSE/RELEASE RESOURCES
    H5Pclose(plist_id);
    H5Fclose(fileHandle);

    delete fileName;

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

reader::~reader() { }
