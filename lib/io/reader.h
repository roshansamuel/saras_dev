#ifndef READER_H
#define READER_H

#include <blitz/array.h>
#include <sys/stat.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <dirent.h>

#include "field.h"
#include "grid.h"
#include "hdf5.h"

class reader {
    private:
        const grid &mesh;

        std::vector<field> &rFields;

#ifdef PLANAR
        blitz::Array<double, 2> fieldData;
#else
        blitz::Array<double, 3> fieldData;
#endif

        std::vector<hid_t> sourceDSpace, targetDSpace;

        std::vector< blitz::TinyVector<int, 3> > localSize;

        double getLastFile();

        void initLimits();

        void copyData(field &outField);

    public:
        reader(const grid &mesh, std::vector<field> &rFields);

        double readData();

        ~reader();
};

/**
 ********************************************************************************************************************************************
 *  \class reader reader.h "lib/reader.h"
 *  \brief Class for all the global variables and functions related to reading input data for the solver.
 *
 *  The computational data for the solver can be read from HDF5 file.
 *  The class allows for both collocated and staggered grid data to be read from separate input files.
 ********************************************************************************************************************************************
 */

#endif
