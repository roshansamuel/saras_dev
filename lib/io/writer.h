#ifndef WRITER_H
#define WRITER_H

#include <blitz/array.h>
#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "field.h"
#include "grid.h"
#include "hdf5.h"
#include "mpi.h"

class writer {
    private:
        const grid &mesh;

        std::vector<field> &wFields;

#ifdef PLANAR
        blitz::Array<double, 2> fieldData;
#else
        blitz::Array<double, 3> fieldData;
#endif

        std::vector<hid_t> sourceDSpace, targetDSpace;

        std::vector< blitz::TinyVector<int, 3> > localSize;

        void outputCheck();

        void initLimits();

        void copyData(field &outField);

    public:
        writer(const grid &mesh, std::vector<field> &wFields);

        void writeData(double time);
        
        ~writer();
};

/**
 ********************************************************************************************************************************************
 *  \class writer writer.h "lib/writer.h"
 *  \brief Class for all the global variables and functions related to writing output data of the solver.
 *
 *  The computational data from the solver is written in HDF5 format in a .h5 file.
 *  The class allows for both collocated and staggered grid data to be written in separate output files.
 ********************************************************************************************************************************************
 */

#endif
