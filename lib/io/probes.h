#ifndef PROBES_H
#define PROBES_H

#include <vector>
#include <fstream>
#include <cstddef>
#include <blitz/array.h>
#include <mpi.h>

#include "field.h"
#include "grid.h"

typedef struct dataStruct {
    int x, y, z;

    // Up to a maximum of 10 field variables can be probed presently
    double probeData[10];

    // Default struct constructor for initializing all elements to 0
    dataStruct(): x(0), y(0), z(0) { for (int i=0; i<10; i++) probeData[i] = 0.0; }
} dataStruct;

class probes {
    private:
        const grid &mesh;

        std::vector<field> &pFields;

        const unsigned int numFields;

        std::vector<blitz::TinyVector<int, 3> > globalProbes, localProbes;

        std::ofstream probeFile;

        MPI_Datatype mpiStruct;

        void getData(dataStruct *outData);
        void gatherData(dataStruct *outData);

        void createMPIStruct();

        void placeProbes();

    public:
        probes(const grid &mesh, std::vector<field> &pFields);

        void probeData(double time);

        ~probes();
};

/**
 ********************************************************************************************************************************************
 *  \class probes probes.h "lib/probes.h"
 *  \brief Handles the writing of data from probes placed in the domain
 *
 *  The class places the probes in probesList provided by user through the parser class.
 *  It also provides an interface to the solver to read data from the probes.
 ********************************************************************************************************************************************
 */

#endif
