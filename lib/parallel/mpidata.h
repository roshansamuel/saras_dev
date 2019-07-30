#ifndef MPIDATA_H
#define MPIDATA_H

#include <blitz/array.h>
#include <mpi.h>

#include "parallel.h"

class mpidata {
    private:
        /** MPI subarray datatype for the slice of data to be sent to the neighbouring sub-domain to left along x-direction */
        MPI_Datatype sendSubarrayX0;
        /** MPI subarray datatype for the slice of data to be sent to the neighbouring sub-domain to right along x-direction */
        MPI_Datatype sendSubarrayX1;
        /** MPI subarray datatype for the slice of data to be sent to the neighbouring sub-domain to left along y-direction (front side) */
        MPI_Datatype sendSubarrayY0;
        /** MPI subarray datatype for the slice of data to be sent to the neighbouring sub-domain to right along y-direction (rear side) */
        MPI_Datatype sendSubarrayY1;

        /** MPI subarray datatype for the slice of data to be received from the neighbouring sub-domain to left along x-direction */
        MPI_Datatype recvSubarrayX0;
        /** MPI subarray datatype for the slice of data to be received from the neighbouring sub-domain to right along x-direction */
        MPI_Datatype recvSubarrayX1;
        /** MPI subarray datatype for the slice of data to be received from the neighbouring sub-domain to left along y-direction (front side) */
        MPI_Datatype recvSubarrayY0;
        /** MPI subarray datatype for the slice of data to be received from the neighbouring sub-domain to right along y-direction (rear side) */
        MPI_Datatype recvSubarrayY1;

        /** An array of MPI_Request data-types necessary for obtaining output from the non-blocking receive MPI_Irecv in the syncData function. */
        blitz::Array<MPI_Request, 1> recvRequest;

        /** An array of MPI_Status data-types necessary for obtaining output from the non-blocking receive MPI_Irecv in the syncData function. */
        blitz::Array<MPI_Status, 1> recvStatus;

        /** Array of values of the data field which needs to be synchronised across processors */
        blitz::Array<double, 3> dataField;

    public:
        /** A const reference to the global variables stored in the parallel class to access rank data */
        const parallel &rankData;

        mpidata(blitz::Array<double, 3> inputArray, const parallel &parallelData);

        void createSubarrays(const blitz::TinyVector<int, 3> globSize,
                             const blitz::TinyVector<int, 3> coreSize,
                             const blitz::TinyVector<int, 3> padWidth,
                             const bool xStag, const bool yStag);

        void syncData();
};

/**
 ********************************************************************************************************************************************
 *  \class mpidata mpidata.h "lib/mpidata.h"
 *  \brief Class to store MPI derived datatypes for individual arrays.
 *
 *  Since the solver uses staggered and collocated grids, the data of its variables are stored in arrays of different limits
 *  depending on whether the variable is staggered or not in each direction.
 *  As a result, the limits of the sub-arrays to be sent across inter-processor boundaries is different for different arrays.
 *  Hence the <B>mpidata</B> class contains MPI_SUBARRAY derived datatypes to be initialized along with different fields in order
 *  to store their sub-arrays for inter-processor communication.
 ********************************************************************************************************************************************
 */

#endif
