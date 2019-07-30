#include "parallel.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the parallel class
 *
 *          The initializing functions of MPI are called in order to get the total number of processes spawned, and
 *          the rank of each process.
 *          The xRank and yRank of each process are calculated and assigned.
 *          Finally, the ranks of neighbouring processes are found and stored in an array for use in MPI communications
 *
 * \param   iDat is a const reference to the global data contained in the parser class
 ********************************************************************************************************************************************
 */
parallel::parallel(const parser &iDat): npX(iDat.npX), npY(iDat.npY) {
    // GET EACH PROCESSES' RANK AND TOTAL NUMBER OF PROCESSES
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);

    // ABORT IF THE NUMBER OF PROCESSORS IN EACH DIRECTION SPECIFIED IN INPUT DOES NOT MATCH WITH AVAILABLE CORES
    if (npX*npY != nProc) {
        if (rank == 0) {
            std::cout << "ERROR: Number of processors specified in input file does not match. Aborting" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }

    // ASSIGN EACH PROCESSES' xRank AND yRank
    assignRanks();

    // GET AND STORE THE RANKS OF ALL NEIGHBOURING PROCESSES FOR FUTURE DATA TRANSFER
    getNeighbours();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to assign the xRank and yRank for each sub-domain according to their global rank
 *
 *          It uses the number of sub-divisions prescribed in each direction, i.e. \ref npX and \ref npY to calculate the
 *          xRank and yRank appropriately.
 ********************************************************************************************************************************************
 */
inline void parallel::assignRanks() {
    xRank = rank % npX;
    yRank = rank / npX;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to get the ranks of each neighbouring sub-domain which shares a face with the given sub-domain
 *
 *          Since the solver uses pencil decomposition, it locates the ranks of a maximum of 4 neighbouring sub-domains.
 ********************************************************************************************************************************************
 */
void parallel::getNeighbours() {
    // EACH PROCESS HAS 4 NEIGHBOURS CORRESPONDING TO THE 4 FACES OF EACH CUBICAL SUB-DOMAIN
    nearRanks.resize(4);

    // EACH PROCESS IS ASSUMED TO HAVE NO NEIGHBOURS INITIALLY
    nearRanks = MPI_PROC_NULL;

    // INITIAL NEIGHBOUR ASSIGNMENTS ARE DONE ASSUMING PERIODIC DOMAIN
    // ALONG X/XI DIRECTION
    nearRanks(0) = findRank(xRank - 1, yRank);
    nearRanks(1) = findRank(xRank + 1, yRank);

    // ALONG Y/ETA DIRECTION
#ifndef PLANAR
    nearRanks(2) = findRank(xRank, yRank - 1);
    nearRanks(3) = findRank(xRank, yRank + 1);
#endif
}
