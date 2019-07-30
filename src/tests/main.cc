#include "unittest.h"
#include "alltests.h"

int rootRank;

int main() {
    // INITIALIZE MPI
    MPI_Init(NULL, NULL);

    int numCases = 0;
    blitz::Array<blitz::TinyVector<int, 4>, 1> testParams;

    // ALL PROCESSES READ THE INPUT PARAMETERS
    parser inputData;

    // INITIALIZE PARALLELIZATION DATA
    parallel mpi(inputData);

    rootRank = mpi.rank;

#ifdef PLANAR
    if (rootRank == 0) {
        std::cout << "\n\033[35m" << std::string(16, ' ') << " RUNNING TESTS FOR 2D CASE\033[0m" << std::endl << std::endl;
    }
    if (mpi.nProc > 1) {
        if (rootRank == 0) {
            std::cout << "\033[35m" << std::string(13, ' ') << " PARALLEL TESTS WITH 4 PROCESSORS\033[0m" << std::endl << std::endl;
        }

        numCases = 1;
        testParams.resize(numCases);
        testParams(0) = 7, 0, 4, 3;
    } else {
        if (rootRank == 0) {
            std::cout << "\033[35m" << std::string(14, ' ') << " SERIAL TESTS WITH 1 PROCESSOR\033[0m" << std::endl << std::endl;
        }

        numCases = 2;
        testParams.resize(numCases);
        testParams(0) = 5, 0, 4, 1;
        //testParams(0) = 4, 0, 5, 1;
        testParams(1) = 5, 0, 5, 2;
    }
#else
    if (rootRank == 0) {
        std::cout << "\n\033[35m" << std::string(16, ' ') << " RUNNING TESTS FOR 3D CASE\033[0m" << std::endl << std::endl;
    }
    if (mpi.nProc > 1) {
        if (rootRank == 0) {
            std::cout << "\033[35m" << std::string(13, ' ') << " PARALLEL TESTS WITH 4 PROCESSORS\033[0m" << std::endl << std::endl;
        }

        if (mpi.npX == 4 and mpi.npY == 1) {
            numCases = 1;
            testParams.resize(numCases);

            testParams(0) = 6, 4, 4, 2;
        } else if (mpi.npX == 1 and mpi.npY == 4) {
            numCases = 1;
            testParams.resize(numCases);

            testParams(0) = 4, 6, 4, 2;
        } else if (mpi.npX == 2 and mpi.npY == 2) {
            numCases = 3;
            testParams.resize(numCases);

            testParams(0) = 5, 5, 4, 1;
            testParams(1) = 6, 7, 4, 3;
            testParams(2) = 7, 4, 5, 2;
        }
    } else {
        if (rootRank == 0) {
            std::cout << "\033[35m" << std::string(14, ' ') << " SERIAL TESTS WITH 1 PROCESSOR\033[0m" << std::endl << std::endl;
        }

        numCases = 2;
        testParams.resize(numCases);
        testParams(0) = 4, 4, 5, 1;
        testParams(1) = 4, 5, 5, 2;
    }
#endif

    for (int i=0; i<numCases; i++) {
        inputData.xInd = testParams(i)(0);
        inputData.yInd = testParams(i)(1);
        inputData.zInd = testParams(i)(2);
        inputData.vcDepth = testParams(i)(3);

        // GRID OBJECT
        grid gridData(inputData, mpi);

        if (rootRank == 0) {
            std::cout << "\033[33m" << std::string(4, ' ') << " Test Case " << i+1 << " of " << numCases << ": "
                      << gridData.globalSize(0) << "x" << gridData.globalSize(1) << "x" << gridData.globalSize(2)
                      << " grid with V-Cycle depth: " << inputData.vcDepth << "\033[0m" << std::endl << std::endl;
        }

        differTest(gridData);
        fieldTest(gridData);
        nlinTest(gridData);
        poissonTest(gridData, inputData);
        hydroTest(gridData, inputData, mpi);
    }

    // FINALIZE AND CLEAN-UP
    MPI_Finalize();

    return 0;
}
