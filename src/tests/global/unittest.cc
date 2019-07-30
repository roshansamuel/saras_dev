#include "unittest.h"

void testError(blitz::Array<double, 3> A, blitz::Array<double, 3> B, int errorMom, double errorTol) {
    int errorCnt;
    double errorVal;

    errorCnt = 0;
    errorVal = 0.0;

    if (errorMom == 1) {
        for (int i=A.lbound(0)+1; i <= A.ubound(0)-1; i++) {
            for (int j=A.lbound(1)+1; j <= A.ubound(1)-1; j++) {
                for (int k=A.lbound(2)+1; k <= A.ubound(2)-1; k++) {
                    errorCnt += 1;
                    errorVal += fabs(A(i, j, k) - B(i, j, k));
                }
            }
        }

        errorVal /= errorCnt;

    } else if (errorMom == 2) {
        for (int i=A.lbound(0)+1; i <= A.ubound(0)-1; i++) {
            for (int j=A.lbound(1)+1; j <= A.ubound(1)-1; j++) {
                for (int k=A.lbound(2)+1; k <= A.ubound(2)-1; k++) {
                    errorCnt += 1;
                    errorVal += pow((A(i, j, k) - B(i, j, k)), 2.0);
                }
            }
        }

        errorVal /= errorCnt;
        errorVal = sqrt(errorVal);
    }

    printResult(errorVal, errorTol);
}

void printResult(double computedValue, double errorTolerance) {
    if (computedValue > errorTolerance) {
        if (rootRank == 0) {
            std::cout << "\033[31m[FAILED]\033[0m" << std::endl << std::endl;
        }
    } else {
        if (rootRank == 0) {
            std::cout << "\033[32m[PASSED]\033[0m" << std::endl << std::endl;
        }
    }
}
