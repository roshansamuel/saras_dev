#ifndef UNITTEST_H
#define UNITTEST_H

#include <blitz/array.h>
#include <iostream>

extern int rootRank;

void testError(blitz::Array<double, 3> A, blitz::Array<double, 3> B, int errorMom, double errorTol);

void printResult(double computedValue, double errorTolerance);

#endif
