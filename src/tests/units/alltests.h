#ifndef ALLTESTS_H
#define ALLTESTS_H

#include "grid.h"
#include "parser.h"
#include "parallel.h"

void fieldTest(grid &gridData);

void differTest(grid &gridData);

void nlinTest(grid &gridData);

void poissonTest(grid &gridData, parser &inputData);

void hydroTest(grid &gridData, parser &inputData, parallel &mpiData);

#endif
