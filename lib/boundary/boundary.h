/********************************************************************************************************************************************
 * Saras
 * 
 * Copyright (C) 2019, Mahendra K. Verma
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     1. Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *     2. Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *     3. Neither the name of the copyright holder nor the
 *        names of its contributors may be used to endorse or promote products
 *        derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ********************************************************************************************************************************************
 */
/*! \file boundary.h
 *
 *  \brief Class declaration of boundary
 *
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <blitz/array.h>

#include "field.h"
#include "grid.h"

class boundary {
    public:
        boundary(const grid &mesh, field &inField, const bool bcType, const int shiftDim);

        void createPatch(int wallNum);

        void imposeBC();

    private:
        const grid &mesh;

        field &dField;

        bool nonHgBC;
        const bool drcBC;
        const int dimVal;

        // The following arrays are allocated memory only if non-homogeneous BCs are being used.
        // Hence attempting to access these arrays in other situations will result in seg-fault.
        blitz::Array<bool, 3> wallMask;
        blitz::Array<double, 3> wallData;

        blitz::Array<double, 1> x, y, z;
        blitz::Array<double, 1> xGlo, yGlo, zGlo;

        void imposeBC_X();
        void imposeBC_Y();
        void imposeBC_Z();

        void setXYZ();
};

/**
 ********************************************************************************************************************************************
 *  \class boundary boundary.h "lib/boundary/boundary.h"
 *  \brief Contains all the global variables related to the imposing of boundary conditions, and functions to impose BCs
 *
 ********************************************************************************************************************************************
 */

#endif
