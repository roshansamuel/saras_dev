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
/*! \file boundary.cc
 *
 *  \brief Definitions for functions of class boundary
 *  \sa boundary.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "boundary.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the boundary class
 *
 *          The class constructor initializes the mesh for computational problem.
 *          Based on the user set values in the parser class, it decides if the simulation is going to be 2D or 3D.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   inField is a reference to scalar field to which the boundary conditions must be applied.
 ********************************************************************************************************************************************
 */
boundary::boundary(const grid &mesh, field &inField, const int bcWall):
                          mesh(mesh), dField(inField), wallNum(bcWall) {
    // By default, rankFlag is true. i.e., the BC will be applied on all sub-domains.
    // This works only for Z-direction (in pencil decomposition), or both Z and Y directions (in slab decomposition).
    // This has to be changed appropriately.
    rankFlag = true;

    // By default, shiftVal is 1. i.e., the BC will be applied on the left wall along a given direction.
    // This has to be changed appropriately for the wall on the other side.
    shiftVal = 1;

    if (wallNum == 0) rankFlag = mesh.rankData.xRank == 0;
    if (wallNum == 1) {
        rankFlag = mesh.rankData.xRank == mesh.rankData.npX - 1;
        shiftVal = -1;
    }
#ifndef PLANAR
    if (wallNum == 2) rankFlag = mesh.rankData.yRank == 0;
    if (wallNum == 3) {
        rankFlag = mesh.rankData.yRank == mesh.rankData.npY - 1;
        shiftVal = -1;
    }
#endif
    if (wallNum == 5) {
        shiftVal = -1;
    }

    shiftDim = (int) wallNum/2;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions on the given field
 *
 *          The function first calls the syncData() function of the field to update the sub-domain pads.
 *          Then the boundary conditions are applied at the full domain boundaries by calling individual functions
 *          to impose the BCs along X, Y and Z directions.
 *
 ********************************************************************************************************************************************
 */
void boundary::imposeBC() { };
