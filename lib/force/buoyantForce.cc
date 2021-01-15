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
/*! \file buoyantForce.cc
 *
 *  \brief Definitions for functions of class buoyantForce
 *  \sa force.h
 *  \author Shashwat Bhattacharya, Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "force.h"

buoyantForce::buoyantForce(const grid &mesh, vfield &U, const sfield &T): force(mesh, U), T(T) {
    switch (mesh.inputParams.rbcType) {
        case 1: Fb = mesh.inputParams.Ra*mesh.inputParams.Pr;
            break;
        case 2: Fb = 1.0;
            break;
        case 3: Fb = mesh.inputParams.Ra;
            break;
        case 4: Fb = mesh.inputParams.Pr;
            break;
    }
}


void buoyantForce::addForcing(plainvf &Hv) {
    //ADD THE BUOYANCY TERM TO THE Vz COMPONENT OF Hv
    V.interTempZ = 0.0;
    for (unsigned int i=0; i < V.Vz.PcIntSlices.size(); i++) {
        //if (mesh.rankData.rank == 0) std::cout << mesh.rankData.rank << "\t" << V.Vz.fCore.ubound() << "\t" << V.Vz.PcIntSlices(i).ubound() << "\t" << T.F.F.ubound() << std::endl;
        //MPI_Barrier(MPI_COMM_WORLD);
        //if (mesh.rankData.rank == 1) std::cout << mesh.rankData.rank << "\t" << V.Vz.fCore.ubound() << "\t" << V.Vz.PcIntSlices(i).ubound() << "\t" << T.F.F.ubound() << std::endl;
        //MPI_Barrier(MPI_COMM_WORLD);
        //if (mesh.rankData.rank == 2) std::cout << mesh.rankData.rank << "\t" << V.Vz.fCore.ubound() << "\t" << V.Vz.PcIntSlices(i).ubound() << "\t" << T.F.F.ubound() << std::endl;
        //MPI_Barrier(MPI_COMM_WORLD);
        //if (mesh.rankData.rank == 3) std::cout << mesh.rankData.rank << "\t" << V.Vz.fCore.ubound() << "\t" << V.Vz.PcIntSlices(i).ubound() << "\t" << T.F.F.ubound() << std::endl;
        //MPI_Barrier(MPI_COMM_WORLD);

        V.interTempZ(V.Vz.fCore) += T.F.F(V.Vz.PcIntSlices(i));
    }
    V.interTempZ /= V.Vz.PcIntSlices.size();

    //if (mesh.rankData.rank == 0) std::cout << mesh.rankData.rank << "\t" << Hv.Vz.ubound() << "\t" << V.interTempZ.ubound() << std::endl;
    //MPI_Barrier(MPI_COMM_WORLD);
    //if (mesh.rankData.rank == 1) std::cout << mesh.rankData.rank << "\t" << Hv.Vz.ubound() << "\t" << V.interTempZ.ubound() << std::endl;
    //MPI_Barrier(MPI_COMM_WORLD);
    //if (mesh.rankData.rank == 2) std::cout << mesh.rankData.rank << "\t" << Hv.Vz.ubound() << "\t" << V.interTempZ.ubound() << std::endl;
    //MPI_Barrier(MPI_COMM_WORLD);
    //if (mesh.rankData.rank == 3) std::cout << mesh.rankData.rank << "\t" << Hv.Vz.ubound() << "\t" << V.interTempZ.ubound() << std::endl;
    //MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Finalize();
    //exit(0);

    Hv.Vz += Fb*V.interTempZ;
    //Hv.Vz += 0.001;

    //if (mesh.rankData.rank == 0) std::cout << V.interTempZ(blitz::Range(25, blitz::toEnd), blitz::Range(25, blitz::toEnd), 5) << std::endl;
    //if (mesh.rankData.rank == 0) std::cout << V.interTempZ(blitz::Range(25, 33), blitz::Range(25, 33), 3) << std::endl;
    //if (mesh.rankData.rank == 3) std::cout << V.interTempZ(blitz::Range(-1, 5), blitz::Range(-1, 5), 3) << std::endl;

    //if (mesh.rankData.rank == 0) std::cout << V.interTempZ(blitz::Range::all(), blitz::Range::all(), 3) << std::endl;
}
