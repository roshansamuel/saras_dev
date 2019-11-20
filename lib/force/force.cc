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
/*! \file force.cc
 *
 *  \brief Definitions for functions of class force
 *  \sa force.h
 *  \author Shashwat Bhattacharya, Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "force.h"

force::force(vfield &U, const parser &solParams, parallel &mpiParam): V(U), mpiData(mpiParam), inputParams(solParams) {
    if (inputParams.probType < 5) {
        Fr = 1.0/inputParams.Ro;
    } else {
        if (inputParams.rbcType == 1) {                     
            Fb = inputParams.Ra*inputParams.Pr;
            Fr = inputParams.Pr*sqrt(inputParams.Ta);
        } else if (inputParams.rbcType == 2) {    
            Fb = 1.0;
            Fr = sqrt(inputParams.Ta*inputParams.Pr/inputParams.Ra); 
        } else if (inputParams.rbcType == 3) {    
            Fb = inputParams.Ra;
            Fr = sqrt(inputParams.Ta);
        } else {                                                            
            Fb = inputParams.Pr;
            Fr = sqrt(inputParams.Ta*inputParams.Pr/inputParams.Ra); 
        }
    }
}


void force::add_VForce(plainvf &Hv) {
    if (inputParams.forceType==0) {
        //No forcing

    } else if (inputParams.forceType==1) {
        //Random forcing
        add_RandomForce(Hv);

    } else if (inputParams.forceType==2) {
        //Coriolis forcing
        add_Coriolis(Hv);

    } else if (inputParams.forceType <= 4) {
        // Throw error message, since forcing 3 and 4 are not allowed for hydrodynamic flows
        if (mpiData.rank == 0) {
            std::cout << "ERROR: Chosen forcing is incompatible with the given problem type. ABORTING" << std::endl;
        }

        MPI_Finalize();
        exit(0);
    } else {
        if (mpiData.rank == 0) {
            std::cout << "ERROR: Invalid forcing parameter. ABORTING" << std::endl;
        }

        MPI_Finalize();
        exit(0);
    }
}


void force::add_VForce(plainvf &Hv, sfield &T) {
    if (inputParams.forceType==0) {
        //No forcing

    } else if (inputParams.forceType==1) {
        //Random forcing
        add_RandomForce(Hv);

    } else if (inputParams.forceType==2) {
        //Coriolis forcing
        add_Coriolis(Hv);

    } else if (inputParams.forceType==3) {
        //Buoyancy forcing
        add_Buoyancy(Hv, T);

    } else if (inputParams.forceType==4) {
        //Buoyancy and Coriolis forcing
        add_Buoyancy(Hv, T);
        add_Coriolis(Hv);

    } else {
        if (mpiData.rank == 0) {
            std::cout << "ERROR: Invalid forcing parameter. ABORTING" << std::endl;
        }

        MPI_Finalize();
        exit(0);
    }
}


void force::add_RandomForce(plainvf &Hv){
    // WARNING: The below implementation involves array definitions and resizings to be done in each time-step
    // This can seriously degrade the performance of the solver
    blitz::Array<double, 3> Force_x, Force_y, Force_z;

    // Currently, random forcing is not implemented. The forcing is 0 for now
    Force_x.resize(V.Vx.fSize);
    Force_x.reindexSelf(V.Vx.flBound);
    Force_x = 0;
    Hv.Vx += Force_x;

#ifndef PLANAR
    Force_y.resize(V.Vy.fSize);
    Force_y.reindexSelf(V.Vy.flBound);
    Force_y = 0;
    Hv.Vy += Force_y;
#endif

    Force_z.resize(V.Vz.fSize);
    Force_z.reindexSelf(V.Vz.flBound);
    Force_z = 0;
    Hv.Vz += Force_z;
}


void force::add_Buoyancy(plainvf &Hv, sfield &T) {
    //ADD THE BUOYANCY TERM TO THE Vz COMPONENT OF Hv
    V.interTempZ = 0.0;
    for (unsigned int i=0; i < V.Vz.PcIntSlices.size(); i++) {
        V.interTempZ(V.Vz.fCore) += T.F.F(V.Vz.PcIntSlices(i));
    }
    V.interTempZ /= V.Vz.PcIntSlices.size();

    Hv.Vz += Fb*V.interTempZ;
}


void force::add_Coriolis(plainvf &Hv){
#ifndef PLANAR
    //ADD THE ROTATING TERM TO THE Vx COMPONENT OF Hv
    V.interTempX = 0.0;
    for (unsigned int i=0; i < V.Vx.VyIntSlices.size(); i++) {
        V.interTempX(V.Vx.fCore) += V.Vy.F(V.Vx.VyIntSlices(i));
    }   
    V.interTempX /= V.Vx.VyIntSlices.size();

    Hv.Vx += Fr*V.interTempX;

    //SUBTRACT THE ROTATING TERM FROM THE Vy COMPONENT of Hv
    V.interTempY = 0.0;
    for (unsigned int i=0; i < V.Vy.VxIntSlices.size(); i++) {
        V.interTempY(V.Vy.fCore) += V.Vx.F(V.Vy.VxIntSlices(i));
    }   
    V.interTempY /= V.Vy.VxIntSlices.size();

    Hv.Vy -= Fr*V.interTempY;
#endif
}


void::force::add_SForce(plainsf &Ht) {
    //Scalar forcing needs to be implemented
    Ht.F += 0;
}
