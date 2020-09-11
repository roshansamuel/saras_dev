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
/*! \file initial.cc
 *
 *  \brief Definitions for functions of class initial
 *  \sa initial.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include <ctime>
#include <cstdlib>
#include "initial.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the initial class
 *
 *          The empty constructer merely initializes the local reference to the global mesh variable.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 ********************************************************************************************************************************************
 */
channelRand::channelRand(const grid &mesh): initial(mesh) { }


/**
 ********************************************************************************************************************************************
 * \brief   Function to impose random initial condition for channel flow on the given input velocity field.
 *
 *          The function generates random values for the entire domain (with independent seeds for individual ranks).
 *          The random values are scaled by multiplying with an inverse parabola and then added to the mean velocity field.
 *          The mean field is assumed to have value 1.0
 *
 ********************************************************************************************************************************************
 */
void channelRand::initializeField(vfield &uField) {
    if (mesh.rankData.rank == 0) std::cout << "Imposing random initial condition for channel flow" << std::endl << std::endl;

    // Seed the random number generator with both time and rank to get different random numbers in different MPI sub-domains
    int randSeed = std::time(0) + mesh.rankData.rank;
    std::srand(randSeed);

    // DEBUG CODE
    //std::srand(mesh.rankData.rank);
    // END DEBUG CODE

    // Von Karman constant for log law
    real vkConst = 0.41;

    // Second constant of the log law
    real cPlus = 5.0;

    // Assuming that the Re specified is in terms of wall shear stress, the y_plus at centreline will be Re_tau.
    // Substituting this value in the log law will give centreline velocity non-dimensionalized by u_tau.
    real vCentre = log(mesh.inputParams.Re)/vkConst + cPlus;

    // Generate inverted parabolic profile with vCentre as velocity at edges
    // Two profiles are generated for cell centered and face centered variables
    blitz::Array<real, 1> proStag, proColl;
    proStag.resize(uField.Vx.fCore.ubound(2) + 1);
    proColl.resize(uField.Vz.fCore.ubound(2) + 1);

    proStag = 0.0;
    for (int k=0; k <= uField.Vx.fCore.ubound(2); k++) {
        real ndZ = mesh.zStaggr(k)/mesh.zLen;
        proStag(k) = vCentre*(4.0*ndZ*(ndZ - 1.0) + 1);
    }

    proColl = 0.0;
    for (int k=0; k <= uField.Vz.fCore.ubound(2); k++) {
        real ndZ = mesh.zColloc(k)/mesh.zLen;
        proColl(k) = vCentre*(4.0*ndZ*(ndZ - 1.0) + 1);
    }

#ifdef PLANAR
    // VELOCITY PERTURBATION FOR PERIODIC CHANNEL FLOW
    // X-VELOCITY
    for (int i=uField.Vx.fCore.lbound(0); i <= uField.Vx.fCore.ubound(0); i++) {
        for (int k=uField.Vx.fCore.lbound(2); k <= uField.Vx.fCore.ubound(2); k++) {
            real randNum = real(std::rand())/RAND_MAX;

            // ALONG X, THERE IS A MEAN FLOW OF vCentre UNITS
            uField.Vx.F(i, 0, k) = vCentre + randNum*proStag(k);
        }
    }

    // Z-VELOCITY
    for (int i=uField.Vz.fCore.lbound(0); i <= uField.Vz.fCore.ubound(0); i++) {
        for (int k=uField.Vz.fCore.lbound(2); k <= uField.Vz.fCore.ubound(2); k++) {
            real randNum = real(std::rand())/RAND_MAX;

            uField.Vz.F(i, 0, k) = randNum*proColl(k);
        }
    }
#else
    // VELOCITY PERTURBATION FOR PERIODIC CHANNEL FLOW
    // X-VELOCITY
    for (int i=uField.Vx.fCore.lbound(0); i <= uField.Vx.fCore.ubound(0); i++) {
        for (int j=uField.Vx.fCore.lbound(1); j <= uField.Vx.fCore.ubound(1); j++) {
            for (int k=uField.Vx.fCore.lbound(2); k <= uField.Vx.fCore.ubound(2); k++) {
                real randNum = real(std::rand())/RAND_MAX;

                // ALONG X, THERE IS A MEAN FLOW OF vCentre UNITS
                uField.Vx.F(i, j, k) = vCentre + randNum*proStag(k);
            }
        }
    }

    // Y-VELOCITY
    for (int i=uField.Vy.fCore.lbound(0); i <= uField.Vy.fCore.ubound(0); i++) {
        for (int j=uField.Vy.fCore.lbound(1); j <= uField.Vy.fCore.ubound(1); j++) {
            for (int k=uField.Vy.fCore.lbound(2); k <= uField.Vy.fCore.ubound(2); k++) {
                real randNum = real(std::rand())/RAND_MAX;

                uField.Vy.F(i, j, k) = randNum*proStag(k);
            }
        }
    }

    // Z-VELOCITY
    for (int i=uField.Vz.fCore.lbound(0); i <= uField.Vz.fCore.ubound(0); i++) {
        for (int j=uField.Vz.fCore.lbound(1); j <= uField.Vz.fCore.ubound(1); j++) {
            for (int k=uField.Vz.fCore.lbound(2); k <= uField.Vz.fCore.ubound(2); k++) {
                real randNum = real(std::rand())/RAND_MAX;

                uField.Vz.F(i, j, k) = randNum*proColl(k);
            }
        }
    }
#endif

    //for (int k=uField.Vx.F.lbound(2); k <= uField.Vx.F.ubound(2); k++) {
    //    real ndZ = mesh.zStaggr(k)/mesh.zLen;
    //    real uP = vCentre*(4.0*ndZ*(ndZ - 1.0) + 1);
    //    if (mesh.rankData.rank == 0) std::cout << k << "\t" << ndZ << "\t" << uP << std::endl;
    //}
    //if (mesh.rankData.rank == 0) std::cout << uField.Vx.F(5, 5, blitz::Range::all()) << std::endl;
    //MPI_Finalize();
    //exit(0);
}
