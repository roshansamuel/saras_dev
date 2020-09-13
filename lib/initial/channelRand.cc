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
 *          The random values are scaled by a factor, and depending on a hardcoded parameter (mainly for debugging),
 *          are either added directly to the mean velocity field or multiplied with an inverse parabola before being
 *          added to the mean velocity field.
 *
 ********************************************************************************************************************************************
 */
void channelRand::initializeField(vfield &uField) {
    if (mesh.rankData.rank == 0) std::cout << "Imposing random initial condition for channel flow" << std::endl << std::endl;

    // Mode of imposing the random fluctuations on mean field:
    // 0 - Add directly to mean field
    // 1 - Scale by an inverse parabola before adding to mean field
    int rfMode = 0;

    // Seed the random number generator with both time and rank to get different random numbers in different MPI sub-domains
    int randSeed = std::time(0) + mesh.rankData.rank;
    std::srand(randSeed);

    // Compute the uniform velocity field
    // Von Karman constant for log law
    real vkConst = 0.41;

    // Second constant of the log law
    real cPlus = 5.0;

    // Assuming that the Re specified is in terms of wall shear stress, the y_plus at centreline will be Re_tau.
    // Substituting this value in the log law will give centreline velocity non-dimensionalized by u_tau.
    real vCentre = log(mesh.inputParams.Re)/vkConst + cPlus;

    // Scaling factor for random velocity fluctuations
    real vScale = vCentre*0.15;

    // Random velocity fluctuation to be computed point-wise
    real vRandom;

    switch (rfMode) {
        case 0:
#ifdef PLANAR
            // X-VELOCITY
            for (int i=uField.Vx.F.lbound(0); i <= uField.Vx.F.ubound(0); i++) {
                for (int k=uField.Vx.F.lbound(2); k <= uField.Vx.F.ubound(2); k++) {
                    vRandom = vScale*(real(std::rand())/RAND_MAX - 0.5);

                    uField.Vx.F(i, 0, k) = vCentre + vRandom;
                }
            }

            // Z-VELOCITY
            for (int i=uField.Vz.F.lbound(0); i <= uField.Vz.F.ubound(0); i++) {
                for (int k=uField.Vz.F.lbound(2); k <= uField.Vz.F.ubound(2); k++) {
                    vRandom = vScale*(real(std::rand())/RAND_MAX - 0.5);

                    uField.Vz.F(i, 0, k) = vRandom;
                }
            }
#else
            // X-VELOCITY
            for (int i=uField.Vx.F.lbound(0); i <= uField.Vx.F.ubound(0); i++) {
                for (int j=uField.Vx.F.lbound(1); j <= uField.Vx.F.ubound(1); j++) {
                    for (int k=uField.Vx.F.lbound(2); k <= uField.Vx.F.ubound(2); k++) {
                        vRandom = vScale*(real(std::rand())/RAND_MAX - 0.5);

                        uField.Vx.F(i, j, k) = vCentre + vRandom;
                    }
                }
            }

            // Y-VELOCITY
            for (int i=uField.Vy.F.lbound(0); i <= uField.Vy.F.ubound(0); i++) {
                for (int j=uField.Vy.F.lbound(1); j <= uField.Vy.F.ubound(1); j++) {
                    for (int k=uField.Vy.F.lbound(2); k <= uField.Vy.F.ubound(2); k++) {
                        vRandom = vScale*(real(std::rand())/RAND_MAX - 0.5);

                        uField.Vy.F(i, j, k) = vRandom;
                    }
                }
            }

            // Z-VELOCITY
            for (int i=uField.Vz.F.lbound(0); i <= uField.Vz.F.ubound(0); i++) {
                for (int j=uField.Vz.F.lbound(1); j <= uField.Vz.F.ubound(1); j++) {
                    for (int k=uField.Vz.F.lbound(2); k <= uField.Vz.F.ubound(2); k++) {
                        vRandom = vScale*(real(std::rand())/RAND_MAX - 0.5);

                        uField.Vz.F(i, j, k) = vRandom;
                    }
                }
            }
#endif
            break;

        case 1:
            real normalZ, randNum;
#ifdef PLANAR
            // X-VELOCITY
            for (int i=uField.Vx.F.lbound(0); i <= uField.Vx.F.ubound(0); i++) {
                for (int k=uField.Vx.F.lbound(2); k <= uField.Vx.F.ubound(2); k++) {
                    randNum = vScale*real(std::rand())/RAND_MAX;
                    normalZ = mesh.zStaggr(k)/mesh.zLen;
                    vRandom = randNum*(4.0*normalZ*(normalZ - 1.0) + 1);

                    uField.Vx.F(i, 0, k) = vCentre + vRandom;
                }
            }

            // Z-VELOCITY
            for (int i=uField.Vz.F.lbound(0); i <= uField.Vz.F.ubound(0); i++) {
                for (int k=uField.Vz.F.lbound(2); k <= uField.Vz.F.ubound(2); k++) {
                    randNum = vScale*real(std::rand())/RAND_MAX;
                    normalZ = mesh.zColloc(k)/mesh.zLen;
                    vRandom = randNum*(4.0*normalZ*(normalZ - 1.0) + 1);

                    uField.Vz.F(i, 0, k) = vRandom;
                }
            }
#else
            // X-VELOCITY
            for (int i=uField.Vx.F.lbound(0); i <= uField.Vx.F.ubound(0); i++) {
                for (int j=uField.Vx.F.lbound(1); j <= uField.Vx.F.ubound(1); j++) {
                    for (int k=uField.Vx.F.lbound(2); k <= uField.Vx.F.ubound(2); k++) {
                        randNum = vScale*real(std::rand())/RAND_MAX;
                        normalZ = mesh.zStaggr(k)/mesh.zLen;
                        vRandom = randNum*(4.0*normalZ*(normalZ - 1.0) + 1);

                        uField.Vx.F(i, j, k) = vCentre + vRandom;
                    }
                }
            }

            // Y-VELOCITY
            for (int i=uField.Vy.F.lbound(0); i <= uField.Vy.F.ubound(0); i++) {
                for (int j=uField.Vy.F.lbound(1); j <= uField.Vy.F.ubound(1); j++) {
                    for (int k=uField.Vy.F.lbound(2); k <= uField.Vy.F.ubound(2); k++) {
                        randNum = vScale*real(std::rand())/RAND_MAX;
                        normalZ = mesh.zStaggr(k)/mesh.zLen;
                        vRandom = randNum*(4.0*normalZ*(normalZ - 1.0) + 1);

                        uField.Vy.F(i, j, k) = vRandom;
                    }
                }
            }

            // Z-VELOCITY
            for (int i=uField.Vz.F.lbound(0); i <= uField.Vz.F.ubound(0); i++) {
                for (int j=uField.Vz.F.lbound(1); j <= uField.Vz.F.ubound(1); j++) {
                    for (int k=uField.Vz.F.lbound(2); k <= uField.Vz.F.ubound(2); k++) {
                        randNum = vScale*real(std::rand())/RAND_MAX;
                        normalZ = mesh.zColloc(k)/mesh.zLen;
                        vRandom = randNum*(4.0*normalZ*(normalZ - 1.0) + 1);

                        uField.Vz.F(i, j, k) = vRandom;
                    }
                }
            }
#endif
            break;
    }
}
