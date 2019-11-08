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
initial::initial(const grid &mesh): mesh(mesh) { }


/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the Taylor Green Vortex initial condition on the given input velocity field
 *
 *          Depending on the preprocessor flag PLANAR, the function applies the Taylor Green vortex equation
 *          to set initial conditions in both 2D and 3D cases.
 *
 ********************************************************************************************************************************************
 */
void initial::initTGV(vfield &uField) {
    if (mesh.rankData.rank == 0) std::cout << "Imposing Taylor-Green vortices initial condition" << std::endl << std::endl;

#ifdef PLANAR
    // X-VELOCITY
    for (int i=uField.Vx.F.lbound(0); i <= uField.Vx.F.ubound(0); i++) {
        for (int k=uField.Vx.F.lbound(2); k <= uField.Vx.F.ubound(2); k++) {
            uField.Vx.F(i, 0, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                   cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
        }
    }

    // Z-VELOCITY
    for (int i=uField.Vz.F.lbound(0); i <= uField.Vz.F.ubound(0); i++) {
        for (int k=uField.Vz.F.lbound(2); k <= uField.Vz.F.ubound(2); k++) {
            uField.Vz.F(i, 0, k) = -cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                    sin(2.0*M_PI*mesh.zColloc(k)/mesh.zLen);
        }
    }

#else
    // X-VELOCITY
    for (int i=uField.Vx.F.lbound(0); i <= uField.Vx.F.ubound(0); i++) {
        for (int j=uField.Vx.F.lbound(1); j <= uField.Vx.F.ubound(1); j++) {
            for (int k=uField.Vx.F.lbound(2); k <= uField.Vx.F.ubound(2); k++) {
                uField.Vx.F(i, j, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                       cos(2.0*M_PI*mesh.yStaggr(j)/mesh.yLen)*
                                       cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
            }
        }
    }

    // Y-VELOCITY
    for (int i=uField.Vy.F.lbound(0); i <= uField.Vy.F.ubound(0); i++) {
        for (int j=uField.Vy.F.lbound(1); j <= uField.Vy.F.ubound(1); j++) {
            for (int k=uField.Vy.F.lbound(2); k <= uField.Vy.F.ubound(2); k++) {
                uField.Vy.F(i, j, k) = -cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                        sin(2.0*M_PI*mesh.yColloc(j)/mesh.yLen)*
                                        cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
            }
        }
    }

    // Z-VELOCITY
    uField.Vz.F = 0.0;
#endif
}


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
void initial::initRandom(vfield &uField) {
    if (mesh.rankData.rank == 0) std::cout << "Imposing random initial condition" << std::endl << std::endl;

    // Seed the random number generator with both time and rank to get different random numbers in different MPI sub-domains
    int randSeed = std::time(0) + mesh.rankData.rank;
    std::srand(randSeed);

    // The below factor was recommended to be set to 26.0 by Anikesh Pal, but that number appears when V is scaled with friction velocity, U_tau
    double vScale = 1.0;

#ifdef PLANAR
    // VELOCITY PERTURBATION FOR PERIODIC CHANNEL FLOW
    // X-VELOCITY
    for (int i=uField.Vx.F.lbound(0); i <= uField.Vx.F.ubound(0); i++) {
        for (int k=uField.Vx.F.lbound(2); k <= uField.Vx.F.ubound(2); k++) {
            double uP, ndZ;
            double randNum = double(std::rand())/RAND_MAX;

            ndZ = mesh.zStaggr(k)/mesh.zLen;
            uP = randNum*vScale*(4.0*ndZ*(ndZ - 1.0) + 1);
            uField.Vx.F(i, 0, k) = 1.0 + uP;
        }
    }

    // Z-VELOCITY
    for (int i=uField.Vz.F.lbound(0); i <= uField.Vz.F.ubound(0); i++) {
        for (int k=uField.Vz.F.lbound(2); k <= uField.Vz.F.ubound(2); k++) {
            double uP, ndZ;
            double randNum = double(std::rand())/RAND_MAX;

            ndZ = mesh.zStaggr(k)/mesh.zLen;
            uP = randNum*vScale*(4.0*ndZ*(ndZ - 1.0) + 1);
            uField.Vz.F(i, 0, k) = uP;
        }
    }

#else
    // VELOCITY PERTURBATION FOR PERIODIC CHANNEL FLOW
    // X-VELOCITY
    for (int i=uField.Vx.F.lbound(0); i <= uField.Vx.F.ubound(0); i++) {
        for (int j=uField.Vx.F.lbound(1); j <= uField.Vx.F.ubound(1); j++) {
            for (int k=uField.Vx.F.lbound(2); k <= uField.Vx.F.ubound(2); k++) {
                double uP, ndZ;
                double randNum = double(std::rand())/RAND_MAX;

                ndZ = mesh.zStaggr(k)/mesh.zLen;
                uP = randNum*vScale*(4.0*ndZ*(ndZ - 1.0) + 1);

                // ALONG X, THERE IS A MEAN FLOW OF 1.0 UNIT VELOCITY
                uField.Vx.F(i, j, k) = 1.0 + uP;
            }
        }
    }

    // Y-VELOCITY
    for (int i=uField.Vy.F.lbound(0); i <= uField.Vy.F.ubound(0); i++) {
        for (int j=uField.Vy.F.lbound(1); j <= uField.Vy.F.ubound(1); j++) {
            for (int k=uField.Vy.F.lbound(2); k <= uField.Vy.F.ubound(2); k++) {
                double uP, ndZ;
                double randNum = double(std::rand())/RAND_MAX;

                ndZ = mesh.zStaggr(k)/mesh.zLen;
                uP = randNum*vScale*(4.0*ndZ*(ndZ - 1.0) + 1);
                uField.Vy.F(i, j, k) = uP;
            }
        }
    }

    // Z-VELOCITY
    for (int i=uField.Vz.F.lbound(0); i <= uField.Vz.F.ubound(0); i++) {
        for (int j=uField.Vz.F.lbound(1); j <= uField.Vz.F.ubound(1); j++) {
            for (int k=uField.Vz.F.lbound(2); k <= uField.Vz.F.ubound(2); k++) {
                double uP, ndZ;
                double randNum = double(std::rand())/RAND_MAX;

                ndZ = mesh.zStaggr(k)/mesh.zLen;
                uP = randNum*vScale*(4.0*ndZ*(ndZ - 1.0) + 1);
                uField.Vz.F(i, j, k) = uP;
            }
        }
    }
#endif
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to generate sinusoidal perturbation for channel flow
 *
 *          The sinusoidal function is multiplied with the equation of a parabola in order to satisfy the channel flow BCs.
 *          This is not divergence-free and has performed poorly in tests so far.
 *
 ********************************************************************************************************************************************
 */
void initial::initSinusoidal(vfield &uField) {
    double kx = 10.0;

    if (mesh.rankData.rank == 0) std::cout << "Imposing sinusoidal perturbation initial condition" << std::endl << std::endl;

#ifdef PLANAR
    // VELOCITY PERTURBATION FOR PERIODIC CHANNEL FLOW
    // X-VELOCITY
    for (int i=uField.Vx.F.lbound(0); i <= uField.Vx.F.ubound(0); i++) {
        for (int k=uField.Vx.F.lbound(2); k <= uField.Vx.F.ubound(2); k++) {
            uField.Vx.F(i, 0, k) = 4.0*(mesh.zStaggr(k) - mesh.zStaggr(k)*mesh.zStaggr(k));
        }
    }

    // Z-VELOCITY
    for (int i=uField.Vz.F.lbound(0); i <= uField.Vz.F.ubound(0); i++) {
        for (int k=uField.Vz.F.lbound(2); k <= uField.Vz.F.ubound(2); k++) {
            uField.Vz.F(i, 0, k) = 0.1*(mesh.zStaggr(k) - mesh.zStaggr(k)*mesh.zStaggr(k))*sin(2.0*M_PI*kx*mesh.xColloc(i)/mesh.xLen);
        }
    }

#else
    // VELOCITY PERTURBATION FOR PERIODIC CHANNEL FLOW
    // X-VELOCITY
    for (int i=uField.Vx.F.lbound(0); i <= uField.Vx.F.ubound(0); i++) {
        for (int j=uField.Vx.F.lbound(1); j <= uField.Vx.F.ubound(1); j++) {
            for (int k=uField.Vx.F.lbound(2); k <= uField.Vx.F.ubound(2); k++) {
                uField.Vx.F(i, j, k) = 4.0*(mesh.zStaggr(k) - mesh.zStaggr(k)*mesh.zStaggr(k));
            }
        }
    }

    // Y-VELOCITY
    uField.Vy.F = 0.0;

    // Z-VELOCITY
    for (int i=uField.Vz.F.lbound(0); i <= uField.Vz.F.ubound(0); i++) {
        for (int j=uField.Vz.F.lbound(1); j <= uField.Vz.F.ubound(1); j++) {
            for (int k=uField.Vz.F.lbound(2); k <= uField.Vz.F.ubound(2); k++) {
                uField.Vz.F(i, j, k) = 0.1*(mesh.zStaggr(k) - mesh.zStaggr(k)*mesh.zStaggr(k))*sin(2.0*M_PI*kx*mesh.xColloc(i)/mesh.xLen);
            }
        }
    }
#endif
}
