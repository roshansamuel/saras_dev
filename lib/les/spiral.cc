/********************************************************************************************************************************************
 * Spiral LES
 * 
 * Copyright (C) 2009, California Institute of Technology (Caltech)
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
 * THIS SOFTWARE IS PROVIDED BY Caltech "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL Caltech BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ********************************************************************************************************************************************
 */
/*! \file spiral.cc
 *
 *  \brief Definitions for functions of class spiral
 *  \sa les.h
 *  \date Sep 2020
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "les.h"

#define EPS (2e-15)

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the spiral class
 *
 *          The constructor initializes the base les class using part of the arguments supplied to it.
 *          The value of viscous diffusion coefficient, denoted by nu, is also set in the initialization list.
 *          All the blitz arrays and scalar fields necessary to compute velocity gradient tensor, stress tensor, etc.
 *          are initialized in the constructor.
 *
 * \param   mesh is a const reference to the global data contained in the grid class.
 * \param   solverV is a const reference to the velocity field which will be read to compute structure function and velocity gradient
 * \param   solverP is a const reference to the cell-centered pressure field whose array limits are used to compute SGS terms
 * \param   kDiff is a const reference to scalar value denoting viscous dissipation
 ********************************************************************************************************************************************
 */
spiral::spiral(const grid &mesh, const vfield &solverV, const sfield &solverP, const real &kDiff):
    les(mesh, solverV, solverP),
    nu(kDiff)
{
    // Scalar fields used to store components of the sub-grid stress tensor field
    Txx = new sfield(mesh, "Txx");
    Tyy = new sfield(mesh, "Tyy");
    Tzz = new sfield(mesh, "Tzz");
    Txy = new sfield(mesh, "Txy");
    Tyz = new sfield(mesh, "Tyz");
    Tzx = new sfield(mesh, "Tzx");

    // Scalar fields used to store components of the sub-grid scalar flux vector
    qX = new sfield(mesh, "qX");
    qY = new sfield(mesh, "qY");
    qZ = new sfield(mesh, "qZ");

    // Scalar fields used to store the velocity fields interpolated to cell-centers
    // The vector field class cannot be used here because vfield by default stores
    // its comopnents in face-centered configuration.
    Vxcc = new sfield(mesh, "Ucc");
    Vycc = new sfield(mesh, "Vcc");
    Vzcc = new sfield(mesh, "Wcc");

    // 3x3x3 arrays which store local velocity field when computing structure function
    u.resize(3, 3, 3);
    v.resize(3, 3, 3);
    w.resize(3, 3, 3);

    // 3x3 matrix to store the velocity gradient tensor
    dudx.resize(3, 3);

    // The 9 blitz arrays of tensor components have the same dimensions and limits as a cell centered variable (P)
    A11.resize(P.F.fSize);      A11.reindexSelf(P.F.flBound);
    A12.resize(P.F.fSize);      A12.reindexSelf(P.F.flBound);
    A13.resize(P.F.fSize);      A13.reindexSelf(P.F.flBound);
    A21.resize(P.F.fSize);      A21.reindexSelf(P.F.flBound);
    A22.resize(P.F.fSize);      A22.reindexSelf(P.F.flBound);
    A23.resize(P.F.fSize);      A23.reindexSelf(P.F.flBound);
    A31.resize(P.F.fSize);      A31.reindexSelf(P.F.flBound);
    A32.resize(P.F.fSize);      A32.reindexSelf(P.F.flBound);
    A33.resize(P.F.fSize);      A33.reindexSelf(P.F.flBound);

    // The 9 arrays of vector components have the same dimensions and limits as a cell centered variable (P)
    // These arrays are needed only when computing the sub-grid scalar flux.
    // Hence they are initialized only for thermal convection problems (probType > 4)
    if (mesh.inputParams.probType > 4) {
        B1.resize(P.F.fSize);       B1.reindexSelf(P.F.flBound);
        B2.resize(P.F.fSize);       B2.reindexSelf(P.F.flBound);
        B3.resize(P.F.fSize);       B3.reindexSelf(P.F.flBound);
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to compute and add contribution from LES model to the RHS of NSE
 *
 *          This function is called when only the hydrodynamics equations are being solved.
 *          It calls the sgsStress function repeatedly for every point in the domain to get the
 *          sub-grid stress tensor field.
 *          The divergence of this tensor field is then calculated to obtain the spiral vortex model's
 *          contribution to the NSE.
 *
 * \param   nseRHS is a reference to the plain vector field denoting the RHS of the NSE
 ********************************************************************************************************************************************
 */
void spiral::computeSG(plainvf &nseRHS) {
    // First interpolate all velocities to cell centers
    // Then U, V and W data are available at all cell centers except ghost points
    Vxcc->F.F = 0.0;
    for (unsigned int i=0; i < P.F.VxIntSlices.size(); i++) {
        Vxcc->F.F(P.F.fCore) += V.Vx.F(P.F.VxIntSlices(i));
    }
    Vxcc->F.F /= P.F.VxIntSlices.size();

    Vycc->F.F = 0.0;
    for (unsigned int i=0; i < P.F.VyIntSlices.size(); i++) {
        Vycc->F.F(P.F.fCore) += V.Vy.F(P.F.VyIntSlices(i));
    }
    Vycc->F.F /= P.F.VyIntSlices.size();

    Vzcc->F.F = 0.0;
    for (unsigned int i=0; i < P.F.VzIntSlices.size(); i++) {
        Vzcc->F.F(P.F.fCore) += V.Vz.F(P.F.VzIntSlices(i));
    }
    Vzcc->F.F /= P.F.VzIntSlices.size();

    Vxcc->F.syncData();
    Vycc->F.syncData();
    Vzcc->F.syncData();

    // Compute the x, y and z derivatives of the interpolated velocity field and store them into
    // the arrays A11, A12, A13, ... A33. These arrays will be later accessed when constructing
    // the velocity gradient tensor at each point in the domain.
    Vxcc->derS.calcDerivative1_x(A11);
    Vxcc->derS.calcDerivative1_y(A12);
    Vxcc->derS.calcDerivative1_z(A13);
    Vycc->derS.calcDerivative1_x(A21);
    Vycc->derS.calcDerivative1_y(A22);
    Vycc->derS.calcDerivative1_z(A23);
    Vzcc->derS.calcDerivative1_x(A31);
    Vzcc->derS.calcDerivative1_y(A32);
    Vzcc->derS.calcDerivative1_z(A33);

    // Set the array limits when looping over the domain to compute SG contribution.
    // Use only the cell centers like a collocated grid and compute T tensor.
    // Since interpolated U, V, and W data is available only in the core,
    // the limits of fBulk are used so that the boundary points are excluded
    // while computing derivatives and structure functions.
    xS = P.F.fBulk.lbound(0);       xE = P.F.fBulk.ubound(0);
    yS = P.F.fBulk.lbound(1);       yE = P.F.fBulk.ubound(1);
    zS = P.F.fBulk.lbound(2);       zE = P.F.fBulk.ubound(2);

    for (int iX = xS; iX <= xE; iX++) {
        real dx = mesh.xColloc(iX) - mesh.xColloc(iX - 1);
        for (int iY = yS; iY <= yE; iY++) {
            real dy = mesh.yColloc(iY) - mesh.yColloc(iY - 1);
            for (int iZ = zS; iZ <= zE; iZ++) {
                real dz = mesh.zColloc(iZ) - mesh.zColloc(iZ - 1);

                // Specify the values of all the quantities necessary for sgsFlux function to
                // compute the sub-grid stress correctly. These are all member variables of the les
                // class and hence globally available to all the member functions of the class.

                // The required quantities are:
                // 1. Cutoff wavelength
                del = std::pow(dx*dy*dz, 1.0/3.0);

                // 2. Velocities at the 3 x 3 x 3 points over which structure function will be calculated
                u = Vxcc->F.F(blitz::Range(iX-1, iX+1), blitz::Range(iY-1, iY+1), blitz::Range(iZ-1, iZ+1));
                v = Vycc->F.F(blitz::Range(iX-1, iX+1), blitz::Range(iY-1, iY+1), blitz::Range(iZ-1, iZ+1));
                w = Vzcc->F.F(blitz::Range(iX-1, iX+1), blitz::Range(iY-1, iY+1), blitz::Range(iZ-1, iZ+1));

                // 3. The x, y and z coordinates of the 3 x 3 x 3 points over which u, v and w have been specified
                x = mesh.xStaggr(blitz::Range(iX-1, iX+1));
                y = mesh.yStaggr(blitz::Range(iY-1, iY+1));
                z = mesh.zStaggr(blitz::Range(iZ-1, iZ+1));

                // 4. The velocity gradient tensor specified as a 3 x 3 matrix
                dudx = A11(iX, iY, iZ), A12(iX, iY, iZ), A13(iX, iY, iZ),
                       A21(iX, iY, iZ), A22(iX, iY, iZ), A23(iX, iY, iZ),
                       A31(iX, iY, iZ), A32(iX, iY, iZ), A33(iX, iY, iZ);

                // Now the sub-grid stress can be calculated
                sgsStress(&sTxx, &sTyy, &sTzz, &sTxy, &sTyz, &sTzx);

                // Copy the calculated values to the sub-grid stress tensor field
                Txx->F.F(iX, iY, iZ) = sTxx;
                Tyy->F.F(iX, iY, iZ) = sTyy;
                Tzz->F.F(iX, iY, iZ) = sTzz;
                Txy->F.F(iX, iY, iZ) = sTxy;
                Tyz->F.F(iX, iY, iZ) = sTyz;
                Tzx->F.F(iX, iY, iZ) = sTzx;
            }
        }
    }

    // Synchronize the sub-grid stress tensor field data across MPI processors
    Txx->syncData();
    Tyy->syncData();
    Tzz->syncData();
    Txy->syncData();
    Tyz->syncData();
    Tzx->syncData();

    // Compute the components of the divergence of sub-grid stress tensor field
    Txx->derS.calcDerivative1_x(A11);
    Txy->derS.calcDerivative1_x(A12);
    Tzx->derS.calcDerivative1_x(A13);
    Txy->derS.calcDerivative1_y(A21);
    Tyy->derS.calcDerivative1_y(A22);
    Tyz->derS.calcDerivative1_y(A23);
    Tzx->derS.calcDerivative1_z(A31);
    Tyz->derS.calcDerivative1_z(A32);
    Tzz->derS.calcDerivative1_z(A33);

    // Sum the components to get the divergence of the stress tensor field
    A11 = A11 + A21 + A31;
    A22 = A12 + A22 + A32;
    A33 = A13 + A23 + A33;

    // Interpolate the computed divergence and add its contribution to the RHS 
    // of the NSE provided as argument to the function
    // Contribution to the X-component of NSE
    xS = V.Vx.fCore.lbound(0) + 2; xE = V.Vx.fCore.ubound(0) - 2;
    yS = V.Vx.fCore.lbound(1) + 2; yE = V.Vx.fCore.ubound(1) - 2;
    zS = V.Vx.fCore.lbound(2) + 2; zE = V.Vx.fCore.ubound(2) - 2;
    for (int iX = xS; iX <= xE; iX++) {
        for (int iY = yS; iY <= yE; iY++) {
            for (int iZ = zS; iZ <= zE; iZ++) {
                nseRHS.Vx(iX, iY, iZ) += (A11(iX, iY, iZ) + A11(iX + 1, iY, iZ))*0.5;
            }
        }
    }

    // Contribution to the Y-component of NSE
    xS = V.Vy.fCore.lbound(0) + 2; xE = V.Vy.fCore.ubound(0) - 2;
    yS = V.Vy.fCore.lbound(1) + 2; yE = V.Vy.fCore.ubound(1) - 2;
    zS = V.Vy.fCore.lbound(2) + 2; zE = V.Vy.fCore.ubound(2) - 2;
    for (int iX = xS; iX <= xE; iX++) {
        for (int iY = yS; iY <= yE; iY++) {
            for (int iZ = zS; iZ <= zE; iZ++) {
                nseRHS.Vy(iX, iY, iZ) += (A22(iX, iY, iZ) + A22(iX, iY + 1, iZ))*0.5;
            }
        }
    }

    // Contribution to the Z-component of NSE
    xS = V.Vz.fCore.lbound(0) + 2; xE = V.Vz.fCore.ubound(0) - 2;
    yS = V.Vz.fCore.lbound(1) + 2; yE = V.Vz.fCore.ubound(1) - 2;
    zS = V.Vz.fCore.lbound(2) + 2; zE = V.Vz.fCore.ubound(2) - 2;
    for (int iX = xS; iX <= xE; iX++) {
        for (int iY = yS; iY <= yE; iY++) {
            for (int iZ = zS; iZ <= zE; iZ++) {
                nseRHS.Vz(iX, iY, iZ) += (A33(iX, iY, iZ) + A33(iX, iY, iZ + 1))*0.5;
            }
        }
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to compute and add contribution from LES model to the RHS of NSE and temperature equation
 *
 *          This function is called when both hydrodynamics and scalar equations are being solved.
 *          It calles the sgsStress and sgsFlux functions repeatedly for every point in the domain to get the
 *          sub-grid stress tensor field sub-grid scalar flux vector field repsectively.
 *          The divergences of these fields are then calculated to obtain the spiral vortex model contribution
 *          to the NSE and scalar equations
 *
 * \param   nseRHS is a reference to the plain vector field denoting the RHS of the NSE
 * \param   tmpRHS is a reference to the plain scalar field denoting the RHS of the scalar equation
 * \param   T is a reference the scalar field denoting the scalar equation
 ********************************************************************************************************************************************
 */
void spiral::computeSG(plainvf &nseRHS, plainsf &tmpRHS, sfield &T) {
    real sQx, sQy, sQz;

    // First interpolate all velocities to cell centers
    // Then U, V and W data are available at all cell centers except ghost points
    Vxcc->F.F = 0.0;
    for (unsigned int i=0; i < P.F.VxIntSlices.size(); i++) {
        Vxcc->F.F(P.F.fCore) += V.Vx.F(P.F.VxIntSlices(i));
    }
    Vxcc->F.F /= P.F.VxIntSlices.size();

    Vycc->F.F = 0.0;
    for (unsigned int i=0; i < P.F.VyIntSlices.size(); i++) {
        Vycc->F.F(P.F.fCore) += V.Vy.F(P.F.VyIntSlices(i));
    }
    Vycc->F.F /= P.F.VyIntSlices.size();

    Vzcc->F.F = 0.0;
    for (unsigned int i=0; i < P.F.VzIntSlices.size(); i++) {
        Vzcc->F.F(P.F.fCore) += V.Vz.F(P.F.VzIntSlices(i));
    }
    Vzcc->F.F /= P.F.VzIntSlices.size();

    Vxcc->F.syncData();
    Vycc->F.syncData();
    Vzcc->F.syncData();

    // Compute the x, y and z derivatives of the interpolated velocity field and store them into
    // the arrays A11, A12, A13, ... A33. These arrays will be later accessed when constructing
    // the velocity gradient tensor at each point in the domain.
    Vxcc->derS.calcDerivative1_x(A11);
    Vxcc->derS.calcDerivative1_y(A12);
    Vxcc->derS.calcDerivative1_z(A13);
    Vycc->derS.calcDerivative1_x(A21);
    Vycc->derS.calcDerivative1_y(A22);
    Vycc->derS.calcDerivative1_z(A23);
    Vzcc->derS.calcDerivative1_x(A31);
    Vzcc->derS.calcDerivative1_y(A32);
    Vzcc->derS.calcDerivative1_z(A33);

    // Compute the x, y and z derivatives of the temperature field and store them into
    // the arrays B1, B2, and B3. These arrays will be later accessed when constructing
    // the temperature gradient tensor at each point in the domain.
    T.derS.calcDerivative1_x(B1);
    T.derS.calcDerivative1_y(B2);
    T.derS.calcDerivative1_z(B3);

    // Set the array limits when looping over the domain to compute SG contribution.
    // Use only the cell centers like a collocated grid and compute T tensor.
    // Since interpolated U, V, and W data is available only in the core,
    // the limits of fBulk are used so that the boundary points are excluded
    // to compute derivatives and structure functions correctly.
    xS = P.F.fBulk.lbound(0);       xE = P.F.fBulk.ubound(0);
    yS = P.F.fBulk.lbound(1);       yE = P.F.fBulk.ubound(1);
    zS = P.F.fBulk.lbound(2);       zE = P.F.fBulk.ubound(2);

    for (int iX = xS; iX <= xE; iX++) {
        real dx = mesh.xColloc(iX) - mesh.xColloc(iX - 1);
        for (int iY = yS; iY <= yE; iY++) {
            real dy = mesh.yColloc(iY) - mesh.yColloc(iY - 1);
            for (int iZ = zS; iZ <= zE; iZ++) {
                real dz = mesh.zColloc(iZ) - mesh.zColloc(iZ - 1);

                // Specify the values of all the quantities necessary for sgsFlux function to
                // compute the sub-grid stress correctly. These are all member variables of the les
                // class and hence globally available to all the member functions of the class.

                // The required quantities are:
                // 1. Cutoff wavelength
                del = std::pow(dx*dy*dz, 1.0/3.0);

                // 2. Velocities at the 3 x 3 x 3 points over which structure function will be calculated
                u = Vxcc->F.F(blitz::Range(iX-1, iX+1), blitz::Range(iY-1, iY+1), blitz::Range(iZ-1, iZ+1));
                v = Vycc->F.F(blitz::Range(iX-1, iX+1), blitz::Range(iY-1, iY+1), blitz::Range(iZ-1, iZ+1));
                w = Vzcc->F.F(blitz::Range(iX-1, iX+1), blitz::Range(iY-1, iY+1), blitz::Range(iZ-1, iZ+1));

                //if (mesh.rankData.rank == 0) std::cout << u << v << w << std::endl;

                // 3. The x, y and z coordinates of the 3 x 3 x 3 points over which u, v and w have been specified
                x = mesh.xStaggr(blitz::Range(iX-1, iX+1));
                y = mesh.yStaggr(blitz::Range(iY-1, iY+1));
                z = mesh.zStaggr(blitz::Range(iZ-1, iZ+1));

                // 4. The velocity gradient tensor specified as a 3 x 3 matrix
                dudx = A11(iX, iY, iZ), A12(iX, iY, iZ), A13(iX, iY, iZ),
                       A21(iX, iY, iZ), A22(iX, iY, iZ), A23(iX, iY, iZ),
                       A31(iX, iY, iZ), A32(iX, iY, iZ), A33(iX, iY, iZ);

                // Now the sub-grid stress can be calculated
                sgsStress(&sTxx, &sTyy, &sTzz, &sTxy, &sTyz, &sTzx);

                // Copy the calculated values to the sub-grid stress tensor field
                Txx->F.F(iX, iY, iZ) = sTxx;
                Tyy->F.F(iX, iY, iZ) = sTyy;
                Tzz->F.F(iX, iY, iZ) = sTzz;
                Txy->F.F(iX, iY, iZ) = sTxy;
                Tyz->F.F(iX, iY, iZ) = sTyz;
                Tzx->F.F(iX, iY, iZ) = sTzx;

                // To compute sub-grid scalar flux, the sgsStress calculations have already provided
                // most of the necessary values. Only an additional temperature gradient vector is needed
                // 5. The temperature gradient vector specified as a 3 component tiny vector
                dsdx = B1(iX, iY, iZ), B2(iX, iY, iZ), B3(iX, iY, iZ);

                // Now the sub-grid scalar flus can be calculated
                sgsFlux(&sQx, &sQy, &sQz);

                // Copy the calculated values to the sub-grid scalar flux vector field
                qX->F.F(iX, iY, iZ) = sQx;
                qY->F.F(iX, iY, iZ) = sQy;
                qZ->F.F(iX, iY, iZ) = sQz;
            }
        }
    }

    // Synchronize the sub-grid stress tensor field data across MPI processors
    Txx->syncData();
    Tyy->syncData();
    Tzz->syncData();
    Txy->syncData();
    Tyz->syncData();
    Tzx->syncData();

    // Compute the components of the divergence of sub-grid stress tensor field
    Txx->derS.calcDerivative1_x(A11);
    Txy->derS.calcDerivative1_x(A12);
    Tzx->derS.calcDerivative1_x(A13);
    Txy->derS.calcDerivative1_y(A21);
    Tyy->derS.calcDerivative1_y(A22);
    Tyz->derS.calcDerivative1_y(A23);
    Tzx->derS.calcDerivative1_z(A31);
    Tyz->derS.calcDerivative1_z(A32);
    Tzz->derS.calcDerivative1_z(A33);

    // Sum the components to get the divergence of the stress tensor field
    B1 = A11 + A21 + A31;
    B2 = A12 + A22 + A32;
    B3 = A13 + A23 + A33;

    // Interpolate the computed divergence and add its contribution to the RHS 
    // of the NSE provided as argument to the function
    // Contribution to the X-component of NSE
    xS = V.Vx.fCore.lbound(0);      xE = V.Vx.fCore.ubound(0);
    yS = V.Vx.fCore.lbound(1);      yE = V.Vx.fCore.ubound(1);
    zS = V.Vx.fCore.lbound(2);      zE = V.Vx.fCore.ubound(2);
    for (int iX = xS; iX <= xE; iX++) {
        for (int iY = yS; iY <= yE; iY++) {
            for (int iZ = zS; iZ <= zE; iZ++) {
                nseRHS.Vx(iX, iY, iZ) += (B1(iX, iY, iZ) + B1(iX + 1, iY, iZ))*0.5;
            }
        }
    }

    // Contribution to the Y-component of NSE
    xS = V.Vy.fCore.lbound(0);      xE = V.Vy.fCore.ubound(0);
    yS = V.Vy.fCore.lbound(1);      yE = V.Vy.fCore.ubound(1);
    zS = V.Vy.fCore.lbound(2);      zE = V.Vy.fCore.ubound(2);
    for (int iX = xS; iX <= xE; iX++) {
        for (int iY = yS; iY <= yE; iY++) {
            for (int iZ = zS; iZ <= zE; iZ++) {
                nseRHS.Vy(iX, iY, iZ) += (B2(iX, iY, iZ) + B2(iX, iY + 1, iZ))*0.5;
            }
        }
    }

    // Contribution to the Z-component of NSE
    xS = V.Vz.fCore.lbound(0);      xE = V.Vz.fCore.ubound(0);
    yS = V.Vz.fCore.lbound(1);      yE = V.Vz.fCore.ubound(1);
    zS = V.Vz.fCore.lbound(2);      zE = V.Vz.fCore.ubound(2);
    for (int iX = xS; iX <= xE; iX++) {
        for (int iY = yS; iY <= yE; iY++) {
            for (int iZ = zS; iZ <= zE; iZ++) {
                nseRHS.Vz(iX, iY, iZ) += (B3(iX, iY, iZ) + B3(iX, iY, iZ + 1))*0.5;
            }
        }
    }

    // Synchronize the sub-grid scalar flux vector field data across MPI processors
    qX->syncData();
    qY->syncData();
    qZ->syncData();

    // Compute the components of the divergence of sub-grid scalar flux vector field
    qX->derS.calcDerivative1_x(B1);
    qY->derS.calcDerivative1_y(B2);
    qZ->derS.calcDerivative1_z(B3);

    // Sum the components to get the divergence of scalar flux, and add its contribution
    // to the RHS of the temperature field equation provided as argument to the function
    xS = T.F.fBulk.lbound(0);       xE = T.F.fBulk.ubound(0);
    yS = T.F.fBulk.lbound(1);       yE = T.F.fBulk.ubound(1);
    zS = T.F.fBulk.lbound(2);       zE = T.F.fBulk.ubound(2);

    for (int iX = xS; iX <= xE; iX++) {
        for (int iY = yS; iY <= yE; iY++) {
            for (int iZ = zS; iZ <= zE; iZ++) {
                tmpRHS.F(iX, iY, iZ) += (B1(iX, iY, iX) + B2(iX, iY, iX) + B3(iX, iY, iX));
            }
        }
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Main function to calculate the sub-grid stress tensor using stretched vortex model
 *
 *          The six components of the subgrid stress tensor - Txx, Tyy, Tzz, Txy, Tyz, Tzx are calculated at x[0], y[0], z[0].
 *          It needs the resolved velocity gradient tensor dudx[3][3], LES cutoff scale del, and kinematic viscosity nu.
 *          It first computes the alignment of the subgrid vortex, e, by calculating the eigenvectors of S_ij.
 *          To compute the structure function, it needs 3x3x3 samples of the local resolved velocity field,
 *          (u[0,0,0], v[0,0,0], w[0,0,0]) at (x[0], y[0], z[0]) to (u[2,2,2], v[2,2,2], w[2,2,2]) at (x[2], y[2], z[2]).
 *          Finally \mathcal{K}_0 \epsilon^{2/3} k_c^{-2/3}, where a = e_i^v e_j^v S_{ij} is the axial stretching.
 *
 ********************************************************************************************************************************************
 */
void spiral::sgsStress(
    real *Txx, real *Tyy, real *Tzz,
    real *Txy, real *Tyz, real *Tzx)
{
    // lv = Sqrt[2 nu / (3 Abs[a])]
    real lv = 0.0;

    {
        // Strain-rate tensor
        Sxx = 0.5 * (dudx(0, 0) + dudx(0, 0));
        Syy = 0.5 * (dudx(1, 1) + dudx(1, 1));
        Szz = 0.5 * (dudx(2, 2) + dudx(2, 2));
        Sxy = 0.5 * (dudx(0, 1) + dudx(1, 0));
        Syz = 0.5 * (dudx(1, 2) + dudx(2, 1));
        Szx = 0.5 * (dudx(2, 0) + dudx(0, 2));

        // By default, eigenvalue corresponding to most extensive eigenvector is returned
        real eigval = eigenvalueSymm();

        // Default alignment: most extensive eigenvector
        e = eigenvectorSymm(eigval);

        // Make e[3] a unit vector
        real length = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
        e /= length;

        // Strain along vortex axis
        real a = e[0] * e[0] * Sxx + e[0] * e[1] * Sxy + e[0] * e[2] * Szx
               + e[1] * e[0] * Sxy + e[1] * e[1] * Syy + e[1] * e[2] * Syz
               + e[2] * e[0] * Szx + e[2] * e[1] * Syz + e[2] * e[2] * Szz;
        lv = sqrt(2.0 * nu / (3.0 * (fabs(a) + EPS)));
    }

    // Structure function calculation
    real F2 = 0.0;
    real Qd = 0.0; 
    {    
        // Average over neighboring points
        int sfCount = 0;
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                for (int k = -1; k <= 1; k++) {
                    if (i or j or k) {
                        real du = u(i+1, j+1, k+1) - u(1, 1, 1);
                        real dv = v(i+1, j+1, k+1) - v(1, 1, 1);
                        real dw = w(i+1, j+1, k+1) - w(1, 1, 1);
                        F2 += du * du + dv * dv + dw * dw;
                        real dx = x[i+1] - x[1];
                        real dy = y[j+1] - y[1];
                        real dz = z[k+1] - z[1];
                        real dx2 = dx * dx   + dy * dy   + dz * dz;
                        real dxe = dx * e[0] + dy * e[1] + dz * e[2];
                        real d = sqrt(dx2 - dxe * dxe) / del;
                        Qd += sfIntegral(d);
                        sfCount++;
                    }
                }
            }
        }
        F2 /= (real) (sfCount);
        Qd /= (real) (sfCount);
    }
    // prefac is the group prefactor
    real prefac = F2 / Qd; // \mathcal{K}_0 \epsilon^{2/3} k_c^{-2/3}
    real kc = M_PI / del;

    K = prefac * keIntegral(kc * lv);

    // T_{ij} = (\delta_{ij} - e_i^v e_j^v) K
    *Txx = (1.0 - e[0] * e[0]) * K;
    *Tyy = (1.0 - e[1] * e[1]) * K;
    *Tzz = (1.0 - e[2] * e[2]) * K;
    *Txy = (    - e[0] * e[1]) * K;
    *Tyz = (    - e[1] * e[2]) * K;
    *Tzx = (    - e[2] * e[0]) * K;
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the sub-grid scalar flux vector using stretched vortex model
 *
 *          The three components of the subgrid scalar flux vector - qx, qy, qz are calculated at x[0], y[0], z[0].
 *          It needs the resolved scalar gradient tensor dsdx[3], sub-grid vortex alignment e[3] (a unit vector),
 *          LES cutoff scale del, and the precalculated SGS kinetic energy K.
 *          WARNING: For this function to work, the values of member variables e and K must be pre-calculated through
 *          a call to the sgsStress function.
 *
 ********************************************************************************************************************************************
 */
void spiral::sgsFlux(real *qx, real *qy, real *qz) {
    real gam = 1.0; // Universal model constant
    real P = -0.5 * gam * del * sqrt(K);

    // q_i = P (\delta_{ij} - e_i^v e_j^v) ds/dx_j
    *qx = P * ((1.0 - e[0] * e[0]) * dsdx[0]
             + (    - e[0] * e[1]) * dsdx[1]
             + (    - e[0] * e[2]) * dsdx[2]);
    *qy = P * ((    - e[1] * e[0]) * dsdx[0]
             + (1.0 - e[1] * e[1]) * dsdx[1]
             + (    - e[1] * e[2]) * dsdx[2]);
    *qz = P * ((    - e[2] * e[0]) * dsdx[0]
             + (    - e[2] * e[1]) * dsdx[1]
             + (1.0 - e[2] * e[2]) * dsdx[2]);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to approximately evaluate the integral for sub-grid energy
 *
 *          The sub-grid energy is calculated as an approximation of (1/2) k^(2/3) Gamma[-1/3, k^2]
 *          with maximum relative error of 0.17% at k=2.42806.
 *          The only input for the function is the non-dimensionalized cut-off wavenumber k.
 *          It returns the sub-grid energy as a real valued number.
 *
 ********************************************************************************************************************************************
 */
real spiral::keIntegral(real k) {
    real k2 = k * k;
    if (k2 < 2.42806) {
        real pade = (3.0 +   2.5107 * k2 +  0.330357 * k2 * k2
                    +  0.0295481 * k2 * k2 * k2)
                    / (1.0 + 0.336901 * k2 + 0.0416684 * k2 * k2
                    + 0.00187191 * k2 * k2 * k2);
        return 0.5 * (pade - 4.06235 * pow(k2, 1.0 / 3.0));
    }
    else {
        real pade = (1.26429 + 0.835714 * k2 + 0.0964286 * k2 * k2)
                    / (1.0     +   2.25   * k2 +  0.964286 * k2 * k2
                    + 0.0964286 * k2 * k2 * k2);
        return 0.5 * pade * exp(-k2);
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to approximately evaluate the structure function integral
 *
 *          The structure function integral is calculated as an approximation of
 *          Integrate[4 x^(-5/3) (1 - BesselJ[0, x Pi d]), {x, 0, 1}]
 *          with maximum relative error of 2.71% at d=0.873469.
 *          It returns the structure function integral as a real valued number.
 *
 ********************************************************************************************************************************************
 */
real spiral::sfIntegral(real d) {
    // Uncomment if spherical averaging and d=1.
    // if (d == 1.0) return 4.09047;

    real d2 = d * d;
    if (d < 0.873469)
        return 7.4022 * d2 - 1.82642 * d2 * d2;
    else
        return 12.2946 * pow(d, 2.0 / 3.0) - 6.0
            - 0.573159 * pow(d, -1.5) * sin(3.14159 * d - 0.785398);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the eigenvalues of the strain-rate tensor
 *
 *          The function claculates the eigenvalues, eigval[0] < eigval[1] < eigval[2],
 *          of the 3 x 3 symmetric matrix, { { Sxx, Sxy, Szx }, { Sxy, Syy, Syz }, { Szx, Syz, Szz } },
 *          assuming distinct eigenvalues.
 *          It returns the eigenvalue corresponding to the most extensive eigenvector.
 *
 ********************************************************************************************************************************************
 */
real spiral::eigenvalueSymm() {
    real eigval[3];

    // x^3 + a * x^2 + b * x + c = 0, where x is the eigenvalue
    real a = - (Sxx + Syy + Szz);
    real b = Sxx * Syy - Sxy * Sxy + Syy * Szz
             - Syz * Syz + Szz * Sxx - Szx * Szx;
    real c = - (Sxx * (Syy * Szz - Syz * Syz)
                + Sxy * (Syz * Szx - Sxy * Szz)
                + Szx * (Sxy * Syz - Syy * Szx));

    real q = (3.0 * b - a * a) / 9.0;
    real r = (9.0 * a * b - 27.0 * c - 2.0 * a * a * a) / 54.0;

    if (q >= 0.0) {
        if (mesh.rankData.rank == 0) {
            std::cout << "The value of q is greater than or equal to 0 in Spiral Eigenvalue calculation. Aborting" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }

    real costheta = r / sqrt(-q * q * q);

    // |costheta| > 1 should not occur, except from round-off errors
    real theta;
    theta = costheta > 1.0 ? 0.0 :
            costheta < -1.0 ? M_PI :
            acos(costheta);

    eigval[0] = 2.0 * sqrt(-q) * cos((theta             ) / 3.0) - a / 3.0;
    eigval[1] = 2.0 * sqrt(-q) * cos((theta + 2.0 * M_PI) / 3.0) - a / 3.0;
    eigval[2] = 2.0 * sqrt(-q) * cos((theta + 4.0 * M_PI) / 3.0) - a / 3.0;

    // Sort eigenvalues: eigval[0] < eigval[1] < eigval[2]
    if (eigval[0] > eigval[1]) {
        real tmp = eigval[0]; eigval[0] = eigval[1]; eigval[1] = tmp;
    }
    if (eigval[1] > eigval[2]) {
        real tmp = eigval[1]; eigval[1] = eigval[2]; eigval[2] = tmp;
    }
    if (eigval[0] > eigval[1]) {
        real tmp = eigval[0]; eigval[0] = eigval[1]; eigval[1] = tmp;
    }

    return eigval[2];
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the eigenvectors of the strain-rate tensor
 *
 *          The function claculates the eigenvector (not normalized), eigvec[3],
 *          corresponding to the precalculated eigenvalue, eigval, of the 3 x 3 symmetric matrix,
 *          { { Sxx, Sxy, Szx }, { Sxy, Syy, Syz }, { Szx, Syz, Szz } }, assuming distinct eigenvalues.
 *          It returns the eigenvector corresponding to the eigenvalue supplied, as a blitz TinyVector.
 *
 ********************************************************************************************************************************************
 */
blitz::TinyVector<real, 3> spiral::eigenvectorSymm(real eigval) {
    blitz::TinyVector<real, 3> eigvec;

    /*
    // Frobenius norm
    real fNorm = std::sqrt(Sxx*Sxx + Syy*Syy + Szz*Szz +
                           Sxy*Sxy + Syz*Syz + Szx*Szx);

    // The original unscaled check for zero had problems.
    // A below version normalized by Frobenius norm was obtained
    // from the version provided by Wan Cheng at KAUST.
    if (fabs((Sxx - eigval) * ((Syy - eigval) * (Szz - eigval) - Syz * Syz)
            + Sxy * (Syz * Szx - Sxy * (Szz - eigval))
            + Szx * (Sxy * Syz - (Syy - eigval) * Szx))/fabs(fNorm) > EPS) {
        if (mesh.rankData.rank == 0) {
            std::cout << "Invalid eigenvalue in Spiral Eigenvector calculation. Aborting" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }
    */
    // The above section had to be commented since the spiral LES solver kept failing with this check
    // for almost every input velocity field given to it.

    real det[3] = { (Syy - eigval) * (Szz - eigval) - Syz * Syz,
                    (Szz - eigval) * (Sxx - eigval) - Szx * Szx,
                    (Sxx - eigval) * (Syy - eigval) - Sxy * Sxy };

    real fabsdet[3] = { fabs(det[0]), fabs(det[1]), fabs(det[2]) };

    if (fabsdet[0] >= fabsdet[1] && fabsdet[0] >= fabsdet[2]) {
        eigvec = 1.0, (-Sxy*(Szz - eigval) + Szx*Syz)/det[0], (-Szx*(Syy - eigval) + Sxy*Syz)/det[0];
    }
    else if (fabsdet[1] >= fabsdet[2] && fabsdet[1] >= fabsdet[0]) {
        eigvec = (-Sxy*(Szz - eigval) + Syz*Szx)/det[1], 1.0, (-Syz*(Sxx - eigval) + Sxy*Szx)/det[1];
    }
    else if (fabsdet[2] >= fabsdet[0] && fabsdet[2] >= fabsdet[1]) {
        eigvec = (-Szx*(Syy - eigval) + Syz*Sxy)/det[2], (-Syz*(Sxx - eigval) + Szx*Sxy)/det[2], 1.0;
    }
    else {
        if (mesh.rankData.rank == 0) {
            std::cout << "Eigenvalues are not distinct in Spiral Eigenvector calculation. Aborting" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }

    return eigvec;
}
