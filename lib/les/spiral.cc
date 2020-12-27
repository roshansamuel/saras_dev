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
 *          The empty constructer merely initializes the local reference to the global mesh variable.
 *
 * \param   mesh is a const reference to the global data contained in the grid class.
 ********************************************************************************************************************************************
 */
spiral::spiral(const grid &mesh): les(mesh) { }

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
void spiral::sgs_stress(
    blitz::Array<real, 3> u,
    blitz::Array<real, 3> v,
    blitz::Array<real, 3> w,
    blitz::Array<real, 2> dudx,
    real *x, real *y, real *z,
    real nu, real del,
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
        real eigval = eigenvalue_symm();

        // Default alignment: most extensive eigenvector
        e = eigenvector_symm(eigval);

        // Make e[3] a unit vector
        real length = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
        e /= length;

        // Strain along vortex axis
        real a = e[0] * e[0] * Sxx + e[0] * e[1] * Sxy + e[0] * e[2] * Szx
                 + e[1] * e[0] * Sxy + e[1] * e[1] * Syy + e[1] * e[2] * Syz
                 + e[2] * e[0] * Szx + e[2] * e[1] * Syz + e[2] * e[2] * Szz;
        lv = sqrt(2.0 * nu / (3.0 * (fabs(a) + EPS)));
    }

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
                        Qd += sf_integral(d);
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

    K = prefac * ke_integral(kc * lv);

    // T_{ij} = (\delta_{ij} - e_i^v e_j^v) K
    *Txx = (1.0 - e[0] * e[0]) * K;
    *Tyy = (1.0 - e[1] * e[1]) * K;
    *Tzz = (1.0 - e[2] * e[2]) * K;
    *Txy = (    - e[0] * e[1]) * K;
    *Tyz = (    - e[1] * e[2]) * K;
    *Tzx = (    - e[2] * e[0]) * K;
}

//
// Calculate the SGS scalar flux (*qx, *qy, *qz) given the resolved scalar
// gradient dsdx[3], vortex alignment e[3] (a unit vector), LES cutoff
// scale del and precalculated SGS kinetic energy K.
// WARNING: For this function to work, the values of member variables e and K
// must be pre-caclculated through a call to the function sgs_stress.
//
void spiral::sgs_flux(
        blitz::TinyVector<real, 3> dsdx,
        real del, real *qx, real *qy, real *qz)
{
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

//
// Approximation of
// (1/2) k^(2/3) Gamma[-1/3, k^2]
// with maximum relative error of 0.17% at k=2.42806.
//
real spiral::ke_integral(real k)
{
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

//
// Approximation of
// Integrate[4 x^(-5/3) (1 - BesselJ[0, x Pi d]), {x, 0, 1}]
// with maximum relative error of 2.71% at d=0.873469.
//
real spiral::sf_integral(real d)
{
    // Uncomment if spherical averaging and d=1.
    // if (d == 1.0) return 4.09047;

    real d2 = d * d;
    if (d < 0.873469)
        return 7.4022 * d2 - 1.82642 * d2 * d2;
    else
        return 12.2946 * pow(d, 2.0 / 3.0) - 6.0
            - 0.573159 * pow(d, -1.5) * sin(3.14159 * d - 0.785398);
}

//
// Calculate the eigenvalues, eigval[0] < eigval[1] < eigval[2],
// of the 3 x 3 symmetric matrix,
// { { Sxx, Sxy, Szx }, { Sxy, Syy, Syz }, { Szx, Syz, Szz } },
// assuming distinct eigenvalues.
//
real spiral::eigenvalue_symm() {
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

//
// Calculate the eigenvector (not normalized), eigvec[3],
// corresponding to the precalculated eigenvalue, eigval,
// of the 3 x 3 symmetric matrix,
// { { Sxx, Sxy, Szx }, { Sxy, Syy, Syz }, { Szx, Syz, Szz } },
// assuming distinct eigenvalues.
//
blitz::TinyVector<real, 3> spiral::eigenvector_symm(real eigval) {
    blitz::TinyVector<real, 3> eigvec;

    // There are problems with this unscaled check for zero.
    // A better one would be to check the zero against the det(S), for instance.
    //real compVal;
    //compVal = fabs((Sxx - eigval) * ((Syy - eigval) * (Szz - eigval) - Syz * Syz)
    //        + Sxy * (Syz * Szx - Sxy * (Szz - eigval))
    //        + Szx * (Sxy * Syz - (Syy - eigval) * Szx));
    //if (compVal < EPS) {
    //    if (mesh.rankData.rank == 0) {
    //        std::cout << "Invalid eigenvalue in Spiral Eigenvector calculation. Aborting" << std::endl;
    //    }
    //    MPI_Finalize();
    //    exit(0);
    //}

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
