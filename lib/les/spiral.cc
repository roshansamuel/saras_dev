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

//
// Calculate the SGS stresses (*Txx, *Tyy, *Tzz, *Txy, *Tyz, *Tzx) at
// (x[0], y[0], z[0]) given n (>=2) samples of the local resolved velocity
// field, (u[0], v[0], w[0]) at (x[0], y[0], z[0])
// to (u[n - 1], v[n - 1], w[n - 1]) at (x[n - 1], y[n - 1], z[n - 1]),
// resolved velocity gradient tensor dudx[3][3], LES cutoff scale del,
// and kinematic viscosity nu. If e[3] == { 0.0, 0.0, 0.0 }, overwrite e[3]
// with default alignment.
// \mathcal{K}_0 \epsilon^{2/3} k_c^{-2/3},
// where a = e_i^v e_j^v S_{ij} is the axial stretching.
//
void spiral::sgs_stress(
    double *u, double *v, double *w,
    double *x, double *y, double *z, int n,
    double dudx[3][3], double e[3], double nu, double del,
    double *Txx, double *Tyy, double *Tzz,
    double *Txy, double *Tyz, double *Tzx)
{
    // lv = Sqrt[2 nu / (3 Abs[a])]
    double lv = 0.0;
    {
        // Strain-rate tensor
        double Sxx = 0.5 * (dudx[0][0] + dudx[0][0]);
        double Syy = 0.5 * (dudx[1][1] + dudx[1][1]);
        double Szz = 0.5 * (dudx[2][2] + dudx[2][2]);
        double Sxy = 0.5 * (dudx[0][1] + dudx[1][0]);
        double Syz = 0.5 * (dudx[1][2] + dudx[2][1]);
        double Szx = 0.5 * (dudx[2][0] + dudx[0][2]);

        if (e[0] == 0.0 && e[1] == 0.0 && e[2] == 0.0) {
            double eigval[3];
            eigenvalue_symm(
                Sxx, Syy, Szz, Sxy, Syz, Szx, eigval);
            // Default alignment: most extensive eigenvector
            eigenvector_symm(
                Sxx, Syy, Szz, Sxy, Syz, Szx, eigval[2], e);
        }

        // Make e[3] a unit vector
        double length = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
        e[0] /= length; e[1] /= length; e[2] /= length;

        // Strain along vortex axis
        double a = e[0] * e[0] * Sxx + e[0] * e[1] * Sxy + e[0] * e[2] * Szx
                 + e[1] * e[0] * Sxy + e[1] * e[1] * Syy + e[1] * e[2] * Syz
                 + e[2] * e[0] * Szx + e[2] * e[1] * Syz + e[2] * e[2] * Szz;
        lv = sqrt(2.0 * nu / (3.0 * (fabs(a) + EPS)));
    }
    double F2 = 0.0;
    double Qd = 0.0; 
    {    
        // Average over neighboring points
        for (int i = 1; i < n; i++) {
            double du = u[i] - u[0];
            double dv = v[i] - v[0];
            double dw = w[i] - w[0];
            F2 += du * du + dv * dv + dw * dw;
            double dx = x[i] - x[0];
            double dy = y[i] - y[0];
            double dz = z[i] - z[0];
            double dx2 = dx * dx   + dy * dy   + dz * dz;
            double dxe = dx * e[0] + dy * e[1] + dz * e[2];
            double d = sqrt(dx2 - dxe * dxe) / del;
            Qd += sf_integral(d);
        }
        F2 /= (double) (n - 1);
        Qd /= (double) (n - 1);
    }
    // prefac is the group prefactor
    double prefac = F2 / Qd; // \mathcal{K}_0 \epsilon^{2/3} k_c^{-2/3}
    double kc = M_PI / del;
    double K = prefac * ke_integral(kc * lv);

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
//
void spiral::sgs_flux(
    double dsdx[3], double e[3], double del, double K,
    double *qx, double *qy, double *qz)
{
    double gam = 1.0; // Universal model constant
    double P = -0.5 * gam * del * sqrt(K);

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
double spiral::ke_integral(double k)
{
    double k2 = k * k;
    if (k2 < 2.42806) {
        double pade = (3.0 +   2.5107 * k2 +  0.330357 * k2 * k2
                    +  0.0295481 * k2 * k2 * k2)
                    / (1.0 + 0.336901 * k2 + 0.0416684 * k2 * k2
                    + 0.00187191 * k2 * k2 * k2);
        return 0.5 * (pade - 4.06235 * pow(k2, 1.0 / 3.0));
    }
    else {
        double pade = (1.26429 + 0.835714 * k2 + 0.0964286 * k2 * k2)
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
double spiral::sf_integral(double d)
{
    // Uncomment if spherical averaging and d=1.
    // if (d == 1.0) return 4.09047;

    double d2 = d * d;
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
void spiral::eigenvalue_symm(
    double Sxx, double Syy, double Szz, double Sxy, double Syz, double Szx,
    double eigval[3])
{
    // x^3 + a * x^2 + b * x + c = 0, where x is the eigenvalue
    double a = - (Sxx + Syy + Szz);
    double b = Sxx * Syy - Sxy * Sxy + Syy * Szz
             - Syz * Syz + Szz * Sxx - Szx * Szx;
    double c = - (Sxx * (Syy * Szz - Syz * Syz)
                + Sxy * (Syz * Szx - Sxy * Szz)
                + Szx * (Sxy * Syz - Syy * Szx));

    double q = (3.0 * b - a * a) / 9.0;
    double r = (9.0 * a * b - 27.0 * c - 2.0 * a * a * a) / 54.0;

    if (q >= 0.0) {
        if (mesh.rankData.rank == 0) {
            std::cout << "The value of q is greater than or equal to 0 in Spiral Eigenvalue calculation. Aborting" << std::endl;
            MPI_Finalize();
            exit(0);
        }
    }

    double costheta = r / sqrt(-q * q * q);

    // |costheta| > 1 should not occur, except from round-off errors
    double theta;
    theta = costheta > 1.0 ? 0.0 :
            costheta < -1.0 ? M_PI :
            acos(costheta);

    eigval[0] = 2.0 * sqrt(-q) * cos((theta             ) / 3.0) - a / 3.0;
    eigval[1] = 2.0 * sqrt(-q) * cos((theta + 2.0 * M_PI) / 3.0) - a / 3.0;
    eigval[2] = 2.0 * sqrt(-q) * cos((theta + 4.0 * M_PI) / 3.0) - a / 3.0;

    // Sort eigenvalues: eigval[0] < eigval[1] < eigval[2]
    if (eigval[0] > eigval[1]) {
        double tmp = eigval[0]; eigval[0] = eigval[1]; eigval[1] = tmp;
    }
    if (eigval[1] > eigval[2]) {
        double tmp = eigval[1]; eigval[1] = eigval[2]; eigval[2] = tmp;
    }
    if (eigval[0] > eigval[1]) {
        double tmp = eigval[0]; eigval[0] = eigval[1]; eigval[1] = tmp;
    }
}

//
// Calculate the eigenvector (not normalized), eigvec[3],
// corresponding to the precalculated eigenvalue, eigval,
// of the 3 x 3 symmetric matrix,
// { { Sxx, Sxy, Szx }, { Sxy, Syy, Syz }, { Szx, Syz, Szz } },
// assuming distinct eigenvalues.
//
void spiral::eigenvector_symm(
    double Sxx, double Syy, double Szz, double Sxy, double Syz, double Szx,
    double eigval, double eigvec[3])
{
    // There are problems with this unscaled check for zero.
    // A better one would be to check the zero against the det(S), for instance.
    if (fabs((Sxx - eigval) * ((Syy - eigval) * (Szz - eigval) - Syz * Syz)
            + Sxy * (Syz * Szx - Sxy * (Szz - eigval))
            + Szx * (Sxy * Syz - (Syy - eigval) * Szx)) > EPS) {
        if (mesh.rankData.rank == 0) {
            std::cout << "Invalid eigenvalue in Spiral Eigenvector calculation. Aborting" << std::endl;
            MPI_Finalize();
            exit(0);
        }
    }

    double det[3] = { (Syy - eigval) * (Szz - eigval) - Syz * Syz,
                      (Szz - eigval) * (Sxx - eigval) - Szx * Szx,
                      (Sxx - eigval) * (Syy - eigval) - Sxy * Sxy };

    double fabsdet[3] = { fabs(det[0]), fabs(det[1]), fabs(det[2]) };

    if (fabsdet[0] >= fabsdet[1] && fabsdet[0] >= fabsdet[2]) {
        eigvec[0] = 1.0;
        eigvec[1] = (-Sxy * (Szz - eigval) + Szx * Syz) / det[0];
        eigvec[2] = (-Szx * (Syy - eigval) + Sxy * Syz) / det[0];
    }
    else if (fabsdet[1] >= fabsdet[2] && fabsdet[1] >= fabsdet[0]) {
        eigvec[0] = (-Sxy * (Szz - eigval) + Syz * Szx) / det[1];
        eigvec[1] = 1.0;
        eigvec[2] = (-Syz * (Sxx - eigval) + Sxy * Szx) / det[1];
    }
    else if (fabsdet[2] >= fabsdet[0] && fabsdet[2] >= fabsdet[1]) {
        eigvec[0] = (-Szx * (Syy - eigval) + Syz * Sxy) / det[2];
        eigvec[1] = (-Syz * (Sxx - eigval) + Szx * Sxy) / det[2];
        eigvec[2] = 1.0;
    }
    else {
        if (mesh.rankData.rank == 0) {
            std::cout << "Eigenvalues are not distinct in Spiral Eigenvector calculation. Aborting" << std::endl;
            MPI_Finalize();
            exit(0);
        }
    }
}
