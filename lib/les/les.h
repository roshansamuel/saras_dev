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
/*! \file les.h
 *
 *  \brief Class declaration of LES Modules
 *
 *  \author Roshan Samuel
 *  \date Sep 2020
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

/********************************************************************************************************************************************
 *
 * The stretched-vortex LES model defined by the derived class, spiral, is adapted
 * from the code provided by Dale I Pullin at Caltech.
 * A detailed description of the module can be found in doc/spiral.pdf
 *
 ********************************************************************************************************************************************
 */

#ifndef LES_H
#define LES_H

#include "grid.h"

class les {
    public:
        les(const grid &mesh);

    protected:
        const grid &mesh;
};

/**
 ********************************************************************************************************************************************
 *  \class les les.h "lib/les/les.h"
 *  \brief Contains all the global variables related to the LES models used by SARAS
 *
 ********************************************************************************************************************************************
 */

class spiral: public les {
    public:
        spiral(const grid &mesh);

        void sgs_stress(
            blitz::Array<double, 3> u,
            blitz::Array<double, 3> v,
            blitz::Array<double, 3> w,
            blitz::Array<double, 2> dudx,
            double *x, double *y, double *z,
            double nu, double del,
            double *Txx, double *Tyy, double *Tzz,
            double *Txy, double *Tyz, double *Tzx);

    private:
        double Sxx, Syy, Szz, Sxy, Syz, Szx;

        double ke_integral(double k);

        double sf_integral(double d);

        double eigenvalue_symm();

        blitz::TinyVector<double, 3> eigenvector_symm(double eigval);
};

/**
 ********************************************************************************************************************************************
 *  \class spiral les.h "lib/les/les.h"
 *  \brief The derived class from les to implement the stretched spiral vortex LES model
 *
 ********************************************************************************************************************************************
 */

#endif
