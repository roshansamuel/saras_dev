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
/*! \file force.h
 *
 *  \brief Class declaration of force
 *
 *  \author Shashwat Bhattacharya, Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef FORCE_H
#define FORCE_H

#include <blitz/array.h>

#include "parallel.h"
#include "plainsf.h"
#include "plainvf.h"
#include "sfield.h"
#include "vfield.h"
#include "parser.h"

class sfield;

class force{
    public:
        force(vfield &U, const parser &solParams, parallel &mpiParam);

        void add_SForce(plainsf &Ht);
        void add_VForce(plainvf &Hv);
        void add_VForce(plainvf &Hv, sfield &T);

    private:
        vfield &V;

        const parallel &mpiData;

        const parser &inputParams;

        double Fb, Fr;

        void add_Coriolis(plainvf &Hv);
        void add_RandomForce(plainvf &Hv);
        void add_Buoyancy(plainvf &Hv, sfield &T);
};

/**
 ********************************************************************************************************************************************
 *  \class force force.h "lib/force/force.h"
 *  \brief Contains all the global variables related to the imposing of forcing, and associated functions
 *
 ********************************************************************************************************************************************
 */

#endif
