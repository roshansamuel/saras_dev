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
/*! \file field.h
 *
 *  \brief Class declaration of field
 *
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef FIELD_H
#define FIELD_H

#include <blitz/array.h>
#include <string>

#include "mpidata.h"
#include "grid.h"

class field {
    private:
        const grid &gridData;

        void setCoreSlice();
        void setBulkSlice();

        void setWallSlices();

        void setInterpolationSlices();

    public:
        blitz::Array<double, 3> F;

        std::string fieldName;

        const bool xStag, yStag, zStag;

        blitz::RectDomain<3> fCore, fBulk;
        blitz::RectDomain<3> fCLft, fCRgt;
        blitz::RectDomain<3> fCFrt, fCBak;
        blitz::RectDomain<3> fCBot, fCTop;

        blitz::Array<blitz::RectDomain<3>, 1> fWalls;

        blitz::Array<blitz::RectDomain<3>, 1> PcIntSlices, QvIntSlices;
        blitz::Array<blitz::RectDomain<3>, 1> VxIntSlices, VyIntSlices, VzIntSlices;
        blitz::Array<blitz::RectDomain<3>, 1> WxIntSlices, WyIntSlices, WzIntSlices;

        blitz::TinyVector<int, 3> fSize;
        blitz::TinyVector<int, 3> flBound, cuBound;

        mpidata *mpiHandle;

        field(const grid &gridData, std::string fieldName, const bool xStag, const bool yStag, const bool zStag);

        void syncData();

        double fieldMax();

        blitz::RectDomain<3> shift(int dim, blitz::RectDomain<3> core, int steps);

        field& operator += (field &a);
        field& operator -= (field &a);

        field& operator += (double a);
        field& operator -= (double a);

        void operator = (field &a);
        void operator = (double a);

        ~field();
};

/**
 ********************************************************************************************************************************************
 *  \class field field.h "lib/field.h"
 *  \brief Field class to store data and perform finite difference operations on the data
 *
 *  The class stores the base data of both scalar and vector fields as blitz arrays.
 *  The data is stored with a uniform grid spacing as in the transformed plane.
 *  The limits of the full domain and its core are also stored in a set of RectDomain and TinyVector objects.
 ********************************************************************************************************************************************
 */

#endif
