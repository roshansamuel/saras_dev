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

        inline void x1Deriv();
        inline void y1Deriv();
        inline void z1Deriv();

        inline void x2Deriv();
        inline void y2Deriv();
        inline void z2Deriv();

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

        field(const grid &gridData, std::string fieldName, const bool xStag, const bool yStag, const bool zStag, const bool allocDerivatives);

        blitz::RectDomain<3> shift(int dim, blitz::RectDomain<3> core, int steps);

        void calcDerivatives1();
        void calcDerivatives2();

        void syncData();

        double fieldMax();

        field& operator += (field &a);
        field& operator -= (field &a);

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
