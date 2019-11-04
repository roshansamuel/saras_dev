#ifndef FIELD_H
#define FIELD_H

#include <blitz/array.h>
#include <string>

#include "mpidata.h"
#include "differ.h"
#include "grid.h"

class field {
    private:
        const grid &gridData;

        blitz::Array<double, 1> x_Metric, y_Metric, z_Metric;
        blitz::Array<double, 1> xxMetric, yyMetric, zzMetric;
        blitz::Array<double, 1> x2Metric, y2Metric, z2Metric;

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

        // The following public arrays for getting derivatives of variables are available *only if allocDerivatives flag is set to true*
        // Attempting to use these arrays of a field with allocDerivatives set to false may give seg-fault!
        blitz::Array<double, 3> d1F_dx1, d1F_dy1, d1F_dz1;
        blitz::Array<double, 3> d2F_dx2, d2F_dy2, d2F_dz2;

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

        differ xDim, yDim, zDim;

        field(const grid &gridData, std::string fieldName, const bool xStag, const bool yStag, const bool zStag);

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
 *  Correspondingly, the finite difference operations are also performed with constant grid spacing.
 *  The limits of the full domain and its core are also stored in a set of RectDomain and TinyVector objects.
 ********************************************************************************************************************************************
 */

#endif
