#ifndef PLAINSF_H
#define PLAINSF_H

#include "sfield.h"
#include "grid.h"

class plainsf {
    private:
        const grid &gridData;

    public:
        blitz::Array<double, 3> F;

        plainsf(const grid &gridData, const sfield &refF);

        inline double fieldMax();

        plainsf& operator += (plainsf &a);
        plainsf& operator -= (plainsf &a);

        plainsf& operator += (sfield &a);
        plainsf& operator -= (sfield &a);

        plainsf& operator *= (double a);

        void operator = (plainsf &a);
        void operator = (sfield &a);

        void operator = (double a);

/**
 ********************************************************************************************************************************************
 * \brief   Function to extract the maximum value from the plain scalar field
 *
 *          The function uses the in-built blitz function to obtain the maximum value in an array.
 *          While performing parallel computation, the function performs an <B>MPI_Allreduce()</B> to get
 *          the global maximum from the entire computational domain.
 *
 * \return  The double precision value of the maximum is returned (it is implicitly assumed that only double precision values are used)
 ********************************************************************************************************************************************
 */
        inline double fxMax() {
            double localMax, globalMax;

            localMax = blitz::max(F);

            MPI_Allreduce(&localMax, &globalMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD);

            return globalMax;
        }

        ~plainsf() { };
};

/**
 ********************************************************************************************************************************************
 *  \class plainsf plainsf.h "lib/plainsf.h"
 *  \brief Plain scalar field class to store simple scalar fields with no differentiation or interpolation
 *
 *  The class stores scalar fields in the form of a Blitz array
 ********************************************************************************************************************************************
 */

#endif
