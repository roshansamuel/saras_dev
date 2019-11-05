#ifndef PLAINSF_H
#define PLAINSF_H

#include "plainvf.h"
#include "sfield.h"
#include "vfield.h"
#include "grid.h"

class plainsf {
    private:
        blitz::firstIndex i;
        blitz::secondIndex j;
        blitz::thirdIndex k;

        const grid &gridData;

    public:
        blitz::Array<double, 3> F;

        blitz::Range xColl, yColl, zColl;

        plainsf(const grid &gridData, const sfield &refF);

        mpidata *mpiHandle;

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
 * \brief   Operator to compute the gradient of the plain scalar field
 *
 *          The gradient operator computes the gradient of the cell centered scalar field, and stores it into a face-centered staggered
 *          plain vector field as defined by the tensor operation:
 *          \f$ \nabla f = \frac{\partial f}{\partial x}i + \frac{\partial f}{\partial y}j + \frac{\partial f}{\partial z}k \f$.
 *
 * \param   gradF is a reference to a plain vector field (plainvf) into which the computed gradient must be written.
 * \param   V is a const reference to a vector field (vfield) whose core slices are used to compute gradient, since plainvf doesn't have them
 ********************************************************************************************************************************************
 */
        inline void gradient(plainvf &gradF, const vfield &V) {
            gradF.Vx(V.Vx.fCore) = gridData.xi_xColloc(xColl)(i)*(F(V.Vx.fCRgt) - F(V.Vx.fCore))/gridData.dXi;
#ifndef PLANAR
            gradF.Vy(V.Vy.fCore) = gridData.et_yColloc(yColl)(j)*(F(V.Vy.fCBak) - F(V.Vy.fCore))/gridData.dEt;
#endif
            gradF.Vz(V.Vz.fCore) = gridData.zt_zColloc(zColl)(k)*(F(V.Vz.fCTop) - F(V.Vz.fCore))/gridData.dZt;
        }

/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronise data across all processors when performing parallel computations
 *
 *          This function calls the \ref mpidata#syncData "syncData" function of mpidata class to perform perform data-transfer and thus update
 *          the sub-domain boundary pads.
 ********************************************************************************************************************************************
 */
        inline void syncData() {
            mpiHandle->syncData();
        }

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
