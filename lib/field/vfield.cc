#include "sfield.h"
#include "vfield.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the vfield class
 *
 *          Three instances of the sfield class to store the data of the three component scalar fields are initialized.
 *          The scalar fields are initialized with appropriate grid staggering to place the components on the cell faces.
 *          The name for the vector field as given by the user is also assigned.
 *
 * \param   gridData is a const reference to the global data contained in the grid class
 * \param   fieldName is a string value set by the user to name and identify the vector field
 ********************************************************************************************************************************************
 */
vfield::vfield(const grid &gridData, std::string fieldName, const bool allocDerivatives):
               gridData(gridData),
               Vx(gridData, "Vx", false, true, true, allocDerivatives),
               Vy(gridData, "Vy", true, false, true, allocDerivatives),
               Vz(gridData, "Vz", true, true, false, allocDerivatives)
{
    this->fieldName = fieldName;

    if (allocDerivatives) {
        interVx2Vx.resize(Vx.fSize);
        interVx2Vx.reindexSelf(Vx.flBound);

        interVx2Vy.resize(Vy.fSize);
        interVx2Vy.reindexSelf(Vy.flBound);

        interVx2Vz.resize(Vz.fSize);
        interVx2Vz.reindexSelf(Vz.flBound);

#ifndef PLANAR
        interVy2Vx.resize(Vx.fSize);
        interVy2Vx.reindexSelf(Vx.flBound);

        interVy2Vy.resize(Vy.fSize);
        interVy2Vy.reindexSelf(Vy.flBound);

        interVy2Vz.resize(Vz.fSize);
        interVy2Vz.reindexSelf(Vz.flBound);
#endif

        interVz2Vx.resize(Vx.fSize);
        interVz2Vx.reindexSelf(Vx.flBound);

        interVz2Vy.resize(Vy.fSize);
        interVz2Vy.reindexSelf(Vy.flBound);

        interVz2Vz.resize(Vz.fSize);
        interVz2Vz.reindexSelf(Vz.flBound);

        // Below array is used only in scalar solvers
        interPc2Vz.resize(Vz.fSize);
        interPc2Vz.reindexSelf(Vz.flBound);
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the diffusion term
 *
 *          It is assumed that the velocity is specified at face-centers, as required by the \ref sfield#computeNLin
 *          "computeNLin" function of sfield.
 *
 * \param   H is a pointer to a vector field (vfield) to which the output of the function is to be written
 ********************************************************************************************************************************************
 */
void vfield::computeDiff(vfield &H) {
    Vx.calcDerivatives2();
    Vy.calcDerivatives2();
    Vz.calcDerivatives2();

#ifdef PLANAR
    H.Vx.F = Vx.d2F_dx2 + Vx.d2F_dz2;
    H.Vy.F = Vy.d2F_dx2 + Vy.d2F_dz2;
    H.Vz.F = Vz.d2F_dx2 + Vz.d2F_dz2;
#else
    H.Vx.F = Vx.d2F_dx2 + Vx.d2F_dy2 + Vx.d2F_dz2;
    H.Vy.F = Vy.d2F_dx2 + Vy.d2F_dy2 + Vy.d2F_dz2;
    H.Vz.F = Vz.d2F_dx2 + Vz.d2F_dy2 + Vz.d2F_dz2;
#endif
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the convective derivative of the vector field
 *
 *          The function computes for the operator \f$ (\mathbf{u}.\nabla)\mathbf{v} \f$ for the vector field, \f$\mathbf{v}\f$.
 *          To do so, the function needs the vector field (vfield) of velocity, \f$\mathbf{u}\f$.
 *          This value is used in three separate calls to the \ref sfield#computeNLin "computeNLin" function of sfield
 *          to compute the derivatives for the three components of the vector field.
 *          It is assumed that the velocity is specified at face-centers, as required by the \ref sfield#computeNLin
 *          "computeNLin" function of sfield.
 *
 * \param   V is a const reference to a vector field (vfield) that specifies the convection velocity at each point
 * \param   H is a pointer to a vector field (vfield) to which the output of the function is to be written
 ********************************************************************************************************************************************
 */
void vfield::computeNLin(const vfield &V, vfield &H) {
    // Compute non-linear term for the Vx component
    Vx.calcDerivatives1();

    interVx2Vx = 0.0;
    for (unsigned int i=0; i < Vx.VxIntSlices.size(); i++) {
        interVx2Vx(Vx.fCore) += V.Vx.F(Vx.VxIntSlices(i));
    }
#ifndef PLANAR
    interVy2Vx = 0.0;
    for (unsigned int i=0; i < Vx.VyIntSlices.size(); i++) {
        interVy2Vx(Vx.fCore) += V.Vy.F(Vx.VyIntSlices(i));
    }
#endif
    interVz2Vx = 0.0;
    for (unsigned int i=0; i < Vx.VzIntSlices.size(); i++) {
        interVz2Vx(Vx.fCore) += V.Vz.F(Vx.VzIntSlices(i));
    }

#ifdef PLANAR
    H.Vx.F(Vx.fCore) -= interVx2Vx(Vx.fCore)*Vx.d1F_dx1(Vx.fCore)/Vx.VxIntSlices.size() +
                        interVz2Vx(Vx.fCore)*Vx.d1F_dz1(Vx.fCore)/Vx.VzIntSlices.size();
#else
    H.Vx.F(Vx.fCore) -= interVx2Vx(Vx.fCore)*Vx.d1F_dx1(Vx.fCore)/Vx.VxIntSlices.size() +
                        interVy2Vx(Vx.fCore)*Vx.d1F_dy1(Vx.fCore)/Vx.VyIntSlices.size() +
                        interVz2Vx(Vx.fCore)*Vx.d1F_dz1(Vx.fCore)/Vx.VzIntSlices.size();
#endif

    // Compute non-linear term for the Vy component
#ifndef PLANAR
    Vy.calcDerivatives1();

    interVx2Vy = 0.0;
    for (unsigned int i=0; i < Vy.VxIntSlices.size(); i++) {
        interVx2Vy(Vy.fCore) += V.Vx.F(Vy.VxIntSlices(i));
    }
    interVy2Vy = 0.0;
    for (unsigned int i=0; i < Vy.VyIntSlices.size(); i++) {
        interVy2Vy(Vy.fCore) += V.Vy.F(Vy.VyIntSlices(i));
    }
    interVz2Vy = 0.0;
    for (unsigned int i=0; i < Vy.VzIntSlices.size(); i++) {
        interVz2Vy(Vy.fCore) += V.Vz.F(Vy.VzIntSlices(i));
    }

    H.Vy.F(Vy.fCore) -= interVx2Vy(Vy.fCore)*Vy.d1F_dx1(Vy.fCore)/Vy.VxIntSlices.size() +
                        interVy2Vy(Vy.fCore)*Vy.d1F_dy1(Vy.fCore)/Vy.VyIntSlices.size() +
                        interVz2Vy(Vy.fCore)*Vy.d1F_dz1(Vy.fCore)/Vy.VzIntSlices.size();
#endif

    // Compute non-linear term for the Vz component
    Vz.calcDerivatives1();

    interVx2Vz = 0.0;
    for (unsigned int i=0; i < Vz.VxIntSlices.size(); i++) {
        interVx2Vz(Vz.fCore) += V.Vx.F(Vz.VxIntSlices(i));
    }
#ifndef PLANAR
    interVy2Vz = 0.0;
    for (unsigned int i=0; i < Vz.VyIntSlices.size(); i++) {
        interVy2Vz(Vz.fCore) += V.Vy.F(Vz.VyIntSlices(i));
    }
#endif
    interVz2Vz = 0.0;
    for (unsigned int i=0; i < Vz.VzIntSlices.size(); i++) {
        interVz2Vz(Vz.fCore) += V.Vz.F(Vz.VzIntSlices(i));
    }

#ifdef PLANAR
    H.Vz.F(Vz.fCore) -= interVx2Vz(Vz.fCore)*Vz.d1F_dx1(Vz.fCore)/Vz.VxIntSlices.size() +
                        interVz2Vz(Vz.fCore)*Vz.d1F_dz1(Vz.fCore)/Vz.VzIntSlices.size();
#else
    H.Vz.F(Vz.fCore) -= interVx2Vz(Vz.fCore)*Vz.d1F_dx1(Vz.fCore)/Vz.VxIntSlices.size() +
                        interVy2Vz(Vz.fCore)*Vz.d1F_dy1(Vz.fCore)/Vz.VyIntSlices.size() +
                        interVz2Vz(Vz.fCore)*Vz.d1F_dz1(Vz.fCore)/Vz.VzIntSlices.size();
#endif
}

/**
 ********************************************************************************************************************************************
 * \brief   Operator to compute the divergence of the vector field
 *
 *          The operator computes the divergence of a face-centered staggered vector field, and stores it into a cell centered
 *          scalar field as defined by the tensor operation:
 *          \f$ \nabla . \mathbf{v} = \frac{\partial \mathbf{v}}{\partial x} +
 *                                    \frac{\partial \mathbf{v}}{\partial y} +
 *                                    \frac{\partial \mathbf{v}}{\partial z} \f$.
 *
 * \param   divV is a pointer to a scalar field (sfield) into which the computed divergence must be written.
 ********************************************************************************************************************************************
 */
void vfield::divergence(sfield &divV) {
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::thirdIndex k;

    blitz::Range xStag, yStag, zStag;

    xStag = blitz::Range(divV.F.fCore.lbound(0), divV.F.fCore.ubound(0));
    yStag = blitz::Range(divV.F.fCore.lbound(1), divV.F.fCore.ubound(1));
    zStag = blitz::Range(divV.F.fCore.lbound(2), divV.F.fCore.ubound(2));

    divV.F.F = 0.0;

#ifdef PLANAR
    divV.F.F(divV.F.fCore) = gridData.xi_xStaggr(xStag)(i)*(Vx.F(divV.F.fCore) - Vx.F(divV.F.fCLft))/gridData.dXi + 
                             gridData.zt_zStaggr(zStag)(k)*(Vz.F(divV.F.fCore) - Vz.F(divV.F.fCBot))/gridData.dZt;
#else
    divV.F.F(divV.F.fCore) = gridData.xi_xStaggr(xStag)(i)*(Vx.F(divV.F.fCore) - Vx.F(divV.F.fCLft))/gridData.dXi + 
                             gridData.et_yStaggr(yStag)(j)*(Vy.F(divV.F.fCore) - Vy.F(divV.F.fCFrt))/gridData.dEt + 
                             gridData.zt_zStaggr(zStag)(k)*(Vz.F(divV.F.fCore) - Vz.F(divV.F.fCBot))/gridData.dZt;
#endif
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronise data across all processors when performing parallel computations
 *
 *          Each of the individual scalar field components have their own subroutine, \ref sfield#syncData "syncData" to send and
 *          receive data across its MPI decomposed sub-domains.
 *          This function calls the \ref sfield#syncData "syncData" function of its components to update the sub-domain boundary pads.
 ********************************************************************************************************************************************
 */
void vfield::syncData() {
    Vx.syncData();
    Vy.syncData();
    Vz.syncData();
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given vector field
 *
 *          The unary operator += adds a given vector field to the entire field stored as vfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another vfield to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
vfield& vfield::operator += (vfield &a) {
    Vx.F += a.Vx.F;
    Vy.F += a.Vy.F;
    Vz.F += a.Vz.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given vector field
 *
 *          The unary operator -= subtracts a given vector field from the entire field stored as vfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another vfield to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
vfield& vfield::operator -= (vfield &a) {
    Vx.F -= a.Vx.F;
    Vy.F -= a.Vy.F;
    Vz.F -= a.Vz.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to multiply a scalar value to the vector field
 *
 *          The unary operator *= multiplies a double precision value to all the fields (Vx, Vy and Vz) stored in vfield and returns
 *          a pointer to itself.
 *
 * \param   a is a double precision number to be multiplied to the vector field
 *
 * \return  A pointer to itself is returned by the vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
vfield& vfield::operator *= (double a) {
    Vx.F *= a;
    Vy.F *= a;
    Vz.F *= a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar value to the vector field
 *
 *          The operator = assigns a double precision value to all the fields (Vx, Vy and Vz) stored in vfield.
 *
 * \param   a is a double precision number to be assigned to the vector field
 ********************************************************************************************************************************************
 */
void vfield::operator = (double a) {
    Vx.F = a;
    Vy.F = a;
    Vz.F = a;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign another vector field to the vector field
 *
 *          The operator = assigns all the three scalar sub-fields of a vfield to all the corresponding fields (Vx, Vy and Vz).
 *
 * \param   a is a vfield to be assigned to the vector field
 ********************************************************************************************************************************************
 */
void vfield::operator = (vfield &a) {
    Vx.F = a.Vx.F;
    Vy.F = a.Vy.F;
    Vz.F = a.Vz.F;
}
