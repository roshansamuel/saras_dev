#include "sfield.h"
#include "vfield.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the sfield class
 *
 *          The instance of the field class to store the data of the scalar field is initialized, and the necessary grid
 *          transformation derivatives along each direction are chosen according to the grid staggering.
 *          The arrays to store the output from various operators like derivatives, convective derivatives, etc. are also
 *          allocated.
 *          Finally, an instance of the <B>mpidata</B> class is initialized to store the sub-arrays to be send/received
 *          across the processors during MPI communication.
 *
 * \param   gridData is a const reference to the global data contained in the grid class
 * \param   fieldName is a string value set by the user to name and identify the scalar field
 ********************************************************************************************************************************************
 */
sfield::sfield(const grid &gridData, std::string fieldName, const bool allocDerivatives):
               gridData(gridData),
               F(gridData, fieldName, true, true, true, allocDerivatives),
               derS(gridData, F)
{
    this->fieldName = fieldName;

    tempMat.resize(F.fSize);
    tempMat.reindexSelf(F.flBound);

    if (allocDerivatives) {
        interVx.resize(F.fSize);
        interVx.reindexSelf(F.flBound);

#ifndef PLANAR
        interVy.resize(F.fSize);
        interVy.reindexSelf(F.flBound);
#endif

        interVz.resize(F.fSize);
        interVz.reindexSelf(F.flBound);
    }

    tempMat = 0.0;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the diffusion term
 *
 *          It is assumed that the velocity is specified at face-centers, as required by the \ref sfield#computeNLin
 *          "computeNLin" function of sfield.
 *
 * \param   H is a pointer to a scalar field (sfield) to which the output of the function is to be written
 ********************************************************************************************************************************************
 */
void sfield::computeDiff(sfield &H) {
    
    derS.calcDerivative2_xx(tempMat);
    H.F.F(F.fCore) = tempMat(F.fCore);
    tempMat=0.0;
    
#ifdef PLANAR
    derS.calcDerivative2_yy(tempMat);
    H.F.F(F.fCore) += tempMat(F.fCore);
    tempMat=0.0;
#endif

    derS.calcDerivative2_zz(tempMat);
    H.F.F(F.fCore) += tempMat(F.fCore);
    tempMat=0.0;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the convective derivative of the scalar field
 *
 *          The function computes for the operator \f$ (\mathbf{u}.\nabla)f \f$ at the grid nodes of the scalar field f.
 *          To do so, the function needs the vector field (vfield) of velocity. It is assumed that the velocity is always
 *          specified at face-centers, and is interpolated accordingly to the scalar field grid points.
 *
 * \param   V is a const reference to a vector field (vfield) that specifies the convection velocity at each point
 ********************************************************************************************************************************************
 */
void sfield::computeNLin(const vfield &V, sfield &H) {
    
    interVx = 0.0;
    for (unsigned int i=0; i < F.VxIntSlices.size(); i++) {
        interVx(F.fCore) += V.Vx.F(F.VxIntSlices(i));
    }
#ifndef PLANAR
    interVy = 0.0;
    for (unsigned int i=0; i < F.VyIntSlices.size(); i++) {
        interVy(F.fCore) += V.Vy.F(F.VyIntSlices(i));
    }
#endif
    interVz = 0.0;
    for (unsigned int i=0; i < F.VzIntSlices.size(); i++) {
        interVz(F.fCore) += V.Vz.F(F.VzIntSlices(i));
    }

    derS.calcDerivative1_x(tempMat);
    H.F.F(F.fCore) -= interVx(F.fCore)*tempMat(F.fCore)/F.VxIntSlices.size();
    tempMat=0.0;

#ifdef PLANAR
    derS.calcDerivative1_y(tempMat);
    H.F.F(F.fCore) -= interVy(F.fCore)*tempMat(F.fCore)/F.VyIntSlices.size();
    tempMat=0.0;
#endif

    derS.calcDerivative1_z(tempMat);
    H.F.F(F.fCore) -= interVz(F.fCore)*tempMat(F.fCore)/F.VzIntSlices.size();
    tempMat=0.0;
}

/**
 ********************************************************************************************************************************************
 * \brief   Operator to compute the gradient of the scalar field
 *
 *          The gradient operator computes the gradient of the cell centered scalar field, and stores it into a face-centered staggered
 *          vector field as defined by the tensor operation:
 *          \f$ \nabla f = \frac{\partial f}{\partial x}i + \frac{\partial f}{\partial y}j + \frac{\partial f}{\partial z}k \f$.
 *
 * \param   gradF is a pointer to a vector field (vfield) into which the computed gradient must be written.
 ********************************************************************************************************************************************
 */
void sfield::gradient(vfield &gradF) {
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::thirdIndex k;

    blitz::Range xColl, yColl, zColl;

    xColl = blitz::Range(gridData.collocCoreDomain.lbound(0), gridData.collocCoreDomain.ubound(0));
    yColl = blitz::Range(gridData.collocCoreDomain.lbound(1), gridData.collocCoreDomain.ubound(1));
    zColl = blitz::Range(gridData.collocCoreDomain.lbound(2), gridData.collocCoreDomain.ubound(2));

    gradF.Vx.F(gradF.Vx.fCore) = gridData.xi_xColloc(xColl)(i)*(F.F(gradF.Vx.fCRgt) - F.F(gradF.Vx.fCore))/gridData.dXi;
#ifndef PLANAR
    gradF.Vy.F(gradF.Vy.fCore) = gridData.et_yColloc(yColl)(j)*(F.F(gradF.Vy.fCBak) - F.F(gradF.Vy.fCore))/gridData.dEt;
#endif
    gradF.Vz.F(gradF.Vz.fCore) = gridData.zt_zColloc(zColl)(k)*(F.F(gradF.Vz.fCTop) - F.F(gradF.Vz.fCore))/gridData.dZt;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronise data across all processors when performing parallel computations
 *
 *          This function calls the \ref mpidata#syncData "syncData" function of mpidata class to perform perform data-transfer and thus update
 *          the sub-domain boundary pads.
 ********************************************************************************************************************************************
 */
void sfield::syncData() {
    F.syncData();
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given scalar field
 *
 *          The unary operator += adds a given scalar field to the entire field stored as sfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another sfield to be added to the member field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
sfield& sfield::operator += (sfield &a) {
    F.F += a.F.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given scalar field
 *
 *          The unary operator -= subtracts a given scalar field from the entire field stored as sfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another sfield to be deducted from the member field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
sfield& sfield::operator -= (sfield &a) {
    F.F -= a.F.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to multiply a scalar value to the scalar field
 *
 *          The unary operator *= multiplies a double precision value to the entire field stored as sfield and returns
 *          a pointer to itself.
 *
 * \param   a is a double precision number to be multiplied to the scalar field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
sfield& sfield::operator *= (double a) {
    F.F *= a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar value to the scalar field
 *
 *          The operator = assigns a double precision value to all the scalar field.
 *
 * \param   a is a double precision number to be assigned to the scalar field
 ********************************************************************************************************************************************
 */
void sfield::operator = (double a) {
    F.F = a;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar field to the scalar field
 *
 *          The operator = copies the contents of the input scalar field to itself.
 *
 * \param   a is the scalar field to be assigned to the scalar field
 ********************************************************************************************************************************************
 */
void sfield::operator = (sfield &a) {
    F.F = a.F.F;
}
