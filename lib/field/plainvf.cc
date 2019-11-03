#include "plainvf.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the plainvf class
 *
 *          Three blitz arrays to store the data of the three component scalar fields are initialized.
 *          The name for the plain vector field as given by the user is also assigned.
 *
 * \param   gridData is a const reference to the global data contained in the grid class
 * \param   fieldName is a string value set by the user to name and identify the plain vector field
 * \param   refV is a const reference to a sample vfield according to whose components the components of plainvf is resized
 ********************************************************************************************************************************************
 */
plainvf::plainvf(const grid &gridData, const vfield &refV): gridData(gridData) {
    Vx.resize(refV.Vx.fSize);
    Vx.reindexSelf(refV.Vx.flBound);

    Vy.resize(refV.Vy.fSize);
    Vy.reindexSelf(refV.Vy.flBound);

    Vz.resize(refV.Vz.fSize);
    Vz.reindexSelf(refV.Vz.flBound);
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given plain vector field
 *
 *          The unary operator += adds a given plain vector field to the entire field stored as plainvf and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another plainvf to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the plain vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainvf& plainvf::operator += (plainvf &a) {
    Vx += a.Vx;
    Vy += a.Vy;
    Vz += a.Vz;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given plain vector field
 *
 *          The unary operator -= subtracts a given plain vector field from the entire field stored as plainvf and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another plainvf to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the plain vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainvf& plainvf::operator -= (plainvf &a) {
    Vx -= a.Vx;
    Vy -= a.Vy;
    Vz -= a.Vz;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given vector field
 *
 *          The unary operator += adds a given vector field to the entire field stored as plainvf and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another vfield to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the plain vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainvf& plainvf::operator += (vfield &a) {
    Vx += a.Vx.F;
    Vy += a.Vy.F;
    Vz += a.Vz.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given vector field
 *
 *          The unary operator -= subtracts a given vector field from the entire field stored as plainvf and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another vfield to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the plain vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainvf& plainvf::operator -= (vfield &a) {
    Vx -= a.Vx.F;
    Vy -= a.Vy.F;
    Vz -= a.Vz.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to multiply a scalar value to the plain vector field
 *
 *          The unary operator *= multiplies a double precision value to all the fields (Vx, Vy and Vz) stored in plainvf and returns
 *          a pointer to itself.
 *
 * \param   a is a double precision number to be multiplied to the plain vector field
 *
 * \return  A pointer to itself is returned by the plain vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainvf& plainvf::operator *= (double a) {
    Vx *= a;
    Vy *= a;
    Vz *= a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign another plain vector field to the plain vector field
 *
 *          The operator = assigns all the three scalar sub-fields of a plainvf to all the corresponding fields (Vx, Vy and Vz).
 *
 * \param   a is a plainvf to be assigned to the plain vector field
 ********************************************************************************************************************************************
 */
void plainvf::operator = (plainvf &a) {
    Vx = a.Vx;
    Vy = a.Vy;
    Vz = a.Vz;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign another vector field to the plain vector field
 *
 *          The operator = assigns all the three scalar sub-fields of a plainvf to all the corresponding fields (Vx, Vy and Vz).
 *
 * \param   a is a vfield to be assigned to the plain vector field
 ********************************************************************************************************************************************
 */
void plainvf::operator = (vfield &a) {
    Vx = a.Vx.F;
    Vy = a.Vy.F;
    Vz = a.Vz.F;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar value to the plain vector field
 *
 *          The operator = assigns a double precision value to all the fields (Vx, Vy and Vz) stored in plainvf.
 *
 * \param   a is a double precision number to be assigned to the plain vector field
 ********************************************************************************************************************************************
 */
void plainvf::operator = (double a) {
    Vx = a;
    Vy = a;
    Vz = a;
}
