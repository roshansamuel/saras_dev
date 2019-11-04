#include "plainsf.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the plainsf class
 *
 *          The instance of the field class to store the data of the scalar field is initialized, and the necessary grid
 *          transformation derivatives along each direction are chosen according to the grid staggering.
 *          The arrays to store the output from various operators like derivatives, convective derivatives, etc. are also
 *          allocated.
 *          Finally, an instance of the <B>mpidata</B> class is initialized to store the sub-arrays to be send/received
 *          across the processors during MPI communication.
 *
 * \param   gridData is a const reference to the global data contained in the grid class
 * \param   fieldName is a string value set by the user to name and identify the plain scalar field
 * \param   refF is a const reference to a sample sfield according to which the plainsf is resized
 ********************************************************************************************************************************************
 */
plainsf::plainsf(const grid &gridData, const sfield &refF): gridData(gridData) {
    F.resize(refF.F.fSize);
    F.reindexSelf(refF.F.flBound);
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given plain scalar field
 *
 *          The unary operator += adds a given plain scalar field to the entire field stored as plainsf and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another plainsf to be added to the member field
 *
 * \return  A pointer to itself is returned by the plain scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainsf& plainsf::operator += (plainsf &a) {
    F += a.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given plain scalar field
 *
 *          The unary operator -= subtracts a given plain scalar field from the entire field stored as plainsf and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another plainsf to be deducted from the member field
 *
 * \return  A pointer to itself is returned by the plain scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainsf& plainsf::operator -= (plainsf &a) {
    F -= a.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given scalar field
 *
 *          The unary operator += adds a given scalar field to the entire field stored as plainsf and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another sfield to be added to the member field
 *
 * \return  A pointer to itself is returned by the plain scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainsf& plainsf::operator += (sfield &a) {
    F += a.F.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given scalar field
 *
 *          The unary operator -= subtracts a given scalar field from the entire field stored as plainsf and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another sfield to be deducted from the member field
 *
 * \return  A pointer to itself is returned by the plain scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainsf& plainsf::operator -= (sfield &a) {
    F -= a.F.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to multiply a scalar value to the scalar field
 *
 *          The unary operator *= multiplies a double precision value to the entire field stored as plainsf and returns
 *          a pointer to itself.
 *
 * \param   a is a double precision number to be multiplied to the scalar field
 *
 * \return  A pointer to itself is returned by the plain scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainsf& plainsf::operator *= (double a) {
    F *= a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign another plain scalar field to the plain scalar field
 *
 *          The operator = copies the contents of the input plain scalar field to itself.
 *
 * \param   a is a plainsf to be assigned to the plain scalar field
 ********************************************************************************************************************************************
 */
void plainsf::operator = (plainsf &a) {
    F = a.F;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign another scalar field to the plain scalar field
 *
 *          The operator = copies the contents of the input scalar field to itself.
 *
 * \param   a is a sfield to be assigned to the scalar field
 ********************************************************************************************************************************************
 */
void plainsf::operator = (sfield &a) {
    F = a.F.F;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar value to the plain scalar field
 *
 *          The operator = assigns a double precision value to all the scalar field.
 *
 * \param   a is a double precision number to be assigned to the plain scalar field
 ********************************************************************************************************************************************
 */
void plainsf::operator = (double a) {
    F = a;
}
