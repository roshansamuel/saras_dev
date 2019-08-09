#include "psfield.h"
#include "pvfield.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the psfield class
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
psfield::psfield(const grid &gridData, std::string fieldName):
               gridData(gridData),
               F(gridData, fieldName, true, true, true)
{
    this->fieldName = fieldName;    
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronise data across all processors when performing parallel computations
 *
 *          This function calls the \ref mpidata#syncData "syncData" function of mpidata class to perform perform data-transfer and thus update
 *          the sub-domain boundary pads.
 ********************************************************************************************************************************************
 */
void psfield::syncData() {
    F.syncData();
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given scalar field
 *
 *          The unary operator += adds a given scalar field to the entire field stored as psfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another psfield to be added to the member field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
psfield& psfield::operator += (psfield &a) {
    F.F += a.F.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given scalar field
 *
 *          The unary operator -= subtracts a given scalar field from the entire field stored as psfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another psfield to be deducted from the member field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
psfield& psfield::operator -= (psfield &a) {
    F.F -= a.F.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to multiply a scalar value to the scalar field
 *
 *          The unary operator *= multiplies a double precision value to the entire field stored as psfield and returns
 *          a pointer to itself.
 *
 * \param   a is a double precision number to be multiplied to the scalar field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
psfield& psfield::operator *= (double a) {
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
void psfield::operator = (double a) {
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
void psfield::operator = (psfield &a) {
    F.F = a.F.F;
}
