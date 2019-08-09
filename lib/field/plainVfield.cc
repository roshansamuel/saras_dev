#include "sfield.h"
#include "plainVfield.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the pvfield class
 *
 *          Three instances of the sfield class to store the data of the three component scalar fields are initialized.
 *          The scalar fields are initialized with appropriate grid staggering to place the components on the cell faces.
 *          The name for the vector field as given by the user is also assigned.
 *
 * \param   gridData is a const reference to the global data contained in the grid class
 * \param   fieldName is a string value set by the user to name and identify the vector field
 ********************************************************************************************************************************************
 */
pvfield::pvfield(const grid &gridData, std::string fieldName):
               gridData(gridData),
               Vx(gridData, "Vx", false, true, true),
               Vy(gridData, "Vy", true, false, true),
               Vz(gridData, "Vz", true, true, false)
{
    this->fieldName = fieldName;
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
void pvfield::syncData() {
    Vx.syncData();
    Vy.syncData();
    Vz.syncData();
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given vector field
 *
 *          The unary operator += adds a given vector field to the entire field stored as pvfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another pvfield to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
pvfield& pvfield::operator += (pvfield &a) {
    Vx.F += a.Vx.F;
    Vy.F += a.Vy.F;
    Vz.F += a.Vz.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given vector field
 *
 *          The unary operator -= subtracts a given vector field from the entire field stored as pvfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another pvfield to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
pvfield& pvfield::operator -= (pvfield &a) {
    Vx.F -= a.Vx.F;
    Vy.F -= a.Vy.F;
    Vz.F -= a.Vz.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to multiply a scalar value to the vector field
 *
 *          The unary operator *= multiplies a double precision value to all the fields (Vx, Vy and Vz) stored in pvfield and returns
 *          a pointer to itself.
 *
 * \param   a is a double precision number to be multiplied to the vector field
 *
 * \return  A pointer to itself is returned by the vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
pvfield& pvfield::operator *= (double a) {
    Vx.F *= a;
    Vy.F *= a;
    Vz.F *= a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar value to the vector field
 *
 *          The operator = assigns a double precision value to all the fields (Vx, Vy and Vz) stored in pvfield.
 *
 * \param   a is a double precision number to be assigned to the vector field
 ********************************************************************************************************************************************
 */
void pvfield::operator = (double a) {
    Vx.F = a;
    Vy.F = a;
    Vz.F = a;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign another vector field to the vector field
 *
 *          The operator = assigns all the three scalar sub-fields of a pvfield to all the corresponding fields (Vx, Vy and Vz).
 *
 * \param   a is a pvfield to be assigned to the vector field
 ********************************************************************************************************************************************
 */
void pvfield::operator = (pvfield &a) {
    Vx.F = a.Vx.F;
    Vy.F = a.Vy.F;
    Vz.F = a.Vz.F;
}
