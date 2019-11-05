#include "derivative.h"
/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the derivative class
 *
 *          The empty constructor of the derivative class only serves to assign values to the two const parameters of the differ class,
 *          namely <B>grid</B> and <B>F</B>.
 *
 * \param   gridData is a const reference to the global data contained in the grid class
 * \param   F is a reference to field on which the derivative operations are to be performed
 ********************************************************************************************************************************************
 */
derivative::derivative(const grid &gridData, const field &F): gridData(gridData), F(F) { 
    tempMat.resize(gridData.staggrFullSize);
    tempMat.reindexSelf(gridData.collocFullDomain.lbound());

    invDelx = 1.0/gridData.dXi;
    invDely = 1.0/gridData.dEt;
    invDelz = 1.0/gridData.dZt; 

    if (F.xStag) {
        x_Metric.reference(gridData.xi_xStaggr);
        xxMetric.reference(gridData.xixxStaggr);
        x2Metric.reference(gridData.xix2Staggr);
    } else {
        x_Metric.reference(gridData.xi_xColloc);
        xxMetric.reference(gridData.xixxColloc);
        x2Metric.reference(gridData.xix2Colloc);
    }

    if (F.yStag) {
        y_Metric.reference(gridData.et_yStaggr);
        yyMetric.reference(gridData.etyyStaggr);
        y2Metric.reference(gridData.ety2Staggr);
    } else {
        y_Metric.reference(gridData.et_yColloc);
        yyMetric.reference(gridData.etyyColloc);
        y2Metric.reference(gridData.ety2Colloc);
    }

    if (F.zStag) {
        z_Metric.reference(gridData.zt_zStaggr);
        zzMetric.reference(gridData.ztzzStaggr);
        z2Metric.reference(gridData.ztz2Staggr);
    } else {
        z_Metric.reference(gridData.zt_zColloc);
        zzMetric.reference(gridData.ztzzColloc);
        z2Metric.reference(gridData.ztz2Colloc);
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the first derivative of the field with respect to x [second order accurate central difference]
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative1_x(blitz::Array<double, 3> outputMat) {
    outputMat(F.fCore) = central12n(F.F, 0);
    outputMat(F.fCore) *= invDelx;
    outputMat = x_Metric(i)*outputMat(i,j,k);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the first derivative of the field with respect to y [second order accurate central difference]
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative1_y(blitz::Array<double, 3> outputMat) {
    outputMat(F.fCore) = central12n(F.F, 1);
    outputMat(F.fCore) *= invDely;
    outputMat = y_Metric(j)*outputMat(i,j,k);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the first derivative of the field with respect to z [second order accurate central difference]
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative1_z(blitz::Array<double, 3> outputMat) {
    outputMat(F.fCore) = central12n(F.F, 2);
    outputMat(F.fCore) *= invDelz;
    outputMat = z_Metric(k)*outputMat(i,j,k);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the second derivative of the field with respect to x [second order accurate central difference]
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative2_xx( blitz::Array<double, 3> outputMat) {
    tempMat(F.fCore) = central12n(F.F, 0);
    tempMat(F.fCore) *= invDelx;
    outputMat(F.fCore) = central22n(F.F, 0);
    outputMat(F.fCore) *= invDelx*invDelx;

    if (gridData.inputParams.iScheme == 1) {
        outputMat = xxMetric(i)*tempMat(i,j,k) + 0.5*x2Metric(i)*outputMat(i,j,k);
    } else {
        outputMat = xxMetric(i)*tempMat(i,j,k) + x2Metric(i)*outputMat(i,j,k);
    }
    
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the second derivatives of the field with respect to y
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative2_yy( blitz::Array<double, 3> outputMat) {
    tempMat(F.fCore) = central12n(F.F, 1);
    tempMat(F.fCore) *= invDely;
    outputMat(F.fCore) = central22n(F.F, 1);
    outputMat(F.fCore) *= invDely*invDely;
    
    if (gridData.inputParams.iScheme == 1) {
        outputMat = yyMetric(j)*tempMat(i,j,k) + 0.5*y2Metric(j)*outputMat(i,j,k);
    } else {
        outputMat = yyMetric(j)*tempMat(i,j,k) + y2Metric(j)*outputMat(i,j,k);
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the second derivatives of the field with respect to z
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative2_zz( blitz::Array<double, 3> outputMat) {
    tempMat(F.fCore) = central12n(F.F, 2);
    tempMat(F.fCore) *= invDelz;
    outputMat(F.fCore) = central22n(F.F, 2);
    outputMat(F.fCore) *= invDelz*invDelz;
    
    if (gridData.inputParams.iScheme == 1) {
        outputMat = zzMetric(k)*tempMat(i,j,k) + 0.5*z2Metric(k)*outputMat(i,j,k);
    } else {
        outputMat = zzMetric(k)*tempMat(i,j,k) + z2Metric(k)*outputMat(i,j,k);
    }
}


/**
*********************************************************************************************************************************************
*********************************************************************************************************************************************
*    FUNCTIONS FOR CROSS DERIVATIVES CAN BE DEVELOPED IN THIS CLASS <see documentation for more details>
*
*
*********************************************************************************************************************************************
*/

derivative::~derivative() { }
