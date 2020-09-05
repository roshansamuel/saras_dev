/********************************************************************************************************************************************
 * Saras
 * 
 * Copyright (C) 2019, Mahendra K. Verma
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     1. Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *     2. Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *     3. Neither the name of the copyright holder nor the
 *        names of its contributors may be used to endorse or promote products
 *        derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ********************************************************************************************************************************************
 */
/*! \file derivative.cc
 *
 *  \brief Definitions for functions of class derivative
 *  \sa derivative.h
 *  \author Ali Asad, Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

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
    // TEMPORARY ARRAY TO STORE DERIVATIVES WHEN CALCULATING 2ND ORDER DERIVATIVES
    tempMat.resize(F.fSize);
    tempMat.reindexSelf(F.flBound);

    // INVERSES OF hx, hy AND hz, WHICH ARE MULTIPLIED TO FINITE-DIFFERENCE STENCILS
    invDelx = 1.0/gridData.dXi;
    invDely = 1.0/gridData.dEt;
    invDelz = 1.0/gridData.dZt; 

    // RANGES OF ARRAY INTO WHICH RESULTS FROM BLITZ STENCIL OPERATORS HAVE TO BE WRITTEN
    fullRange = blitz::Range::all();
    xRange = blitz::Range(0, F.fCore.ubound(0), 1);
    yRange = blitz::Range(0, F.fCore.ubound(1), 1);
    zRange = blitz::Range(0, F.fCore.ubound(2), 1);

    setWallRectDomains();

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

    tempMat = 0.0;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the first derivative of the field with respect to x [second order accurate central difference]
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative1_x(blitz::Array<real, 3> outputMat) {
    if (gridData.inputParams.dScheme == 1) {
        outputMat(xRange, fullRange, fullRange) = central12n(F.F, 0);
    } else if (gridData.inputParams.dScheme == 2) {
        outputMat(xRange, fullRange, fullRange) = central14n(F.F, 0);

        // 2ND ORDER CENTRAL DIFFERENCE AT BOUNDARIES
        if (gridData.rankData.xRank == 0) {
            outputMat(x0Mid) = 0.5*(outputMat(x0Rgt) - outputMat(x0Lft));
        }
        if (gridData.rankData.xRank == gridData.rankData.npX - 1) {
            outputMat(x1Mid) = 0.5*(outputMat(x1Rgt) - outputMat(x1Lft));
        }
    }
    outputMat *= invDelx;

    outputMat = x_Metric(i)*outputMat(i, j, k);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the first derivative of the field with respect to y [second order accurate central difference]
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative1_y(blitz::Array<real, 3> outputMat) {
    if (gridData.inputParams.dScheme == 1) {
        outputMat(fullRange, yRange, fullRange) = central12n(F.F, 1);
    } else if (gridData.inputParams.dScheme == 2) {
        outputMat(fullRange, yRange, fullRange) = central14n(F.F, 1);

        // 2ND ORDER CENTRAL DIFFERENCE AT BOUNDARIES
        if (gridData.rankData.yRank == 0) {
            outputMat(y0Mid) = 0.5*(outputMat(y0Rgt) - outputMat(y0Lft));
        }
        if (gridData.rankData.yRank == gridData.rankData.npY - 1) {
            outputMat(y1Mid) = 0.5*(outputMat(y1Rgt) - outputMat(y1Lft));
        }
    }
    outputMat *= invDely;

    outputMat = y_Metric(j)*outputMat(i, j, k);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the first derivative of the field with respect to z [second order accurate central difference]
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative1_z(blitz::Array<real, 3> outputMat) {
    if (gridData.inputParams.dScheme == 1) {
        outputMat(fullRange, fullRange, zRange) = central12n(F.F, 2);
    } else if (gridData.inputParams.dScheme == 2) {
        outputMat(fullRange, fullRange, zRange) = central14n(F.F, 2);

        // 2ND ORDER CENTRAL DIFFERENCE AT BOUNDARIES
        outputMat(z0Mid) = 0.5*(outputMat(z0Rgt) - outputMat(z0Lft));
        outputMat(z1Mid) = 0.5*(outputMat(z1Rgt) - outputMat(z1Lft));
    }
    outputMat *= invDelz;

    outputMat = z_Metric(k)*outputMat(i, j, k);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the second derivative of the field with respect to x [second order accurate central difference]
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative2xx(blitz::Array<real, 3> outputMat) {
    tempMat(xRange, fullRange, fullRange) = central12n(F.F, 0);
    tempMat *= invDelx;

    outputMat(xRange, fullRange, fullRange) = central22n(F.F, 0);
    outputMat *= invDelx*invDelx;
    if (gridData.inputParams.iScheme == 1) {
        outputMat = xxMetric(i)*tempMat(i, j, k) + 0.5*x2Metric(i)*outputMat(i, j, k);
    } else {
        outputMat = xxMetric(i)*tempMat(i, j, k) + x2Metric(i)*outputMat(i, j, k);
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
void derivative::calcDerivative2yy(blitz::Array<real, 3> outputMat) {
    tempMat(fullRange, yRange, fullRange) = central12n(F.F, 1);
    tempMat *= invDely;

    outputMat(fullRange, yRange, fullRange) = central22n(F.F, 1);
    outputMat *= invDely*invDely;
    if (gridData.inputParams.iScheme == 1) {
        outputMat = yyMetric(j)*tempMat(i, j, k) + 0.5*y2Metric(j)*outputMat(i, j, k);
    } else {
        outputMat = yyMetric(j)*tempMat(i, j, k) + y2Metric(j)*outputMat(i, j, k);
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
void derivative::calcDerivative2zz(blitz::Array<real, 3> outputMat) {
    tempMat(fullRange, fullRange, zRange) = central12n(F.F, 2);
    tempMat *= invDelz;

    outputMat(fullRange, fullRange, zRange) = central22n(F.F, 2);
    outputMat *= invDelz*invDelz;
    if (gridData.inputParams.iScheme == 1) {
        outputMat = zzMetric(k)*tempMat(i, j, k) + 0.5*z2Metric(k)*outputMat(i, j, k);
    } else {
        outputMat = zzMetric(k)*tempMat(i, j, k) + z2Metric(k)*outputMat(i, j, k);
    }
}

/**
*********************************************************************************************************************************************
*    FUNCTIONS FOR CROSS DERIVATIVES CAN BE DEVELOPED IN THIS CLASS (SEE DOCUMENTATION FOR MORE DETAILS)
*
*********************************************************************************************************************************************
*/

/**
 ********************************************************************************************************************************************
 * \brief   Function to set the RectDomain objects to be used for computing derivatives at the domain boundaries
 *
 *          When using 4th order stencil, the derivatives at the boundaries need to be recomputed with second order schemes
 *          in order to avoid spurious input from ghost points (this is because the BCs are applied only at the ghost points
 *          next to the boundary and not the ones beyond).
 *          These RectDomain objects can be used to compute these derivatives quickly after using the 4th order stencils.
 *          
 ********************************************************************************************************************************************
 */
void derivative::setWallRectDomains() {
    blitz::TinyVector<int, 3> lb, ub;

    lb = -gridData.padWidths;       lb(0) = 0;
    ub = F.F.ubound();              ub(0) = 0;
    x0Mid = blitz::RectDomain<3>(lb, ub);
    x0Lft = x0Mid;      x0Lft.lbound()(0) -= 1;      x0Lft.ubound()(0) -= 1;
    x0Rgt = x0Mid;      x0Rgt.lbound()(0) += 1;      x0Rgt.ubound()(0) += 1;

    lb = -gridData.padWidths;       lb(0) = F.fCore.ubound(0);
    ub = F.F.ubound();              ub(0) = F.fCore.ubound(0);
    x1Mid = blitz::RectDomain<3>(lb, ub);
    x1Lft = x1Mid;      x1Lft.lbound()(0) -= 1;      x1Lft.ubound()(0) -= 1;
    x1Rgt = x1Mid;      x1Rgt.lbound()(0) += 1;      x1Rgt.ubound()(0) += 1;


    lb = -gridData.padWidths;       lb(1) = 0;
    ub = F.F.ubound();              ub(1) = 0;
    y0Mid = blitz::RectDomain<3>(lb, ub);
    y0Lft = y0Mid;      y0Lft.lbound()(1) -= 1;      y0Lft.ubound()(1) -= 1;
    y0Rgt = y0Mid;      y0Rgt.lbound()(1) += 1;      y0Rgt.ubound()(1) += 1;

    lb = -gridData.padWidths;       lb(1) = F.fCore.ubound(1);
    ub = F.F.ubound();              ub(1) = F.fCore.ubound(1);
    y1Mid = blitz::RectDomain<3>(lb, ub);
    y1Lft = y1Mid;      y1Lft.lbound()(1) -= 1;      y1Lft.ubound()(1) -= 1;
    y1Rgt = y1Mid;      y1Rgt.lbound()(1) += 1;      y1Rgt.ubound()(1) += 1;


    lb = -gridData.padWidths;       lb(2) = 0;
    ub = F.F.ubound();              ub(2) = 0;
    z0Mid = blitz::RectDomain<3>(lb, ub);
    z0Lft = z0Mid;      z0Lft.lbound()(2) -= 1;      z0Lft.ubound()(2) -= 1;
    z0Rgt = z0Mid;      z0Rgt.lbound()(2) += 1;      z0Rgt.ubound()(2) += 1;

    lb = -gridData.padWidths;       lb(2) = F.fCore.ubound(2);
    ub = F.F.ubound();              ub(2) = F.fCore.ubound(2);
    z1Mid = blitz::RectDomain<3>(lb, ub);
    z1Lft = z1Mid;      z1Lft.lbound()(2) -= 1;      z1Lft.ubound()(2) -= 1;
    z1Rgt = z1Mid;      z1Rgt.lbound()(2) += 1;      z1Rgt.ubound()(2) += 1;
}
