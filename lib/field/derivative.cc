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
 *  \author Roshan Samuel, Ali Asad
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
 *          The constructor assigns values to the two const parameters of the derivative class,
 *          namely <B>grid</B> and <B>F</B>.
 *          It resizes tmpArray, the blitz array used to hold temporary data while calculating
 *          derivatives, computes the factors to be used with blitz stencils, and assigns
 *          the appropriate references to grid derivatives for performing finite-difference
 *          operations on non-uniform grids.
 *
 * \param   gridData is a const reference to the global data in the grid class
 * \param   F is a reference to the field on which finite-difference operations will be performed
 ********************************************************************************************************************************************
 */

derivative::derivative(const grid &gridData, const field &F): gridData(gridData), F(F) {
    // TEMPORARY ARRAY TO STORE DERIVATIVES WHEN CALCULATING 2ND ORDER DERIVATIVES
    tmpArray.resize(F.fSize);
    tmpArray.reindexSelf(F.flBound);

    // INVERSES OF hx, hy AND hz, WHICH ARE MULTIPLIED TO FINITE-DIFFERENCE STENCILS
    ihx = 1.0/gridData.dXi;         ihx2 = pow(ihx, 2.0);
    ihy = 1.0/gridData.dEt;         ihy2 = pow(ihy, 2.0);
    ihz = 1.0/gridData.dZt;         ihz2 = pow(ihz, 2.0);

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

    tmpArray = 0.0;

    xfr = (gridData.rankData.xRank == 0)? true: false;
    yfr = (gridData.rankData.yRank == 0)? true: false;
    xlr = (gridData.rankData.xRank == gridData.rankData.npX - 1)? true: false;
    ylr = (gridData.rankData.yRank == gridData.rankData.npY - 1)? true: false;
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the first derivative of the field with respect to x
 *
 *          This function must be called using an output array whose shape and size
 *          should be same as that of the field.
 *          It uses blitz stencils to calculate derivatives using central differencing.
 *          
 * \param   outArray is the blitz array into which result will be written.
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative1_x(blitz::Array<real, 3> outArray) {
    if (gridData.inputParams.dScheme == 1) {
        outArray(xRange, fullRange, fullRange) = central12n(F.F, 0);

    } else if (gridData.inputParams.dScheme == 2) {
        outArray(xRange, fullRange, fullRange) = central14n(F.F, 0);

        // 2ND ORDER CENTRAL DIFFERENCE AT BOUNDARIES
        if (xfr) outArray(x0Mid) = 0.5*(F.F(x0Rgt) - F.F(x0Lft));
        if (xlr) outArray(x1Mid) = 0.5*(F.F(x1Rgt) - F.F(x1Lft));
    }

    outArray *= ihx;

    outArray = x_Metric(i)*outArray(i, j, k);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the first derivative of the field with respect to y
 *
 *          This function must be called using an output array whose shape and size
 *          should be same as that of the field.
 *          It uses blitz stencils to calculate derivatives using central differencing.
 *          
 * \param   outArray is the blitz array into which result will be written.
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative1_y(blitz::Array<real, 3> outArray) {
    if (gridData.inputParams.dScheme == 1) {
        outArray(fullRange, yRange, fullRange) = central12n(F.F, 1);

    } else if (gridData.inputParams.dScheme == 2) {
        outArray(fullRange, yRange, fullRange) = central14n(F.F, 1);

        // 2ND ORDER CENTRAL DIFFERENCE AT BOUNDARIES
        if (yfr) outArray(y0Mid) = 0.5*(F.F(y0Rgt) - F.F(y0Lft));
        if (ylr) outArray(y1Mid) = 0.5*(F.F(y1Rgt) - F.F(y1Lft));
    }
    outArray *= ihy;

    outArray = y_Metric(j)*outArray(i, j, k);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the first derivative of the field with respect to z
 *
 *          This function must be called using an output array whose shape and size
 *          should be same as that of the field.
 *          It uses blitz stencils to calculate derivatives using central differencing.
 *          
 * \param   outArray is the blitz array into which result will be written.
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative1_z(blitz::Array<real, 3> outArray) {
    if (gridData.inputParams.dScheme == 1) {
        outArray(fullRange, fullRange, zRange) = central12n(F.F, 2);

    } else if (gridData.inputParams.dScheme == 2) {
        outArray(fullRange, fullRange, zRange) = central14n(F.F, 2);

        // 2ND ORDER CENTRAL DIFFERENCE AT BOUNDARIES
        outArray(z0Mid) = 0.5*(F.F(z0Rgt) - F.F(z0Lft));
        outArray(z1Mid) = 0.5*(F.F(z1Rgt) - F.F(z1Lft));
    }
    outArray *= ihz;

    outArray = z_Metric(k)*outArray(i, j, k);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the second derivative of the field with respect to x
 *
 *          This function must be called using an output array whose shape and size
 *          should be same as that of the field.
 *          It uses blitz stencils to calculate derivatives using central differencing.
 *          
 * \param   outArray is the blitz array into which result will be written.
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative2xx(blitz::Array<real, 3> outArray) {
    tmpArray(xRange, fullRange, fullRange) = central12n(F.F, 0);
    tmpArray *= ihx;

    outArray(xRange, fullRange, fullRange) = central22n(F.F, 0);
    outArray *= ihx2;

    outArray = xxMetric(i)*tmpArray(i, j, k) + x2Metric(i)*outArray(i, j, k);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the second derivatives of the field with respect to y
 *
 *          This function must be called using an output array whose shape and size
 *          should be same as that of the field.
 *          It uses blitz stencils to calculate derivatives using central differencing.
 *          
 * \param   outArray is the blitz array into which result will be written.
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative2yy(blitz::Array<real, 3> outArray) {
    tmpArray(fullRange, yRange, fullRange) = central12n(F.F, 1);
    tmpArray *= ihy;

    outArray(fullRange, yRange, fullRange) = central22n(F.F, 1);
    outArray *= ihy2;

    outArray = yyMetric(j)*tmpArray(i, j, k) + y2Metric(j)*outArray(i, j, k);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the second derivatives of the field with respect to z
 *
 *          This function must be called using an output array whose shape and size
 *          should be same as that of the field.
 *          It uses blitz stencils to calculate derivatives using central differencing.
 *          
 * \param   outArray is the blitz array into which result will be written.
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative2zz(blitz::Array<real, 3> outArray) {
    tmpArray(fullRange, fullRange, zRange) = central12n(F.F, 2);
    tmpArray *= ihz;

    outArray(fullRange, fullRange, zRange) = central22n(F.F, 2);
    outArray *= ihz2;

    outArray = zzMetric(k)*tmpArray(i, j, k) + z2Metric(k)*outArray(i, j, k);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to set the RectDomain objects used for computing derivatives at domain boundaries
 *
 *          When using 4th order stencil, derivatives at the boundaries need to be recomputed with second order schemes.
 *          This is to avoid spurious input from ghost points (because the BCs are applied only at the ghost points
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
