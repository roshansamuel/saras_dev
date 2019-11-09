#include "field.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the field class
 *
 *          The field class decides the limits necessary for a 3D array to store the data as per the specified grid staggering details.
 *          It initializes and stores necessary RectDomain objects for getting the core slice and various offset slices for performing
 *          finite difference operations.
 *          Correspondingly, three instances of the \ref differ class are also initialized to compute derivatives along the three directions.
 *          The upper and lower bounds necessary for the array are also calculated depending on the directions along which the mesh is
 *          staggered and those along which it is collocated.
 *          Finally, a blitz array to store the data of the field is resized according to the limits and initialized to 0.
 *
 * \param   gridData is a const reference to the global data contained in the grid class
 * \param   xStag is a const boolean value that is <B>true</B> when the grid is staggered along the x-direction and <B>false</B> when it is not
 * \param   yStag is a const boolean value that is <B>true</B> when the grid is staggered along the y-direction and <B>false</B> when it is not
 * \param   zStag is a const boolean value that is <B>true</B> when the grid is staggered along the z-direction and <B>false</B> when it is not
 ********************************************************************************************************************************************
 */
field::field(const grid &gridData, std::string fieldName, const bool xStag, const bool yStag, const bool zStag):
             gridData(gridData),
             xStag(xStag), yStag(yStag), zStag(zStag),
             xDim(0, gridData.dXi), yDim(1, gridData.dEt), zDim(2, gridData.dZt)
{
    this->fieldName = fieldName;

    fSize = gridData.collocFullSize;
    flBound = gridData.collocFullDomain.lbound();

    if (xStag) {
        fSize(0) = gridData.staggrFullSize(0);
        flBound(0) = gridData.staggrFullDomain.lbound()(0);

        x_Metric.reference(gridData.xi_xStaggr);
        xxMetric.reference(gridData.xixxStaggr);
        x2Metric.reference(gridData.xix2Staggr);
    } else {
        x_Metric.reference(gridData.xi_xColloc);
        xxMetric.reference(gridData.xixxColloc);
        x2Metric.reference(gridData.xix2Colloc);
    }

    if (yStag) {
        fSize(1) = gridData.staggrFullSize(1);
        flBound(1) = gridData.staggrFullDomain.lbound()(1);

        y_Metric.reference(gridData.et_yStaggr);
        yyMetric.reference(gridData.etyyStaggr);
        y2Metric.reference(gridData.ety2Staggr);
    } else {
        y_Metric.reference(gridData.et_yColloc);
        yyMetric.reference(gridData.etyyColloc);
        y2Metric.reference(gridData.ety2Colloc);
    }

    if (zStag) {
        fSize(2) = gridData.staggrFullSize(2);
        flBound(2) = gridData.staggrFullDomain.lbound()(2);

        z_Metric.reference(gridData.zt_zStaggr);
        zzMetric.reference(gridData.ztzzStaggr);
        z2Metric.reference(gridData.ztz2Staggr);
    } else {
        z_Metric.reference(gridData.zt_zColloc);
        zzMetric.reference(gridData.ztzzColloc);
        z2Metric.reference(gridData.ztz2Colloc);
    }

    F.resize(fSize);
    F.reindexSelf(flBound);

    mpiHandle = new mpidata(F, gridData.rankData);

    d1F_dx1.resize(fSize);
    d1F_dx1.reindexSelf(flBound);

    d1F_dy1.resize(fSize);
    d1F_dy1.reindexSelf(flBound);

    d1F_dz1.resize(fSize);
    d1F_dz1.reindexSelf(flBound);

    d2F_dx2.resize(fSize);
    d2F_dx2.reindexSelf(flBound);

    d2F_dy2.resize(fSize);
    d2F_dy2.reindexSelf(flBound);

    d2F_dz2.resize(fSize);
    d2F_dz2.reindexSelf(flBound);

    setCoreSlice();
    setBulkSlice();

    setWallSlices();

    setInterpolationSlices();

    mpiHandle->createSubarrays(fSize, cuBound + 1, gridData.padWidths, xStag, yStag);

    F = 0.0;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the core slice and its offset views
 *
 *          The core and full slices of the field differentiates the domain into the computational sub-domain and
 *          the ghost point/pad point regions which play only an auxilliary role in computing the derivatives.
 *          The core slice defined here also includes the walls of the domain.
 *          The core slice is also offset in all the directions for ease in computation of numerical derivatives.
 *
 ********************************************************************************************************************************************
 */
void field::setCoreSlice() {
    cuBound = gridData.collocCoreDomain.ubound();

    if (xStag) {
        cuBound(0) = gridData.staggrCoreDomain.ubound()(0);
    }

    if (yStag) {
        cuBound(1) = gridData.staggrCoreDomain.ubound()(1);
    }

    if (zStag) {
        cuBound(2) = gridData.staggrCoreDomain.ubound()(2);
    }

    fCore = blitz::RectDomain<3>(blitz::TinyVector<int, 3>(0, 0, 0), cuBound);

    fCLft = xDim.shift(fCore, -1);
    fCRgt = xDim.shift(fCore,  1);

    fCFrt = yDim.shift(fCore, -1);
    fCBak = yDim.shift(fCore,  1);

    fCBot = zDim.shift(fCore, -1);
    fCTop = zDim.shift(fCore,  1);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the bulk slice and its offset views
 *
 *          The bulk and wall slices of the field differentiates the domain into the bulk of the fluid and
 *          the walls of the domain.
 *          The bulk slice is also offset in all the directions for ease in computation of numerical derivatives.
 *
 ********************************************************************************************************************************************
 */
void field::setBulkSlice() {
    blitz::TinyVector<int, 3> blBound;
    blitz::TinyVector<int, 3> buBound;

    blBound = gridData.collocCoreDomain.lbound();
    buBound = gridData.collocCoreDomain.ubound();

    if (xStag) {
        blBound(0) = gridData.staggrCoreDomain.lbound()(0);
        buBound(0) = gridData.staggrCoreDomain.ubound()(0);
    }

    if (yStag) {
        blBound(1) = gridData.staggrCoreDomain.lbound()(1);
        buBound(1) = gridData.staggrCoreDomain.ubound()(1);
    }

    if (zStag) {
        blBound(2) = gridData.staggrCoreDomain.lbound()(2);
        buBound(2) = gridData.staggrCoreDomain.ubound()(2);
    }

    // Bulk and core slices are differentiated only in the boundary sub-domains,
    // and that differentiation is imposed in the following lines
    // At all interior sub-domains after performing MPI domain decomposition,
    // the bulk and core slices are identical

    if (xStag and gridData.subarrayStarts(0) == 0) blBound(0) += 1;

    if (xStag and gridData.subarrayEnds(0) == gridData.globalSize(0) - 1) buBound(0) -= 1;

    if (yStag and gridData.subarrayStarts(1) == 0) blBound(1) += 1;

    if (yStag and gridData.subarrayEnds(1) == gridData.globalSize(1) - 1) buBound(1) -= 1;

    if (zStag) {
        blBound(2) += 1;
        buBound(2) -= 1;
    }

#ifdef PLANAR
    blBound(1) = 0;
    buBound(1) = 0;
#endif

    // THE FOLLOWING LINES ENSURE THAT BOUNDARY PADS GET CORRECTLY UPDATED IN PERIODIC BC
    if (xStag and gridData.inputParams.xPer) {
        if (gridData.rankData.xRank == 0) {
            blBound(0) -= 1;
        }
        if (gridData.rankData.xRank == gridData.rankData.npX - 1) {
            buBound(0) += 1;
        }
    }

    if (yStag and gridData.inputParams.yPer) {
        if (gridData.rankData.yRank == 0) {
            blBound(1) -= 1;
        }

        if (gridData.rankData.yRank == gridData.rankData.npY - 1) {
            buBound(1) += 1;
        }
    }

    fBulk = blitz::RectDomain<3>(blBound, buBound);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the wall slices for the sub-domains
 *
 *          The wall slices of the sub-domain are for imposing the full-domain boundary conditions and hence are of importance
 *          only for the near-boundary sub-domains in parallel computations.
 *          Moreover, only the collocated grid has points on the boundary.
 *          The staggered grid points lie on either side of the domain boundaries.
 *          As a result the wall slices are defined only for those fields for which at least one of \ref xStag, \ref yStag or
 *          \ref zStag is false
 ********************************************************************************************************************************************
 */
void field::setWallSlices() {
    blitz::Array<blitz::TinyVector<int, 3>, 1> wlBound;
    blitz::Array<blitz::TinyVector<int, 3>, 1> wuBound;

    // 6 slices are stored in fWalls corresponding to the 6 faces of the 3D box
    fWalls.resize(6);

    wlBound.resize(6);
    wuBound.resize(6);

    // Wall slices are the locations where the BC (both Neumann and Dirichlet) is imposed.
    // In the places where these slices are being used, they should be on the LHS of equation.
    for (int i=0; i<6; i++) {
        wlBound(i) = F.lbound();
        wuBound(i) = F.ubound();
    }

    // The bulk slice corresponds to the part of the fluid within which all variables are computed at each time step.
    // Correspondingly, the boundary conditions are imposed on the layer just outside the bulk

    // UPPER BOUNDS OF LEFT WALL
    wlBound(0)(0) = wuBound(0)(0) = fBulk.lbound(0) - 1;

    // LOWER BOUNDS OF RIGHT WALL
    wuBound(1)(0) = wlBound(1)(0) = fBulk.ubound(0) + 1;

    // UPPER BOUNDS OF FRONT WALL
    wlBound(2)(1) = wuBound(2)(1) = fBulk.lbound(1) - 1;

    // LOWER BOUNDS OF BACK WALL
    wuBound(3)(1) = wlBound(3)(1) = fBulk.ubound(1) + 1;

    // UPPER BOUNDS OF BOTTOM WALL
    wlBound(4)(2) = wuBound(4)(2) = fBulk.lbound(2) - 1;

    // LOWER BOUNDS OF TOP WALL
    wuBound(5)(2) = wlBound(5)(2) = fBulk.ubound(2) + 1;

    for (int i=0; i<6; i++) {
        fWalls(i) = blitz::RectDomain<3>(wlBound(i), wuBound(i));
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the slices for interpolation while computing convective derivative
 *
 *          This function must be called before using the values contained in the arrays d2F_dx2, d2F_dy2 and d2F_dz2.
 ********************************************************************************************************************************************
 */
void field::setInterpolationSlices() {
    // INTERPOLATION SLICES FOR INTERPOLATING VALUES OF Vx FROM THE vfield
    // In all the below slices, we are considering interpolations between the following 8 variables
    //
    // Vx, Vy, Vz - these are face centered variables sitting on X, Y and Z planes respectively (like velocity)
    // Wx, Wy, Wz - these are edge centered variables sitting on X, Y and Z axes respectively (like vorticity)
    // Pc - this a cell centered variable (like temperature)
    // Qv - this a vertex centered variable (I don't know if anything sits here. But hey! completeness!)
    if (not xStag) {
        if (yStag) {
            if (zStag) {
                /* coll stag stag */
                /**X-Face centered configuration - Vx **/

                /* Interpolation of data - Vx ---> Vx
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => no change
                 **/
                VxIntSlices.resize(1);
                VxIntSlices(0) = fCore;

                /* Interpolation of data - Vy ---> Vx
                 * Interpolation types:
                 *      - X direction => staggered to collocated
                 *      - Y direction => collocated to staggered
                 *      - Z direction => no change
                 **/
                VyIntSlices.resize(4);
                VyIntSlices(0) = fCore;
                VyIntSlices(1) = yDim.shift(VyIntSlices(0), -1);
                VyIntSlices(2) = xDim.shift(VyIntSlices(0), 1);
                VyIntSlices(3) = yDim.shift(VyIntSlices(2), -1);

                /* Interpolation of data - Vz ---> Vx
                 * Interpolation types:
                 *      - X direction => staggered to collocated
                 *      - Y direction => no change
                 *      - Z direction => collocated to staggered
                 **/
                VzIntSlices.resize(4);
                VzIntSlices(0) = fCore;
                VzIntSlices(1) = zDim.shift(VzIntSlices(0), -1);
                VzIntSlices(2) = xDim.shift(VzIntSlices(0), 1);
                VzIntSlices(3) = zDim.shift(VzIntSlices(2), -1);

                /* Interpolation of data - Pc ---> Vx
                 * Interpolation types:
                 *      - X direction => staggered to collocated
                 *      - Y direction => no change
                 *      - Z direction => no change
                 **/
                PcIntSlices.resize(2);
                PcIntSlices(0) = fCore;
                PcIntSlices(1) = xDim.shift(PcIntSlices(0), 1);

            } else {
                /* coll stag coll */
                /**Y-Edge centered configuration - Wy **/

                /* Interpolation of data - Vx ---> Wy
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => staggered to collocated
                 **/
            }
        } else {
            if (zStag) {
                /* coll coll stag */
                /**Z-Edge centered configuration - Wz **/

                /* Interpolation of data - Vx ---> Wz
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => staggered to collocated
                 *      - Z direction => no change
                 **/
            } else {
                /* coll coll coll */
                /**Vertex centered configuration - Qv **/

                /* Interpolation of data - Vx ---> Qv
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => staggered to collocated
                 *      - Z direction => staggered to collocated
                 **/
            }
        }
    } else {
        if (yStag) {
            if (zStag) {
                /* stag stag stag */
                /**Cell centered configuration - Pc **/

                /* Interpolation of data - Vx ---> Pc
                 * Interpolation types:
                 *      - X direction => collocated to staggered
                 *      - Y direction => no change
                 *      - Z direction => no change
                 **/
                VxIntSlices.resize(2);
                VxIntSlices(0) = fCore;
                VxIntSlices(1) = xDim.shift(VxIntSlices(0), -1);

                /* Interpolation of data - Vy ---> Pc
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => collocated to staggered
                 *      - Z direction => no change
                 **/
                VyIntSlices.resize(2);
                VyIntSlices(0) = fCore;
                VyIntSlices(1) = yDim.shift(VyIntSlices(0), -1);

                /* Interpolation of data - Vz ---> Pc
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => collocated to staggered
                 **/
                VzIntSlices.resize(2);
                VzIntSlices(0) = fCore;
                VzIntSlices(1) = zDim.shift(VzIntSlices(0), -1);

                /* Interpolation of data - Pc ---> Pc
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => no change
                 **/
                PcIntSlices.resize(1);
                PcIntSlices(0) = fCore;

            } else {
                /* stag stag coll */
                /**Z-Face centered configuration - Vz **/

                /* Interpolation of data - Vx ---> Vz
                 * Interpolation types:
                 *      - X direction => collocated to staggered
                 *      - Y direction => no change
                 *      - Z direction => staggered to collocated
                 **/
                VxIntSlices.resize(4);
                VxIntSlices(0) = fCore;
                VxIntSlices(1) = xDim.shift(VxIntSlices(0), -1);
                VxIntSlices(2) = zDim.shift(VxIntSlices(0), 1);
                VxIntSlices(3) = xDim.shift(VxIntSlices(2), -1);

                /* Interpolation of data - Vy ---> Vz
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => collocated to staggered
                 *      - Z direction => staggered to collocated
                 **/
                VyIntSlices.resize(4);
                VyIntSlices(0) = fCore;
                VyIntSlices(1) = yDim.shift(VyIntSlices(0), -1);
                VyIntSlices(2) = zDim.shift(VyIntSlices(0), 1);
                VyIntSlices(3) = yDim.shift(VyIntSlices(2), -1);

                /* Interpolation of data - Vz ---> Vz
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => no change
                 **/
                VzIntSlices.resize(1);
                VzIntSlices(0) = fCore;

                /* Interpolation of data - Pc ---> Vz
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => staggered to collocated
                 **/
                PcIntSlices.resize(2);
                PcIntSlices(0) = fCore;
                PcIntSlices(1) = zDim.shift(PcIntSlices(0), 1);

            }
        } else {
            if (zStag) {
                /* stag coll stag */
                /**Y-Face centered configuration - Vy **/

                /* Interpolation of data - Vx ---> Vy
                 * Interpolation types:
                 *      - X direction => collocated to staggered
                 *      - Y direction => staggered to collocated
                 *      - Z direction => no change
                 **/
                VxIntSlices.resize(4);
                VxIntSlices(0) = fCore;
                VxIntSlices(1) = xDim.shift(VxIntSlices(0), -1);
                VxIntSlices(2) = yDim.shift(VxIntSlices(0), 1);
                VxIntSlices(3) = xDim.shift(VxIntSlices(2), -1);

                /* Interpolation of data - Vy ---> Vy
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => no change
                 **/
                VyIntSlices.resize(1);
                VyIntSlices(0) = fCore;

                /* Interpolation of data - Vz ---> Vy
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => staggered to collocated
                 *      - Z direction => collocated to staggered
                 **/
                VzIntSlices.resize(4);
                VzIntSlices(0) = fCore;
                VzIntSlices(1) = zDim.shift(VzIntSlices(0), -1);
                VzIntSlices(2) = yDim.shift(VzIntSlices(0), 1);
                VzIntSlices(3) = zDim.shift(VzIntSlices(2), -1);

                /* Interpolation of data - Pc ---> Vy
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => staggered to collocated
                 *      - Z direction => no change
                 **/
                PcIntSlices.resize(2);
                PcIntSlices(0) = fCore;
                PcIntSlices(1) = yDim.shift(PcIntSlices(0), 1);

            } else {
                /* stag coll coll */
                /**X-Edge centered configuration - Wx **/

                /* Interpolation of data - Vx ---> Wx
                 * Interpolation types:
                 *      - X direction => collocated to staggered
                 *      - Y direction => staggered to collocated
                 *      - Z direction => staggered to collocated
                 **/
            }
        }
    }


    // RESET INTERPOLATION SLICES FOR PLANAR GRID
#ifdef PLANAR
    if (fieldName == "Vy") {
        VxIntSlices.resize(1);
        VyIntSlices.resize(1);
        VzIntSlices.resize(1);

        VxIntSlices(0) = fCore;
        VyIntSlices(0) = fCore;
        VzIntSlices(0) = fCore;
    } else if (fieldName == "Vx") {
        VyIntSlices.resize(1);

        VyIntSlices(0) = fCore;
    } else if (fieldName == "Vz") {
        VyIntSlices.resize(1);

        VyIntSlices(0) = fCore;
    }
#endif
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the first derivative along X direction
 *
 *          The first derivative of the field values are calculated using functions from the \ref differ class.
 *          The corresponding \ref differ object has already been initialized in the constructor.
 *
 ********************************************************************************************************************************************
 */
inline void field::x1Deriv() {
    xDim.D1D(F, d1F_dx1, fCore);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the first derivative along Y direction
 *
 *          The first derivative of the field values are calculated using functions from the \ref differ class.
 *          The corresponding \ref differ object has already been initialized in the constructor.
 *
 ********************************************************************************************************************************************
 */
inline void field::y1Deriv() {
    yDim.D1D(F, d1F_dy1, fCore);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the first derivative along Z direction
 *
 *          The first derivative of the field values are calculated using functions from the \ref differ class.
 *          The corresponding \ref differ object has already been initialized in the constructor.
 *
 ********************************************************************************************************************************************
 */
inline void field::z1Deriv() {
    zDim.D1D(F, d1F_dz1, fCore);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the second derivative along X direction
 *
 *          The first derivative of the field values are calculated using functions from the \ref differ class.
 *          The corresponding \ref differ object has already been initialized in the constructor.
 *
 ********************************************************************************************************************************************
 */
inline void field::x2Deriv() {
    xDim.D2D(F, d2F_dx2, fCore);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the second derivative along Y direction
 *
 *          The first derivative of the field values are calculated using functions from the \ref differ class.
 *          The corresponding \ref differ object has already been initialized in the constructor.
 *
 ********************************************************************************************************************************************
 */
inline void field::y2Deriv() {
    yDim.D2D(F, d2F_dy2, fCore);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the second derivative along Z direction
 *
 *          The first derivative of the field values are calculated using functions from the \ref differ class.
 *          The corresponding \ref differ object has already been initialized in the constructor.
 *
 ********************************************************************************************************************************************
 */
inline void field::z2Deriv() {
    zDim.D2D(F, d2F_dz2, fCore);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the first derivatives of the scalar field along all the three directions
 *
 *          This function must be called before using the values contained in the arrays d1F_dx1, d1F_dy1 and d1F_dz1.
 *          The values of these arrays are updated by this function in each call.
 ********************************************************************************************************************************************
 */
void field::calcDerivatives1() {
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::thirdIndex k;

    x1Deriv();
    d1F_dx1 = x_Metric(i)*d1F_dx1(i,j,k);

#ifndef PLANAR
    y1Deriv();
    d1F_dy1 = y_Metric(j)*d1F_dy1(i,j,k);
#endif

    z1Deriv();
    d1F_dz1 = z_Metric(k)*d1F_dz1(i,j,k);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the second derivatives of the scalar field
 *
 *          This function must be called before using the values contained in the arrays d2F_dx2, d2F_dy2 and d2F_dz2.
 *          The values of these arrays are updated by this function in each call.
 *          <B>This function changes the values contained in the arrays d1F_dx1, d1F_dy1 and d1F_dz1</B>.
 *          As a result, these arrays must be appropriately updated with call to calcDerivatives1 before using them again.
 ********************************************************************************************************************************************
 */
void field::calcDerivatives2() {
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::thirdIndex k;

    x1Deriv();
    x2Deriv();
    if (gridData.inputParams.iScheme == 1) {
        d2F_dx2 = xxMetric(i)*d1F_dx1(i,j,k) + 0.5*x2Metric(i)*d2F_dx2(i,j,k);
    } else {
        d2F_dx2 = xxMetric(i)*d1F_dx1(i,j,k) + x2Metric(i)*d2F_dx2(i,j,k);
    }

#ifndef PLANAR
    y1Deriv();
    y2Deriv();
    if (gridData.inputParams.iScheme == 1) {
        d2F_dy2 = yyMetric(j)*d1F_dy1(i,j,k) + 0.5*y2Metric(j)*d2F_dy2(i,j,k);
    } else {
        d2F_dy2 = yyMetric(j)*d1F_dy1(i,j,k) + y2Metric(j)*d2F_dy2(i,j,k);
    }
#endif

    z1Deriv();
    z2Deriv();
    if (gridData.inputParams.iScheme == 1) {
        d2F_dz2 = zzMetric(k)*d1F_dz1(i,j,k) + 0.5*z2Metric(k)*d2F_dz2(i,j,k);
    } else {
        d2F_dz2 = zzMetric(k)*d1F_dz1(i,j,k) + z2Metric(k)*d2F_dz2(i,j,k);
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronise data across all processors when performing parallel computations
 *
 *          This function calls the \ref mpidata#syncData "syncData" function of mpidata class to perform perform data-transfer and thus update
 *          the sub-domain boundary pads.
 ********************************************************************************************************************************************
 */
void field::syncData() {
    mpiHandle->syncData();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to extract the maximum value from the field
 *
 *          The function uses the in-built blitz function to obtain the maximum value in an array.
 *          While performing parallel computation, the function performs an <B>MPI_Allreduce()</B> to get
 *          the global maximum from the entire computational domain.
 *
 * \return  The double precision value of the maximum is returned (it is implicitly assumed that only double precision values are used)
 ********************************************************************************************************************************************
 */
double field::fieldMax() {
    double localMax, globalMax;

    localMax = blitz::max(F);

    /***************************************************************************************************************
     * DID YOU KNOW?                                                                                               *
     * In the line above, most compilers will not complain even if you omitted the namespace specification blitz:: *
     * This behaviour wasted an hour of my development time (including the effort of making this nice box).        *
     * Check Ref. [4] in README for explanation.                                                                   *
     ***************************************************************************************************************/

    MPI_Allreduce(&localMax, &globalMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD);

    return globalMax;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given field
 *
 *          The unary operator += adds a given field to the field stored by the class and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another to be added to the member field
 *
 * \return  A pointer to itself is returned by the field class to which the operator belongs
 ********************************************************************************************************************************************
 */
field& field::operator += (field &a) {
    F += a.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given field
 *
 *          The unary operator -= subtracts a given field from the field stored by the class and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another field to be deducted from the member field
 *
 * \return  A pointer to itself is returned by the field class to which the operator belongs
 ********************************************************************************************************************************************
 */
field& field::operator -= (field &a) {
    F -= a.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given scalar value
 *
 *          The unary operator += adds a given constant scalar value to the field stored by the class and returns
 *          a pointer to itself.
 *
 * \param   a is a double precision number to be added to the field
 *
 * \return  A pointer to itself is returned by the field class to which the operator belongs
 ********************************************************************************************************************************************
 */
field& field::operator += (double a) {
    F += a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given scalar value
 *
 *          The unary operator -= subtracts a given constant scalar value from the field stored by the class and returns
 *          a pointer to itself.
 *
 * \param   a is a double precision number to be subtracted from the field
 *
 * \return  A pointer to itself is returned by the field class to which the operator belongs
 ********************************************************************************************************************************************
 */
field& field::operator -= (double a) {
    F -= a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar value to the field
 *
 *          The operator = assigns a double precision value to the entire field.
 *
 * \param   a is a double precision number to be assigned to the field
 ********************************************************************************************************************************************
 */
void field::operator = (double a) {
    F = a;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a field to the field
 *
 *          The operator = copies the contents of the input field to itself.
 *
 * \param   a is the field to be assigned to the field
 ********************************************************************************************************************************************
 */
void field::operator = (field &a) {
    F = a.F;
}

field::~field() { }
