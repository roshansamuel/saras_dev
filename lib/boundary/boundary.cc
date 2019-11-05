#include "boundary.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the boundary class
 *
 *          The class constructor initializes the mesh for computational problem.
 *          Based on the user set values in the parser class, it decides if the simulation is going to be 2D or 3D.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   inField is a reference to scalar field to which the boundary conditions must be applied.
 ********************************************************************************************************************************************
 */
boundary::boundary(const grid &mesh, sfield &inField): mesh(mesh),
                                                       dField(inField) {
    nonHgBC = false;
    setXYZ();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions on the given input scalar field
 *
 *          The function first calls the syncData() function of the sfield to update the sub-domain pads.
 *          Then the boundary conditions are applied at the full domain boundaries.
 *
 ********************************************************************************************************************************************
 */
void boundary::imposeBC() {
    dField.syncData();

    // IMPOSE BC FOR Vz ALONG LEFT AND RIGHT WALLS. FOR PERIODIC CASE, DATA TRANSFER AUTOMATICALLY IMPOSES BC
    if (not mesh.inputParams.xPer) {
        // ADIABATIC BC FOR RBC, SST AND RRBC
        if (mesh.inputParams.probType == 5 || mesh.inputParams.probType == 6 || mesh.inputParams.probType == 8) {
            // T LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND T IS STAGGERED ALONG X
            if (mesh.rankData.xRank == 0) {
                dField.F.F(dField.F.fWalls(0)) = dField.F.F(dField.F.shift(0, dField.F.fWalls(0), 1));
            }
            // T LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND T IS STAGGERED ALONG X
            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                dField.F.F(dField.F.fWalls(1)) = dField.F.F(dField.F.shift(0, dField.F.fWalls(1), -1));
            }

        // CONDUCTING BC FOR VERTICAL CONVECTION
        } else if (mesh.inputParams.probType == 7) {
            // T LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND T IS STAGGERED ALONG X
            if (mesh.rankData.xRank == 0) {
                dField.F.F(dField.F.fWalls(0)) = 1.0;
            }
            // T LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND T IS STAGGERED ALONG X
            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                dField.F.F(dField.F.fWalls(1)) = 0.0;
            }
        }
    }

#ifndef PLANAR
    // IMPOSE BC FOR T ALONG FRONT AND BACK WALLS. FOR PERIODIC CASE, DATA TRANSFER AUTOMATICALLY IMPOSES BC
    if (not mesh.inputParams.yPer) {
        // ADIABATIC BCS
        // T LIES ON THE FRONT WALL AS THE WALL IS ON STAGGERED POINT AND T IS STAGGERED ALONG Y
        if (mesh.rankData.yRank == 0) {
            dField.F.F(dField.F.fWalls(2)) = dField.F.F(dField.F.shift(1, dField.F.fWalls(2), 1));
        }
        // T LIES ON THE BACK WALL AS THE WALL IS ON STAGGERED POINT AND T IS STAGGERED ALONG Y
        if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
            dField.F.F(dField.F.fWalls(3)) = dField.F.F(dField.F.shift(1, dField.F.fWalls(3), -1));
        }
    }
#endif

    // IMPOSE BC FOR T ALONG TOP AND BOTTOM WALLS
    if (mesh.inputParams.zPer) {
        // PERIODIC BCS
        dField.F.F(dField.F.fWalls(4)) = dField.F.F(dField.F.shift(2, dField.F.fWalls(5), -1));
        dField.F.F(dField.F.fWalls(5)) = dField.F.F(dField.F.shift(2, dField.F.fWalls(4), 1));
    } else {
        // NON PERIODIC BCS
        // HOT PLATE AT BOTTOM AND COLD PLATE AT TOP FOR RBC AND RRBC
        if (mesh.inputParams.probType == 5 || mesh.inputParams.probType == 8) {
            if (nonHgBC) {
                // First impose Neumann BC everywhere for adiabatic wall
                dField.F.F(dField.F.fWalls(4)) = dField.F.F(dField.F.shift(2, dField.F.fWalls(4), 1));

                // Now in the area of the circular patch, make all values 0 in order to apply conducting BC
                dField.F.F(dField.F.fWalls(4)) = dField.F.F(dField.F.fWalls(4))*wallMask;

                // Finally apply conducting BC in the circular patch alone
                dField.F.F(dField.F.fWalls(4)) = dField.F.F(dField.F.fWalls(4)) +
                        wallData*(2.0*wallData - dField.F.F(dField.F.shift(2, dField.F.fWalls(4), 1)));

            } else {
                // T LIES ON THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND T IS STAGGERED ALONG Z
                dField.F.F(dField.F.fWalls(4)) = 1.0;
            }

            // T LIES ON THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND T IS STAGGERED ALONG Z
            dField.F.F(dField.F.fWalls(5)) = 0.0;

        // COLD PLATE AT BOTTOM AND HOT PLATE AT TOP FOR SST
        } else if (mesh.inputParams.probType == 6) {
            // T LIES ON THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND T IS STAGGERED ALONG Z
            dField.F.F(dField.F.fWalls(4)) = 0.0;
            // T LIES ON THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND T IS STAGGERED ALONG Z
            dField.F.F(dField.F.fWalls(5)) = 1.0;

        // ADIABATIC BC FOR VERTICAL CONVECTION
        } else if (mesh.inputParams.probType == 7) {
            // T LIES ON EITHER SIDE OF THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND T IS COLLOCATED ALONG Z
            dField.F.F(dField.F.fWalls(4)) = dField.F.F(dField.F.shift(2, dField.F.fWalls(4), 1));
            // T LIES ON EITHER SIDE OF THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND T IS COLLOCATED ALONG Z
            dField.F.F(dField.F.fWalls(5)) = dField.F.F(dField.F.shift(2, dField.F.fWalls(5), -1));
        }
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to create a heating patch for non-homogeneous boundary conditions
 *
 *          The current implementation if for a conducting circular heating plate surrounded by adiabatic insulated walls.
 *          This may replaced in future with an external implementation of non-homogeneous BC.
 *
 ********************************************************************************************************************************************
 */
void boundary::createPatch(int wallNum) {
    // FOLLOWING CRUDE IMPLEMENTATION IS A TEST FOR NON-HOMOGENEOUS BCs
    // IT WORKS FOR ONLY SINGLE CORE FOR NOW
    //
    // wallNum is the number of the wall on which the BC is implemented

    blitz::TinyVector<int, 3> centPt;

    centPt = mesh.globalSize/2;
    centPt(2) = dField.F.fWalls(wallNum).lbound()(2);
    wallData.resize(blitz::TinyVector<int, 3>(dField.F.fWalls(wallNum).ubound() - dField.F.fWalls(wallNum).lbound() + 1));
    wallData.reindexSelf(dField.F.fWalls(wallNum).lbound());
    wallData = 0.0;

    wallMask.resize(blitz::TinyVector<int, 3>(dField.F.fWalls(wallNum).ubound() - dField.F.fWalls(wallNum).lbound() + 1));
    wallMask.reindexSelf(dField.F.fWalls(wallNum).lbound());
    wallMask = true;

    for (int iX = dField.F.fWalls(wallNum).lbound(0); iX <= dField.F.fWalls(wallNum).ubound(0); iX++) {
        for (int iY = dField.F.fWalls(wallNum).lbound(1); iY <= dField.F.fWalls(wallNum).ubound(1); iY++) {
            for (int iZ = dField.F.fWalls(wallNum).lbound(2); iZ <= dField.F.fWalls(wallNum).ubound(2); iZ++) {
                double ptRadius = sqrt(pow((x(iX) - xGlo(centPt(0))), 2) + pow((y(iY) - yGlo(centPt(1))), 2) + pow((z(iZ) - zGlo(centPt(2))), 2));
                if (ptRadius <= mesh.inputParams.patchRadius) {
                    wallData(iX, iY, iZ) = 1.0;
                    wallMask(iX, iY, iZ) = false;
                }
            }
        }
    }

    nonHgBC = true;
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to set the x, y and z arrays by reference to the values in grid class
 *
 *          Depending on whether the variable is collocated or staggered along a given direction,
 *          the function sets the reference to each grid array, x, y and z, to either the
 *          array of locations of collocated points or staggered points.
 ********************************************************************************************************************************************
 */
void boundary::setXYZ() {
    if (dField.F.xStag) {
        x.reference(mesh.xStaggr);
        xGlo.reference(mesh.xStaggrGlobal);
    } else {
        x.reference(mesh.xColloc);
        xGlo.reference(mesh.xCollocGlobal);
    }

    if (dField.F.yStag) {
        y.reference(mesh.yStaggr);
        yGlo.reference(mesh.yStaggrGlobal);
    } else {
        y.reference(mesh.yColloc);
        yGlo.reference(mesh.yCollocGlobal);
    }

    if (dField.F.zStag) {
        z.reference(mesh.zStaggr);
        zGlo.reference(mesh.zStaggrGlobal);
    } else {
        z.reference(mesh.zColloc);
        zGlo.reference(mesh.zCollocGlobal);
    }
}
