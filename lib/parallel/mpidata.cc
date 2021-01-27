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
/*! \file mpidata.cc
 *
 *  \brief Definitions for functions of class mpidata
 *  \sa mpidata.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "mpidata.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the mpidata class
 *
 *          The short constructor of mpidata class merely resizes the array of MPI_Status and MPI_Request datatypes.
 *          The former is used in non-blocking communication of MPI_Irecv, while the later is used in the MPI_Waitall
 *          function to complete the non-blocking communication call.
 *
 * \param   inputArray is the blitz array whose sub-arrays have to be created and synchronised across processors
 * \param   parallelData is a const reference to the global data contained in the parallel class
 ********************************************************************************************************************************************
 */
mpidata::mpidata(blitz::Array<real, 3> inputArray, const parallel &parallelData): dataField(inputArray), rankData(parallelData) {
    recvStatus.resize(4);
    recvRequest.resize(4);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the subarray MPI_Datatypes
 *
 *          Must be called only after the grid class has been initialized.
 *          The subarray data-types cannot be created within the constructor of the parallel class as it needs the grid parameters for
 *          setting the limits of the subarrays.
 *          For this, the grid class will have to be included in the parallel class.
 *
 *          However, the grid object cannot be passed to the parallel class as the grid class already includes the parallel object
 *          within itself, and a reverse include will raise cyclic dependency error.
 *          As a result, the mpidata class offers an additional layer over the parallel class for grid specific data transfer functions.
 *
 * \param   globSize stores the global size of a sub-domain - including core and pads
 * \param   coreSize stores the size of the core of the sub-domain and is similar to the collocCoreSize variable in the grid class
 * \param   padWidth contains the widths of pads along the 3 directions, namely padWidths TinyVector from the grid class
 * \param   xStag specifies whether the array to which the instance of \ref mpidata class is associated with has its data points staggered in x-direction or not
 * \param   yStag specifies whether the array to which the instance of \ref mpidata class is associated with has its data points staggered in y-direction or not
 ********************************************************************************************************************************************
 */
void mpidata::createSubarrays(const blitz::TinyVector<int, 3> globSize,
                              const blitz::TinyVector<int, 3> coreSize,
                              const blitz::TinyVector<int, 3> padWidth,
                              const bool xStag, const bool yStag,
                              const bool xPer, const bool yPer) {
    int count;

    /** The <B>loclSize</B> variable holds the local size of the sub-array slice to be sent/received within the sub-domain. */
    blitz::TinyVector<int, 3> loclSize;

    /** The <B>saStarts</B> variable holds the starting coordinates of the sub-array slice. */
    blitz::TinyVector<int, 3> saStarts;

    /** The <B>globCopy</B> variable holds a copy of the global size of the sub-domain. This keeps the original array safe*/
    blitz::TinyVector<int, 3> globCopy;

    cSize = coreSize;
    pSize = padWidth;
    globCopy = globSize;

    xsFlag = xStag;
    ysFlag = yStag;

    xsPer = xPer;
    ysPer = yPer;

    // CREATING SUBARRAYS FOR TRANSFER ACROSS THE 4 FACES OF EACH SUB-DOMAIN

    /************************************************************** NOTE ************************************************************\
     * MPI Subarrays assume that the starting index of arrays is 0, 0, 0                                                            *
     * But the arrays used here through Blitz start with the index (-padX, -padY, -padZ)                                            *
     * Hence the saStarts variable is shifted accordingly                                                                           *
    \********************************************************************************************************************************/

    //*****************************************************! ALONG XI-DIRECTION !***************************************************//
    // SEND SUB-ARRAY ON LEFT SIDE
    saStarts = padWidth;
    loclSize = coreSize;            loclSize(0) = padWidth(0);

    // STAGGERED GRID SHARE A POINT ACROSS SUB-DOMAIN BOUNDARIES.
    // THE DATA ON THIS POINT IS AVERAGED FROM THE VALUES IN BOTH SUB-DOMAINS.
    // HENCE, THE SHARED POINT AND ONE INTERIOR POINT ARE SENT.
    if (xStag) {
        loclSize(0) += 1;
    }

    MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &sendSubarrayX0);
    MPI_Type_commit(&sendSubarrayX0);

    //if (rankData.rank == 0) std::cout << loclSize << std::endl;
    //MPI_Finalize();
    //exit(0);

    // RECEIVE SUB-ARRAY ON LEFT SIDE
    // FOR FACE-CENTERED (STAGGERED) DATA, THE RECEIVED DATA IS WRITTEN
    // INTO A TEMPORARY ARRAY, AND HENCE USES A DIFFERENT MPI DATATYPE
    if (xStag) {
        count = loclSize(0)*loclSize(1)*loclSize(2);

        MPI_Type_contiguous(count, MPI_FP_REAL, &recvSubarrayX0);

        recvDataX0.resize(loclSize);

        padX0 = blitz::RectDomain<3>(blitz::TinyVector<int, 3>(-pSize(0), 0, 0),
                                     blitz::TinyVector<int, 3>(-1, cSize(1) - 1, cSize(2) - 1));
        wallX0 = blitz::RectDomain<3>(blitz::TinyVector<int, 3>(0, 0, 0),
                                      blitz::TinyVector<int, 3>(0, cSize(1) - 1, cSize(2) - 1));
    } else {
        saStarts = padWidth;            saStarts(0) = 0;
        loclSize = coreSize;            loclSize(0) = padWidth(0);

        MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &recvSubarrayX0);
    }
    MPI_Type_commit(&recvSubarrayX0);


    // SEND SUB-ARRAY ON RIGHT SIDE
    saStarts = padWidth;            saStarts(0) = coreSize(0);
    loclSize = coreSize;            loclSize(0) = padWidth(0);

    // STAGGERED GRID SHARE A POINT ACROSS SUB-DOMAIN BOUNDARIES.
    // THE DATA ON THIS POINT IS AVERAGED FROM THE VALUES IN BOTH SUB-DOMAINS.
    // HENCE, THE SHARED POINT AND ONE INTERIOR POINT ARE SENT.
    if (xStag) {
        saStarts(0) -= 1;
        loclSize(0) += 1;
    }

    MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &sendSubarrayX1);
    MPI_Type_commit(&sendSubarrayX1);

    // RECEIVE SUB-ARRAY ON RIGHT SIDE
    // FOR FACE-CENTERED (STAGGERED) DATA, THE RECEIVED DATA IS WRITTEN
    // INTO A TEMPORARY ARRAY, AND HENCE USES A DIFFERENT MPI DATATYPE
    if (xStag) {
        count = loclSize(0)*loclSize(1)*loclSize(2);

        MPI_Type_contiguous(count, MPI_FP_REAL, &recvSubarrayX1);

        recvDataX1.resize(loclSize);

        padX1 = blitz::RectDomain<3>(blitz::TinyVector<int, 3>(cSize(0), 0, 0),
                                     blitz::TinyVector<int, 3>(cSize(0) + pSize(0) - 1, cSize(1) - 1, cSize(2) - 1));

        wallX1 = blitz::RectDomain<3>(blitz::TinyVector<int, 3>(cSize(0) - 1, 0, 0),
                                      blitz::TinyVector<int, 3>(cSize(0) - 1, cSize(1) - 1, cSize(2) - 1));
    } else {
        saStarts = padWidth;            saStarts(0) = coreSize(0) + padWidth(0);
        loclSize = coreSize;            loclSize(0) = padWidth(0);

        MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &recvSubarrayX1);
    }
    MPI_Type_commit(&recvSubarrayX1);


    //****************************************************! ALONG ETA-DIRECTION !***************************************************//
    // SEND SUB-ARRAY ON FRONT SIDE
    saStarts = padWidth;
    loclSize = coreSize;            loclSize(1) = padWidth(1);

    // STAGGERED GRID SHARE A POINT ACROSS SUB-DOMAIN BOUNDARIES.
    // THE DATA ON THIS POINT IS AVERAGED FROM THE VALUES IN BOTH SUB-DOMAINS.
    // HENCE, THE SHARED POINT AND ONE INTERIOR POINT ARE SENT.
    if (yStag) {
        loclSize(1) += 1;
    }

    MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &sendSubarrayY0);
    MPI_Type_commit(&sendSubarrayY0);

    // RECEIVE SUB-ARRAY ON FRONT SIDE
    // FOR FACE-CENTERED (STAGGERED) DATA, THE RECEIVED DATA IS WRITTEN
    // INTO A TEMPORARY ARRAY, AND HENCE USES A DIFFERENT MPI DATATYPE
    if (yStag) {
        count = loclSize(0)*loclSize(1)*loclSize(2);

        MPI_Type_contiguous(count, MPI_FP_REAL, &recvSubarrayY0);

        recvDataY0.resize(loclSize);

        padY0 = blitz::RectDomain<3>(blitz::TinyVector<int, 3>(0, -pSize(1), 0),
                                     blitz::TinyVector<int, 3>(cSize(0) - 1, -1, cSize(2) - 1));

        wallY0 = blitz::RectDomain<3>(blitz::TinyVector<int, 3>(0, 0, 0),
                                      blitz::TinyVector<int, 3>(cSize(0) - 1, 0, cSize(2) - 1));
    } else {
        saStarts = padWidth;            saStarts(1) = 0;
        loclSize = coreSize;            loclSize(1) = padWidth(1);

        MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &recvSubarrayY0);
    }
    MPI_Type_commit(&recvSubarrayY0);


    // SEND SUB-ARRAY ON REAR SIDE
    saStarts = padWidth;            saStarts(1) = coreSize(1);
    loclSize = coreSize;            loclSize(1) = padWidth(1);

    // STAGGERED GRID SHARE A POINT ACROSS SUB-DOMAIN BOUNDARIES.
    // THE DATA ON THIS POINT IS AVERAGED FROM THE VALUES IN BOTH SUB-DOMAINS.
    // HENCE, THE SHARED POINT AND ONE INTERIOR POINT ARE SENT.
    if (yStag) {
        saStarts(1) -= 1;
        loclSize(1) += 1;
    }

    MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &sendSubarrayY1);
    MPI_Type_commit(&sendSubarrayY1);

    // RECEIVE SUB-ARRAY ON REAR SIDE
    // FOR FACE-CENTERED (STAGGERED) DATA, THE RECEIVED DATA IS WRITTEN
    // INTO A TEMPORARY ARRAY, AND HENCE USES A DIFFERENT MPI DATATYPE
    if (yStag) {
        count = loclSize(0)*loclSize(1)*loclSize(2);

        MPI_Type_contiguous(count, MPI_FP_REAL, &recvSubarrayY1);

        recvDataY1.resize(loclSize);

        padY1 = blitz::RectDomain<3>(blitz::TinyVector<int, 3>(0, cSize(1), 0),
                                     blitz::TinyVector<int, 3>(cSize(0) - 1, cSize(1) + pSize(1) - 1, cSize(2) - 1));

        wallY1 = blitz::RectDomain<3>(blitz::TinyVector<int, 3>(0, cSize(1) - 1, 0),
                                      blitz::TinyVector<int, 3>(cSize(0) - 1, cSize(1) - 1, cSize(2) - 1));
    } else {
        saStarts = padWidth;            saStarts(1) = coreSize(1) + padWidth(1);
        loclSize = coreSize;            loclSize(1) = padWidth(1);

        MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &recvSubarrayY1);
    }
    MPI_Type_commit(&recvSubarrayY1);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to send data across all sub-domain faces
 *
 *          This is the core function of the mpidata class.
 *          The end slices of each sub-domain recieves data from their corresponding neighbouring sub-domains,
 *          while the interior slices of each sub-domain sends data to their corresponding neighbouring sub-domains.
 *
 *          All the data slices are send as subarray MPI derived data-types created in the \ref createSubarrays function.
 *          As a result, \ref syncData must be called only after the subarrays have been created.
 *
 *          The data transfer is implemented here with a mixture of blocking and non-blocking communication calls.
 *          The receives are non-blocking, while the sends are blocking. This combination prevents inter-processor deadlock.
 ********************************************************************************************************************************************
 */
void mpidata::syncData() {
    blitz::Range all;

    all = blitz::Range::all();
    recvRequest = MPI_REQUEST_NULL;

    //if (xsFlag) {
    //    if (rankData.rank == 0) std::cout << "0 Before \t" << dataField(blitz::Range(60, blitz::toEnd), 10, 10) << std::endl;
    //    if (rankData.rank == 1) std::cout << "1 Before \t" << dataField(blitz::Range(blitz::fromStart, 6), 10, 10) << std::endl;
    //}
    if (xsFlag) {
        MPI_Irecv(recvDataX0.dataFirst(), 1, recvSubarrayX0, rankData.faceRanks(0), 1, MPI_COMM_WORLD, &recvRequest(0));
        MPI_Irecv(recvDataX1.dataFirst(), 1, recvSubarrayX1, rankData.faceRanks(1), 2, MPI_COMM_WORLD, &recvRequest(1));
    } else {
        MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayX0, rankData.faceRanks(0), 1, MPI_COMM_WORLD, &recvRequest(0));
        MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayX1, rankData.faceRanks(1), 2, MPI_COMM_WORLD, &recvRequest(1));
    }

    if (ysFlag) {
        MPI_Irecv(recvDataY0.dataFirst(), 1, recvSubarrayY0, rankData.faceRanks(2), 3, MPI_COMM_WORLD, &recvRequest(2));
        MPI_Irecv(recvDataY1.dataFirst(), 1, recvSubarrayY1, rankData.faceRanks(3), 4, MPI_COMM_WORLD, &recvRequest(3));
    } else {
        MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayY0, rankData.faceRanks(2), 3, MPI_COMM_WORLD, &recvRequest(2));
        MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayY1, rankData.faceRanks(3), 4, MPI_COMM_WORLD, &recvRequest(3));
    }

    MPI_Send(dataField.dataFirst(), 1, sendSubarrayX0, rankData.faceRanks(0), 2, MPI_COMM_WORLD);
    MPI_Send(dataField.dataFirst(), 1, sendSubarrayX1, rankData.faceRanks(1), 1, MPI_COMM_WORLD);
    MPI_Send(dataField.dataFirst(), 1, sendSubarrayY0, rankData.faceRanks(2), 4, MPI_COMM_WORLD);
    MPI_Send(dataField.dataFirst(), 1, sendSubarrayY1, rankData.faceRanks(3), 3, MPI_COMM_WORLD);

    MPI_Waitall(4, recvRequest.dataFirst(), recvStatus.dataFirst());

    if (xsFlag) {
        //if (rankData.rank == 0) dataField = 0.0;
        //if (rankData.rank == 0) std::cout << "0\t" << dataField(all, 6, 6) << std::endl;
        //if (rankData.rank == 1) std::cout << "1\t" << dataField(all, 6, 6) << std::endl;

        dataField(padX0) = recvDataX0(blitz::Range(0, pSize(0) - 1), all, all);
        dataField(padX1) = recvDataX1(blitz::Range(1, pSize(0)), all, all);

        //if (rankData.rank == 0) std::cout << pSize(0) << std::endl;
        //MPI_Finalize();
        //exit(0);

        if (xsPer) {
            dataField(wallX0) = (dataField(wallX0) + recvDataX0(blitz::Range(pSize(0), pSize(0)), all, all))*0.5;
            dataField(wallX1) = (dataField(wallX1) + recvDataX1(blitz::Range(0, 0), all, all))*0.5;
        } else {
            if (rankData.xRank > 0) {
                dataField(wallX0) = (dataField(wallX0) + recvDataX0(blitz::Range(pSize(0), pSize(0)), all, all))*0.5;
            }
            if (rankData.xRank < rankData.npX - 1) {
                dataField(wallX1) = (dataField(wallX1) + recvDataX1(blitz::Range(0, 0), all, all))*0.5;
            }
        }

        //if (rankData.rank == 0) std::cout << padX0.ubound() << padX0.lbound() << std::endl;
        //if (rankData.rank == 0) std::cout << dataField.ubound() << dataField.lbound() << std::endl;
        //if (rankData.rank == 0) std::cout << recvDataX0.ubound() << recvDataX0.lbound() << std::endl;
        //if (rankData.rank == 0) std::cout << recvDataX1.ubound() << recvDataX1.lbound() << std::endl;

        //if (rankData.rank == 0) std::cout << recvDataX0(all, 6, 6) << std::endl;
        //if (rankData.rank == 0) std::cout << recvDataX1(all, 6, 6) << std::endl;
        //if (rankData.rank == 0) std::cout << dataField(all, 6, 6) << std::endl;
        //if (rankData.rank == 2) std::cout << recvDataX1(6, all, 6) << std::endl;
        //MPI_Finalize();
        //exit(0);
    }

    if (ysFlag) {
        dataField(padY0) = recvDataY0(all, blitz::Range(0, pSize(1) - 1), all);
        dataField(padY1) = recvDataY1(all, blitz::Range(1, pSize(1)), all);

        if (ysPer) {
            dataField(wallY0) = (dataField(wallY0) + recvDataY0(all, blitz::Range(pSize(1), pSize(1)), all))*0.5;
            dataField(wallY1) = (dataField(wallY1) + recvDataY1(all, blitz::Range(0, 0), all))*0.5;
        } else {
            if (rankData.yRank > 0) {
                dataField(wallY0) = (dataField(wallY0) + recvDataY0(all, blitz::Range(pSize(1), pSize(1)), all))*0.5;
            }
            if (rankData.yRank < rankData.npY - 1) {
                dataField(wallY1) = (dataField(wallY1) + recvDataY1(all, blitz::Range(0, 0), all))*0.5;
            }
        }
    }
    //if (xsFlag) {
    //    if (rankData.rank == 0) std::cout << "0 After \t" << dataField(blitz::Range(60, blitz::toEnd), 10, 10) << std::endl;
    //    if (rankData.rank == 1) std::cout << "1 After \t" << dataField(blitz::Range(blitz::fromStart, 6), 10, 10) << std::endl;
    //}
}
