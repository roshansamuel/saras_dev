#include "force.h"

force::force(vfield &U, const parser &solParams, parallel &mpiParam):
    V(U),
    inputParams(solParams),
    mpiData(mpiParam)
{
    if (inputParams.probType < 5){
        Fr = 1.0/inputParams.Ro;
    }
    else{
        if (inputParams.rbcType == 1) {                     
            Fb = inputParams.Ra*inputParams.Pr;
            Fr = inputParams.Pr*sqrt(inputParams.Ta);
        } else if (inputParams.rbcType == 2) {    
            Fb = 1.0;
            Fr = sqrt(inputParams.Ta*inputParams.Pr/inputParams.Ra); 
        } else if (inputParams.rbcType == 3) {    
            Fb = inputParams.Ra;
            Fr = sqrt(inputParams.Ta);
        } else {                                                            
            Fb = inputParams.Pr;
            Fr = sqrt(inputParams.Ta*inputParams.Pr/inputParams.Ra); 
        }
    }
}


void force::add_VForce(plainvf &Hv)
{
    if (inputParams.Force==0){ //No forcing
    }    
    else if (inputParams.Force==1){ //Random forcing
        add_RandomForce(Hv);
    } 

    else if (inputParams.Force==2){ //Coriolis forcing
        add_Coriolis(Hv);
    }

    else if (inputParams.Force <= 4) { // Throw error message, since forcing 3 and 4 are not allowed for hydrodynamic flows
        if (mpiData.rank==0){
            std::cout<<"This forcing is not allowed for the given problem type. Aborting..."<<std::endl;
            exit(0);
        }
    }

    else {
        if (mpiData.rank==0){
            std::cout<<"Invalid Forcing. Program aborting..."<<std::endl;
            exit(0);
        }
    }

    
}

void force::add_VForce(plainvf &Hv, sfield &T)
{
    if (inputParams.Force==0){ //No forcing
    }

    else if (inputParams.Force==1){ //Random forcing
        add_RandomForce(Hv);
    } 

    else if (inputParams.Force==2){ //Coriolis forcing
        add_Coriolis(Hv);
    }

    else if (inputParams.Force==3){ //Buoyancy forcing
        add_Buoyancy(Hv, T);
    }
    
    else if (inputParams.Force==4){ //Buoyancy and Coriolis forcing
        add_Buoyancy(Hv, T);
        add_Coriolis(Hv);
    }

    else {
        if (mpiData.rank==0){
            std::cout<<"Invalid Forcing. Program aborting..."<<std::endl;
            exit(0);
        }
    }
}

void force::add_RandomForce(plainvf &Hv){
    blitz::Array<double, 3> Force_x, Force_y, Force_z;

    Force_x.resize(V.Vx.fSize);
    Force_x.reindexSelf(V.Vx.flBound);
    Force_x = 0; //Needs to be developed; currently set to zero
    Hv.Vx += Force_x;

#ifndef PLANAR
    Force_y.resize(V.Vy.fSize);
    Force_y.reindexSelf(V.Vy.flBound);
    Force_y = 0; //Needs to be developed; currently set to zero
    Hv.Vy += Force_y;
#endif
    Force_z.resize(V.Vz.fSize);
    Force_z.reindexSelf(V.Vz.flBound);
    Force_z = 0; //Needs to be developed; currently set to zero
    Hv.Vz += Force_z;
}


void force::add_Buoyancy(plainvf &Hv, sfield &T) {
    V.interPc2Vz = 0.0;
    for (unsigned int i=0; i < V.Vz.PcIntSlices.size(); i++) {
        V.interPc2Vz(V.Vz.fCore) += T.F.F(V.Vz.PcIntSlices(i));
    }   
    V.interPc2Vz /= V.Vz.PcIntSlices.size();
    Hv.Vz += Fb*V.interPc2Vz;
}

void force::add_Coriolis(plainvf &Hv){
#ifndef PLANAR
    //ADD THE ROTATING TERM IN THE Vx COMPONENT OF Hv
    V.interVy2Vx = 0.0;
    for (unsigned int i=0; i < V.Vx.VyIntSlices.size(); i++) {
        V.interVy2Vx(V.Vx.fCore) += V.Vy.F(V.Vx.VyIntSlices(i));
    }   
    V.interVy2Vx /= V.Vx.VyIntSlices.size();

    Hv.Vx += Fr*V.interVy2Vx;

    //SUBTRACT THE ROTATING TERM IN THE Vy COMPONENT of Hv
    V.interVx2Vy = 0.0;
    for (unsigned int i=0; i < V.Vy.VxIntSlices.size(); i++) {
        V.interVx2Vy(V.Vy.fCore) += V.Vx.F(V.Vy.VxIntSlices(i));
    }   
    V.interVx2Vy /= V.Vy.VxIntSlices.size();

    Hv.Vy -= Fr*V.interVx2Vy;
#endif
}

void::force::add_SForce(plainsf &Ht) {
    Ht.F += 0; //Scalar forcing needs to be implemented
}

force::~force() { }
