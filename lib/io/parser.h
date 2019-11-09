#ifndef PARSER_H
#define PARSER_H

#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <blitz/array.h>
#include <yaml-cpp/yaml.h>

class parser {
    private:
        std::string domainType;
        std::string probeCoords;

        void parseYAML();
        void checkData();

        void testProbes();
        void parseProbes();

        void setGrids();
        void setPeriodicity();

    public:
        int ioCnt;
        int rbcType;
        int nThreads;
        int npY, npX;
        int forceType;
        int xInd, yInd, zInd;
        int vcDepth, vcCount;
        int preSmooth, postSmooth;

        int iScheme;
        int probType;
        int xGrid, yGrid, zGrid;

        bool useCFL;
        bool nonHgBC;
        bool readProbes;
        bool restartFlag;
        bool xPer, yPer, zPer;

        double Re;
        double Ra;
        double Pr;
        double Ta;
        double Ro;

        double fwInt;
        double rsInt;
        double prInt;
        double tolerance;
        double Lx, Ly, Lz;
        double tStp, tMax;
        double patchRadius;
        double betaX, betaY, betaZ;
        double courantNumber;

        std::string dScheme;
        std::string meshType;

        std::vector<int> interSmooth;
        std::vector<blitz::TinyVector<int, 3> > probesList;

        parser();

        void writeParams();
};

/**
 ********************************************************************************************************************************************
 *  \class parser parser.h "lib/parser.h"
 *  \brief  Contains all the global variables set by the user through the yaml file
 *
 *  The class parses the paramters.yaml file and stores all the simulation paramters in publicly accessible constants.
 *  The class also has a function to check the consistency of the user set paramters and throw exceptions.
 *  The class is best initialized as a constant to prevent inadvertent tampering of the global variables it contains.
 ********************************************************************************************************************************************
 */

#endif
