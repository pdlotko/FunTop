//  computePersistenceOfFunction.cpp
#include <iostream>
#include <list>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <stdexcept>


using namespace std;

#include "../../include/cell.hpp"
#include "../../include/cellComplex.hpp"
#include "configuration.hpp"
#include "../../include/partiallyConstantApproximationOfFunctionOnTopDimensionalCells.h"

//this program compute partiall constatn approximation of a vector-valued 
//function by using only top dimensional rectangles.

//In this example, we consider a function f: R^2 -> R^2 defined on a bax
//[0,2]x[0,2]. In this example:
//f(x,y) = [ x+y , x^2+y^2 ]
//Epsilon is set to 0.5.

//N.B. you may need to run:
//export LD_LIBRARY_PATH=/path_to_your_cxsc/cxsc-2-5-4/lib:$LD_LIBRARY_PATH
//in case you get an error while loading shared libraries!!

int main(int argc,char **argv)
{
    std::vector< std::pair<double,double> > initialRectangle;
    initialRectangle.push_back( std::make_pair(0,2) );
    initialRectangle.push_back( std::make_pair(0,2) );
    
    double epsilon = 0.1;

    std::vector< HessType (*)( const HTvector& x ) > followingCoordinatesOfFunction;
    followingCoordinatesOfFunction.push_back(xPlusY);
    followingCoordinatesOfFunction.push_back(xSquared_plus_ySquared);

    try
    {
        clock_t begin = clock();
        topDimensionalCellDecompositionOfDomain decomposition( initialRectangle , followingCoordinatesOfFunction , epsilon);
        decomposition.writeResultToFile((char*)"PC_approximation_of_vector_valued_function");
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << "Time elapsed : " << elapsed_secs << " seconds \n";
    }
    catch ( const std::runtime_error& e )
    {
        cerr << "We are back in main(). Caught exception : " << e.what() << endl;
        return 0;
    }

    cout << "We are done, please press enter!\n";
    std::cin.ignore();
    return 0;
}
