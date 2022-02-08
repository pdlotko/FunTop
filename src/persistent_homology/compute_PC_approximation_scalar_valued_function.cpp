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

//this program compute partiall constant approximation of a function on
//the top dimensional recatngles and return the result (be careful, the 
//output file may be large!

//In this example, the Rastrigin3d function on the domain [-2,2]x[-2,2]x[-2,2] is considered. 
//The accuracy (epsilon parameter) is set to 2 

//N.B. you may need to run:
//export LD_LIBRARY_PATH=/path_to_your_cxsc/cxsc-2-5-4/lib:$LD_LIBRARY_PATH
//in case you get an error while loading shared libraries!!

int main(int argc,char **argv)
{
    std::vector< std::pair<double,double> > initialRectangle;
    initialRectangle.push_back( std::make_pair(-2,2) );
    initialRectangle.push_back( std::make_pair(-2,2) );
    initialRectangle.push_back( std::make_pair(-2,2) );
    
    double epsilon = 2;
 
    try
    {
        clock_t begin = clock();
        topDimensionalCellDecompositionOfDomain decomposition( initialRectangle , Rastrigin3d , epsilon);
        std::cout << "NUmber of validated cells : " << decomposition.number_of_validated_cells() << std::endl;
        
        decomposition.writeResultToFile((char*)"PC_approximation");
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << "Time elapsed : " << elapsed_secs << " seconds \n";
    }
    catch ( const std::runtime_error& e )
    {
        cerr << "We are back in main(). Caught exception : " << e.what() << endl;
        return 0;
    }
    
    

    return 0;
}
