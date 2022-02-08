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
#include "../../include/computePersistenceOfFunction.h"



//this is the main function in this implementation. It computes persistence of a function.
//In this example, the Ackley's function on the domain [0,1]x[0,1].
//The accuracy (epsilon parameter) is set to 0.02 

//N.B. you may need to run:
//export LD_LIBRARY_PATH=/path_to_your_cxsc/cxsc-2-5-4/lib:$LD_LIBRARY_PATH
//in case you get an error while loading shared libraries!!

int main(int argc,char **argv)
{
	//open log file
	ofstream log;
	log.open("LOG");
	
	//The function will be approximated on [0,1]x[0,1]. Please select here
	//the right box for your purpose.
    std::vector< std::pair<double,double> > range;
    range.push_back( std::make_pair(0,1) );
    range.push_back( std::make_pair(0,1) );
    
    //Select the accuracy of the approximation
    double epsilon = 0.1;
   
    std::vector< std::vector< std::pair<cxsc::real,cxsc::real> > > intervals;
    try
    {
		clock_t begin = clock();
		//At the moment we select Ackley's function for the computations.
		//Please consider other options, or implenent your own function in 
		//configuration.hpp
		intervals = computePersistenceOfAFunction<double>(range, Ackleys, epsilon );
		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		std::cout << "Time elapsed : " << elapsed_secs << " seconds \n";
		log <<  "Time elapsed : " << elapsed_secs << " seconds \n";
    }
    catch ( const std::runtime_error& e )
    {
        cerr << "We are back in main(). Caught exception : " << e.what() << endl;
        log <<  "We are back in main(). Caught exception : " << e.what() << endl;
        return 0;
    }

    cout << "Writing persistence to a file \n";
    log << "Writing persistence to a file \n";
    for ( size_t dim = 0 ; dim != intervals.size() ; ++dim )
    {
		ofstream out;
		ostringstream name;
        name << "BarcodesInDimension_" << dim << ".txt";
        log << "BarcodesInDimension_" << dim << ".txt";
        std::string nameStr  = name.str();
        const char* filename1 = nameStr.c_str();
        out.open( filename1 );
		for ( size_t i = 0 ;  i != intervals[dim].size() ; ++i )
        {
			out << intervals[dim][i].first << " " << intervals[dim][i].second << endl;
		}
		out.close();
    }
    log.close();
    return 0;
}
