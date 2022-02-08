// computePersistenceOfFunction.h
//
#include <queue>
#include <stdexcept>
#include <set>
#include <fstream>
//Phat include
#include <phat/compute_persistence_pairs.h>
#include <phat/representations/vector_vector.h>
#include <phat/representations/vector_set.h>
#include <phat/representations/vector_list.h>
#include <phat/representations/sparse_pivot_column.h>
#include <phat/representations/full_pivot_column.h>
#include <phat/representations/bit_tree_pivot_column.h>
#include <phat/algorithms/twist_reduction.h>
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/chunk_reduction.h>




template <typename T>
bool compareCellsAccordingToFiltration( cell<T>* c1 , cell<T>* c2 )
{
    return ( c1->val() < c2->val() );
}

/*
template <typename T>
class compareCellsAccordingToFiltration
{
public:
  compareCellsAccordingToFiltration(){}
  bool operator() (const cell<T>* c1 , const cell<T>* c2) const
  {
     return ( c1->val() < c2->val() );
  }
};
*/




//this function returns the persistence intervals. The last parameter indicate if the program should output the intervals which are shorter than the given value of error. Due to the stability theorem for persistence,
//we cannot guarantee that those intervals are not the noise. That is why it is set default to false.
bool computePersistenceOfAFunctionDBG = false;
template <typename T>
std::vector< std::vector< std::pair<cxsc::real,cxsc::real> > > computePersistenceOfAFunction( std::vector< std::pair<T,T> > point , HessType (*f)( const HTvector& x ) , double epsilon , cxsc::real min_ = -INT_MAX,  cxsc::real max_ = INT_MAX, bool showIntervalsShorterThanError = false )throw (std::runtime_error)
{
    unsigned dimensionOfDomain = point.size();
    //for random subdivsions in case they are needed.
    srand( time(0) );

    cellComplex<T>* cmplx = new cellComplex<T>( point );
    if (computePersistenceOfAFunctionDBG)
    {
        cerr << "Initial complex, i.e. the complex build on a given initial rectangle : " << endl <<  *cmplx << endl<< endl<< endl<< endl;
    }

    cell<T>* initialCube = *cmplx->elemen()[ point.size() ].begin();
    if (computePersistenceOfAFunctionDBG)cerr << *initialCube << endl;
    extern unsigned initialSubdivision;
    if ( initialSubdivision )

    //CAUTION! Do not use initialCube after this procedure. The initialCube is deleted there! That caused me hard to find segmentation faults already
    cmplx->divideTheCubeByGivenFactor(initialCube,initialSubdivision);

    std::queue< cell<T>* > queueOfCells;
    for ( size_t i = 0 ; i != cmplx->elemen()[cmplx->elemen().size()-1].size() ; ++i )
    {
        queueOfCells.push( cmplx->elemen()[cmplx->elemen().size()-1][i] );
    }

    while ( !queueOfCells.empty() )
    {
        cell<T>* currentCell = queueOfCells.front();
        queueOfCells.pop();

        if ( currentCell->isVerified() )continue;

        if ( computePersistenceOfAFunctionDBG ) cerr << "currentCell : " << *currentCell << endl;

        HessType fx;
        interval f_x;
        HTvector xx = currentCell->getHTvectorRepresentationOfCell();
        fx = f(xx);
        f_x = fValue(fx);

        if ( computePersistenceOfAFunctionDBG )
        {
            cerr <<  "diam(f_x) : " << diam(f_x) << endl;
            cerr << "2*epsilon : " << 2*epsilon << endl;
        }

        if ( (min_ != -INT_MAX) || (max_ != INT_MAX) )
        {
            if ( Sup(f_x) < min_ )
            {
                currentCell->isVerified() = true;
                currentCell->val() = min_;
            }
            if ( Inf(f_x) > max_ )
            {
                currentCell->isVerified() = true;
                currentCell->val() = max_;
            }
        }
        //in case we have verified the cell above, continue.
        if ( currentCell->isVerified() )continue;



        if ( diam(f_x) <= 2*epsilon )
        {
            if ( computePersistenceOfAFunctionDBG )
            {
                cerr << "Difeerences of values of the function of a cell is not greater than : " << 2*epsilon << ", we are cone with this cell." << endl;
                getchar();
            }
            //we are done for this cell
            currentCell->isVerified() = true;
            currentCell->val() = 0.5*(Sup(f_x)+Inf(f_x));
        }
        else
        {
            if ( computePersistenceOfAFunctionDBG ) cerr << "We need to subdivide this cell \n";
            //subdivide if we did not reach max subdivision depth:
            if ( currentCell->subdivDepth()+1 == subdivisionDepth )
            {
                if ( computePersistenceOfAFunctionDBG )
                {
                    cerr << "1 Maximal subdivision depth has been reached, program terminates \n";
                    if(computePersistenceOfAFunctionDBG)getchar();
                }
                throw std::runtime_error("Too many subdivisions\n");
            }


            std::vector< cell<double>* > newCells = cmplx->divideCubeInAllDirectionsItHasNonzeroLenght(currentCell);
            //and add the resulting cells of the top dimension to the queue.
            for ( size_t i = 0 ; i != newCells.size() ; ++i )
            {
                if ( newCells[i]->dim() == dimensionOfDomain )
                {
                    queueOfCells.push(newCells[i]);
                }
            }
        }
    }

    if (computePersistenceOfAFunctionDBG )cerr << "Exiting from while, imposing lower star filtration \n";


    //now, we filter all the cells according to lower star filtration:
    for ( unsigned dim = dimensionOfDomain ; dim > 0 ; --dim )
    {
        for ( size_t nr = 0 ; nr != cmplx->elemen()[dim].size() ; ++nr )
        {
            for ( typename cell<T>::BdIterator bd = cmplx->elemen()[dim][nr]->bdBegin() ; bd != cmplx->elemen()[dim][nr]->bdEnd() ; ++bd )
            {
                if ( (*bd)->val() > cmplx->elemen()[dim][nr]->val() )
                {
                    (*bd)->val() = cmplx->elemen()[dim][nr]->val();
                }
            }
        }
    }

    //if (computePersistenceOfAFunctionDBG )
    cerr << "Done with imposing lower star filtration \n";



    phat::boundary_matrix< phat::vector_vector > boundary_matrix;

    int numberOfCells = 0;
    for ( int dim = 0 ; dim != cmplx->elemen().size() ; ++dim )
	numberOfCells += cmplx->elemen()[dim].size();


    //put all the cells into a vector:
    std::vector< cell<T>* > cells(numberOfCells);
    unsigned counter = 0;
    for ( int dim = 0 ; dim != cmplx->elemen().size() ; dim++ )
    {
        for ( size_t i = 0 ; i != cmplx->elemen()[dim].size() ; ++i )
        {
            cells[counter] = cmplx->elemen()[dim][i];
            ++counter;
        }
    }

    if (computePersistenceOfAFunctionDBG )cerr << "cells.size() : " << cells.size() << " , numberOfCells : " << numberOfCells << endl;


    if (computePersistenceOfAFunctionDBG )cerr << "Begin sorting\n";

    //sort the vector according to the filtration values.
    std::sort( cells.begin() , cells.end() , compareCellsAccordingToFiltration<T> );

    if (computePersistenceOfAFunctionDBG )cerr << "Done with sorting\n";



    cerr << "numberOfCells : " << numberOfCells << endl;

    // set the number of columns
    boundary_matrix.set_num_cols( numberOfCells );
    // set the dimension of the cell that a column represents:






    cell<T>** numberOfGeneratorToDim = new cell<T>*[numberOfCells+1];
    int nrOfGenerator = 0;
    for ( size_t i = 0 ; i != cells.size() ; ++i )
    {
            boundary_matrix.set_dim( nrOfGenerator, cells[i]->dim() );
            cells[i]->number() = nrOfGenerator;
            numberOfGeneratorToDim[nrOfGenerator] = cells[i];
            ++nrOfGenerator;
    }

    // set the respective columns -- the columns entries have to be sorted
    std::vector< phat::index > temp_col;
    for ( size_t i = 0 ; i != cells.size() ; ++i )
    {
        cell<T>* aa = cells[i];
        int numberElInBd = 0;
        for ( typename cell<T>::BdIterator bd = aa->bdBegin() ; bd != aa->bdEnd() ; ++bd )
        {
            temp_col.push_back( (*bd)->numberInCmplx() );
            numberElInBd++;
        }

        std::sort( temp_col.begin() , temp_col.end() );
        boundary_matrix.set_col( cells[i]->numberInCmplx() , temp_col );
        temp_col.clear();
    }

    // define the object to hold the resulting persistence pairs
    phat::persistence_pairs pairs;


    //TODO -- pick the best option!
    // choose an algorithm (choice affects performance) and compute the persistence pair
    // (modifies boundary_matrix)
    //phat::compute_persistence_pairs< phat::chunk_reduction >( pairs, boundary_matrix );
    cerr << "Entering phat \n";
    //phat::compute_persistence_pairs_dualized< phat::chunk_reduction >( pairs, boundary_matrix );
    phat::compute_persistence_pairs< phat::chunk_reduction >( pairs, boundary_matrix );
    cerr << "Exiting phat \n";

    std::vector< std::vector< std::pair<cxsc::real,cxsc::real> > > intervals;
    //first, initialize this structure:
    for ( size_t dim = 0 ; dim != cmplx->elemen().size() ; ++dim )
    {
	  std::vector< std::pair<cxsc::real,cxsc::real> > intInThisDimension;
	  intervals.push_back(intInThisDimension);
    }


    bool* pairedCollumns = new bool[numberOfCells];
    for ( size_t i = 0 ; i != numberOfCells ; ++i )
    {
	  pairedCollumns[i] = false;
    }
    for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
    {
        if (computePersistenceOfAFunctionDBG )std::cout << "Birth: " << pairs.get_pair( idx ).first << ", Death: " << pairs.get_pair( idx ).second << std::endl;
            if (computePersistenceOfAFunctionDBG )std::cout << "Birth: " << pairs.get_pair( idx ).first << ", Death: " << pairs.get_pair( idx ).second << std::endl;
        cell<T>* first = numberOfGeneratorToDim[pairs.get_pair( idx ).first];
        cell<T>* second = numberOfGeneratorToDim[pairs.get_pair( idx ).second];
	 pairedCollumns[pairs.get_pair( idx ).first] = true;
	 pairedCollumns[pairs.get_pair( idx ).second] = true;

	 if ( first->val() != second->val() )
	 {
		if (computePersistenceOfAFunctionDBG )cerr << "Pair, if : " << *first << "of a dimension : " << first->dim() << " and value : " << first->val() << "\n and \n" << *second << " of a dimension : " << second->dim() << " and value : " << second->val()  << endl;
		if ( first->val() < second->val() )
		{
			if (!showIntervalsShorterThanError)
			{
				if ( second->val()-first->val() > epsilon )
				{
					//cerr << "second->val()-first->val() : " << second->val()-first->val() << endl;
	 				intervals[ first->dim() ].push_back( std::make_pair( first->val() , second->val() ) );
				}
			}
			else
			{
				intervals[ first->dim() ].push_back( std::make_pair( first->val() , second->val() ) );
			}
		}
		else
		{
			if (!showIntervalsShorterThanError)
			{
				if ( first->val()-second->val() > epsilon )
				{
					cerr << "first->val()-second->val() : " << first->val()-second->val() << endl;
	 				intervals[ first->dim() ].push_back( std::make_pair( second->val() , first->val() ) );
				}
			}
			else
			{
				intervals[ first->dim() ].push_back( std::make_pair( second->val() , first->val() ) );
			}
		}
	 }
    }

    if (computePersistenceOfAFunctionDBG )
    {
    	cerr << "boundary_matrix.get_num_cols() : " << boundary_matrix.get_num_cols() << endl;
    	cerr << "numberOfCells  : " << numberOfCells  << endl;
    }


    for ( size_t i = 0 ; i != numberOfCells ; ++i )
    {
	  if ( !pairedCollumns[i] )
	  {
              if (computePersistenceOfAFunctionDBG )
		{
			cerr << "The collumn : " << i << " is empty and it is not an element of any pair!! \n";
			cerr << *numberOfGeneratorToDim[i] << endl;
		}
		//TODO -- set infinity into some better way!
		intervals[ numberOfGeneratorToDim[i]->dim() ].push_back( std::make_pair( numberOfGeneratorToDim[i]->val() , 1000000 ) );
	  }
    }



    for ( size_t dim = 0 ; dim != cmplx->elemen().size() ; ++dim )
    {
	  if (computePersistenceOfAFunctionDBG )cerr << "intervals["<<dim<<"].size() : " << intervals[dim].size() << endl;
    }





    cxsc::real level = 0.1;

    for ( size_t dim = 0 ; dim != cmplx->elemen().size() ; ++dim )
    {
          for ( size_t nr  = 0 ; nr != cmplx->elemen()[dim].size() ; ++nr )
	   {
		if ( cmplx->elemen()[dim][nr]->val() < level  )
		{
//cerr << "Deleteing cell in dimension : " << dim << endl;
			cmplx->elemen()[dim][nr]->del();
		}
	   }
    }





    return intervals;
}//computeNodalDomainOfAFunction


























































template <typename T>
std::vector< std::vector< std::pair<cxsc::real,cxsc::real> > > computePersistenceOfAFunctionUseInformationAboutDerivatives( std::vector< std::pair<T,T> > point , HessType (*f)( const HTvector& x ) , double epsilon , cxsc::real min_ = -INT_MAX,  cxsc::real max_ = INT_MAX, bool showIntervalsShorterThanError = false )throw (std::runtime_error)
{
    bool computePersistenceOfAFunctionDBG = false;

    unsigned dimensionOfDomain = point.size();
    //for random subdivsions in case they are needed.
    srand( time(0) );

    cellComplex<T>* cmplx = new cellComplex<T>( point );
    if (computePersistenceOfAFunctionDBG)
    {
        cerr << "Initial complex, i.e. the complex build on a given initial rectangle : " << endl <<  *cmplx << endl<< endl<< endl<< endl;
    }

    cell<T>* initialCube = *cmplx->elemen()[ point.size() ].begin();
    if (computePersistenceOfAFunctionDBG)cerr << *initialCube << endl;
    extern unsigned initialSubdivision;
    if ( initialSubdivision )

    //CAUTION! Do not use initialCube after this procedure. The initialCube is deleted there! That caused me hard to find segmentation faults already
    cmplx->divideTheCubeByGivenFactor(initialCube,initialSubdivision);

    std::queue< cell<T>* > queueOfCells;
    for ( size_t i = 0 ; i != cmplx->elemen()[cmplx->elemen().size()-1].size() ; ++i )
    {
        queueOfCells.push( cmplx->elemen()[cmplx->elemen().size()-1][i] );
    }

    bool doWeHaveToCareAbout123456789Cells = false;

    while ( !queueOfCells.empty() )
    {
        cell<T>* currentCell = queueOfCells.front();
        queueOfCells.pop();

        if ( currentCell->isVerified() )continue;

        if ( computePersistenceOfAFunctionDBG ) cerr << "currentCell : " << *currentCell << endl;

        HessType fx;
        interval f_x;
        HTvector xx = currentCell->getHTvectorRepresentationOfCell();
        fx = f(xx);
        f_x = fValue(fx);

        ivector Gfx(dimensionOfCoDomain);
        Gfx = gradValue(fx);

        bool doesDerivativesInAllDirectionsContainZero = true;
        for ( size_t dir = 0 ; dir != currentCell->coords().size() ; ++dir )
        {
            if ( currentCell->coords()[dir].first == currentCell->coords()[dir].second )continue;
            interval directionalDeriv = Gfx[dir+1];
            if ( !(0.0 <= directionalDeriv) )
            {
                doesDerivativesInAllDirectionsContainZero = false;
                break;
            }
        }

        if ( doesDerivativesInAllDirectionsContainZero == false )
        {
            //function is monotone in this direction, there is no need to validate this cell, since homology do not change here.
            currentCell->isVerified() = true;
            currentCell->val() = -123456789;//to be corrected
            doWeHaveToCareAbout123456789Cells = true;
            //in case we have verified the cell above, continue.
            continue;
        }


        if ( computePersistenceOfAFunctionDBG )
        {
            cerr <<  "diam(f_x) : " << diam(f_x) << endl;
            cerr << "2*epsilon : " << 2*epsilon << endl;
        }

        if ( (min_ != -INT_MAX) || (max_ != INT_MAX) )
        {
            if ( Sup(f_x) < min_ )
            {
                currentCell->isVerified() = true;
                currentCell->val() = min_;
                //in case we have verified the cell above, continue.
                continue;
            }
            if ( Inf(f_x) > max_ )
            {
                currentCell->isVerified() = true;
                currentCell->val() = max_;
                //in case we have verified the cell above, continue.
                continue;
            }
        }



        //if non of the conditions above holds, we continue.
        if ( diam(f_x) <= 2*epsilon )
        {
            if ( computePersistenceOfAFunctionDBG )
            {
                cerr << "Difeerences of values of the function of a cell is not greater than : " << 2*epsilon << ", we are cone with this cell." << endl;
                getchar();
            }
            //we are done for this cell
            currentCell->isVerified() = true;
            currentCell->val() = 0.5*(Sup(f_x)+Inf(f_x));
        }
        else
        {
            if ( computePersistenceOfAFunctionDBG ) cerr << "We need to subdivide this cell \n";
            //subdivide if we did not reach max subdivision depth:
            if ( currentCell->subdivDepth()+1 == subdivisionDepth )
            {
                if ( computePersistenceOfAFunctionDBG )
                {
                    cerr << "1 Maximal subdivision depth has been reached, program terminates \n";
                    if(computePersistenceOfAFunctionDBG)getchar();
                }
                throw std::runtime_error("Too many subdivisions\n");
            }


            std::vector< cell<double>* > newCells = cmplx->divideCubeInAllDirectionsItHasNonzeroLenght(currentCell);
            //and add the resulting cells of the top dimension to the queue.
            for ( size_t i = 0 ; i != newCells.size() ; ++i )
            {
                if ( newCells[i]->dim() == dimensionOfDomain )
                {
                    queueOfCells.push(newCells[i]);
                }
            }
        }
    }

    if ( doWeHaveToCareAbout123456789Cells )
    {
        cerr << "while doWeHaveToCareAbout123456789Cells \n";
        while ( true )
        {

            bool do_we_need_to_continue = false;
            size_t dim = cmplx->elemen().size()-1;
            for ( size_t nr = 0 ; nr != cmplx->elemen()[ dim ].size() ; ++nr )
            {
                if ( cmplx->elemen()[ dim ][nr]->val() == -123456789 )
                {
                    /*
                    cxsc::real max_value_of_neigh = -INT_MAX;
                    //find maximal neighbor:
                    for ( typename cell<T>::BdIterator bd = cmplx->elemen()[dim][nr]->bdBegin() ; bd != cmplx->elemen()[dim][nr]->bdEnd() ; ++bd )
                    {
                        for ( typename cell<T>::CbdIterator cbd = (*bd)->cbdBegin() ; cbd != (*bd)->cbdEnd() ; ++cbd )
                        {
                            if ( (*cbd)->val() > max_value_of_neigh )max_value_of_neigh = (*cbd)->val();
                        }
                    }
                    cmplx->elemen()[ dim ][nr]->val() = max_value_of_neigh;
                    if ( max_value_of_neigh == -123456789 )
                    {
                        do_we_need_to_continue = true;
                    }
                    */
                    HessType fx;
                    interval f_x;
                    HTvector xx = cmplx->elemen()[ dim ][nr]->getHTvectorRepresentationOfCell();
                    fx = f(xx);
                    f_x = fValue(fx);

                    ivector Gfx(dimensionOfCoDomain);
                    Gfx = gradValue(fx);

                    int direction_in_which_function_is_monotone = -1;
                    for ( size_t dir = 0 ; dir != cmplx->elemen()[ dim ][nr]->coords().size() ; ++dir )
                    {
                        if ( cmplx->elemen()[ dim ][nr]->coords()[dir].first == cmplx->elemen()[ dim ][nr]->coords()[dir].second )continue;
                        interval directionalDeriv = Gfx[dir+1];
                        if ( !(0.0 <= directionalDeriv) )
                        {
                            direction_in_which_function_is_monotone = dir;
                            break;
                        }
                    }

                    cerr << "We have 123456789 cell. direction_in_which_function_is_monotone : " << direction_in_which_function_is_monotone << endl;

                    if ( direction_in_which_function_is_monotone == -1 )throw("This should not happened \n");

                    //now we need to find neighbor in the direction direction_in_which_function_is_monotone
                    cxsc::real maxx = -INT_MAX;
                    for ( typename cell<T>::BdIterator bd = cmplx->elemen()[dim][nr]->bdBegin() ; bd != cmplx->elemen()[dim][nr]->bdEnd() ; ++bd )
                    {
                        for ( typename cell<T>::CbdIterator cbd = (*bd)->cbdBegin() ; cbd != (*bd)->cbdEnd() ; ++cbd )
                        {
                            if ( (*cbd)->coords()[direction_in_which_function_is_monotone].first == cmplx->elemen()[dim][nr]->coords()[direction_in_which_function_is_monotone].second
                                 ||
                                 (*cbd)->coords()[direction_in_which_function_is_monotone].second == cmplx->elemen()[dim][nr]->coords()[direction_in_which_function_is_monotone].first
                               )
                            {
                                if ( (*cbd)->val() > maxx )maxx = (*cbd)->val();
                            }
                        }
                    }

                    cerr << "New value : " << maxx << endl;

                    cmplx->elemen()[ dim ][nr]->val() = maxx;
                    if ( maxx == -123456789 )
                    {
                        do_we_need_to_continue = true;
                    }
                }
            }
            if ( !do_we_need_to_continue )break;
        }
    }



    if (computePersistenceOfAFunctionDBG )cerr << "Exiting from while, imposing lower star filtration \n";


    //now, we filter all the cells according to lower star filtration:
    for ( unsigned dim = dimensionOfDomain ; dim > 0 ; --dim )
    {
        for ( size_t nr = 0 ; nr != cmplx->elemen()[dim].size() ; ++nr )
        {
            for ( typename cell<T>::BdIterator bd = cmplx->elemen()[dim][nr]->bdBegin() ; bd != cmplx->elemen()[dim][nr]->bdEnd() ; ++bd )
            {
                if ( (*bd)->val() > cmplx->elemen()[dim][nr]->val() )
                {
                    (*bd)->val() = cmplx->elemen()[dim][nr]->val();
                }
            }
        }
    }
    //if (computePersistenceOfAFunctionDBG )
    cerr << "Done with imposing lower star filtration \n";



    phat::boundary_matrix< phat::vector_vector > boundary_matrix;

    int numberOfCells = 0;
    for ( int dim = 0 ; dim != cmplx->elemen().size() ; ++dim )
	numberOfCells += cmplx->elemen()[dim].size();


    //put all the cells into a vector:
    std::vector< cell<T>* > cells(numberOfCells);
    unsigned counter = 0;
    for ( int dim = 0 ; dim != cmplx->elemen().size() ; dim++ )
    {
        for ( size_t i = 0 ; i != cmplx->elemen()[dim].size() ; ++i )
        {
            cells[counter] = cmplx->elemen()[dim][i];
            ++counter;
        }
    }

    if (computePersistenceOfAFunctionDBG )cerr << "cells.size() : " << cells.size() << " , numberOfCells : " << numberOfCells << endl;


    if (computePersistenceOfAFunctionDBG )cerr << "Begin sorting\n";

    //sort the vector according to the filtration values.
    std::sort( cells.begin() , cells.end() , compareCellsAccordingToFiltration<T> );

    if (computePersistenceOfAFunctionDBG )cerr << "Done with sorting\n";



    cerr << "numberOfCells : " << numberOfCells << endl;

    // set the number of columns
    boundary_matrix.set_num_cols( numberOfCells );
    // set the dimension of the cell that a column represents:






    cell<T>** numberOfGeneratorToDim = new cell<T>*[numberOfCells+1];
    int nrOfGenerator = 0;
    for ( size_t i = 0 ; i != cells.size() ; ++i )
    {
            boundary_matrix.set_dim( nrOfGenerator, cells[i]->dim() );
            cells[i]->number() = nrOfGenerator;
            numberOfGeneratorToDim[nrOfGenerator] = cells[i];
            ++nrOfGenerator;
    }

    // set the respective columns -- the columns entries have to be sorted
    std::vector< phat::index > temp_col;
    for ( size_t i = 0 ; i != cells.size() ; ++i )
    {
        cell<T>* aa = cells[i];
        int numberElInBd = 0;
        for ( typename cell<T>::BdIterator bd = aa->bdBegin() ; bd != aa->bdEnd() ; ++bd )
        {
            temp_col.push_back( (*bd)->numberInCmplx() );
            numberElInBd++;
        }

        std::sort( temp_col.begin() , temp_col.end() );
        boundary_matrix.set_col( cells[i]->numberInCmplx() , temp_col );
        temp_col.clear();
    }

    // define the object to hold the resulting persistence pairs
    phat::persistence_pairs pairs;


    //TODO -- pick the best option!
    // choose an algorithm (choice affects performance) and compute the persistence pair
    // (modifies boundary_matrix)
    //phat::compute_persistence_pairs< phat::chunk_reduction >( pairs, boundary_matrix );
    cerr << "Entering phat \n";
    //phat::compute_persistence_pairs_dualized< phat::chunk_reduction >( pairs, boundary_matrix );
    phat::compute_persistence_pairs< phat::chunk_reduction >( pairs, boundary_matrix );
    cerr << "Exiting phat \n";

    std::vector< std::vector< std::pair<cxsc::real,cxsc::real> > > intervals;
    //first, initialize this structure:
    for ( size_t dim = 0 ; dim != cmplx->elemen().size() ; ++dim )
    {
	  std::vector< std::pair<cxsc::real,cxsc::real> > intInThisDimension;
	  intervals.push_back(intInThisDimension);
    }


    bool* pairedCollumns = new bool[numberOfCells];
    for ( size_t i = 0 ; i != numberOfCells ; ++i )
    {
	  pairedCollumns[i] = false;
    }
    for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
    {
        if (computePersistenceOfAFunctionDBG )std::cout << "Birth: " << pairs.get_pair( idx ).first << ", Death: " << pairs.get_pair( idx ).second << std::endl;
            if (computePersistenceOfAFunctionDBG )std::cout << "Birth: " << pairs.get_pair( idx ).first << ", Death: " << pairs.get_pair( idx ).second << std::endl;
        cell<T>* first = numberOfGeneratorToDim[pairs.get_pair( idx ).first];
        cell<T>* second = numberOfGeneratorToDim[pairs.get_pair( idx ).second];
	 pairedCollumns[pairs.get_pair( idx ).first] = true;
	 pairedCollumns[pairs.get_pair( idx ).second] = true;

	 if ( first->val() != second->val() )
	 {
		if (computePersistenceOfAFunctionDBG )cerr << "Pair, if : " << *first << "of a dimension : " << first->dim() << " and value : " << first->val() << "\n and \n" << *second << " of a dimension : " << second->dim() << " and value : " << second->val()  << endl;
		if ( first->val() < second->val() )
		{
			if (!showIntervalsShorterThanError)
			{
				if ( second->val()-first->val() > epsilon )
				{
					//cerr << "second->val()-first->val() : " << second->val()-first->val() << endl;
	 				intervals[ first->dim() ].push_back( std::make_pair( first->val() , second->val() ) );
				}
			}
			else
			{
				intervals[ first->dim() ].push_back( std::make_pair( first->val() , second->val() ) );
			}
		}
		else
		{
			if (!showIntervalsShorterThanError)
			{
				if ( first->val()-second->val() > epsilon )
				{
					cerr << "first->val()-second->val() : " << first->val()-second->val() << endl;
	 				intervals[ first->dim() ].push_back( std::make_pair( second->val() , first->val() ) );
				}
			}
			else
			{
				intervals[ first->dim() ].push_back( std::make_pair( second->val() , first->val() ) );
			}
		}
	 }
    }

    if (computePersistenceOfAFunctionDBG )
    {
    	cerr << "boundary_matrix.get_num_cols() : " << boundary_matrix.get_num_cols() << endl;
    	cerr << "numberOfCells  : " << numberOfCells  << endl;
    }


    for ( size_t i = 0 ; i != numberOfCells ; ++i )
    {
	  if ( !pairedCollumns[i] )
	  {
              if (computePersistenceOfAFunctionDBG )
		{
			cerr << "The collumn : " << i << " is empty and it is not an element of any pair!! \n";
			cerr << *numberOfGeneratorToDim[i] << endl;
		}
		//TODO -- set infinity into some better way!
		intervals[ numberOfGeneratorToDim[i]->dim() ].push_back( std::make_pair( numberOfGeneratorToDim[i]->val() , 1000000 ) );
	  }
    }



    for ( size_t dim = 0 ; dim != cmplx->elemen().size() ; ++dim )
    {
	  if (computePersistenceOfAFunctionDBG )cerr << "intervals["<<dim<<"].size() : " << intervals[dim].size() << endl;
    }





    cxsc::real level = 0.1;

    for ( size_t dim = 0 ; dim != cmplx->elemen().size() ; ++dim )
    {
          for ( size_t nr  = 0 ; nr != cmplx->elemen()[dim].size() ; ++nr )
	   {
		if ( cmplx->elemen()[dim][nr]->val() < level  )
		{
//cerr << "Deleteing cell in dimension : " << dim << endl;
			cmplx->elemen()[dim][nr]->del();
		}
	   }
    }





    return intervals;
}//computeNodalDomainOfAFunction
