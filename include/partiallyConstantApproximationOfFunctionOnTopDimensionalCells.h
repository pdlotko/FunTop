#include <iostream>
#include <vector>
#include <queue>

using namespace std;

class topDimensionalCellDecompositionOfDomain
{
public:
//scalar valued function constructors
    topDimensionalCellDecompositionOfDomain( std::vector< std::vector< std::pair<double,double> > > initialRectangles , HessType (*f)( const HTvector& x ) , double epsilon  );
    topDimensionalCellDecompositionOfDomain( std::vector< std::pair<double,double> > initialRectangle , HessType (*f)( const HTvector& x ) , double epsilon  );

//vector valued function constructor
    topDimensionalCellDecompositionOfDomain( std::vector< std::pair<double,double> > initialRectangle , std::vector< HessType (*)( const HTvector& x ) > functionsForFollowingCoordinates, double epsilon  );

//writing to a file subroutine (what else should be in this class?)
    void writeResultToFile( char* filename );
    
    size_t number_of_validated_cells(){return this->validatedCells.size();}
    size_t number_of_not_validated_cells(){return this->notValidatedCells.size();}
private:
    void construct( std::vector< std::vector< std::pair<double,double> > > initialRectangles , HessType (*f)( const HTvector& x ) , double epsilon  );
    std::vector< std::pair<double,double> > initialRectangle;
    std::vector< std::pair< std::vector< std::pair<double,double> > , std::vector<cxsc::real> > > validatedCells;
    std::queue< std::vector< std::pair<double,double> > > notValidatedCells;
};



interval computeValueOfFunctionOnRectangle( std::vector< std::pair<double,double> > topRectangle , HessType (*f)( const HTvector& x ) )
{
    HessType fx;
    interval f_x;

    cxsc::ivector x( topRectangle.size() );
    for (size_t i=1; i != topRectangle.size()+1; ++i)
    {
        x[i] = cxsc::interval(topRectangle[i-1].first,topRectangle[i-1].second);
    }
    HTvector xx = HessVar(x);


    fx = f(xx);
    f_x = fValue(fx);
    return f_x;
}


//vector valued function constructor
topDimensionalCellDecompositionOfDomain::topDimensionalCellDecompositionOfDomain( std::vector< std::pair<double,double> > initialRectangle , std::vector< HessType (*)( const HTvector& x ) > functionsForFollowingCoordinates , double epsilon  )
{
    std::vector< std::vector< std::pair<double,double> > > initialRectangles;
    initialRectangles.push_back( initialRectangle );

    for ( size_t dim = 0 ; dim != functionsForFollowingCoordinates.size() ; ++dim )
    {
        if ( !this->validatedCells.empty() )
        {
            //we enter here only if we have already been in this loop.
            initialRectangles.clear();
            for ( size_t i = 0 ; i != this->validatedCells.size() ; ++i )
            {
                initialRectangles.push_back( this->validatedCells[i].first );
            }
            this->validatedCells.clear();
        }

        this->construct( initialRectangles , functionsForFollowingCoordinates[dim] , epsilon );
    }

    //now, once we are done here, we need to compute the value for each validated rectangle of each function:
    for ( size_t i = 0 ; i != this->validatedCells.size() ; ++i )
    {
        std::vector<cxsc::real> values;
        for ( size_t dim = 0 ; dim != functionsForFollowingCoordinates.size() ; ++dim )
        {
            interval f_x = computeValueOfFunctionOnRectangle( this->validatedCells[i].first , functionsForFollowingCoordinates[dim] );
            values.push_back( 0.5*(Sup(f_x)+Inf(f_x)) );
        }
        this->validatedCells[i].second = values;
    }
}


//scalar valued function constructors
topDimensionalCellDecompositionOfDomain::topDimensionalCellDecompositionOfDomain( std::vector< std::vector< std::pair<double,double> > > initialRectangles , HessType (*f)( const HTvector& x ) , double epsilon  )
{
    this->construct( initialRectangles , f , epsilon );
}


topDimensionalCellDecompositionOfDomain::topDimensionalCellDecompositionOfDomain( std::vector< std::pair<double,double> > initialRectangle , HessType (*f)( const HTvector& x ) , double epsilon  )
{
    std::vector< std::vector< std::pair<double,double> > > initialRectangles;
    initialRectangles.push_back( initialRectangle );
    this->construct( initialRectangles , f , epsilon );
}

void topDimensionalCellDecompositionOfDomain::construct( std::vector< std::vector< std::pair<double,double> > > initialRectangles , HessType (*f)( const HTvector& x ) , double epsilon  )
{
    bool dbg = false;
    for ( size_t i = 0 ; i != initialRectangles.size() ; ++i )
    {
        this->notValidatedCells.push( initialRectangles[i] );
    }

    while ( !this->notValidatedCells.empty() )
    {
        std::vector< std::pair<double,double> > topRectangle = this->notValidatedCells.front();
        this->notValidatedCells.pop();

        if ( dbg )
        {
            cerr << "Considering non validated cell : ";
            for ( size_t dim = 0 ; dim != topRectangle.size() ; ++dim )
            {
                cerr << "[ " << topRectangle[dim].first << " , " << topRectangle[dim].second << "] ";
            }
            cerr << endl;
        }

        //now trying to validate topRectangle:
         interval f_x = computeValueOfFunctionOnRectangle( topRectangle , f );

         if ( dbg )
         {
             cerr << "Diameter of the result : " << diam(f_x) << endl;
         }

         if ( diam(f_x) <= 2*epsilon )
         {
             if ( dbg ){cerr << "This cell is validated!\n";}
             std::vector< cxsc::real > value;
             value.push_back( 0.5*(Sup(f_x)+Inf(f_x)) );
             this->validatedCells.push_back( std::make_pair( topRectangle , value ) );
         }
         else
         {
             if ( dbg ){cerr << "This cell will be subdivided!\n";}
             //subdivide this rectangle in all directions:
             std::vector< std::vector< std::pair<double,double> > > toSubdivide;
             toSubdivide.push_back( topRectangle );
             for ( size_t direction = 0 ; direction != topRectangle.size() ; ++direction )
             {
                 if ( dbg )
                 {
                     cerr << "Subdividing all cells in direction : " << direction << endl;

                     cerr << "Here are the cells we are to subdivide : \n";
                     for ( size_t j = 0 ; j != toSubdivide.size() ; ++j )
                     {
                          for ( size_t dim = 0 ; dim != toSubdivide[j].size() ; ++dim )
                          {
                              cerr << "[ " << toSubdivide[j][dim].first << " , " << toSubdivide[j][dim].second << " ]";
                          }
                          cerr << endl;
                     }
                     cerr << endl;
                 }
                 std::vector< std::vector< std::pair<double,double> > > toSubdivideNew;
                 //subdivide every rectangle from topRectangle into direction 'direction' if only it has nonzero thickness there:
                 for ( size_t i = 0 ; i != toSubdivide.size() ; ++i )
                 {
                     if ( toSubdivide[i][direction].first != toSubdivide[i][direction].second )
                     {
                         std::vector< std::pair<double,double> > first( toSubdivide[i].size() );
                         std::vector< std::pair<double,double> > second( toSubdivide[i].size() );
                         for ( size_t localDimension = 0 ; localDimension != toSubdivide[i].size() ; ++localDimension )
                         {
                             if ( localDimension != direction )
                             {
                                 first[localDimension] = std::make_pair( toSubdivide[i][localDimension].first , toSubdivide[i][localDimension].second );
                                 second[localDimension] = std::make_pair( toSubdivide[i][localDimension].first , toSubdivide[i][localDimension].second );
                             }
                             else
                             {
                                 first[localDimension] = std::make_pair( toSubdivide[i][localDimension].first , 0.5*(toSubdivide[i][localDimension].first+toSubdivide[i][localDimension].second) );
                                 second[localDimension] = std::make_pair( 0.5*(toSubdivide[i][localDimension].first+toSubdivide[i][localDimension].second) , toSubdivide[i][localDimension].second );
                             }
                         }
                         toSubdivideNew.push_back( first );
                         toSubdivideNew.push_back( second );
                     }
                     else
                     {
                         toSubdivideNew.push_back( toSubdivide[i] );
                     }
                 }
                 if ( dbg )
                 {
                     cerr << "Elements after subdivision : \n";
                     for ( size_t j = 0 ; j != toSubdivideNew.size() ; ++j )
                     {
                          for ( size_t dim = 0 ; dim != toSubdivideNew[j].size() ; ++dim )
                          {
                              cerr << "[ " << toSubdivideNew[j][dim].first << " , " << toSubdivideNew[j][dim].second << " ]";
                          }
                          cerr << endl;
                     }
                     getchar();
                 }
                 toSubdivide = toSubdivideNew;
             }

            if ( dbg )
            {
                cerr << "Here are the elements after subdivision in this direction: \n";
            }
             for ( size_t i = 0 ; i != toSubdivide.size() ; ++i )
             {
                 this->notValidatedCells.push( toSubdivide[i] );
                 if ( dbg )
                 {
                     for ( size_t dim = 0 ; dim != toSubdivide[i].size() ; ++dim )
                     {
                         cerr << "[ " << toSubdivide[i][dim].first << " , " << toSubdivide[i][dim].second << " ]";
                     }
                     cerr << endl;
                 }
             }

         }
         if ( dbg )std::cin.ignore();
    }
}

void topDimensionalCellDecompositionOfDomain::writeResultToFile( char* filename )
{
    ofstream out;
    out.open( filename );

    for ( size_t i = 0 ; i != this->validatedCells.size() ; ++i )
    {
        //std::vector< std::pair< std::vector< std::pair<double,double> > > > validatedCells
        for ( size_t localDim = 0 ; localDim != this->validatedCells[i].first.size() ; ++localDim )
        {
            out << "[" << this->validatedCells[i].first[localDim].first << "," << this->validatedCells[i].first[localDim].second << "]";
            if ( localDim != this->validatedCells[i].first.size()-1 )out << " x ";
        }
        out << " -> [";
        for ( size_t numberOfFunction = 0 ; numberOfFunction != this->validatedCells[i].second.size() ; ++numberOfFunction )
        {
            out << this->validatedCells[i].second[numberOfFunction] << " ";
        }
        out << "] " << endl;
    }
    out.close();
}
