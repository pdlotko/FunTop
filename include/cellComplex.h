// cellComplex.h

#include<iostream>
#include<list>
#include<vector>

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



enum subdivisionType {dyadic, randm, goldenRatio};

template <typename T>
class cellComplex
{
public:
       cellComplex();

       cellComplex( std::vector< std::pair<T,T> > cube);

       void writeToFile(char* filename);


       typename std::vector< cell<T>* > divideCube( cell<T>* cube , T x , int cord , typename std::vector< cell<T>* >& vec , bool addDepth = true );


       typename std::vector< cell<T>* > divideCubeAlongLongestAxis( cell<T>* cube , typename std::vector< cell<T>* >& vec )
       {
          int dominantDirection = -1;
          double longestSide    = -1;

          for ( size_t side = 0 ; side != cube->coords().size() ; ++side )
          {
            if ( cube->coords()[side].second - cube->coords()[side].first > longestSide )
            {
              longestSide = cube->coords()[side].second - cube->coords()[side].first;
              dominantDirection = side;
            }
          }
          typename std::vector< cell<T>* > result;
          //result = divideCube( cube , dominantDirection , vec );
          divideCube( cube , dominantDirection , vec );
          return vec;
       }
       

       typename std::vector< cell<T>* > divideCube( cell<T>* cube , int cord , typename std::vector< cell<T>* >& vec )
       {


           //in case of random subdivision the values are taken from the following
           double beginInterval = 0.3;
           double endInterval = 0.7;

           typename std::vector< cell<T>* > result;
           extern subdivisionType usedSubdivisionType;
           if ( usedSubdivisionType == dyadic )
           {
                result = divideCube( cube , 0.5*(cube->coords()[cord].first+cube->coords()[cord].second) , cord , vec , true );
           }
           if ( usedSubdivisionType == randm )
           {
                double point = (rand()/(double)RAND_MAX)*(endInterval - beginInterval) + beginInterval;
                point *= (cube->coords()[cord].second-cube->coords()[cord].first);
                point += cube->coords()[cord].first;
                result = divideCube( cube , point , cord , vec , true );
           }
           if ( usedSubdivisionType == goldenRatio )
           {
                double goldenRatioConst = 0;
                int nr = rand()%2;

                if ( rand()%2 == 0 )
                {
                    goldenRatioConst = 0.618034;
                }
                else
                {
                    goldenRatioConst = 0.381966;
                }

                result = divideCube( cube , (goldenRatioConst*(cube->coords()[cord].second-cube->coords()[cord].first)+cube->coords()[cord].first) , cord , vec , true );
           }
           return result;
       }

       typename std::vector< cell<T>* > divideCubeInAllDirectionsItHasNonzeroLenght( cell<T>* cube );

       ~cellComplex()
       {
            for ( unsigned int i = 0 ; i != this->elements.size() ; ++i )
            {
                typename std::vector< cell<T>* >::iterator it;
                for ( it = this->elements[i].begin() ; it != this->elements[i].end() ; ++it )
                {
                    delete (*it);
                }
            }
       }

       std::vector< unsigned > computeBettiNumbersOfNonDeletedPart();


        template <class K>
        friend std::ostream& operator<<(std::ostream& out, cellComplex<K>& cmplx)
        {
             bool displayInformationAboutBoundaryAndCoboundary = false;
             for ( unsigned int i = 0 ; i < cmplx.elements.size() ; ++i )
             {
                 out << "\n dimension : " << i <<"\n";
                 for ( unsigned nr = 0 ; nr != cmplx.elements[i].size() ; ++nr  )
                 {
                     out << "\n" ;
                     out << "element: "<< *cmplx.elements[i][nr] << " , " << cmplx.elements[i][nr]->subdivDepth();
                     out <<"\n";
                     if ( displayInformationAboutBoundaryAndCoboundary )
                     {
                         typename cell<K>::BdIterator bd;
                         out << "Boundary : \n";
                         for ( bd = cmplx.elements[i][nr]->bdBegin() ; bd != cmplx.elements[i][nr]->bdEnd() ; ++bd )
                         {
                             out << **bd << " , ";
                         }
                         out << "\n Coboundary: \n";
                         typename cell<K>::CbdIterator cbd;
                         for ( cbd = cmplx.elements[i][nr]->cbdBegin() ; cbd != cmplx.elements[i][nr]->cbdEnd() ; ++cbd )
                         {
                             out << **cbd << " , ";
                         }
                         out << "\n";
                     }
                 }
             }
             return out;
        }

        /** Return this.elements, a std::vector containing (vectors of)
         * cells of this cellComplex.
         * Note: elements.size() is the dimension of this cellComplex,
         * and the components of this.elements are vectors of cells of the
         * corresponding dimension.
         * */
        typename std::vector< std::vector< cell<T>* > >& elemen(){return this->elements;}



        //this procedure is used to divide the cube 'cub' to n*n cubes having the same size
        typename std::vector< cell<T>* > divideTheCubeByGivenFactor( cell<T>* cub , int n );

        int dim(){return this->elements.size();};

       void checkComplex()
       {
           for ( size_t dim = 0  ; dim != this->elements.size() ; ++dim )
           {
               for ( size_t nr = 0 ; nr != this->elements[dim].size() ; ++nr )
               {
                   if ( this->elements[dim][nr]->pos != nr )
                   {
                       cerr << "Nonconsistent enumeration in the dimension : " << dim << " at a complex : " << *this->elements[dim][nr] << endl;
                       cerr << "Should be : " << nr << " while there is " << this->elements[dim][nr]->pos  << endl;
                       cin.ignore();
                   }
               }
           }
       }

       //this method is used in imposig periodic boundary conditions to check if 'element' should be mapped to some other element or not
       bool isAtLastOneCoordinateMaximal( cell<T>* element )
       {
           for ( size_t i = 0 ; i != element->coef.size() ; ++i )
           {
                if ( element->coef[i].first != element->coef[i].second )continue;
                //in this case we know that (element->coef[i].first == element->coef.second). We should check if this is a maximal one
                if ( element->coef[i].first == this->initialCube[i].second )
                {
                    return true;
                }
           }
           return false;
       }

       cell<T>* findPairedCell( cell<T>* element , std::vector< cell<T>* >& bCellsInThisDimension );


private:
       void tie( cell<T>* coface, cell<T>* face );

       void substitute( cell<T>* source , cell<T>* target )
       {
             typename std::list< cell<T>* >::iterator it;
             for ( it = source->bd().begin() ; it != source->bd().end() ; ++it )
             {
                 target->bound.push_back( (*it) );
             }

             for ( it = source->cbd().begin() ; it != source->cbd().end() ; ++it )
             {
                 target->coBound.push_back( (*it) );
             }
             target->delet = source->delet;

             delete source;
       }

protected:
       //TODO : checkif using set here is optimal. Probably after constructing it it would bebetter to use vector or so....
       std::vector< std::pair<T,T> > initialCube;
       typename std::vector< std::vector< cell<T>* > > elements; //!< elements has as many components as this cellComplex has dimensions.  Each component is a vector of cells making up that dimension's structure.
};//simplicalComplex
