// cell.h

/**************************************************************/
#include<iostream>
#include<list>
#include<vector>
#include<set>
#include <cassert>
#include <climits>
#include <iomanip>
#include "interval.hpp"  // predefined interval arithmetic
#include "hess_ari.hpp"


template <typename T>
class cell
{
public:
    //constructors:
    cell();
    cell( const cell& original);
    cell( typename std::vector< std::pair<T,T> > coef)
    {
         for ( unsigned int i = 0 ; i != coef.size() ; ++i )
         {
             this->coef.push_back( std::make_pair(coef[i].first , coef[i].second) );
         }
         this->delet = true;
         this->subdivisionDepth = 0;
         this->verified = false;
         ++this->numberOfAlCells;
         this->numberOfCell = this->numberOfAlCells;
	  this->value = INT_MAX;
    }

    //this procedure checkes if non-deleted part of the boundary of cube is conractible to deleted part of the boundary
    bool isNonDeletedPartOfBoundaryContractibleToDeletedOne();


    ivector getivectorRepresentationOfCell() {

      cxsc::ivector x(this->coef.size());
      for (size_t i=1; i != this->coef.size()+1; ++i) {
        x[i] = cxsc::interval(this->coef[i-1].first,this->coef[i-1].second);
      }
      return x;

    }
    HTvector getHTvectorRepresentationOfCell() {

      cxsc::ivector x(this->coef.size());
      for (size_t i=1; i != this->coef.size()+1; ++i) {
        x[i] = cxsc::interval(this->coef[i-1].first,this->coef[i-1].second);
      }
      return HessVar(x);

    }




    /** Return integer dimension of cell.
     *
     * \code
      For a two-dimensional cell: (0,1) *----* (1,1)
                                        |    |
                                        |    |
                                        |    |
                                  (0,0) *----* (1,0)

      cell->dim() would be 2.
      \endcode
     *
     * \section ex1 Example
     * \snippet /Users/tds3/ranval3d/code/tests/testCellComplex.cpp  dim_test
     */
    int dim()const
    {

        //counting the number of non degenerated intervals in coef:
        assert ( this->coef.size() != 0 );
        int dim = 0;
        for ( unsigned int i = 0 ; i != this->coef.size() ; ++i )
        {
            if ( this->coef[i].first != this->coef[i].second )
            {
                 ++dim;
            }
        }
        return dim;

    }

    //TODO -- check if using vectors here at the end would not be a better idea as for complexity (i.e. profile the code)
    std::list< cell* >& bd() {return this->bound;};
    std::list< cell* > bd() const {return this->bound;};
    std::list< cell* >& cbd(){return this->coBound;};
    std::list< cell* > cbd()const{return this->coBound;};

    bool& isVerified(){return this->verified;}

    typedef typename std::list< cell* >::iterator BdIterator;
    typedef typename std::list< cell* >::iterator CbdIterator;

    const BdIterator bdBegin(){return this->bound.begin();};
    const BdIterator bdEnd(){return this->bound.end();};
    const CbdIterator cbdBegin(){return this->coBound.begin();};
    const CbdIterator cbdEnd(){return this->coBound.end();};


    unsigned int bdSize();
    unsigned int cbdSize();
    inline void del() {
      this->delet = true; 
      this->verified = true;
    }

    inline void undel() {
      this->delet = false; 
      this->verified = true;
    }

    void undelWithAllBdElem();
    void delWithAllBdElem();

    inline bool& deleted(){return this->delet;};

    template < typename A >
    friend std::ostream& operator<<(std::ostream& out, cell<A>& sim)
    {
         out << std::fixed;
         //out << "[ "  ;
         for ( unsigned int i = 0 ; i != sim.coef.size() ; ++i )
         {
              out << std::setprecision (30) << "[" << sim.coef[i].first << "," << sim.coef[i].second << "]";
              if ( i != sim.coef.size()-1 ) out << "x";
         }
         //out << " ] validat : " << sim.verified << " , deleted : " << sim.delet << endl;
         //out << " ]";
         return out;
    };

    template < typename A >
    inline friend bool operator == ( const cell<A>& t1, const cell<A>& t2 )
    {
           if ( t1.coef.size() == t2.coef.size() )
           {
                for ( unsigned int i = 0 ; i != t1.coef.size() ; ++i )
                {
                    if ( t1.coef[i].first != t2.coef[i].first )
                    {
                         return false;
                    }
                    if ( t1.coef[i].second != t2.coef[i].second )
                    {
                         return false;
                    }
                }
                return true;
           }
           return false;
    }

    template <typename A>
    friend class cellComplex;

    template <typename A>
    friend int computeIncidence( cell<A>* first , cell<A>* second );

    template <typename A>
    friend bool checkIfIsASubset( cell<A>* coface , cell<A>* face );

    typename std::vector< cell<T>* > computeWholeBoundaryOfCube();


    /** Return vector containing coordinates describing the extent of
     * the cell.
     *
     * \code
      For a two-dimensional cell: (0,2) *----* (1,2)
                                        |    |
                                        |    |
                                        |    |
                                  (0,0) *----* (1,0)

      cell->coords() would be [ [0,1]x[0,2] ].
      \endcode
     *
     * \section ex1 Example
     * \snippet /Users/tds3/ranval3d/code/tests/testCellComplex.cpp  coords_test
     */
    typename std::vector< std::pair<T,T> > coords(){ return this->coef; }


    size_t& position(){return this->pos;}

    int subdivDepth()const{return this->subdivisionDepth;}

    unsigned& number(){return this->numberOfCell;}
    unsigned numberInCmplx()const{return this->numberOfCell;}
    cxsc::real& val(){return value;}

protected:
    std::list<cell*> bound;
    std::list<cell*> coBound;
    std::vector< std::pair<T,T> > coef;
    bool delet;
    size_t pos;
    int subdivisionDepth;
    bool verified;
    static unsigned numberOfAlCells;
    unsigned numberOfCell;
    cxsc::real value; // we need it to compute persistence induced by the map.
};//simplex

template <typename T>
unsigned cell<T>::numberOfAlCells = 0;

template <typename T>
class cellComparison
{
public:
  cellComparison(){}
  bool operator() (const cell<T>* lhs , const cell<T>* rhs) const
  {
     if ( lhs->dim() < rhs->dim() )return false;  // cell of greatest dimension is on top of priority_queue
     if ( lhs->dim() > rhs->dim() )return true;
     //if ( lhs->dim() < rhs->dim() )return true;     // cell of least dimension is on top of priority_queue
     //if ( lhs->dim() > rhs->dim() )return false;
     return ( lhs->numberInCmplx() < rhs->numberInCmplx() );
  }
};

template <typename T>
bool compareCells(const cell<T>* lhs , const cell<T>* rhs)
{
    cellComparison<T> compare;
    return compare(lhs,rhs);
}




template <typename A>
bool checkIfIsASubset( cell<A>* coface , cell<A>* face )
{
     if ( coface->coef.size() != face->coef.size() )
     {
          return false;
     }
     for ( unsigned int i = 0 ; i != coface->coef.size() ; ++i )
     {
         if ( !(
             (coface->coef[i].first <= face->coef[i].first)
             &&
             (coface->coef[i].second >= face->coef[i].second)
             )
            )
            {
                return false;
            }
     }
     return true;
}//checkIfIsASubset

template <typename A>
int computeIncidence( cell<A>* coface , cell<A>* face )
{
    if ( checkIfIsASubset(coface, face)==false ){return 0;}

    int sumOfDim = 0;
    unsigned int i = 0;
    while ( i != coface->coef.size() )
    {
        if ( (coface->coef[i].first != coface->coef[i].second)
             &&
             ( face->coef[i].first != face->coef[i].second )
           )
        {
             ++sumOfDim;
        }
        else
        {
            //std::cerr << "Roznia sie na wymairze : " << i << "\n";

            //if coface is non degenerated at this coord:
            if ( coface->coef[i].first != coface->coef[i].second )
            {
                 int minusOneToPowSum = 1;
                 if ( sumOfDim%2 == 1 ){minusOneToPowSum = -1;}
                 //in this case we are sure, that face is degenerated at this cord
                 if ( coface->coef[i].first == face->coef[i].first )
                 {
                      return (-1)*minusOneToPowSum;
                 }
                 if ( coface->coef[i].second == face->coef[i].first )
                 {
                      return minusOneToPowSum;
                 }
            }
        }
        ++i;
    }
    return 0;
}
