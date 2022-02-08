// cell.hpp


#include"cell.h"
#include <algorithm>
#include <climits>
using namespace std;


/** Return boundary size.
 */
template <typename T>
inline unsigned int cell<T>::bdSize()
{
     unsigned int num = 0;
     BdIterator holder;
     for (BdIterator it = this->bound.begin();
     it!=this->bound.end();++it)
     {
         if (!(*it)->delet){num++;}
     }
     return num;
}//bdSize

/** Return coboundary size.
 */
template <typename T>
inline unsigned int cell<T>::cbdSize()
{
     unsigned int num = 0;
     CbdIterator holder;
     for (CbdIterator it = this->coBound.begin();
     it!=this->coBound.end();++it)
     {
           if (!(*it)->delet){num++;}
     }
     return num;
}//cbdSize


/** Default constructor.
 * */
template <typename T>
cell<T>::cell()
{
     this->delet = true;
     this->subdivisionDepth = 0;
     this->verified = false;
     ++this->numberOfAlCells;
     this->numberOfCell = this->numberOfAlCells;
     this->value = INT_MAX;
}

/** Full constructor.
 * */
template <typename T>
cell<T>::cell( const cell& original)
{
     typename std::list< cell<T>* >::iterator it;
     for ( it = original.bd().begin() ; it != original.bd().end() ; ++it )
     {
         this->bound.push_back( (*it) );
     }

     for ( it = original.cbd().begin() ; it != original.cbd().end() ; ++it )
     {
         this->coBound.push_back( (*it) );
     }
     this->delet = original.delet;
     this->toRemove = original.toRemove;
     this->id = original.id;
     this->subdivisionDepth = original.subdivisionDepth;
     this->verified = original.verified;
     ++this->numberOfAlCells;
     this->numberOfCell = this->numberOfAlCells;
     this->vale = original.value;
}//copy constructor

bool isNonDeletedPartOfBoundaryContractibleToDeletedOneDBG = false;
template <typename T>
bool cell<T>::isNonDeletedPartOfBoundaryContractibleToDeletedOne()
{

    std::vector< cell<T>* > deletedElementsInCube;
    std::set< cell<T>* > delElms;
    std::vector< cell<T>* > elems;
    elems.push_back(this);
    delElms.insert(this);
    deletedElementsInCube.push_back(this);
    while ( !elems.empty() )
    {
        std::vector< cell<T>* > NewElems;
        for ( size_t i = 0 ; i != elems.size() ; ++i )
        {
            if ( elems[i]->deleted() )
            {
                deletedElementsInCube.push_back(elems[i]);
                delElms.insert(elems[i]);
            }
            for ( typename cell<T>::BdIterator bd = elems[i]->bdBegin() ; bd != elems[i]->bdEnd() ; ++bd )
            {
                NewElems.push_back( *bd );
            }
        }
        std::sort(NewElems.begin(), NewElems.end());
        NewElems.erase( std::unique(NewElems.begin(), NewElems.end()),  NewElems.end());
        elems.swap(NewElems);
    }

    while ( true )
    {
        typename std::set< cell<T>* >::iterator it = delElms.begin();
        cell<T>* freeFace = 0;
        while ( it != delElms.end() )
        {
            //check if *it is a free coface:
            int numberOfDeletedElementsInBoundary = 0;
            for ( cell<T>::BdIterator bd = (*it)->bdBegin() ; bd != (*it)->bdEnd() ; ++bd )
            {
                if ( (*bd)->deleted() )
                {
                    freeFace = *bd;
                    ++numberOfDeletedElementsInBoundary;
                }
            }
            if ( numberOfDeletedElementsInBoundary == 1 )
            {
                break;
            }
            else
            {
                ++it;
            }
        }
        if ( it != delElms.end() )
        {
            //cerr << "Collpase of : " << **it << " and : " << *freeFace << endl;
            (*it)->undel();
            freeFace->undel();
            delElms.erase(it);
            delElms.erase( delElms.find(freeFace) );
        }
        else
        {
            break;
        }
    }

    //puting back all previously deleted elements as deleted
    for ( size_t i = 0 ; i != deletedElementsInCube.size() ; ++i )
    {
        deletedElementsInCube[i]->del();
    }

    return ( delElms.size() == 0 );

    /*
    std::set< cell<T>* , cellComparison<T> > nonDeletedElementsInCube;
    std::vector< cell<T>* > elems;
    elems.push_back( this );
    while ( !elems.empty() )
    {
        std::vector< cell<T>* > NewElems;

        if ( isNonDeletedPartOfBoundaryContractibleToDeletedOneDBG )
        {
            cerr << "elems: " << endl;
            for ( size_t i = 0 ; i != elems.size() ; ++i )
            {
                cerr << *elems[i] << " " << elems[i]->deleted() << endl;
            }
            getchar();
        }

        for ( size_t i = 0 ; i != elems.size() ; ++i )
        {
            if ( elems[i]->deleted() ){nonDeletedElementsInCube.insert(elems[i]);}
            for ( typename cell<T>::BdIterator bd = elems[i]->bdBegin() ; bd != elems[i]->bdEnd() ; ++bd )
            {
                NewElems.push_back( *bd );
            }
        }
        std::sort(NewElems.begin(), NewElems.end());
        NewElems.erase( std::unique(NewElems.begin(), NewElems.end()),  NewElems.end());
        elems.swap(NewElems);
    }

    if (isNonDeletedPartOfBoundaryContractibleToDeletedOneDBG ){std::cerr << "nonDeletedElementsInCube.size() : " << nonDeletedElementsInCube.size() << "\n";}

    std::vector<cell<T>*> elementsDeletedInTheProcess;

    //There are a few possible strategy to check that. I will use coreductions which is the simples one here:
    //This code can be optymized for higher dimensions:
    while (true)
    {
       size_t initialSize = elementsDeletedInTheProcess.size();

	//for every non deleted element
	for ( typename std::set< cell<T>* , cellComparison<T> >::iterator it = nonDeletedElementsInCube.begin() ; it != nonDeletedElementsInCube.end() ; ++it )
	{
		//if for (*it) there exist unique nonDeleted element in boundary of (*it)
		int numbeOfNonDelInBd = 0;
		cell<T>* nonDel = 0;
		for ( typename cell<T>::BdIterator bd = (*it)->bdBegin() ; bd != (*it)->bdEnd() ; ++bd )
		{
			if ( (*bd)->deleted() == false )
			{
				++numbeOfNonDelInBd ;
				nonDel = *bd;
			}
		}
		if ( numbeOfNonDelInBd  == 1 )
		{
			elementsDeletedInTheProcess.push_back(*it);
			elementsDeletedInTheProcess.push_back(nonDel);
			(*it)->del();
			nonDel->del();
			if ( isNonDeletedPartOfBoundaryContractibleToDeletedOneDBG  )
			{
				cerr << "Removing pair : " << **it << " and " << *nonDel << endl;
			}
			if ( (*it)==this )
            {
                for ( size_t dir = 0 ; dir != this->coef.size() ; ++dir )
                {
                    if ( this->coef[dir] != nonDel->coef[dir] )
                    {
                        direction = dir;
                        break;
                    }
                }
            }
		}
	}
	std::vector< typename std::set<cell<T>* , cellComparison<T> >::iterator > elementsToDelete;
	for ( typename std::set< cell<T>* , cellComparison<T> >::iterator it = nonDeletedElementsInCube.begin() ; it != nonDeletedElementsInCube.end() ; ++it )
	{
		if ( (*it)->deleted() )elementsToDelete.push_back( it );
	}
	for ( size_t i = 0 ; i != elementsToDelete.size() ; ++i )
	{
		nonDeletedElementsInCube.erase( elementsToDelete[i] );
	}

	size_t finalSize = elementsDeletedInTheProcess.size();
	if ( initialSize  == finalSize  )break;
    }

    for ( size_t i = 0 ; i != elementsDeletedInTheProcess.size() ; ++i )
    {
	  elementsDeletedInTheProcess[i]->deleted() = false;
    }

    bool isContractible = ( nonDeletedElementsInCube.size() == 0 );
    return isContractible;*/
}


template <typename T>
void cell<T>::undelWithAllBdElem()
{
    std::vector< cell<T>* > elems;
    elems.push_back( this );
    while ( !elems.empty() )
    {
       std::vector< cell<T>* > NewElems;
	for ( size_t i = 0 ; i != elems.size() ; ++i )
	{
	    elems[i]->undel();
	    for ( typename cell<T>::BdIterator bd = elems[i]->bdBegin() ; bd != elems[i]->bdEnd() ; ++bd )
	    {
		  NewElems.push_back( *bd );
	    }
	}
	std::sort(NewElems.begin(), NewElems.end());
	NewElems.erase( std::unique(NewElems.begin(), NewElems.end()),  NewElems.end());
	elems.swap(NewElems);
    }
}

template <typename T>
void cell<T>::delWithAllBdElem()
{
    std::vector< cell<T>* > elems;
    elems.push_back( this );
    while ( !elems.empty() )
    {
       std::vector< cell<T>* > NewElems;
	for ( size_t i = 0 ; i != elems.size() ; ++i )
	{
	    elems[i]->del();
	    for ( typename cell<T>::BdIterator bd = elems[i]->bdBegin() ; bd != elems[i]->bdEnd() ; ++bd )
	    {
		  NewElems.push_back( *bd );
	    }
	}
	std::sort(NewElems.begin(), NewElems.end());
	NewElems.erase( std::unique(NewElems.begin(), NewElems.end()),  NewElems.end());
	elems.swap(NewElems);
    }
}



template <typename T>
std::vector< cell<T>* > cell<T>::computeWholeBoundaryOfCube()
{
    std::vector< cell<T>* > result;
    std::vector< cell<T>* > boundaryAtThisLevel;
    boundaryAtThisLevel.push_back(this);
    while ( !boundaryAtThisLevel.empty() )
    {
        std::vector< cell<T>* > newBoundaryAtThisLevel;
        for ( size_t i = 0 ; i != boundaryAtThisLevel.size() ; ++i )
        {
            for ( typename cell<T>::BdIterator bd = boundaryAtThisLevel[i]->bdBegin() ; bd != boundaryAtThisLevel[i]->bdEnd() ; ++bd )
            {
                newBoundaryAtThisLevel.push_back(*bd);
            }
        }
        std::sort( newBoundaryAtThisLevel.begin() , newBoundaryAtThisLevel.end() );
        newBoundaryAtThisLevel.erase(std::unique(newBoundaryAtThisLevel.begin(), newBoundaryAtThisLevel.end()), newBoundaryAtThisLevel.end());
        boundaryAtThisLevel = newBoundaryAtThisLevel;
        result.insert( result.end() , newBoundaryAtThisLevel.begin() , newBoundaryAtThisLevel.end() );
    }
    return result;
}
