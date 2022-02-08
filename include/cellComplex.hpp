// cellComplex.hpp

#include"cellComplex.h"
#include<algorithm>
#include<map>
#include<assert.h>
#include<stdexcept>



//TODO -- now we are computing positive nodal domain in a sense that the complex obtained as positive nodal domain is closed. It is trivial to obtain from this complex
//a closed complex that represents negative nodal domain. Add such a functionality!!



/**
 * */
template <typename T>
cellComplex<T>::cellComplex( std::vector< std::pair<T,T> > cube )
{
    for ( size_t i = 0 ; i != cube.size() ; ++i )
    {
        if (cube[i].first > cube[i].second)
        {
            throw std::runtime_error("Wrong coordinates of the initial cube\n");
        }
    }

    for ( size_t i = 0 ; i != cube.size() ; ++i )
    {
        this->initialCube.push_back( std::make_pair(cube[i].first , cube[i].second) );
    }
    //first make sure that cube is sorted in the following way:
    //cube[i].first <= cube[i].second
    for ( unsigned int i = 0 ; i != cube.size() ; ++i )
    {
        if ( cube[i].first > cube[i].second )
        {
             T buf = cube[i].first;
             cube[i].first = cube[i].second;
             cube[i].second = buf;
        }
    }

    //we start from a given cube:
    cell<T>* initialCube = new cell<T>( cube );

    //compute the dimension of the cube:
    int dim = 0;
    for ( unsigned int i = 0 ; i != cube.size() ; ++i )
    {
        //if the interval is not denenerated, then increase the dimension:
        if ( (cube[i].first-cube[i].second) != 0 )
        {
             ++dim;
        }
    }


    std::vector<std::list<cell<T>*> > temporaryStructureToStoreInitialCube;
    for ( int i = 0 ; i <= dim ; ++i )
    {
        std::list< cell<T>* > zbi;
        temporaryStructureToStoreInitialCube.push_back( zbi );
    }

    temporaryStructureToStoreInitialCube[dim].push_back(initialCube);
    initialCube->position() = temporaryStructureToStoreInitialCube[dim].size()-1;

    //cout << "temp[" << dim << "]: " << *initialCube << endl;

    //for every dimension:
    for ( int i = dim-1 ; i >= 0 ; --i )
    {
        //for every cube in this dimension:
        for ( typename std::list<cell<T>*>::iterator li = temporaryStructureToStoreInitialCube[i+1].begin() ; li != temporaryStructureToStoreInitialCube[i+1].end() ; ++li )
        {
            //cout << "i=" << i << ",  li: "  << endl;
            //cout << "(" << (*li)->coef[0].first << ", " << (*li)->coef[0].second << ") x (" <<
            //               (*li)->coef[1].first << ", " << (*li)->coef[1].second << ")" << endl << endl;
            
            //take a cube *li and generate all its boundary elements:
            for ( unsigned int num = 0 ; num != (*li)->coef.size() ; ++num )
            {
                //if (*li)->coef[num] is non degenerated:
                if ( (*li)->coef[num].first != (*li)->coef[num].second )
                {
                     //create 2 cells of dimension (i-1)
                     std::vector< std::pair<T,T> > first, second;
                     for ( unsigned int j = 0 ; j != (*li)->coef.size() ; ++j )
                     {
                         if ( num != j )
                         {
                             first.push_back( std::make_pair( (*li)->coef[j].first , (*li)->coef[j].second ) );
                             second.push_back( std::make_pair( (*li)->coef[j].first , (*li)->coef[j].second ) );
                         }
                         else
                         {
                             first.push_back( std::make_pair( (*li)->coef[j].first , (*li)->coef[j].first ) );
                             second.push_back( std::make_pair( (*li)->coef[j].second , (*li)->coef[j].second ) );
                         }
                     }

                     cell<T>* firstC = new cell<T>( first );
                     cell<T>* secondC = new cell<T>( second );

                     //set the boundary and coboundary elements:
                     (*li)->bound.push_back( firstC );
                     (*li)->bound.push_back( secondC );
                     firstC->coBound.push_back( *li );
                     secondC->coBound.push_back( *li );

                     //put firstC and secondC to the list vectorOfLists[i-1]
                     temporaryStructureToStoreInitialCube[i].push_back( firstC );
                     firstC->position() = temporaryStructureToStoreInitialCube[i].size()-1;
                     temporaryStructureToStoreInitialCube[i].push_back( secondC );
                     secondC->position() = temporaryStructureToStoreInitialCube[i].size()-1;
                }
            }
        }

        //remove all the repetitions in vectorOfLists[i-1] list and merge coboundary lists
        //tehere are just two repetitions of each cell in the list.
        typename std::list<cell<T>*>::iterator li = temporaryStructureToStoreInitialCube[i].begin();
        while ( li != temporaryStructureToStoreInitialCube[i].end() )
        {
              typename std::list<cell<T>*>::iterator li1 = li;
              ++li1;
              

              while ( li1 != temporaryStructureToStoreInitialCube[i].end() )
              {
                   if ( (**li) == (**li1) )
                   {
                        //cout << "li: " << *li << ", li1: " << *li1 << endl;

                        //cout << "(" << (*li1)->coef[0].first << ", " << (*li1)->coef[0].second << ") x (" <<
                        //               (*li1)->coef[1].first << ", " << (*li1)->coef[1].second << ")" << endl;
                        //cout << "--------------------------------------------------------------" << endl;
                        //cout << "found dupe! " << endl;
                        
                        
                        //join the lists of coboundary elements of *li1 to li
                        for ( typename cell<T>::CbdIterator cbd = (*li1)->cbdBegin() ; cbd != (*li1)->cbdEnd() ; ++cbd )
                        {
                            (*li)->coBound.push_back( *cbd );
                        }
                        typename std::list<cell<T>*>::iterator toErase = li1;
                        ++li1;

                        //tutaj jeszcze trzeba w kobrzegu toErase ustawic odpowiednio wskazniki!!!
                        for ( typename cell<T>::CbdIterator cbd = (*toErase)->cbdBegin() ; cbd != (*toErase)->cbdEnd() ; ++cbd )
                        {
                            //find toErase in (*cbd)->bound list:
                            typename cell<T>::BdIterator bd = (*cbd)->bdBegin();
                            while ( bd != (*cbd)->bdEnd() )
                            {
                                  if ( *bd == *toErase )
                                  {
                                        (*cbd)->bound.erase(bd);
                                        break;
                                  }
                                  ++bd;
                            }
                            (*cbd)->bound.push_back( *li );
                        }
                        temporaryStructureToStoreInitialCube[i].erase( toErase );
                   }
                   else
                   {
                        ++li1;
                   }
              }
              ++li;
        }//while ( li != this->elements[i].end() )
    }


    //rewriting temporaryStructureToStoreInitialCube to the vector<vector<cell> > in the structure of cellComplex.
    for ( unsigned dim = 0 ; dim != temporaryStructureToStoreInitialCube.size() ; ++dim )
    {
        //cerr << "dim : " << dim << endl;
        std::vector< cell<T>* > aa( temporaryStructureToStoreInitialCube[dim].size() );
        int i = 0;
        for ( typename std::list<cell<T>*>::iterator it = temporaryStructureToStoreInitialCube[dim].begin() ; it != temporaryStructureToStoreInitialCube[dim].end() ; ++it)
        {
            aa[i] = *it;
            aa[i]->position() = i;
            ++i;
            //cerr << "it : " << **it << endl;
        }
        this->elements.push_back(aa);
    }
}//cellComplex<T>::cellComplex( std::vector< std::pair<T,T> > cube )



bool divideTheCubeByGivenFactorDBG = false;
/** Divides cube by integer factor n
 * @param cub a cube
 * */
template <typename T>
typename std::vector< cell<T>* > cellComplex<T>::divideTheCubeByGivenFactor( cell<T>* cub , int n )
{
     std::vector< cell<T>* > cubes2subdiv;
     cubes2subdiv.push_back( cub );
     if ( (n==0)||(n==1) )return cubes2subdiv;

     std::vector< std::pair<T,T> > coefCords;
     for ( unsigned int i = 0 ; i != cub->coef.size() ; ++i )
     {
         coefCords.push_back( std::make_pair( cub->coef[i].first , cub->coef[i].second ) );
     }
     //in each dimension divide to n cubes:
     for ( unsigned int i = 0 ; i != coefCords.size() ; ++i )
     {
         if ( divideTheCubeByGivenFactorDBG )
         {
             std::cout << "i: " << i << "\n";
         }

         //if only cube is not degenerated in this dimension:
         if ( coefCords[i].first != coefCords[i].second )
         {
              std::list< cell<T>* > newCubes;

              if ( divideTheCubeByGivenFactorDBG )
              {
                  std::cout << "cubes2subdiv : \n";
                  for ( typename std::vector< cell<T>* >::iterator it = cubes2subdiv.begin() ; it != cubes2subdiv.end() ; ++it )
                  {
                      std::cout << **it << ", ";
                  }
                  std::cout << "\n\n\n";
              }

              for ( int div = 1 ; div != n ; ++div )
              {
                  if ( divideTheCubeByGivenFactorDBG )
                  {
                      std::cout << "div : " << div << "\n";
                  }

                  std::list< cell<T>* > aa;

                  //divide every cube in cubes2subdiv in cub->coef[div].first+(cub->coef[div].second - cub->coef[div].first)*div/n
                  for ( typename std::vector< cell<T>* >::iterator it = cubes2subdiv.begin() ; it != cubes2subdiv.end() ; ++it )
                  {
                      //and add the new elements to aa and newCubes list respectivelly
                      if ( divideTheCubeByGivenFactorDBG )
                      {
                          std::cout << "points of the division : " << (coefCords[i].first+(coefCords[i].second - coefCords[i].first)*div/n) << "\n";
                      }
                      std::vector< cell<T>* > newCubesOfLowerDimension;
                      std::vector< cell<T>* > vect = this->divideCube( (*it) , coefCords[i].first+(coefCords[i].second - coefCords[i].first)*div/n , i , newCubesOfLowerDimension, false );
                      if ( div != (n-1) )
                      {
                          newCubes.push_back(vect[0]);
                          aa.push_back( vect[1] );
                      }
                      else
                      {
                          newCubes.push_back(vect[0]);newCubes.push_back(vect[1]);
                      }

                      if ( divideTheCubeByGivenFactorDBG )
                      {
                          std::cout << "The remianing cubes : " << *vect[0] << "\n";
                          std::cout << "I am adding to aa the following stuff  : " << *vect[1] << "\n";
                      }
                  }
                  cubes2subdiv.clear();
                  for ( typename std::list< cell<T>* >::iterator bla = aa.begin() ; bla != aa.end() ; ++bla )
                  {
                      cubes2subdiv.push_back(*bla);
                  }
              }
              cubes2subdiv.clear();
              //add elements from newCubes to cubes2subdiv lists
              for ( typename std::list< cell<T>* >::iterator it = newCubes.begin() ; it != newCubes.end() ; ++it )
              {
                  cubes2subdiv.push_back( *it );
              }
         }
     }

     //in this part we will impose a boundary conditions. For this, the inital cube have to be divided by a factor at least 2.
     extern bool imposeBoundaryConditions;
     bool debugImposingBdConditions = false;
     if ( (n>1)&&(imposeBoundaryConditions) )
     {
         cout << "Imposing periodic boundary conditions... \n";
         std::vector<std::vector< cell<T>* > > boundaryCells(this->elements.size()-1);
         for ( size_t i = 0 ; i != this->elements[ this->elements.size()-2 ].size() ; ++i )
         {
             if(debugImposingBdConditions){cout << "this->elements[ this->elements.size()-2 ][i] : " << this->elements[ this->elements.size()-2 ][i] <<  endl;}
             if ( this->elements[ this->elements.size()-2 ][i]->coBound.size() == 1 )
             {
                   boundaryCells[ this->elements.size()-2 ].push_back(this->elements[ this->elements.size()-2 ][i]);
                   if(debugImposingBdConditions){cout << "is in boundary \n";}
             }
         }
         if(debugImposingBdConditions){cout << "boundaryCells[ this->elements.size()-2 ].size() : " << boundaryCells[ this->elements.size()-2 ].size() << endl;}
         for ( size_t dim = this->elements.size()-2 ; dim != 0 ; --dim )
         {
              if(debugImposingBdConditions){cerr << "dim : " << dim << endl;}
              for ( size_t nr = 0 ; nr != boundaryCells[dim].size() ; ++nr )
              {
                  for ( typename cell<T>::BdIterator bd = boundaryCells[dim][nr]->bdBegin() ; bd != boundaryCells[dim][nr]->bdEnd() ; ++bd )
                  {
                      boundaryCells[dim-1].push_back( *bd );
                  }
              }
              if(debugImposingBdConditions){cerr << "boundaryCells["<<dim-1<<"].size() : " << boundaryCells[dim-1].size() << endl;}
              //removing repetition
              std::sort( boundaryCells[dim-1].begin() , boundaryCells[dim-1].end() );
              boundaryCells[dim-1].erase(std::unique(boundaryCells[dim-1].begin(), boundaryCells[dim-1].end()), boundaryCells[dim-1].end());
              if(debugImposingBdConditions){cerr << "boundaryCells["<<dim-1<<"].size() : " << boundaryCells[dim-1].size() << endl;}
         }


         //this is a map which for a given cell returns the value of a corresponding cell for periodic boundary conditions
         std::map< cell<T>* , cell<T>* , cellComparison<T> > correspondenceMapping;
         for ( int dim = (int)boundaryCells.size()-1 ; dim >= 0 ; --dim )
         {
             if(debugImposingBdConditions){cerr << "dim : " << dim << endl;}
             for ( size_t nr = 0 ; nr != boundaryCells[dim].size() ; ++nr )
             {
                 cell<T>* paired = this->findPairedCell( boundaryCells[dim][nr] , boundaryCells[dim] );
                 if ( paired )
                 {
                     if(debugImposingBdConditions){cerr << * boundaryCells[dim][nr] << " is paired with  " << *paired << endl;}
                     correspondenceMapping[boundaryCells[dim][nr]] = paired;
                 }
             }
         }

         //first we need to substitute in boundary and coboundary lists all elements that are paired
         for ( int dim = 0 ; dim != boundaryCells.size() ; ++dim )
         {
             for ( size_t nr = 0 ; nr != boundaryCells[dim].size() ; ++nr )
             {
                 for ( typename cell<T>::BdIterator bd = boundaryCells[dim][nr]->bdBegin() ; bd != boundaryCells[dim][nr]->bdEnd() ; ++bd )
                 {
                     for ( typename cell<T>::CbdIterator cbd = (*bd)->cbdBegin() ; cbd != (*bd)->cbdEnd() ; ++cbd )
                     {
                         if ( correspondenceMapping.find(*cbd) != correspondenceMapping.end() )
                         {
                             *cbd = correspondenceMapping[*cbd];
                         }
                     }
                 }
                 for ( typename cell<T>::CbdIterator cbd = boundaryCells[dim][nr]->cbdBegin() ; cbd != boundaryCells[dim][nr]->cbdEnd() ; ++cbd )
                 {
                     for ( typename cell<T>::BdIterator bd = (*cbd)->bdBegin() ; bd != (*cbd)->bdEnd() ; ++bd )
                     {
                         if ( correspondenceMapping.find( *bd ) != correspondenceMapping.end() )
                         {
                             if(debugImposingBdConditions){cerr << "Changing boundary of : "<< *boundaryCells[dim][nr] << endl;cerr << "used to be : " << **bd  << endl;}
                             (*bd) = correspondenceMapping.find(*bd)->second;
                             if(debugImposingBdConditions){cerr << "it is now : " << **bd << endl;cin.ignore();}
                         }
                     }
                 }
             }
         }

         if(debugImposingBdConditions){cerr << "Boundary and coboundary elements has been updated \n";}

          cell<T>* dummyCell = new cell<T>();
         //now we need to remove all the elements A that has correspondenceMapping.find(A) != correspondenceMapping.end;
         for ( typename std::map< cell<T>* , cell<T>* , cellComparison<T> >::iterator it = correspondenceMapping.begin() ; it != correspondenceMapping.end() ; ++it )
         {
             this->elements[it->first->dim()][it->first->position()] = dummyCell;
             delete it->first;
         }
         for ( size_t dim = 0 ; dim != this->elements.size() ; ++dim )
         {
             if(debugImposingBdConditions){cerr << "this->elements[" << dim<< "].size = " << this->elements[dim].size() << "\n";}
             //std::remove( this->elements[dim].begin() , this->elements[dim].end() , dummyCell );
             this->elements[dim].erase(std::remove(this->elements[dim].begin(), this->elements[dim].end(), dummyCell), this->elements[dim].end());
             if(debugImposingBdConditions){cerr << "this->elements[" << dim<< "].size = " << this->elements[dim].size() << "\n";cin.ignore();}
         }
         delete dummyCell;
         cout << "Periodic boundary conditions has been imposed.\n";
     }
     return cubes2subdiv;
}

template <typename T>
inline bool containedInBdry( const std::vector< std::pair<T,T> > coface , const std::vector< std::pair<T,T> > face )
{
     if ( coface.size() == face.size() )
     {
          for ( unsigned int i = 0 ; i != coface.size() ; ++i )
          {
              //if face[i] \subset coface[i]
              if ( (face[i].first < coface[i].first)||(face[i].first>coface[i].second) )
              {
                   return false;
              }
              if ( (face[i].second < coface[i].first)||(face[i].second>coface[i].second) )
              {
                   return false;
              }
          }
          return true;
     }
     return false;
}//containedInBdry

template <typename T>
inline void cellComplex<T>::tie( cell<T>* coface, cell<T>* face )
{
     coface->bound.push_back( face );
     face->coBound.push_back( coface );
}

bool debugeDivCube = false;
//bool debugeDivCube = true;
template <typename T>
typename std::vector< cell<T>* > cellComplex<T>::divideCube( cell<T>* cube , T x , int coef , typename std::vector< cell<T>* >& vec , bool howToComputeDepth )
{
      std::vector< cell<T>* > result;

      //if cube has zero thickness in the direction of coef, then there is nothink to subdivide:
      if ( cube->coef[coef].first == cube->coef[coef].second )
      {
          cerr << "if ( cube->coef[coef].first == cube->coef[coef].second ), return \n";
           return result;
      }

      //if x \not \in [cube->coef[coef].first , cube->coef[coef].second], then there is nothink to subdivide:
      if ( (x <= cube->coef[coef].first) || ( x>=cube->coef[coef].second ) )
      {
          cerr << "(x <= cube->coef[coef].first) || ( x>=cube->coef[coef].second ), return \n";
          cerr << "x : " << x << " , cube->coef[coef].first : " << cube->coef[coef].first << " ,cube->coef[coef].second : " << cube->coef[coef].second << endl;
          getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();getchar();
           return result;
      }

      //Here I should create 3 new cells. A and B should have the same dimension as cube, and C should
      //have dim(cube)-1. Then after all the changes in the structure of boundary and coboundary I should
      //put A, B and C to the complex and remove cube from the complex.
      //However removing cube may be costly, therefore the following trick is used instead:
      //cube will not be removed, it will be changed to A. Then only B and C needs to be created and added
      //to the complex.

      //first we save the original boundary and coboundary of cube:
      std::vector< cell<T>* > boundaryOfCube;
      std::vector< cell<T>* > coboundaryOfCube;

      for ( typename cell<T>::BdIterator bd = cube->bdBegin() ; bd != cube->bdEnd() ; ++bd )
      {
          boundaryOfCube.push_back( *bd );
      }
      for ( typename cell<T>::CbdIterator cbd = cube->cbdBegin() ; cbd != cube->cbdEnd() ; ++cbd )
      {
          coboundaryOfCube.push_back( *cbd );
      }


      //now we create A, B and C:
      std::vector< std::pair<T,T> > CCoef;
      for ( unsigned int i = 0 ; i != cube->coef.size() ; ++i )
      {
          if ( i != coef )
          {
              CCoef.push_back( std::make_pair( cube->coef[i].first , cube->coef[i].second ) );
          }
          else
          {
              CCoef.push_back( std::make_pair( x , x ) );
          }
      }

      std::vector< std::pair<T,T> > BCoef;
      for ( unsigned int i = 0 ; i != cube->coef.size() ; ++i )
      {
          if ( i != coef )
          {
              BCoef.push_back( std::make_pair( cube->coef[i].first , cube->coef[i].second ) );
          }
          else
          {
              BCoef.push_back( std::make_pair( x , cube->coef[i].second ) );
          }
      }

      std::vector< std::pair<T,T> > ACoef;
      for ( unsigned int i = 0 ; i != cube->coef.size() ; ++i )
      {
          if ( i != coef )
          {
              ACoef.push_back( std::make_pair( cube->coef[i].first , cube->coef[i].second ) );
          }
          else
          {
              ACoef.push_back( std::make_pair( cube->coef[i].first , x ) );
          }
      }



      //create A, B and C:
      cell<T>* A = new cell<T>(ACoef);
      cell<T>* B = new cell<T>(BCoef);
      cell<T>* C = new cell<T>(CCoef);
      result.push_back(A);
      result.push_back(B);
      result.push_back(C);

      vec.push_back(A);
      vec.push_back(B);
      vec.push_back(C);


      if ( howToComputeDepth )
      {
     	   A->subdivisionDepth = cube->subdivisionDepth+1;
      	 B->subdivisionDepth = cube->subdivisionDepth+1;
      	 C->subdivisionDepth = cube->subdivisionDepth+1;
      }
      else
      {
         A->subdivisionDepth = cube->subdivisionDepth;
      	 B->subdivisionDepth = cube->subdivisionDepth;
      	 C->subdivisionDepth = cube->subdivisionDepth;
      }

      extern unsigned subdivisionDepth;
      if ( cube->subdivisionDepth+1 == subdivisionDepth )
      {
          cerr << "Maximal subdivision depth has been reached, program terminated \n";
          throw "Maximal subdivision depth has been reached, program terminated \n";
      }

      if ( debugeDivCube )
      {
          std::cout << "A, B i C has been created \n";
          std::cout << "A  : " << *A << "\n";
          std::cout << "B  : " << *B << "\n";
          std::cout << "C  : " << *C << "\n";
          getchar();
      }

      //remove cube from the complex:
      int cubeDim = cube->dim();

      //add A, B and C to complex:
      //TODO - temporary change!!!
      //this->elements[ cubeDim ].push_back( A );
      //A->position() = this->elements[ cubeDim ].size()-1;

      size_t pos = cube->position();
      this->elemen()[cube->dim()][pos] = A;
      A->position() = pos;



      this->elements[ cubeDim ].push_back( B );
      B->position() = this->elements[ cubeDim ].size()-1;
      this->elements[ cubeDim-1 ].push_back( C );
      C->position() = this->elements[ cubeDim-1 ].size()-1;

      if ( debugeDivCube )
      {
          cerr << "cubeDim : " << cubeDim <<  endl;
	      cout << "Put : " << *A << " to the list of dimension : " << cubeDim << " at a position " << A->position() <<  endl;
	      cout << "Put : " << *B << " to the list of dimension : " << cubeDim << " at a position " << B->position() << endl;
	      cout << "Put : " << *C << " to the list of dimension : " << cubeDim-1 << " at a position " << C->position() << endl;
	      getchar();
      }


      //Put C to the boundary of cube and B. Put cube and B to coboundary of C:

      A->bound.push_back( C );
      B->bound.push_back( C );
      C->coBound.push_back( A );
      C->coBound.push_back( B );
      /*
      cout << "\n\n\n " << *this  << "\n\n\n";
      std::cin.ignore();
      */


      //let us now iterate through the coboundary of cube and for each coboundary element E substitute
      //cube to A and B in E.boundary, and add E to A.coBound and B.cobound:
      for ( typename std::vector<cell<T>*>::iterator gamma = coboundaryOfCube.begin() ; gamma != coboundaryOfCube.end() ; ++gamma )
      {
          typename cell<T>::BdIterator bd = (*gamma)->bdBegin();
          while ( bd != (*gamma)->bdEnd() )
          {
                if ( *bd == cube )
                {
                     break;
                }
                ++bd;
          }
          (*gamma)->bound.erase( bd );
          //put A and B to boundary of (*gamma):
          (*gamma)->bound.push_back(A);
          (*gamma)->bound.push_back(B);
          //put (*gamma) to coboundary of A and B:
          A->coBound.push_back( *gamma );
          B->coBound.push_back( *gamma );
      }

      //now let us iterate through the boundaryOfCube and let us remove cube from the coboundary of beta:
      for ( typename std::vector<cell<T>*>::iterator beta = boundaryOfCube.begin() ; beta != boundaryOfCube.end() ; ++beta )
      {
          typename cell<T>::CbdIterator cbd = (*beta)->cbdBegin();
          while ( cbd != (*beta)->cbdEnd() )
          {
                if ( (*cbd) == cube )
                {
                     break;
                }
                ++cbd;
          }
          (*beta)->coBound.erase( cbd );
      }

    /*
      cout << "tutaj : \n\n";
      cout << *this << "\n\n";
      cin.ignore();
    */

      //let us now fill the bdryElementsToSubdiv list and set the boundary elements of A and B:
      for ( typename std::vector<cell<T>*>::iterator beta = boundaryOfCube.begin() ; beta != boundaryOfCube.end() ; ++beta )
      {
          //cout << "Wywolanie petli for dla : " << **beta << endl;
          //cout << *this << "\n\n";
          //cin.ignore();

          //first check if there is no element in the boundary of (*beta) that can be direcly
          //used as a boundary of C:
          for ( typename cell<T>::BdIterator bd = (*beta)->bdBegin()  ; bd != (*beta)->bdEnd() ; ++bd )
          {
              //if (*bd) is not already in the boundary of C:
              bool isBdAlreadyInBoundaryOfC = false;
              for ( typename cell<T>::BdIterator bdOfC = C->bdBegin() ; bdOfC != C->bdEnd() ; ++bdOfC )
              {
                  if ( (*bdOfC) == (*bd) )
                  {
                       isBdAlreadyInBoundaryOfC = true;
                  }
              }
              if ( isBdAlreadyInBoundaryOfC ){continue;}

              if ( containedInBdry( C->coef , (*bd)->coef ) )
              {
                   if ( debugeDivCube )
                   {
                       std::cout << "tie(" <<  *C  << " , " <<  (**bd) << ") : \n";
                   }
                   tie( C , (*bd) );
              }
          }


          //cout << "jestesmy za forem " << endl;
          //cout << *this << "\n\n";
          //cin.ignore();


          //if *beta is in the boundary of only cube A or only B, then add it to suitable boundary lists
          //and modify (*beta)->coBound.
          if ( containedInBdry( A->coef , (*beta)->coef ) )
          {
               if ( debugeDivCube )
               {
                   std::cout << "tie ( " << *A << " , "<<**beta << " ) \n";
               }
               tie( A , *beta );
          }
          else
          {
              if ( containedInBdry( B->coef , (*beta)->coef ) )
              {
                   if ( debugeDivCube )
                   {
                        std::cout << "tie ( " << *B << " , "<<**beta << " ) \n";
                   }
                   tie( B , (*beta) );
              }
              else
              {
                   if ( debugeDivCube ){cerr << "beta : " << **beta << endl;}
                  //cout << "\n\n Przed wywolaniem : " << endl;
                  //cout << *this << "\n\n\n";
                  //cin.ignore();
                  //cout << "Wywolanie rekurencyjne funkcji dla parametrow : " << **beta << " oraz : " << x << "\n";
                  //std::cin.ignore();

                  if ( debugeDivCube )
                  {
                     cerr << "Recurrent call of divideCube for : " << **beta << " , " << x << " , " << coef << "\n";
                     getchar();
                  }

                  typename std::vector< cell<T>* > faceSubdivision = divideCube( (*beta) , x , coef , vec , howToComputeDepth );



                  //cout << "faceSubdivision[0] : " << *faceSubdivision[0] << endl;
                  //cout << "faceSubdivision[1] : " << *faceSubdivision[1] << endl;
                  //cout << "faceSubdivision[2] : " << *faceSubdivision[2] << endl;
                  //std::cin.ignore();
                  //std::cin.ignore();

                  if ( faceSubdivision.size() )
                  {
                       int currentDim = A->dim();
                       for ( unsigned int i = 0 ; i != faceSubdivision.size() ; ++i )
                       {
                           if ( faceSubdivision[i]->dim() == (currentDim-1) )
                           {
                               if ( containedInBdry( A->coef , faceSubdivision[i]->coef ) )
                               {
                                   if ( debugeDivCube )
                                   {
                                       std::cout << " tie ( " << *A << " , " << *faceSubdivision[i] << "\n";
                                   }
                                   tie( A , faceSubdivision[i] );
                               }
                               else
                               {
                                   //containedInBdry( B->coef , faceSubdivision[i]->coef ) == true
                                   if ( debugeDivCube )
                                   {
                                       std::cout << " tie ( " << *B << " , " << *faceSubdivision[i] << "\n";
                                   }
                                   tie( B , faceSubdivision[i] );
                               }
                           }
                           else
                           {
                               //faceSubdivision[i]->dim() == (currentDim-2)
                               if ( debugeDivCube )
                               {
                                   std::cout << " tie ( " << *C << " , " << *faceSubdivision[i] << "\n";
                               }
                               tie( C , faceSubdivision[i] );
                           }
                       }
                  }
              }
          }
      }

      if ( debugeDivCube )
      {
        cerr << "End of the procedure divideCube \n";
        getchar();
      }

      //we are checking if cube is not somewere in the complex:
      for ( size_t dim = 0 ; dim != this->elements.size() ; ++dim )
      {
          for ( size_t nr = 0 ; nr != this->elements[dim].size() ; ++nr )
          {
              if ( this->elements[dim][nr] == cube )
              {
                  cerr << "We hae a prolem, cube is still in the complex!!! \n\n\n";
                  getchar();getchar();getchar();getchar();
              }
          }
      }
      delete cube;
      return result;
}//divideCube


bool divideCubeInAllDirectionsItHasNonzeroLenghtDBG = false;
template <typename T>
typename std::vector< cell<T>* > cellComplex<T>::divideCubeInAllDirectionsItHasNonzeroLenght( cell<T>* cub )
{
     typename std::vector< std::pair<T,T> > coefs;
     for ( size_t i = 0 ; i != cub->coef.size() ; ++i )
     {
         coefs.push_back( std::make_pair( cub->coef[i].first , cub->coef[i].second ) );
     }

     std::vector<cell<T>*> cubesToSubdivide;
     cubesToSubdivide.push_back(cub);
     unsigned dimensionOfTheInitialCube = cub->dim();

     for ( size_t dim = 0 ; dim != coefs.size() ; ++dim )
     {
         if ( divideCubeInAllDirectionsItHasNonzeroLenghtDBG ){cerr << "cubesToSubdivide.size() : " << cubesToSubdivide.size() << endl;}
         if ( coefs[dim].first == coefs[dim].second )continue;
         //in this case, we know that coefs[dim].first != coefs[dim].second. We will subdivide all the cubes in the collection in this direction:
         if ( divideCubeInAllDirectionsItHasNonzeroLenghtDBG ){cerr << "dim : " << dim << endl;}

         std::vector<cell<T>*> cubesToSubdivideNew;
         for ( size_t i = 0 ; i != cubesToSubdivide.size() ; ++i )
         {
             if ( divideCubeInAllDirectionsItHasNonzeroLenghtDBG ){cerr << "Dividing cube : " << *cubesToSubdivide[i] << endl;}
             std::vector< cell<T>* > vec;
             std::vector< cell<T>* > cubesAfterDivision = divideCube( cubesToSubdivide[i] , dim , vec );
             for ( size_t j = 0 ; j != cubesAfterDivision.size() ; ++j )
             {
                 if ( cubesAfterDivision[j]->dim() == dimensionOfTheInitialCube )
                 {
                     cubesToSubdivideNew.push_back( cubesAfterDivision[j] );
                     if ( divideCubeInAllDirectionsItHasNonzeroLenghtDBG ){cerr << "New cube resulting from subdiv : " << *cubesAfterDivision[j] << endl;}
                 }
             }
         }
         cubesToSubdivide.clear();
         cubesToSubdivide.insert( cubesToSubdivide.begin() , cubesToSubdivideNew.begin() , cubesToSubdivideNew.end() );
     }


     if ( divideCubeInAllDirectionsItHasNonzeroLenghtDBG ){cerr << "Writing down the results \n";}
     std::vector< cell<T>* > result;
     for ( size_t i = 0 ; i != cubesToSubdivide.size() ; ++i )
     {
         std::vector< cell<T>* > bd = cubesToSubdivide[i]->computeWholeBoundaryOfCube();
         for ( size_t j = 0 ; j != bd.size() ; ++j )
         {
             if ( !bd[j]->isVerified() )
             {
                 result.push_back( bd[j] );
             }
         }
     }
     result.insert(result.end() , cubesToSubdivide.begin() , cubesToSubdivide.end());
     sort( result.begin(), result.end() );
     result.erase( unique( result.begin(), result.end() ), result.end() );
     return result;

}//divideCubeInAllDirectionsItHasNonzeroLenght


/*
template <typename T>
typename std::vector< cell<T>* > cellComplex<T>::divideCubeInAllDirectionsItHasNonzeroLenght( cell<T>* cub )
{
     std::vector< cell<T>* > cubes2subdiv;
     cubes2subdiv.push_back( cub );

     std::vector< std::pair<T,T> > coefCords;
     for ( unsigned int i = 0 ; i != cub->coef.size() ; ++i )
     {
         coefCords.push_back( std::make_pair( cub->coef[i].first , cub->coef[i].second ) );
     }
     //in each dimension divide to n cubes:
     for ( unsigned int i = 0 ; i != coefCords.size() ; ++i )
     {
         if ( divideTheCubeByGivenFactorDBG )
         {
             std::cout << "i: " << i << "\n";
         }

         //if only cube is not degenerated in this dimension:
         if ( coefCords[i].first != coefCords[i].second )
         {
             cerr << "Her 1 \n";
              std::list< cell<T>* > newCubes;

              if ( divideTheCubeByGivenFactorDBG )
              {
                  std::cout << "cubes2subdiv : \n";
                  for ( typename std::vector< cell<T>* >::iterator it = cubes2subdiv.begin() ; it != cubes2subdiv.end() ; ++it )
                  {
                      std::cout << **it << ", ";
                  }
                  std::cout << "\n\n\n";
              }

              cerr << "Her 2 \n";

              for ( int div = 1 ; div != 2 ; ++div )
              {
                  if ( divideTheCubeByGivenFactorDBG )
                  {
                      std::cout << "div : " << div << "\n";
                  }

                  std::list< cell<T>* > aa;

                  cerr << "Her 2.5 \n";

                  //divide every cube in cubes2subdiv in cub->coef[div].first+(cub->coef[div].second - cub->coef[div].first)*div/n
                  for ( typename std::vector< cell<T>* >::iterator it = cubes2subdiv.begin() ; it != cubes2subdiv.end() ; ++it )
                  {
                      cerr << "**it : " << **it << endl;

                      //and add the new elements to aa and newCubes list respectivelly
                      if ( divideTheCubeByGivenFactorDBG )
                      {
                          std::cout << "points of the division : " << (coefCords[i].first+(coefCords[i].second - coefCords[i].first)*div/2) << "\n";
                      }
                      std::vector< cell<T>* > newCubesOfLowerDimension;
                      std::vector< cell<T>* > vect = this->divideCube( (*it) , i , newCubesOfLowerDimension );
                      if ( div != 1 )
                      {
                          newCubes.push_back(vect[0]);
                          aa.push_back( vect[1] );
                      }
                      else
                      {
                          newCubes.push_back(vect[0]);newCubes.push_back(vect[1]);
                      }

                      if ( divideTheCubeByGivenFactorDBG )
                      {
                          std::cout << "The remianing cubes : " << *vect[0] << "\n";
                          std::cout << "I am adding to aa the following stuff  : " << *vect[1] << "\n";
                      }
                  }

                  cerr << "Her 3 \n";

                  cubes2subdiv.clear();
                  for ( typename std::list< cell<T>* >::iterator bla = aa.begin() ; bla != aa.end() ; ++bla )
                  {
                      cubes2subdiv.push_back(*bla);
                  }

                  cerr << "Her 4 \n";
              }
              cubes2subdiv.clear();
              //add elements from newCubes to cubes2subdiv lists
              for ( typename std::list< cell<T>* >::iterator it = newCubes.begin() ; it != newCubes.end() ; ++it )
              {
                  cubes2subdiv.push_back( *it );
              }
              cerr << "Her 5 \n";
         }
     }

    std::vector< cell<T>* > result;
    for ( size_t i = 0 ; i != cubes2subdiv.size() ; ++i )
    {
        std::vector< cell<T>* > boundary = cubes2subdiv[i]->computeWholeBoundaryOfCube();
        for ( size_t bd = 0 ; bd != boundary.size() ; ++bd )
        {
            if ( boundary[bd]->isVerified() == false )
            {
                result.push_back( boundary[bd] );
            }
        }
    }
    result.insert(result.end(),cubes2subdiv.begin(),cubes2subdiv.end());

    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());

    return result;

}//divideCubeInAllDirectionsItHasNonzeroLenght

*/



template <typename T>
cell<T>* cellComplex<T>::findPairedCell( cell<T>* element , std::vector< cell<T>* >& bCellsInThisDimension )
{
   bool findPairedCellDebug = false;
   std::set<size_t> maxCord;
   bool isAtLastOneCoordinateMax = false;
   for ( size_t i = 0 ; i != element->coef.size() ; ++i )
   {
        if ( element->coef[i].first != element->coef[i].second )continue;
        //in this case we know that (element->coef[i].first == element->coef.second). We should check if this is a maximal one
        if ( element->coef[i].first == this->initialCube[i].second )
        {
            isAtLastOneCoordinateMax = true;
            maxCord.insert(i);
        }
   }
   std::vector< std::pair<T,T> > coordOfCellWeLookFor;
   if ( isAtLastOneCoordinateMax )
   {
        //maxCord contains maximal coords. We should change them to minimal cord in a cube:
        for ( size_t i = 0 ; i != element->coef.size() ; ++i )
        {
            if ( maxCord.find(i) != maxCord.end() )
            {
                //change this coord to minimal in this direction
                coordOfCellWeLookFor.push_back( std::make_pair( this->initialCube[i].first , this->initialCube[i].first ) );
            }
            else
            {
                coordOfCellWeLookFor.push_back( std::make_pair( element->coef[i].first , element->coef[i].second ) );
            }
        }
   }

   if ( findPairedCellDebug )
   {
       cerr << "The initial cube : " << endl;
       cout << *element;
       cerr << "\n The cube we are looking for : \n";
       if ( coordOfCellWeLookFor.size() )
       {
           for ( size_t i = 0 ; i != element->coef.size() ; ++i )
           {
               cout << "[" << coordOfCellWeLookFor[i].first << "," << coordOfCellWeLookFor[i].second << "] , ";
           }
       }
       else
       {
           cout << "This cube is not mapped anywere \n";
       }
       cin.ignore();
   }

   if ( coordOfCellWeLookFor.size() == 0 )
   {
       return 0;
   }

   for ( size_t nr = 0 ; nr != bCellsInThisDimension.size() ; ++nr )
   {
       bool isItThisCell = true;
       for ( size_t i = 0 ; i != bCellsInThisDimension[nr]->coef.size() ; ++i )
       {
            if ( (bCellsInThisDimension[nr]->coef[i].first != coordOfCellWeLookFor[i].first )
                  ||
                 (bCellsInThisDimension[nr]->coef[i].second != coordOfCellWeLookFor[i].second )
               )
            {
                isItThisCell = false;
                break;
            }
       }
       if ( isItThisCell )
       {
           return bCellsInThisDimension[nr];
       }
   }
   return 0;
}


template <typename T>
void cellComplex<T>::writeToFile(char* filename)
{
    ofstream out;
    out.open(filename);
    for ( int dim = 0 ; dim != this->elemen().size() ; ++dim )
    {
        out << "dim " << dim << endl;
        for ( size_t nr = 0 ; nr != this->elemen()[dim].size() ; ++nr )
        {
            out << *this->elemen()[dim][nr] << " " << this->elemen()[dim][nr]->deleted() << endl;
        }
    }
    out.close();
}



template <typename T>
std::vector< unsigned > cellComplex<T>::computeBettiNumbersOfNonDeletedPart()
{
    phat::boundary_matrix< phat::vector_vector > boundary_matrix;

    std::vector<unsigned> bettiNumbers(this->elemen().size());
    for ( size_t i = 0  ; i != this->elemen().size() ; ++i )
    {
        bettiNumbers[i] = 0;
    }

    int numberOfCells = 0;
    for ( int dim = 0 ; dim != this->elemen().size() ; ++dim )
    {
        for ( size_t nr = 0 ; nr != this->elemen()[dim].size() ; ++nr )
        {
            if ( !this->elemen()[dim][nr]->deleted() )
            {
                ++numberOfCells;
                bettiNumbers[dim]++;
            }
        }
    }

    //cerr << "numberOfCells  : " << numberOfCells  << endl;
    //for ( size_t i = 0 ; i != bettiNumbers.size() ; ++i )
    //{
    //    cerr << "B_" << i << " = " << bettiNumbers[i] << endl;
    //}
    //getchar();

    //put all the cells into a vector:
    std::vector< cell<T>* > cells(numberOfCells);
    unsigned counter = 0;
    for ( int dim = 0 ; dim != this->elemen().size() ; dim++ )
    {
        //cerr << "dim : " << dim << endl;
        for ( size_t i = 0 ; i != this->elemen()[dim].size() ; ++i )
        {
            if ( this->elemen()[dim][i]->deleted() )continue;
            cells[counter] = this->elemen()[dim][i];
            ++counter;
        }
    }

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
            if ( !(*bd)->deleted() )//I do not think we need this if statement...
            {
                temp_col.push_back( (*bd)->numberInCmplx() );
                numberElInBd++;
            }
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
    phat::compute_persistence_pairs_dualized< phat::chunk_reduction >( pairs, boundary_matrix );

    for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
    {
        cell<T>* first = numberOfGeneratorToDim[pairs.get_pair( idx ).first];
        cell<T>* second = numberOfGeneratorToDim[pairs.get_pair( idx ).second];
        bettiNumbers[ first->dim() ]--;
        bettiNumbers[ second->dim() ]--;
    }
    return bettiNumbers;
}//computeBettiNumbersOfNonDeletedPart
