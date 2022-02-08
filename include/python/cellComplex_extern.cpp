// cellComplex_extern.cpp

#include <iostream>
#include <stdexcept>

#include "../cell.hpp"
#include "../cellComplex.hpp"
//#include "../../examples/configuration.hpp"
#include "configuration.hpp"
//#include "../computePersistenceOfFunction.h"
#include "../validate.h"
#include "hess_ari.hpp"


extern "C" {
  
  /*//////////////////////////////////////////////////////////////
   *  the following 'get' function returns a pointer to a mathematical
   *  function defined in configuration.hpp */
  void * get_userDefinedFunction() { 
    /* R^n --> R,  scalar valued function!*/
    return (void*) userDefinedFunction; 
  }
  void * get_vector_field() { 
    /* R^n --> R,  scalar valued function!*/
    return (void*) vector_field; 
  }
  void * get_exit_set_function() { 
    /* R^n --> R,  scalar valued function!*/
    return (void*) exit_set; 
  }
  void * get_tangency_test_function() { 
    /* R^n --> R,  scalar valued function!*/
    return (void*) tangency_test; 
  }
  // void * get_energy_density_function() { 
  //   /* R^n --> R,  scalar valued function!*/
  //   return (void*) energy_density; 
  // }
  void * get_block_parameterization() { 
    /* R^n --> R^(n+1),  vector valued function!*/
    return (void*) block_parameterization; 
  }
  
  double * eval_vector_field(double x, double y, double z) {
    /* R^3 --> R^3,  vector valued function!*/
    
    ivector ivec(3);
    interval ix = interval(x); 
    interval iy = interval(y);
    interval iz = interval(z);
    ivec[1] = ix;
    ivec[2] = iy;
    ivec[3] = iz;
    ivector range_encl(3);
    
    HTvector fx(3);
    HTvector htvec = HessVar(ivec);
    fx = vector_field(htvec);
    range_encl = fValue(fx);

    double* return_value = new double[3];
    rvector mid_range_encl = mid(range_encl);
    return_value[0] = _double(mid_range_encl[1]);
    return_value[1] = _double(mid_range_encl[2]);
    return_value[2] = _double(mid_range_encl[3]);
    
    return return_value;
  }

  double * eval_autograd(double p[][2], size_t dim, void *f) {
 
    cout << " inside double * eval_autograd(), externC" << endl;
    std::vector< std::pair<double,double> > point;
    for (size_t i=0; i<dim; ++i) {
      point.push_back( std::make_pair(p[i][0],p[i][1]));
    }
    typedef HessType (*fptr)(const HTvector& x);
    fptr func = (fptr)f;

    cellComplex<double>* cmplx = new cellComplex<double>( point );
    cell<double>* currentCell  = *cmplx->elemen()[ point.size() ].begin();
    
    ivector xx;
    xx = currentCell->getivectorRepresentationOfCell();

    HessType fx;
    HTvector htvec = HessVar(xx);
    fx = func(htvec);
    cout << "current cell:  [" 
         << p[0][0] << ", " << p[0][1] << "] x [" 
         << p[1][0] << ", " << p[1][1] << "]" << endl; 
    cout << "f at on current cell: " << fValue(fx) << endl;
    ivector rangeEncl_Gfx;
    rangeEncl_Gfx = gradValue(fx);

    //interval func_at_xx;
    //fgEvalH(func,xx,func_at_xx,rangeEncl_Gfx);

    cout << "Gfx_1: " << rangeEncl_Gfx[1] << endl;
    cout << "Gfx_2: " << rangeEncl_Gfx[2] << endl;
    
    double* gfx = new double[4];
    gfx[0] = _double(Inf(rangeEncl_Gfx[1]));
    gfx[1] = _double(Sup(rangeEncl_Gfx[1]));

    gfx[2] = _double(Inf(rangeEncl_Gfx[2]));
    gfx[3] = _double(Sup(rangeEncl_Gfx[2]));


    return gfx;
  }


  double eval_exit_set(double x, double y) {
    /* R^2 --> R,  scalar valued function!*/
    
    ivector ivec(2);
    interval ix = interval(x);
    interval iy = interval(y);
    ivec[1] = ix;
    ivec[2] = iy;
    interval range_encl;
    
    HessType fx;
    HTvector htvec = HessVar(ivec);
    fx = exit_set(htvec);
    range_encl = fValue(fx);

    return _double(Mid(range_encl));
  }
  
  double eval_tangency_test(double x, double y) {
    /* R^2 --> R,  scalar valued function!*/
    
    ivector ivec(2);
    interval ix = interval(x);
    interval iy = interval(y);
    ivec[1] = ix;
    ivec[2] = iy;
    interval range_encl;
    
    HessType fx;
    HTvector htvec = HessVar(ivec);
    fx = tangency_test(htvec);
    range_encl = fValue(fx);

    return _double(Mid(range_encl));
  }
  
  // double eval_energy_density(double x, double y) {
  //   /* R^2 --> R,  scalar valued function!*/
  //   
  //   ivector ivec(2);
  //   interval ix = interval(x);
  //   interval iy = interval(y);
  //   ivec[1] = ix;
  //   ivec[2] = iy;
  //   interval range_encl;
  //   
  //   HessType fx;
  //   HTvector htvec = HessVar(ivec);
  //   fx = energy_density(htvec);
  //   range_encl = fValue(fx);

  //   return _double(Mid(range_encl));
  // }
  
  
  double * eval_block_parameterization(double x, double y) {
    /* R^2 --> R^3,  vector valued function!*/
    
    ivector ivec(2);
    interval ix = interval(x); 
    interval iy = interval(y);
    ivec[1] = ix;
    ivec[2] = iy;
    ivector range_encl(3);
    
    HTvector fx(3);
    HTvector htvec = HessVar(ivec);
    fx = block_parameterization(htvec);
    range_encl = fValue(fx);


    double* return_value = new double[3];
    rvector mid_range_encl = mid(range_encl);
    return_value[0] = _double(mid_range_encl[1]);
    return_value[1] = _double(mid_range_encl[2]);
    return_value[2] = _double(mid_range_encl[3]);

    //cout << "tx from C: " << return_value[0] << ", ";
    return return_value;
  }


  double pointEvalUserDefinedFunction(int dim, double x[][1]) {
    
    ivector ivec(dim);
    interval ix = interval(x[0][0]);
    interval iy = interval(x[1][0]);
    
    ivec[1] = ix;
    ivec[2] = iy;
    
    if (dim==3) {
      interval iz = interval(x[2][0]);
      ivec[3] = iz;
    }

    interval range_encl;
    
    HessType fx;
    HTvector htvec = HessVar(ivec);
    fx = userDefinedFunction(htvec);
    range_encl = fValue(fx);

    return _double(Mid(range_encl));
  }

  double pointEvalUserDefinedFunction_2d(double x, double y) {
    
    ivector ivec(2);
    interval ix = interval(x);
    interval iy = interval(y);
    ivec[1] = ix;
    ivec[2] = iy;
    interval range_encl;
    
    HessType fx;
    HTvector htvec = HessVar(ivec);
    fx = userDefinedFunction(htvec);
    range_encl = fValue(fx);

    return _double(Mid(range_encl));
  }
  
  double pointEvalUserDefinedFunction_3d(double x, double y, double z) {
    
    ivector ivec(3);
    interval ix = interval(x);
    interval iy = interval(y);
    interval iz = interval(z);
    ivec[1] = ix;
    ivec[2] = iy;
    ivec[3] = iz;
    interval range_encl;
    
    HessType fx;
    HTvector htvec = HessVar(ivec);
    fx = userDefinedFunction(htvec);
    range_encl = fValue(fx);

    return _double(Mid(range_encl));
  }
  
  void * new_cellComplex(double p[][2], size_t dim) {
    std::vector< std::pair<double,double> > point;
    for (size_t i=0; i<dim; ++i) {
      point.push_back( std::make_pair(p[i][0],p[i][1]));
    }
    cellComplex<double>* cmplx = new cellComplex<double>(point);
    return cmplx;
  }

  int dim_cellComplex(void *ptr) {
    cellComplex<double>* ref = static_cast<cellComplex<double>* >(ptr);
    return ref->dim();
  }

  void * divideTheCubeByGivenFactor(void *cmplxPtr, void *cellPtr, int n) {
    cellComplex<double>* cmplx = static_cast<cellComplex<double>* >(cmplxPtr);
    cell<double>* c = static_cast<cell<double>* >(cellPtr);
    cmplx->divideTheCubeByGivenFactor(c,n);
    return cmplx;
  }
  
  struct twoCells 
  {
    cell<double>* cell_one;
    cell<double>* cell_two;
  };

  struct twoCells divideCubeAlongLongestAxis(void *cmplxPtr, void *cellPtr) {
    cellComplex<double>* cmplx = static_cast<cellComplex<double>* >(cmplxPtr);
    cell<double>* c = static_cast<cell<double>* >(cellPtr);
    
    std::vector< cell<double>* > vec;
    cmplx->divideCubeAlongLongestAxis(c,vec);
    
    twoCells twoCells_struct = {.cell_one = vec[0], .cell_two = vec[1]} ;
    return twoCells_struct;
  }

  struct twoCells divideCube(void *cmplxPtr, void *cellPtr, int n) {
    cellComplex<double>* cmplx = static_cast<cellComplex<double>* >(cmplxPtr);
    cell<double>* c            = static_cast<cell<double>* >(cellPtr);
    
    std::vector< cell<double>* > vec;
    cmplx->divideCube(c,n,vec);
    
    twoCells twoCells_struct = {.cell_one = vec[0], .cell_two = vec[1]} ;
    return twoCells_struct;
  }

  std::vector<void *>* get_elementsAtSpecifiedDim(void *ptr, int dim) {
    cellComplex<double>* cmplx = static_cast<cellComplex<double>* >(ptr);
    std::vector<cell<double>* > e = cmplx->elemen()[dim];
    std::vector<void *>* v;
    for (int i=e.size()-1; i>=0; --i) {
      v->push_back(e[i]);
    }
    return v;
  }

  void * get_elementAtSpecifiedDimAndLoc(void *ptr, int dim, int nr) {
    cellComplex<double>* cmplx = static_cast<cellComplex<double>* >(ptr);
    cell<double>* c = cmplx->elemen()[dim][nr];
    return c;
  }

  int size_elementsAtSpecifiedDim(void *ptr, int dim) {
    cellComplex<double>* cmplx = static_cast<cellComplex<double>* >(ptr);
    return cmplx->elemen()[dim].size();
  }
  
  int size_wholeBoundaryOfCube(void *cellPtr, int dim) {
    cell<double>* c = static_cast<cell<double>* >(cellPtr);
    return c->computeWholeBoundaryOfCube().size();
  }
  
  void * get_boundaryElementAtPosition(void *cellPtr, int pos) {
    cell<double>* c = static_cast<cell<double>* >(cellPtr);
    cell<double>* bd_element = c->computeWholeBoundaryOfCube()[pos];
    return bd_element;
  }



  /// !!!!!
  /// !!!!!
  /// !!!!!
  /// !!!!!
  void * new_cell(double p[][2], size_t dim) {

    std::vector< std::pair<double,double> > point;
    for (size_t i=0; i<dim; ++i) {
      point.push_back( std::make_pair(p[i][0],p[i][1]));
    }
    cell<double>* c = new cell<double>(point);
    return c;
  }

  double * get_coordsOfCellAtLoc(void *cellPtr, int i) {
    cell<double>* ref = static_cast<cell<double>* >(cellPtr);
    double* c = new double[2];
    c[0] = ref->coords()[i].first;
    c[1] = ref->coords()[i].second;
    return c;
  }

  int get_dimOfContainingComplex(void *cellPtr) {
    cell<double>* c = static_cast<cell<double>* >(cellPtr);
    return c->coords().size();
  }

  int dim_cell(void *cellPtr) {
    cell<double>* c = static_cast<cell<double>* >(cellPtr);
    return c->dim();
  }

  void * delet(void *cellPtr) {
    cell<double>* c = static_cast<cell<double>* >(cellPtr);
    c->del();
    return c;
  }

  void * undelet(void *ptr) {
    cell<double>* ref = static_cast<cell<double>* >(ptr);
    ref->undel();
    return ref;
  }

  bool isDeleted(void *ptr) {
    cell<double>* ref = static_cast<cell<double>* >(ptr);
    return ref->deleted();
  }

  std::vector< cell<double>* > computeWholeBoundaryOfCube(void *cellPtr) {
    cell<double>* c = static_cast<cell<double>* >(cellPtr);
    std::vector< cell<double>*> whole_bd;
    
    cout << "\nin extern C: \n";
    
    whole_bd = c->computeWholeBoundaryOfCube();
    
    cout << "\nin extern C: \n" <<  whole_bd.front();
    
    return whole_bd;
  }


  void* computePersistence(double p[][2], size_t dim) {
    return *p;
    //std::vector< std::pair<double,double> > point;
    //for (size_t i=0; i<dim; ++i) {
    //  point.push_back( std::make_pair(p[i][0],p[i][1]));
    //}
    ////typedef HessType (*fptr)(const HTvector& x);
    ////fptr* func = getCircleFunc();
    //cellComplex<double>* validatedCmplx_withPersist;
    ////validatedCmplx = computeNodalDomainOfAFunction<double>(point, *func);
    ////validatedCmplx = computeNodalDomainOfAFunction<double>(point, halfPlane);
    //cout << "we are going to validate and compute persistence  now" << endl;
    ////validatedCmplx = computeNodalDomainOfAFunction<double>(point, xSquared_plus_ySquared);
    ////validatedCmplx = computeNodalDomainOfAFunction<double>(point, circle);
    ////validatedCmplx = computeNodalDomainOfAFunction<double>(point, doubleWellPotential);
    ////validatedCmplx = computeNodalDomainOfAFunction<double>(point, doubleWellPotential3d);
    //validatedCmplx_withPersist = computePersistenceOfAFunction<double>(point,R2FunctionWithMinimumOnCircleOfRadiusTwo,0.09);
    //cout << "we have validated cmplx and computed persistence" << endl;
    //return validatedCmplx_withPersist;
  }

  void* validateDomainByValidatingTopDimensionalCells(double p[][2], size_t dim, void *f) {
    std::vector< std::pair<double,double> > point;
    for (size_t i=0; i<dim; ++i) {
      point.push_back( std::make_pair(p[i][0],p[i][1]));
    }
    typedef HessType (*fptr)(const HTvector& x);
    fptr func = (fptr)f;
   
    cellComplex<double>* validatedCmplx;
    
    cout << "we are going to validate (via top dim) now" << endl;
    validatedCmplx = validateDomainByValidatingTopDimensionalCells<double>(point,func);
    cout << "we have validated cmplx" << endl;
    
    return validatedCmplx;
 
  }

  bool is_positive_on_entire_cell(double p[][2], size_t dim, void *f) {
    std::vector< std::pair<double,double> > point;
    for (size_t i=0; i<dim; ++i) {
      point.push_back( std::make_pair(p[i][0],p[i][1]));
    }
    typedef HessType (*fptr)(const HTvector& x);
    fptr func = (fptr)f;

    bool is_positive;
    is_positive = is_positive_on_entire_cell<double>(point,func);
    
    return is_positive;
    
  }

  void* computeNodalDomain(double p[][2], size_t dim, void *f) {
    std::vector< std::pair<double,double> > point;
    for (size_t i=0; i<dim; ++i) {
      point.push_back( std::make_pair(p[i][0],p[i][1]));
    }
    typedef HessType (*fptr)(const HTvector& x);
    fptr func = (fptr)f;
   
    cellComplex<double>* validatedCmplx;
    
    cout << "we are going to validate now" << endl;
    validatedCmplx = computeNodalDomainOfAFunction<double>(point,func);
    cout << "we have validated cmplx" << endl;
    
    return validatedCmplx;
  }



}
