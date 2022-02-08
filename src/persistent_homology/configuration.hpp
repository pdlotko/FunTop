// configuration.hpp

#include<cmath>
//parameters to be set by the user
#include "interval.hpp"  // predefined interval arithmetic
#include "hess_ari.hpp"  // data structures allowing auto-differentiation

using namespace cxsc;

//////////
//////////
//////////  PARAMETERS ARE SET AT THE BOTTOM OF THIS FILE!!!!!!
//////////
//////////
//////////
//////////


unsigned dimensionOfCoDomain;

double PI = 3.14159265359;
double e = 2.7182818284590452353602874713527;

//functions for optymization taken from https://en.wikipedia.org/wiki/Test_functions_for_optimization

//Ackley's function, https://en.wikipedia.org/wiki/Test_functions_for_optimization
HessType Ackleys( const HTvector& x)
  // suggested domain: x = [[-5,5],[-5,5]]
{
    assert( x.Dim() == 2 );
    return -20*(exp(-0.2*power((0.5*power((x[1]+x[2]),2)),2))) 
    - 
    exp(0.5*(cos(2*PI*x[1]) + cos(2*PI*x[2]))) + e + 20;    
}

//https://en.wikipedia.org/wiki/Test_functions_for_optimization
HessType Rosenbrock3d( const HTvector& x)
  // suggested domain: x = [[-5,5],[-5,5]]
{
    assert( x.Dim() == 3 );
    return 100*power((x[3]-power(x[2],2)),2)+power((x[2]-1),2)
           +
           100*power((x[2]-power(x[1],2)),2)+power((x[1]-1),2);
}

//https://en.wikipedia.org/wiki/Test_functions_for_optimization
HessType Rosenbrock4d( const HTvector& x)
  // suggested domain: x = [[-5,5],[-5,5],[-5,5]]
{
    assert( x.Dim() == 4 );
    return 100*power((x[4]-power(x[3],2)),2)+power((x[3]-1),2)
           +
           100*power((x[3]-power(x[2],2)),2)+power((x[2]-1),2)
           +
           100*power((x[2]-power(x[1],2)),2)+power((x[1]-1),2);
}

//https://en.wikipedia.org/wiki/Test_functions_for_optimization
HessType Rosenbrock5d( const HTvector& x)
  // suggested domain: x = [[-5,5],[-5,5],[-5,5],[-5,5]]
{
    assert( x.Dim() == 5 );    
    return 100*power((x[5]-power(x[4],2)),2)+power((x[4]-1),2)
           +
		   100*power((x[4]-power(x[3],2)),2)+power((x[3]-1),2)
           +
           100*power((x[3]-power(x[2],2)),2)+power((x[2]-1),2)
           +
           100*power((x[2]-power(x[1],2)),2)+power((x[1]-1),2);
}

//https://en.wikipedia.org/wiki/Test_functions_for_optimization
HessType Beales( const HTvector& x)
  // suggested domain: x = [[-4.5,4.5],[-4.5,4.5]]
{
    assert( x.Dim() == 2 );
    return power((1.5-x[1]+x[1]*x[2]),2)
		   +
		   power((2.25-x[1]+x[1]*power(x[2],2)),2)
		   +
		   power((2.625-x[1]+x[1]*power(x[2],3)),2)
			;
}

//https://en.wikipedia.org/wiki/Test_functions_for_optimization
HessType Goldstein_Price( const HTvector& x)
  // suggested domain: x = [[-2,2],[-2,2]]
{
    assert( x.Dim() == 2 );
    return 
    (1+power(x[1]+x[2]+1,2)*(19-14*x[1]+3*power(x[1],2)-14*x[2]+6*x[1]*x[2]+3*power(x[2],2)))
    *
    (30+power(2*x[1]-3*x[2],2)*(18-32*x[1]+12*power(x[1],2)+48*x[2]-36*x[1]*x[2]+27*power(x[2],2)));
}

//https://en.wikipedia.org/wiki/Test_functions_for_optimization
HessType McCormick( const HTvector& x)
  // suggested domain: x = [[-1,5,4],[-3,4]]
{
    assert( x.Dim() == 2 );
    return 
    sin( x[1]+x[2] )+power(x[1]-x[2],2)-1.5*x[1]+2.5*x[2]+1;
}


//http://www.geatbx.com/docu/fcnindex-01.html
HessType Griewangks_8_dim2( const HTvector& x)
  // suggested domain: x = [[-600,600],[-600,600]]
{
    assert( x.Dim() == 2 );
    return 
    power(x[1],2)/4000 + power(x[2],2)/4000 -
    cos( x[1] / sqrt(1) )
    *
    cos( x[2] / sqrt(2) )
    +1;
}

//http://www.geatbx.com/docu/fcnindex-01.html
HessType Griewangks_8_dim3( const HTvector& x)
  // suggested domain: x = [[-600,600],[-600,600],[-600,600]]
{
    assert( x.Dim() == 3 );
    return 
    power(x[1],2)/4000 + power(x[2],2)/4000 + power(x[3],2)/4000   -
    cos( x[1] / sqrt(1) )
    *
    cos( x[2] / sqrt(2) )
    *
    cos( x[3] / sqrt(3) )
    +1;
}

//http://www.geatbx.com/docu/fcnindex-01.html
HessType Griewangks_8_dim4( const HTvector& x)
  // suggested domain: x = [[-600,600],[-600,600],[-600,600],[-600,600]]
{
    assert( x.Dim() == 4 );
    return 
    power(x[1],2)/4000 + power(x[2],2)/4000 + power(x[3],2)/4000 + power(x[4],2)/4000   -
    cos( x[1] / sqrt(1) )
    *
    cos( x[2] / sqrt(2) )
    *
    cos( x[3] / sqrt(3) )
    *
    cos( x[4] / sqrt(4) )
    +1;
}

//http://www.geatbx.com/docu/fcnindex-01.html
HessType Michalewiczs_12_dim2( const HTvector& x)
  // suggested domain: x = [[0,PI],[0,PI]]
{
    assert( x.Dim() == 2 );
    return 
    -sin(x[1])*power( sin(1*power( x[1],2 )/PI) , 2*x.Dim() )
    -
    sin(x[2])*power( sin(2*power( x[2],2 )/PI) , 2*x.Dim() )
    ;
}
//http://www.geatbx.com/docu/fcnindex-01.html
HessType Michalewiczs_12_dim3( const HTvector& x)
  // suggested domain: x = [[0,PI],[0,PI],[0,PI]]
{
    assert( x.Dim() == 3 );
    return 
    -sin(x[1])*power( sin(1*power( x[1],2 )/PI) , 2*x.Dim() )
    -
    sin(x[2])*power( sin(2*power( x[2],2 )/PI) , 2*x.Dim() )
    -
    sin(x[3])*power( sin(3*power( x[3],2 )/PI) , 2*x.Dim() )
    ;
}

//http://www.geatbx.com/docu/fcnindex-01.html
HessType Michalewiczs_12_dim4( const HTvector& x)
  // suggested domain: x = [[0,PI],[0,PI],[0,PI]]
{
    assert( x.Dim() == 4 );
    return 
    -sin(x[1])*power( sin(1*power( x[1],2 )/PI) , 2*x.Dim() )
    -
    sin(x[2])*power( sin(2*power( x[2],2 )/PI) , 2*x.Dim() )
    -
    sin(x[3])*power( sin(3*power( x[3],2 )/PI) , 2*x.Dim() )
    -
    sin(x[4])*power( sin(4*power( x[4],2 )/PI) , 2*x.Dim() )
    ;
}

HessType Rastrigin2d( const HTvector& x)
  // suggested domain: x = [[-5.12,5.12],[-5.12,5.12]]
{
    assert( x.Dim() == 2 );
    return 
    10*2+
    (power(x[1],2)-10*cos(2*PI*x[1]))
    +
    (power(x[2],2)-10*cos(2*PI*x[2]));
}

HessType Rastrigin3d( const HTvector& x)
  // suggested domain: x = [[-5.12,5.12],[-5.12,5.12],[-5.12,5.12]]
{
    assert( x.Dim() == 3 );
    return 
    10*3+     
    (power(x[1],2)-10*cos(2*PI*x[1]))
    +
    (power(x[2],2)-10*cos(2*PI*x[2]))
    +
    (power(x[3],2)-10*cos(2*PI*x[3]));
}


//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1882.htm
HessType Shubert( const HTvector& x)
  // suggested domain: x = [[-5.12,5.12],[-5.12,5.12]]
{
    assert( x.Dim() == 2 );
    return 
    (
        1*cos( (1+1)*x[1] + 1 )  + 
        2*cos( (2+1)*x[1] + 2 )  +
        3*cos( (3+1)*x[1] + 3 )  +
        4*cos( (4+1)*x[1] + 4 )  +
        5*cos( (5+1)*x[1] + 5 )
    )
    *
    (
		1*cos( (1+1)*x[2] + 1 )  + 
        2*cos( (2+1)*x[2] + 2 )  +
        3*cos( (3+1)*x[2] + 3 )  +
        4*cos( (4+1)*x[2] + 4 )  +
        5*cos( (5+1)*x[2] + 5 )
    )
    ;
}

HessType R3FunctionWithMinimumOnSphereOfRadiusOne(const HTvector& x)
{
  return power( power(x[1],2) + power(x[2],2) + power(x[3],2) - 1,2);
}



//function to validate given as a computer program. This function goes from R^dimensionOfCoDomain to R.
//x^2+y^2
//template <typename T>
//T function( std::vector<T> x )
  // suggested domain: x = [[],[],[]]
//{
//    dimensionOfCoDomain = 2;
//    assert( x.size() == dimensionOfCoDomain );
//    return power(x[1],2) + power(x[2],2) - 1;
//}

HessType R2FunctionWithMinimumOnCircleOfRadiusTwo(const HTvector& x)
{
  return power(power(x[1],2) + power(x[2],2) - 4,2);
}





HessType random_trig_poly_a(const HTvector& x) {
  /*{{{*/
  // suggested domain: x = [[0,1],[0,1]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );

  HessType val(x.Dim());
  val = 0;

  double PI    = 3.1415926535897932385;

  HessType x1 = x[1];
  HessType x2 = x[2];

  int N = 5;
double A[5][5]; double B[5][5]; double C[5][5]; double D[5][5];
A[0][0] = -0.4336; A[1][0] =  3.0349; A[2][0] = -0.1241; A[3][0] = -1.2075; A[4][0] =  0.7269;
A[0][1] =  0.3426; A[1][1] =  0.7254; A[2][1] =  1.4897; A[3][1] =  0.7172; A[4][1] = -0.3034;
A[0][2] =  3.5784; A[1][2] = -0.0631; A[2][2] =  1.4090; A[3][2] =  1.6302; A[4][2] =  0.2939;
A[0][3] =  2.7694; A[1][3] =  0.7147; A[2][3] =  1.4172; A[3][3] =  0.4889; A[4][3] = -0.7873;
A[0][4] = -1.3499; A[1][4] = -0.2050; A[2][4] =  0.6715; A[3][4] =  1.0347; A[4][4] =  0.8884;

B[0][0] = -1.1471; B[1][0] =  0.3252; B[2][0] = -0.2414; B[3][0] = -0.1649; B[4][0] =  0.0774;
B[0][1] = -1.0689; B[1][1] = -0.7549; B[2][1] =  0.3192; B[3][1] =  0.6277; B[4][1] = -1.2141;
B[0][2] = -0.8095; B[1][2] =  1.3703; B[2][2] =  0.3129; B[3][2] =  1.0933; B[4][2] = -1.1135;
B[0][3] = -2.9443; B[1][3] = -1.7115; B[2][3] = -0.8649; B[3][3] =  1.1093; B[4][3] = -0.0068;
B[0][4] =  1.4384; B[1][4] = -0.1022; B[2][4] = -0.0301; B[3][4] = -0.8637; B[4][4] =  1.5326;
		   
C[0][0] = -0.7697; C[1][0] = -1.4916; C[2][0] = -1.4224; C[3][0] =  0.8351; C[4][0] = -0.0825;
C[0][1] =  0.3714; C[1][1] = -0.7423; C[2][1] =  0.4882; C[3][1] = -0.2437; C[4][1] = -1.9330;
C[0][2] = -0.2256; C[1][2] = -1.0616; C[2][2] = -0.1774; C[3][2] =  0.2157; C[4][2] = -0.4390;
C[0][3] =  1.1174; C[1][3] =  2.3505; C[2][3] = -0.1961; C[3][3] = -1.1658; C[4][3] = -1.7947;
C[0][4] = -1.0891; C[1][4] = -0.6156; C[2][4] =  1.4193; C[3][4] = -1.1480; C[4][4] =  0.8404;

D[0][0] =  0.0326; D[1][0] =  0.7481; D[2][0] =  0.2916; D[3][0] =  0.1049; D[4][0] = -0.8880;
D[0][1] =  0.5525; D[1][1] = -0.1924; D[2][1] =  0.1978; D[3][1] =  0.7223; D[4][1] =  0.1001;
D[0][2] =  1.1006; D[1][2] =  0.8886; D[2][2] =  1.5877; D[3][2] =  2.5855; D[4][2] = -0.5445;
D[0][3] =  1.5442; D[1][3] = -0.7648; D[2][3] = -0.8045; D[3][3] = -0.6669; D[4][3] =  0.3035;
D[0][4] =  0.0859; D[1][4] = -1.4023; D[2][4] =  0.6966; D[3][4] =  0.1873; D[4][4] = -0.6003;

  for (int i=0; i!=N; ++i) {
    for (int j=0; j!=N; ++j) { 
      val =  val + A[i][j]*sin(2.0*PI*(i+1)*x1)*sin(2.0*PI*(j+1)*x2) + B[i][j]*sin(2.0*PI*(i+1)*x1)*cos(2.0*PI*(j+1)*x2) +
                   C[i][j]*cos(2.0*PI*(i+1)*x1)*sin(2.0*PI*(j+1)*x2) + D[i][j]*cos(2.0*PI*(i+1)*x1)*cos(2.0*PI*(j+1)*x2);
    }
  }
  return val;
  /*}}}*/
}











HessType xPlusY(const HTvector& x) {
  // suggested domain: x = [[-2.5,2.5],[-2.5,2.5]]
  // suggested domain: x = [[0.1],[0,1]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );
  const int p1 = 2;
  const int p2 = 2;
  return x[1] + x[2];
}




HessType xSquared_plus_ySquared(const HTvector& x) {
  // suggested domain: x = [[-2.5,2.5],[-2.5,2.5]]
  // suggested domain: x = [[0.1],[0,1]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );
  const int p1 = 2;
  const int p2 = 2;
  return power(x[1]-0.5,p1) + power(x[2]-0.5,p2) - 0.2;
}

HessType halfPlane(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );
  return x[1]+x[2]-0.8;
}

HessType circle(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );
  return power(x[1],2)+power(x[2],2)-0.2;
}

HessType doubleWellPotential(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );
  double alpha = -0.75;
  //return power((power(x[1],2)-1),2) + 10*power(x[2],2) + alpha;
  return power(power(x[1],2)-1,2) + power(x[2],2) + alpha;
}

HessType shiftedDoubleWell(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );
  double alpha = -0.15;
  return power(x[1]-0.25,2) + power(x[2]-0.25,2) + alpha;
}



HessType doubleWellPotential3d(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  dimensionOfCoDomain = 3;
  assert( x.Dim() == dimensionOfCoDomain );
  //return power((power(x[1],2)-1),2) + 10*power(x[2],2) + alpha;
  return power(power(x[1],2)-1,2) + power(x[2],2) + x[3];
}


HessType otherFunction(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  dimensionOfCoDomain = 3;
  assert( x.Dim() == dimensionOfCoDomain );
  return power((power(x[1],2)-1),2) + 10*power(x[2],2) + x[3];
}

HessType constantFunction(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  dimensionOfCoDomain = x.Dim();
  HessType c(x.Dim());
  c = 1; // puts c = ([1,1], [1,1], ...) length = x.Dim()
  return c;
}

HessType TenMinusUnitSphere(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  dimensionOfCoDomain = 3;
  assert( x.Dim() == dimensionOfCoDomain );
    return 10- power(x[1],2) + power(x[2],2) + power(x[3],2);
}

//HessType thing wrong here, always evaluates to zero... (TS)
//  // suggested domain: x = [[],[],[]]
//  dimensionOfCoDomain = 2;
//  assert( x.Dim() == dimensionOfCoDomain );
//  return (-4/25) * (sqrt(((-1)*power(sin(x[2]),2) * ((-49)+(-15)*cos(2*x[2])+
//     30*cos(2*x[1])*power(sin(x[2]),2))))) *((-81)*power(cos(x[1]),2)*cos(x[2])*
//     power(sin(x[2]),3)+sin(x[1])*power(sin(x[2]),2)*((-50)+(-30)*cos(x[2])+
//     81*power(cos(x[1]),2)*power(sin(x[2]),2))+25*sin(2*x[2]));
//}

HessType twtest3_torus(const HTvector& x) {
  // suggested domain: x = [[-2,2],[-2,2],[-1,1]]
  dimensionOfCoDomain = 3;
  assert( x.Dim() == dimensionOfCoDomain );
  return power(power((power(x[1],2) + power(x[2],2)),2) - power(x[1],2) + power(x[2],2) ,2) + power(x[3],2) - 1.0/100.0;
  //return power((x[1]*x[1] + x[2]*x[2])*(x[1]*x[1] + x[2]*x[2]) - x[1]*x[1] + x[2]*x[2],2) + x[3]*x[3] - 1.0/100.0;
  //return -power(x[1]*(x[1]-1.0)*(x[1]-1.0)*(x[1]-2.0) + x[2]*x[2], 2) - x[3]*x[3] + 1.0/100.0;
  //return -(x[1]*(x[1]-1.0)*(x[1]-1.0)*(x[1]-2.0) + x[2]*x[2]) * (x[1]*(x[1]-1.0)*(x[1]-1.0)*(x[1]-2.0) + x[2]*x[2]) - x[3]*x[3] + 1.0/20.0;
}

HessType thin_torus(const HTvector& x) {
  // suggested domain: x = [[-1.5,1.5],[-1.5,1.5],[-0.5,0.5]]
  dimensionOfCoDomain = 3;
  assert( x.Dim() == dimensionOfCoDomain );
  return power(power((power(x[1],2) + power(0.5*x[2],2)),2) - power(x[1],2) + power(0.5*x[2],2) ,2) + power(x[3],2) - 1.0/100.0;
  //return power((x[1]*x[1] + x[2]*x[2])*(x[1]*x[1] + x[2]*x[2]) - x[1]*x[1] + x[2]*x[2],2) + x[3]*x[3] - 1.0/100.0;
  //return -power(x[1]*(x[1]-1.0)*(x[1]-1.0)*(x[1]-2.0) + x[2]*x[2], 2) - x[3]*x[3] + 1.0/100.0;
  //return -(x[1]*(x[1]-1.0)*(x[1]-1.0)*(x[1]-2.0) + x[2]*x[2]) * (x[1]*(x[1]-1.0)*(x[1]-1.0)*(x[1]-2.0) + x[2]*x[2]) - x[3]*x[3] + 1.0/20.0;
}

HessType twtest4(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  return power(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+0.6*0.6-0.2*0.2, 2) - 4.0*0.6*0.6*(x[1]*x[1] + x[2]*x[2]);
}


HessType sine_plus_linear_function(const HTvector& x) {
  // suggested domain: x = [[0,1],[0,1]]
  dimensionOfCoDomain = 2;
  double PI    = 3.1415926535897932385;
  assert( x.Dim() == dimensionOfCoDomain );
  double k = 5.0;
  double alpha = 2.0;
  double beta = -1.01;
  return sin(k*PI*x[1]) + alpha*x[2] + beta;
}

HessType cosTimeSine(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );
  return cos(x[1])*sin(x[2]);
}

HessType cosTimesArcSine(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );
  return power(cos(x[1]),2)*x[1]*asin(x[2]);
}


HessType rotatedDoubleWell(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );

  double PI    = 3.1415926535897932385;
  double theta = 0.25*PI;
  double C     = -0.01;

  HessType x1 = x[1] - 0.3*std::sqrt(3.0);
  HessType y1 = x[2] - 0.4*std::sqrt(2.0);

  HessType x2 = 5.0 * ( std::cos(theta)*x1 + std::sin(theta)*y1);
  HessType y2 = 5.0 * (-std::sin(theta)*x1 + std::cos(theta)*y1);

  return C + 0.25*(2.0*power(x2,2) - power(x2,4) - 2.0*power(y2,2));
}


HessType rotatedDoubleWell_simple(const HTvector& x) {
  // suggested domain: x = [[0,1],[0,1]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );

  double C     = -0.1;

  //x1 = 5*(x - 0.3 * sqrt(3));
  //y1 = 5*(y - 0.4 * sqrt(2))
  return C + 0.25*(2.0*power(5.0*(x[1]-0.3*std::sqrt(3.0)),2) - power(5.0*(x[1]-0.3*std::sqrt(3.0)),4) - 2.0*power(5.0*(x[2] - 0.4*std::sqrt(2.0)),2));
}

HessType matlab_peaks(const HTvector& x) {
  // suggested domain: x = [[-3,3],[-3,3]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );
  //z =  3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ... 
  //     - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) ... 
  //        - 1/3*exp(-(x+1).^2 - y.^2) 

  return  3.0*power(1.0-x[1],2)*exp(-power(x[1],2)) - power(x[2]+1.0,2)- 10.0*(x[1]/5.0 - power(x[1],3) - power(x[2],5))*exp(-power(x[1],2)-power(x[2],2))- 0.33*exp(-power(x[1]+1,2) - power(x[2],2)); 
}

HessType cosSquared(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );
  return power(cos(x[1]),2)*power(x[2],2) - 0.5;
}


HessType RatschekRokne6pt12(const HTvector& x) {
  // suggested domain: x = [1.5,1.5]x[1.5,1.5]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );

  double eps=0.01;
  return -(power(x[1],2) + power(x[2],2) - 1)*(power(x[1],2) + power(x[2],2) - 1 - eps);


}

HessType RatschekRokne6pt13(const HTvector& x) {
  // suggested domain: x = [[],[],[]]
  dimensionOfCoDomain = 2;
  assert( x.Dim() == dimensionOfCoDomain );

  return (x[1] + x[2])*(x[1] - x[2]);
}


double** randomNumbers;
const unsigned degreeOfPolynimial = 2;

void initializeRandomNumbers()
{
    srand( time(0) );
    randomNumbers = new double*[degreeOfPolynimial];
    for ( unsigned i = 0 ; i != degreeOfPolynimial ; ++i )
    {
        randomNumbers[i] = new double[4];
        randomNumbers[i][0] = rand()/(double)RAND_MAX;
        randomNumbers[i][1] = rand()/(double)RAND_MAX;
        randomNumbers[i][2] = rand()/(double)RAND_MAX;
        randomNumbers[i][3] = rand()/(double)RAND_MAX;
    }
}


HessType randomTrigonometricPolynomial(const HTvector& x)
{
    double PI    = 3.1415926535897932385;
    HessType result(x.Dim());
    result = 0;
    initializeRandomNumbers();
    for ( unsigned k = 0 ; k != degreeOfPolynimial ; ++k )
    {
        for ( unsigned l = 0 ; l != degreeOfPolynimial ; ++l )
        {
            //cerr << "k : " <<k << " , l : " << l << endl;
            result = result + randomNumbers[k][0]* randomNumbers[l][0] * cos( 2*PI*k*x[1] )*cos( 2*PI*l*x[2] )+
                              randomNumbers[k][1]* randomNumbers[l][1] * cos( 2*PI*k*x[1] )*sin( 2*PI*l*x[2] )+
                              randomNumbers[k][2]* randomNumbers[l][2] * sin( 2*PI*k*x[1] )*cos( 2*PI*l*x[2] )+
                              randomNumbers[k][3]* randomNumbers[l][3] * sin( 2*PI*k*x[1] )*sin( 2*PI*l*x[2] );
        }
    }
    return result;
}


// // // // // // // // // // // // // // // //
/* start boczko exit set tube definitions */
//
/*{{{*/
//HTvector boczko(const ivector& x) {
HTvector boczko(const HTvector& x) {
  /*{{{*/
  /* boczko: R^3 --> R^3, differential equation in Boczko, Example 4.2
   * Equation 4.3 in Wanner Stephens isoblockval
   */
  HTvector rhs(3);
  rhs[1] = -power(x[1],2) - x[1] - x[3];
  rhs[2] = 2*x[2] + 6*x[1]*x[2] - power(x[2],2) + 3*x[1] + x[3];
  rhs[3] = 2*x[3] - x[1]*x[3] + 5*x[2]*x[3];
  return rhs;
 /*}}}*/
}

HTvector rk(const HTvector& phiTheta) {
  /*{{{*/
  /* rk: R^2 --> R^3, requires k and its 1st and 2nd derivatives
   */
  HTvector rk(3);
  double r = 0.25;

  int dim = phiTheta.Dim();
  HessType phi   = phiTheta[1];
  HessType theta = phiTheta[2];
  HessType k1(dim),k2(dim),k3(dim);
  HessType k1t(dim),k2t(dim),k3t(dim);
  HessType k1tt(dim),k2tt(dim),k3tt(dim);

  double kh, kc;
  kh = 0.5; kc = 0.2;
  k1 = theta-1;
  k2 = theta-1;
  k3 = 4*kh*( theta-power(theta,2)) + kc;

  k1t = 1;
  k2t = 1;
  k3t = 4*kh*(1 - 2*theta);

  k1tt = 0;
  k2tt = 0;
  k3tt = -8*kh;

  HessType bb1t(dim), bb2t(dim), bb3t(dim), nn1t(dim), nn2t(dim), nn3t(dim);

  bb1t = k2t*k3tt - k3t*k2tt;
  bb2t = k3t*k1tt - k1t*k3tt;
  bb3t = k1t*k2tt - k2t*k1tt;
  nn1t = bb2t*k3t - bb3t*k2t;
  nn2t = bb3t*k1t - bb1t*k3t;
  nn3t = bb1t*k2t - bb2t*k1t;

  HessType bbnorm, nnnorm;
  bbnorm = sqrt(power(bb1t,2) + power(bb2t,2) + power(bb3t,2));
  nnnorm = sqrt(power(nn1t,2) + power(nn2t,2) + power(nn3t,2));

  HessType bb1, bb2, bb3;
  HessType nn1, nn2, nn3;
  bb1 = bb1t/bbnorm;
  bb2 = bb2t/bbnorm;
  bb3 = bb3t/bbnorm;
  nn1 = nn1t/nnnorm;
  nn2 = nn2t/nnnorm;
  nn3 = nn3t/nnnorm;

  rk[1] = k1 - r*cos(phi)*nn1 - r*sin(phi)*bb1;
  rk[2] = k2 - r*cos(phi)*nn2 - r*sin(phi)*bb2;
  rk[3] = k3 - r*cos(phi)*nn3 - r*sin(phi)*bb3;

  return rk;
 /*}}}*/
}

HessType boczko_exitSetTube(const HTvector& phiTheta) {
  /*{{{*/
  /* boczko_exitSetTube: R^2 --> R, represents dot product between
   * boczko and outward unit normal of parametrized surface defined by
   * rk and k, dk, ddk
   * Intended to be evaluate on phi \in [0,2pi], \theta \in [0,pi]
   */
  dimensionOfCoDomain = 2;
  assert( phiTheta.Dim() == dimensionOfCoDomain );

  int nvecfactor = -1;

  HessType rk1,rk2,rk3;
  HTvector htvec_rk = rk(phiTheta);
  //HTvector rk = fValue(htvec_rk);
  rk1 = htvec_rk[1];
  rk2 = htvec_rk[2];
  rk3 = htvec_rk[3];

  imatrix jac_rk = JacValue(htvec_rk);
  interval rk1_phi = jac_rk[1][1];
  interval rk2_phi = jac_rk[2][1];
  interval rk3_phi = jac_rk[3][1];

  interval rk1_theta = jac_rk[1][2];
  interval rk2_theta = jac_rk[2][2];
  interval rk3_theta = jac_rk[3][2];

  interval n1,n2,n3;
  interval n1t,n2t,n3t;
  real nvnorm;
  n1t = rk2_theta*rk3_phi - rk3_theta*rk2_phi;
  n2t = rk3_theta*rk1_phi - rk1_theta*rk3_phi;
  n3t = rk1_theta*rk2_phi - rk2_theta*rk1_phi;
  // hack!! ??? grabbing midpoint of an interval here...
  nvnorm = Mid(sqrt(power(n1t,2) + power(n2t,2) + power(n3t,2)));
  n1 = nvecfactor*n1t/nvnorm;
  n2 = nvecfactor*n2t/nvnorm;
  n3 = nvecfactor*n3t/nvnorm;

  HessType sf1, sf2, sf3;
  HTvector sf = boczko(htvec_rk);
  sf1 = sf[1];
  sf2 = sf[2];
  sf3 = sf[3];

  return sf1*n1 + sf2*n2 + sf3*n3;
 /*}}}*/
}
/*}}}*/
//
/* end of boczko exit set tube definitions */
// // // // // // // // // // // // // // // //


//// // // // // // // // // // // // // // // //
///* start eberlein gradient flow, Example 4.1, page 1867 from Wanner Stephens SIADS 2014  */
///*{{{*/
//HTvector vector_field_eberlein_gradient_flow(const HTvector& x) {
//  /*{{{*/
//  /* eberlein gradient flow R^3 --> R^3, differential equation in isoblockval.tex
//   */
//  HTvector rhs(3);
//  
//  rhs[1] = 2*x[1]*(x[3] - x[2]);
//  rhs[2] = 1 + x[3] - power(x[1],2);
//  rhs[3] = -1 + x[2] + power(x[1],2);
//  
//  return rhs;
// /*}}}*/
//}
//
//HTvector candidate_block_parameterization_eberlein_gradient_flow(const HTvector& phiTheta) {
//  // use (phi,theta) \in [0,2pi]x[0,pi]
//  dimensionOfCoDomain = 2;
//  assert( phiTheta.Dim() == dimensionOfCoDomain );
//  
//  // convert phiTheta to phi, theta
//  HessType phi; HessType theta;
//  phi   = phiTheta[1];
//  theta = phiTheta[2];
//  
//  /* 
//   * Define parametrization of candidate isolating block, i.e. a closed 2D surface in 3D: (phi,theta) |--> (x,y,z) 
//   * */
//  HessType x,y,z;
//  real sqrt5 = sqrt(5);
//  x = (3.0/sqrt5)*sin(theta)*cos(phi) - (4.0/sqrt5)*sin(theta)*sin(phi);
//  y = (6.0/sqrt5)*sin(theta)*cos(phi) - (2.0/sqrt5)*sin(theta)*sin(phi);
//  z = 2*cos(theta);
//
//  HTvector xyz(3);
//  xyz[1] = x;
//  xyz[2] = y;
//  xyz[3] = z;
//  return xyz;
//}
//
//
//HessType exit_set_eberlein3d(const HTvector& phiTheta) {
//  /*{{{*/
//  // Example from Wanner Stephens, isoblockval, Example 4.3, page 1874
//  // use (phi,theta) \in [0,2pi]x[0,pi]
//  dimensionOfCoDomain = 2;
//  assert( phiTheta.Dim() == dimensionOfCoDomain );
//
//  HessType x,y,z;
//  HTvector xyz = candidate_block_parameterization_eberlein3d(phiTheta);
//  x = xyz[1];
//  y = xyz[2];
//  z = xyz[3];
//  
//  return (1.0/180)*( 12*power(x,2) - 80*power(x,4) + 34*x*y + 20*power(x,3)*y
//                   - 25*power(y,2) + 450*power(z,2) );
//  /*}}}*/
//}
//
//HessType tangency_test_eberlein3d(const HTvector& phiTheta) {
//  /*{{{*/
//  // Example from Wanner Stephens, isoblockval, Example 4.3, page 1874
//  // use (phi,theta) \in [0,2pi]x[0,pi]
//  dimensionOfCoDomain = 2;
//  assert( phiTheta.Dim() == dimensionOfCoDomain );
//
//  HessType x,y,z;
//  HTvector xyz = candidate_block_parameterization_eberlein3d(phiTheta);
//  x = xyz[1];
//  y = xyz[2];
//  z = xyz[3];
//
//  // looking for w = <Df(p)f(p), grad Phi(p)> + <hess Phi(p)f(p), f(p)>
//  HessType f1,f2,f3,D11(2),D12(2),D13(2),D21(2),D22(2),D23(2),D31(2),D32(2),D33(2),gradPHI1(2),gradPHI2(2),gradPHI3(2),Dff_dot_gradPHI,
//           hessPHI11(2),hessPHI12(2),hessPHI13(2),hessPHI21(2),hessPHI22(2),hessPHI23(2),hessPHI31(2),hessPHI32(2),hessPHI33(2),
//           hessPHIf_dot_f, w;
//  //
//  // f(p)
//  real lmbda = 1.0/10.0;
//  HTvector rhs = vector_field_eberlein3d(xyz,lmbda);
//  f1 = rhs[1];
//  f2 = rhs[2];
//  f3 = rhs[3];
//  //
//  //
//  // Df(p)
//  D11 = -3.0*power(x,2) + lmbda; D12 = 1.0;   D13 = 0;
//  D21 = -1.0;                    D22 = lmbda; D23 = 0;
//  D31 =  0.0;                    D32 = 0.0;   D33 = 5.0;
//  //
//  //
//  // grad Phi(p)
//  gradPHI1 = (1.0/9.0)*(4.0*x - y);
//  gradPHI2 = (1.0/9.0)*(-x + (5.0/2.0)*y);
//  gradPHI3 =  0.5*z;
//  //
//  //
//  Dff_dot_gradPHI = D11*f1*gradPHI1 + D12*f2*gradPHI1 + D13*f3*gradPHI1 +
//                    D21*f1*gradPHI2 + D22*f2*gradPHI2 + D23*f3*gradPHI2 +
//                    D31*f1*gradPHI3 + D32*f2*gradPHI3 + D33*f3*gradPHI3 ;
//  //
//  //
//  // hess Phi(p)f(p)
//  hessPHI11 =  4.0/9.0; hessPHI12 = -1.0/9.0; hessPHI13 = 0.0;
//  hessPHI21 = -1.0/9.0; hessPHI22 = 5.0/18.0; hessPHI23 = 0.0;
//  hessPHI31 =   0.0;    hessPHI32 = 0.0;      hessPHI33 = 0.5;
//  //
//  //
//  // <Df(p)f(p), grad Phi(p)>
//  //
//  // <hess Phi(p)f(p), f(p)>
//  hessPHIf_dot_f = hessPHI11*power(f1,2) + hessPHI12*f2*f1       + hessPHI13*f3*f1      +
//                   hessPHI21*f1*f2       + hessPHI22*power(f2,2) + hessPHI23*f3*f2      +
//                   hessPHI31*f1*f3       + hessPHI32*f2*f3       + hessPHI33*power(f3,2);
//  
//  w = Dff_dot_gradPHI + hessPHIf_dot_f;
//
//  return w;
//
//  /*}}}*/
//}
//
//////
///* end of eberlien3d exit set tube definitions */
///*}}}*/
//// // // // // // // // // // // // // // // //

// // // // // // // // // // // // // // // //
/* start eberlein3d, Example 4.3 page 1874 from Wanner Stephens SIADS 2014 */
/*{{{*/
HTvector vector_field_eberlein3d(const HTvector& X) {
  /*{{{*/
  /* eberlien3d R^3 --> R^3, differential equation in isoblockval.tex
   */
  unsigned dim_of_arg = 2; 
  assert( X.Dim() == dim_of_arg );
  
  real lmbda = 1.0/10.0;
  
  HTvector rhs(3);
  rhs[1] = lmbda*X[1] + X[2] - power(X[1],3);
  rhs[2] = -X[1] + lmbda*X[2];
  rhs[3] = 5*X[3];
  
  return rhs;
 /*}}}*/
}

HessType candidate_block_zero_set_expression_eberlein3d(const HTvector& X) {
  unsigned dim_of_arg = 2; 
  assert( X.Dim() == dim_of_arg );

  return 2.0*power(X[1],2)/9.0 - X[1]*X[2]/9.0 + 5.0*power(X[2],2)/36.0 + power(X[3],2)/4.0 - 1.0;
}


HTvector candidate_block_parameterization_eberlein3d(const HTvector& phiTheta) {
  // use (phi,theta) \in [0,2pi]x[0,pi]
  unsigned dim_of_arg = 2; 
  assert( phiTheta.Dim() == dim_of_arg );
  
  // convert phiTheta to phi, theta
  HessType phi; HessType theta;
  phi   = phiTheta[1];
  theta = phiTheta[2];
  
  /* 
   * Define parametrization of candidate isolating block, i.e. a closed 2D surface in 3D: (phi,theta) |--> (x,y,z) 
   * */
  HTvector X(3);
  real sqrt5 = sqrt(5);
  X[1] = (3.0/sqrt5)*sin(theta)*cos(phi) - (4.0/sqrt5)*sin(theta)*sin(phi);
  X[2] = (6.0/sqrt5)*sin(theta)*cos(phi) - (2.0/sqrt5)*sin(theta)*sin(phi);
  X[3] = 2*cos(theta);

  return X;
}


HessType exit_set_eberlein3d(const HTvector& phiTheta) {
  /*{{{*/
  // Example from Wanner Stephens, isoblockval, Example 4.3, page 1874 with lambda = 1/10
  // use (phi,theta) \in [0,2pi]x[0,pi]
  dimensionOfCoDomain = 2;
  assert( phiTheta.Dim() == dimensionOfCoDomain );

  HessType x,y,z;
  HTvector xyz = candidate_block_parameterization_eberlein3d(phiTheta);
  x = xyz[1];
  y = xyz[2];
  z = xyz[3];
  
  return (1.0/180)*( 12*power(x,2) - 80*power(x,4) + 34*x*y + 20*power(x,3)*y
                   - 25*power(y,2) + 450*power(z,2) );
  /*}}}*/
}

HessType tangency_test_eberlein3d(const HTvector& phiTheta) {
  /*{{{*/
  // Example from Wanner Stephens, isoblockval, Example 4.3, page 1874
  // use (phi,theta) \in [0,2pi]x[0,pi]
  dimensionOfCoDomain = 2;
  assert( phiTheta.Dim() == dimensionOfCoDomain );

  HessType x,y,z;
  HTvector xyz = candidate_block_parameterization_eberlein3d(phiTheta);
  x = xyz[1];
  y = xyz[2];
  z = xyz[3];

  // looking for w = <Df(p)f(p), grad Phi(p)> + <hess Phi(p)f(p), f(p)>
  HessType f1,f2,f3,D11(2),D12(2),D13(2),D21(2),D22(2),D23(2),D31(2),D32(2),D33(2),gradPHI1(2),gradPHI2(2),gradPHI3(2),Dff_dot_gradPHI,
           hessPHI11(2),hessPHI12(2),hessPHI13(2),hessPHI21(2),hessPHI22(2),hessPHI23(2),hessPHI31(2),hessPHI32(2),hessPHI33(2),
           hessPHIf_dot_f, w;
  //
  // f(p)
  real lmbda = 1.0/10.0;
  real params[1] = {lmbda};
  HTvector rhs = vector_field_eberlein3d(xyz);
  f1 = rhs[1];
  f2 = rhs[2];
  f3 = rhs[3];
  //
  //
  // Df(p)
  D11 = -3.0*power(x,2) + lmbda; D12 = 1.0;   D13 = 0;
  D21 = -1.0;                    D22 = lmbda; D23 = 0;
  D31 =  0.0;                    D32 = 0.0;   D33 = 5.0;
  //
  //
  // grad Phi(p)
  gradPHI1 = (1.0/9.0)*(4.0*x - y);
  gradPHI2 = (1.0/9.0)*(-x + (5.0/2.0)*y);
  gradPHI3 =  0.5*z;
  //
  //
  Dff_dot_gradPHI = D11*f1*gradPHI1 + D12*f2*gradPHI1 + D13*f3*gradPHI1 +
                    D21*f1*gradPHI2 + D22*f2*gradPHI2 + D23*f3*gradPHI2 +
                    D31*f1*gradPHI3 + D32*f2*gradPHI3 + D33*f3*gradPHI3 ;
  //
  //
  // hess Phi(p)f(p)
  hessPHI11 =  4.0/9.0; hessPHI12 = -1.0/9.0; hessPHI13 = 0.0;
  hessPHI21 = -1.0/9.0; hessPHI22 = 5.0/18.0; hessPHI23 = 0.0;
  hessPHI31 =   0.0;    hessPHI32 = 0.0;      hessPHI33 = 0.5;
  //
  //
  // <Df(p)f(p), grad Phi(p)>
  //
  // <hess Phi(p)f(p), f(p)>
  hessPHIf_dot_f = hessPHI11*power(f1,2) + hessPHI12*f2*f1       + hessPHI13*f3*f1      +
                   hessPHI21*f1*f2       + hessPHI22*power(f2,2) + hessPHI23*f3*f2      +
                   hessPHI31*f1*f3       + hessPHI32*f2*f3       + hessPHI33*power(f3,2);
  
  w = Dff_dot_gradPHI + hessPHIf_dot_f;

  return w;

  /*}}}*/
}

////
/* end of eberlien3d exit set tube definitions */
/*}}}*/
// // // // // // // // // // // // // // // //


//define the type of subdivision that is going to be used
//dyadic, randm, goldenRatio
//subdivisionType usedSubdivisionType = dyadic;
//subdivisionType usedSubdivisionType = randm;
subdivisionType usedSubdivisionType = goldenRatio;

//maximal subdivision depth
unsigned subdivisionDepth = 150;


//This variable indicate if the user wants to use preliminary checking based on Discrete Morse Theory before checking derivatives. This test check if the rectangle collapses
//to its positive part. If it does not, then we are sure that monotonicity conditions cannot be satisfied.
//bool useDiscreteMorseTheoryToDoPreliminaryTestsOnContractability = true;

//Poisson boundary conditions
//bool imposeBoundaryConditions = true;
bool imposeBoundaryConditions = false;

unsigned initialSubdivision = 0;

// specify user defined functions
typedef HessType (*fptr)(const HTvector& x);
typedef HTvector (*v_fptr)(const HTvector& x);
