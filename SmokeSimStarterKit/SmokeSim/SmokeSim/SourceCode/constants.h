// Written by Peter Kutz.

#ifndef CONSTANTS_H
#define CONSTANTS_H



#define LERP(a,b,t) (1-t)*a + t*b
#define CUBIC(qi,di,dip1,delta_p,x) qi + di*x + (3.0f*delta_p - 2.0f*di - dip1)*x*x + (-2.0f*delta_p+di+dip1)*x*x*x

// Don't try to modify the values of these here.
// Modify the values of these in constants.cpp instead.
extern const int theMillisecondsPerFrame;
extern const int theDim[3];
extern const double theCellSize;



#endif