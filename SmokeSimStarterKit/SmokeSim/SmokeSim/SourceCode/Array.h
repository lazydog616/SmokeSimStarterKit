#pragma once
#ifndef ARRAY__H
#define ARRAY_H

#include "util.h"
#include "grid_data.h"
#define UINT_MAX 0xffffffff;
inline unsigned int randhash(unsigned int seed)
{
   unsigned int i=(seed^0xA3C59AC3u)*2654435769u;
   i^=(i>>16);
   i*=2654435769u;
   i^=(i>>16);
   i*=2654435769u;
   return i;
}
inline float randhashf(unsigned int seed)
{ return randhash(seed)/(float)UINT_MAX; }

inline void get_barycentric(double x, int& i, double& f, int i_low, int i_high)
{
   double s=std::floor(x);
   i=(int)s;
   if(i<i_low){
      i=i_low;
      f=0;
   }else if(i>i_high-2){
      i=i_high-2;
      f=1;
   }else
      f=(x-s);
}
inline double lerp(const double& value0, const double& value1, double f)
{ return (1-f)*value0 + f*value1; }
inline double bilerp(const double& v00, const double& v10, 
                const double& v01, const double& v11, 
                double fx, double fy)
{ 
   return lerp(lerp(v00, v10, fx),
               lerp(v01, v11, fx), 
               fy);
}
inline double trilerp(const double& v000, const double& v100,
                 const double& v010, const double& v110,
                 const double& v001, const double& v101,  
                 const double& v011, const double& v111,
                 double fx, double fy, double fz) 
{
   return lerp(bilerp(v000, v100, v010, v110, fx, fy),
               bilerp(v001, v101, v011, v111, fx, fy),
               fz);
}


double interpolate_value(const vec3& point, const GridData& grid) {
   int i,j,k;
   double fi,fj,fk;

   get_barycentric(point[0], i, fi, 0, theDim[MACGrid::X]);
   get_barycentric(point[1], j, fj, 0, theDim[MACGrid::Y]);
   get_barycentric(point[2], k, fk, 0, theDim[MACGrid::Z]);

   return trilerp(
         grid(i,j,k), grid(i+1,j,k), grid(i,j+1,k), grid(i+1,j+1,k), 
         grid(i,j,k+1), grid(i+1,j,k+1), grid(i,j+1,k+1), grid(i+1,j+1,k+1), 
         fi,fj,fk);
}


double interpolate_gradient(vec3& gradient, const vec3& point, const GridData& grid) {
   int i,j,k;
   double fx,fy,fz;
   
   get_barycentric(point[0], i, fx, 0, theDim[MACGrid::X]);
   get_barycentric(point[1], j, fy, 0, theDim[MACGrid::Y]);
   get_barycentric(point[2], k, fz, 0, theDim[MACGrid::Z]);
   
   double v000 = grid(i,j,k);
   double v001 = grid(i,j,k+1);
   double v010 = grid(i,j+1,k);
   double v011 = grid(i,j+1,k+1);
   double v100 = grid(i+1,j,k);
   double v101 = grid(i+1,j,k+1);
   double v110 = grid(i+1,j+1,k);
   double v111 = grid(i+1,j+1,k+1);

   double ddx00 = (v100 - v000);
   double ddx10 = (v110 - v010);
   double ddx01 = (v101 - v001);
   double ddx11 = (v111 - v011);
   double dv_dx = bilerp(ddx00,ddx10,ddx01,ddx11, fy,fz);

   double ddy00 = (v010 - v000);
   double ddy10 = (v110 - v100);
   double ddy01 = (v011 - v001);
   double ddy11 = (v111 - v101);
   double dv_dy = bilerp(ddy00,ddy10,ddy01,ddy11, fx,fz);

   double ddz00 = (v001 - v000);
   double ddz10 = (v101 - v100);
   double ddz01 = (v011 - v010);
   double ddz11 = (v111 - v110);
   double dv_dz = bilerp(ddz00,ddz10,ddz01,ddz11, fx,fy);

   gradient[0] = dv_dx;
   gradient[1] = dv_dy;
   gradient[2] = dv_dz;
   
   //return value for good measure.
   return trilerp(
      v000, v100,
      v010, v110, 
      v001, v101,
      v011, v111,
      fx, fy, fz);
}

#endif