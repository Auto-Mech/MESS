/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2013, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

#ifndef QUATERNION_HH
#define QUATERNION_HH

#include "array.hh"
#include "error.hh"

/********************************************************************************************/
// normalized quaternion-to-matrix transformation
//
template <typename T>
void nquat2mat (const double* quat, T& mat)
{
  double    q0, q1, q2, q3, qq0, qq1, qq2, qq3;
  double    a01, a02, a03, a12, a13, a23;
  
  q0 = quat[0];  q1 = quat[1]; q2 = quat[2];  q3 = quat[3];

  qq0 = q0 * q0;  qq1 = q1 * q1;  qq2 = q2 * q2;  qq3 = q3 * q3;
 
  a01 = 2.0 * q0 * q1;  a02 = 2.0 * q0 * q2;  a03 = 2.0 * q0 * q3;
  a12 = 2.0 * q1 * q2;  a13 = 2.0 * q1 * q3;
  a23 = 2.0 * q2 * q3;

  mat (0, 1) = a12 + a03;
  mat (1, 0) = a12 - a03;
  mat (0, 2) = a13 - a02;
  mat (2, 0) = a13 + a02;
  mat (1, 2) = a23 + a01;
  mat (2, 1) = a23 - a01;
 
  mat (0, 0) = (qq0 + qq1 - qq2 - qq3);
  mat (1, 1) = (qq0 - qq1 + qq2 - qq3);
  mat (2, 2) = (qq0 - qq1 - qq2 + qq3);
}

// quaternion-to-matrix transformation
//
template <typename T>
void quat2mat (const double* quat, T& mat)
{
  double    q0, q1, q2, q3, qq0, qq1, qq2, qq3, qq;
  double    a01, a02, a03, a12, a13, a23;
  
  q0 = quat[0];  q1 = quat[1]; q2 = quat[2];  q3 = quat[3];

  qq0 = q0 * q0;  qq1 = q1 * q1;  qq2 = q2 * q2;  qq3 = q3 * q3;
 
  qq = qq0 + qq1 + qq2 + qq3;

  a01 = 2.0 * q0 * q1 / qq;  a02 = 2.0 * q0 * q2 / qq;  a03 = 2.0 * q0 * q3 / qq;
  a12 = 2.0 * q1 * q2 / qq;  a13 = 2.0 * q1 * q3 / qq;
  a23 = 2.0 * q2 * q3 / qq;

  mat (0, 1) = a12 + a03;
  mat (1, 0) = a12 - a03;
  mat (0, 2) = a13 - a02;
  mat (2, 0) = a13 + a02;
  mat (1, 2) = a23 + a01;
  mat (2, 1) = a23 - a01;
 
  mat (0, 0) = (qq0 + qq1 - qq2 - qq3)/qq;
  mat (1, 1) = (qq0 - qq1 + qq2 - qq3)/qq;
  mat (2, 2) = (qq0 - qq1 - qq2 + qq3)/qq;
}

void axis_rotation_quaternion (int, double, double*);

void quat_product (const double*, const double*, double*);

void euler2quat      (const double*, char,double*);
//void euler2mat       (const double*, char, TMatrix<double,3>&);

void polar2cart (const double*, double*);

void polar2av   (const double*, const double*, double*);
void polar2lv   (const double*, const double*, double*);

void euler2mf_av     (const double*, char, const double*, double*);
void euler2lf_av     (const double*, char, const double*, double*);

void euler_d2euler_m (const double*, char, const double*, const double*, double*);
void euler_m2euler_d (const double*, char, const double*, const double*, double*);

#endif
