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

#include "quaternion.hh"

/*************************************************************************
 *************************** Quaternions *********************************
 *************************************************************************/

void axis_rotation_quaternion (int axis, double angle, double* q)
{
  q[0] =        std::cos(angle / 2.);
  q[axis + 1] = std::sin(angle / 2.);
}

// q1*q2
//
void quat_product (const double* q1, const double* q2, double* q) 
{
  const double* v1 = q1 + 1;
  const double* v2 = q2 + 1;
  double*       v  = q  + 1;

  vector_product(v1, v2, v);

  for(int i = 0; i < 3; ++i)
    //
    v[i] += q1[0] * v2[i] + q2[0] * v1[i];

  q[0] = q1[0] * q2[0] - vdot(v1, v2, 3);
}

void euler2quat (const double* euler_ang, char conv, double* q)
{// transformation from euler angles in appropriate convention to quaternion

   const double& theta = euler_ang [0];
   double phi, psi;
   switch(conv)
     {
     case 'X':
       phi   = euler_ang [1];
       psi   = euler_ang [2];
       break;
     case 'Y':
       phi   = euler_ang [1] + M_PI_2;
       psi   = euler_ang [2] - M_PI_2;
       break;
     default:
       error ("euler2quat: wrong convention");
     }

   const double two = 2.0;

   const double sum = (phi+psi)/two, dif = (phi-psi)/two;

   const double cos_theta_2 = cos(theta/two), sin_theta_2 = sin(theta/two);

   q[0] = cos_theta_2 * cos(sum);
   q[1] = sin_theta_2 * cos(dif);
   q[2] = sin_theta_2 * sin(dif);
   q[3] = cos_theta_2 * sin(sum);
}

void euler2mat (const double* euler_ang, char conv, TMatrix<double,3>& mat)
{// rotational matrix in terms of euler angles in appropriate convention
 
   const double& theta = euler_ang [0];
   double phi, psi;
   switch(conv)
     {
     case 'X':
       phi   = euler_ang [1];
       psi   = euler_ang [2];
       break;
     case 'Y':
       phi   = euler_ang [1] + M_PI_2;
       psi   = euler_ang [2] - M_PI_2;
       break;
     default:
       error ("euler2mat: wrong convention");
     }

   const double cos_theta = cos(theta), cos_phi = cos(phi), cos_psi = cos(psi);
   const double sin_theta = sin(theta), sin_phi = sin(phi), sin_psi = sin(psi);

   mat(0,0) =  cos_phi * cos_psi - sin_phi * cos_theta * sin_psi;
   mat(0,1) =  sin_phi * cos_psi + cos_phi * cos_theta * sin_psi;
   mat(0,2) =  sin_theta * sin_psi;
   mat(1,0) = -cos_phi * sin_psi - sin_phi * cos_theta * cos_psi;
   mat(1,1) = -sin_phi * sin_psi + cos_phi * cos_theta * cos_psi;
   mat(1,2) =  sin_theta * cos_psi;
   mat(2,0) =  sin_phi * sin_theta;
   mat(2,1) = -cos_phi * sin_theta;
   mat(2,2) =  cos_theta;
}

void polar2cart (const double* polar_ang, double* n)
{// polar_angles-to-directional_cosines transformation

   const double& theta = polar_ang[0];
   const double& phi   = polar_ang[1];

   const double sin_theta = sin(theta);

   n[0] = sin_theta * cos(phi);
   n[1] = sin_theta * sin(phi);
   n[2] = cos(theta);
}


void euler2mf_av(const double* euler_ang, char conv, const double* euler_der, double* av)
{// transformation of euler angles derivatives to angular velocity in molecular frame

  const double& theta = euler_ang [0];
  double psi;
  switch(conv) {

  case 'X':

    psi   = euler_ang [2];
    break;
  case 'Y':

    psi   = euler_ang [2] - M_PI_2;
    break;
  default:
    error ("euler2mf_av: wrong convention");
  }

  const double& theta_d = euler_der [0];
  const double& phi_d   = euler_der [1];
  const double& psi_d   = euler_der [2];

  const double sin_theta = sin(theta), sin_psi = sin(psi), cos_psi = cos(psi);

  av[0] = phi_d * sin_theta * sin_psi + theta_d * cos_psi;
  av[1] = phi_d * sin_theta * cos_psi - theta_d * sin_psi;
  av[2] = phi_d * cos(theta) + psi_d;
}

void euler2lf_av(const double* euler_ang, char conv, const double* euler_der, double* av)
{// transformation of euler angles derivatives to angular velocity in laboratory frame

  const double& theta = euler_ang [0];
  double phi;
  switch(conv) {

  case 'X':
    phi   = euler_ang [1];

    break;
  case 'Y':
    phi   = euler_ang [1] + M_PI_2;

    break;
  default:
    error ("euler2lf_av: wrong convention");
  }

  const double& theta_d = euler_der [0];
  const double& phi_d   = euler_der [1];
  const double& psi_d   = euler_der [2];

  const double sin_theta = sin(theta), sin_phi = sin(phi), cos_phi = cos(phi);

  av[0] = theta_d * cos_phi + psi_d * sin_theta * sin_phi;
  av[1] = theta_d * sin_phi - psi_d * sin_theta * cos_phi;
  av[2] = phi_d + psi_d * cos(theta);
}

void polar2av(const double* polar_ang, const double* polar_der, double* av)
{// transformation from polar angles derivatives to angular velocity

   const double& theta = polar_ang [0];
   const double& phi   = polar_ang [1];

   const double& theta_d = polar_der [0];
   const double& phi_d   = polar_der [1];

   const double sin_phi = sin(phi), cos_phi = cos(phi), sin_theta = sin(theta);
   const double p = sin_theta * cos(theta);
   
   av[0] = - theta_d * sin_phi - phi_d * p * cos_phi;
   av[1] =   theta_d * cos_phi - phi_d * p * sin_phi;
   av[2] =                       phi_d * sin_theta * sin_theta;
}

void polar2lv(const double* polar_ang, const double* polar_der, double* lv)
{// transformation from polar angles derivatives to linear velocity on the unity sphere

   const double& theta = polar_ang [0];
   const double& phi   = polar_ang [1];

   const double& theta_d = polar_der [0];
   const double& phi_d   = polar_der [1];

   const double sin_phi = sin(phi), cos_phi = cos(phi), 
     sin_theta = sin(theta), cos_theta = cos(theta);
   
   lv[0] =   theta_d * cos_theta * cos_phi - phi_d * sin_theta * sin_phi;
   lv[1] =   theta_d * cos_theta * sin_phi + phi_d * sin_theta * cos_phi;
   lv[2] = - theta_d * sin_theta;
}

void euler_d2euler_m(const double* euler_ang, char conv, const double* iner_mom, 
		     const double* euler_d, double* euler_m)
{// transformation from euler angles derivatives to corresponding generalized momenta

  const double& theta = euler_ang [0];
  double psi;
  switch(conv) {

  case 'X':

    psi   = euler_ang [2];
    break;
  case 'Y':

    psi   = euler_ang [2] - M_PI_2;
    break;
  default:
    error ("euler_d2euler_m: wrong convention");
  }

  const double two = 2.0;

  const double cos_2psi = cos(two * psi), sin_2psi = sin(two * psi);
  const double cos_theta = cos(theta), sin_theta = sin(theta);
  const double im_sum = (iner_mom[0] + iner_mom[1]) / two;
  const double im_dif = (iner_mom[0] - iner_mom[1]) / two;

  TMatrix<double, 3> trans;

  trans(0,0) = im_sum + im_dif*cos_2psi;
  trans(1,1) = sin_theta * sin_theta * (im_sum - im_dif * cos_2psi) +  
    iner_mom[2] * cos_theta * cos_theta;
  trans(2,2) = iner_mom[2];

  double temp = im_dif * sin_theta * sin_2psi;
  trans(0,1) = temp;
  trans(1,0) = temp;

  temp = iner_mom[2] * cos_theta;
  trans(1,2) = temp;
  trans(2,1) = temp;

  trans(0,2) = 0.0;
  trans(2,0) = 0.0;

  matrix_vector_product (trans, euler_d, euler_m);
}

void euler_m2euler_d(const double* euler_ang, char conv, const double* iner_mom, 
		     const double* euler_m, double* euler_d)
{// transformation from euler angles momenta to derivatives

   const double& theta = euler_ang [0];
   double psi;
   switch(conv)   {

     case 'X':

        psi   = euler_ang [2];
       break;
     case 'Y':

        psi   = euler_ang [2] - M_PI_2;
       break;
     default:
       error ("euler_m2euler_d: wrong convention");
     }

   const double two = 2.0;

   const double cos_2psi = cos(two * psi), sin_2psi = sin(two * psi);
   const double cos_theta = cos(theta), sin_theta = sin(theta);
   const double im_sum = (iner_mom[0] + iner_mom[1]) / two;
   const double im_dif = (iner_mom[0] - iner_mom[1]) / two;

   TMatrix<double, 3> trans;

   trans(0,0) = im_sum + im_dif * cos_2psi;
   trans(1,1) = sin_theta * sin_theta * (im_sum - im_dif * cos_2psi) +  
      iner_mom[2] * cos_theta * cos_theta;
   trans(2,2) = iner_mom[2];

   double temp = im_dif * sin_theta * sin_2psi;
   trans(0,1) = temp;
   trans(1,0) = temp;

   temp = iner_mom[2] * cos_theta;
   trans(1,2) = temp;
   trans(2,1) = temp;

   trans(0,2) = 0.0;
   trans(2,0) = 0.0;

   trans.invert();

   matrix_vector_product (trans, euler_m, euler_d);
}

