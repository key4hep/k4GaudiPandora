#include "ClusterShapes.h"

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector.h>

#include <algorithm>
#include <iostream>

// #################################################
// #####                                       #####
// #####  Additional Structures and Functions  #####
// #####                                       #####
// #################################################

struct data {
  int n;
  float* x;
  float* y;
  float* z;
  float* ex;
  float* ey;
  float* ez;
};

int signum(float x) {
  // computes the signum of x. Needed for the 3. parametrisation

  if (x >= 0)
    return 1; // x == 0 is taken as positive
  else
    return -1;
}

int functParametrisation1(const gsl_vector* par, void* d, gsl_vector* f) {
  //     For helix fitting
  // calculate fit function f0[i] =
  // ( (x0 + R*cos(b*z[i] + phi0)) - x[i] ) for i = 0 to n-1
  //                    and f1[i] =
  // ( (y0 + R*sin(b*z[i] + phi0)) - y[i] ) for i = n to dim*n - 1
  // That means, minimise the two functions f0[i] and f1[i]

  float x0 = gsl_vector_get(par, 0);
  float y0 = gsl_vector_get(par, 1);
  float R = gsl_vector_get(par, 2);
  float b = gsl_vector_get(par, 3);
  float phi0 = gsl_vector_get(par, 4);
  int n = ((struct data*)d)->n;
  float* x = ((struct data*)d)->x;
  float* y = ((struct data*)d)->y;
  float* z = ((struct data*)d)->z;
  // float* ex = ((struct data*)d)->ex;
  // float* ey = ((struct data*)d)->ey;
  // float* ez = ((struct data*)d)->ez;
  float fi = 0.0;

  // first dimension
  for (int i = 0; i < n; i++) {
    fi = (x0 + R * cos(b * z[i] + phi0)) - x[i];
    //    float err = sqrt(ex[i]*ex[i]+R*b*R*b*sin(b*z[i] + phi0)*R*b*R*b*sin(b*z[i] + phi0)*ez[i]*ez[i]);
    //    fi = fi/err;
    gsl_vector_set(f, i, fi);
  }
  // second dimension
  for (int i = 0; i < n; i++) {
    fi = (y0 + R * sin(b * z[i] + phi0)) - y[i];
    //    float err = sqrt(ey[i]*ey[i]+R*b*R*b*cos(b*z[i] + phi0)*R*b*R*b*cos(b*z[i] + phi0)*ez[i]*ez[i]);
    //    fi = fi/err;
    gsl_vector_set(f, i + n, fi);
  }

  return GSL_SUCCESS;
}

int dfunctParametrisation1(const gsl_vector* par, void* d, gsl_matrix* J) {
  //     For helix fitting
  float R = gsl_vector_get(par, 2);
  float b = gsl_vector_get(par, 3);
  float phi0 = gsl_vector_get(par, 4);

  int n = ((struct data*)d)->n;
  float* z = ((struct data*)d)->z;
  // float* ex = ((struct data*)d)->ex;
  // float* ey = ((struct data*)d)->ey;
  // float* ez = ((struct data*)d)->ez;

  // calculate Jacobi's matrix J[i][j] = dfi/dparj

  // part of Jacobi's matrix corresponding to first dimension
  for (int i = 0; i < n; i++) {
    //    float err = sqrt(ex[i]*ex[i]+R*b*R*b*sin(b*z[i] + phi0)*R*b*R*b*sin(b*z[i] + phi0)*ez[i]*ez[i]);
    gsl_matrix_set(J, i, 0, 1);
    gsl_matrix_set(J, i, 1, 0);
    gsl_matrix_set(J, i, 2, cos(b * z[i] + phi0));
    gsl_matrix_set(J, i, 3, -z[i] * R * sin(b * z[i] + phi0));
    gsl_matrix_set(J, i, 4, -R * sin(b * z[i] + phi0));
  }

  // part of Jacobi's matrix corresponding to second dimension
  for (int i = 0; i < n; i++) {
    //    float err = sqrt(ey[i]*ey[i]+R*b*R*b*cos(b*z[i] + phi0)*R*b*R*b*cos(b*z[i] + phi0)*ez[i]*ez[i]);
    gsl_matrix_set(J, i + n, 0, 0);
    gsl_matrix_set(J, i + n, 1, 1);
    gsl_matrix_set(J, i + n, 2, sin(b * z[i] + phi0));
    gsl_matrix_set(J, i + n, 3, z[i] * R * cos(b * z[i] + phi0));
    gsl_matrix_set(J, i + n, 4, R * cos(b * z[i] + phi0));
  }

  return GSL_SUCCESS;
}

int fdfParametrisation1(const gsl_vector* par, void* d, gsl_vector* f, gsl_matrix* J) {
  //     For helix fitting
  functParametrisation1(par, d, f);
  dfunctParametrisation1(par, d, J);

  return GSL_SUCCESS;
}

int functParametrisation2(const gsl_vector* par, void* d, gsl_vector* f) {
  //     For helix fitting
  // calculate fit function f0[i] =
  // ( (x0 + R*cos(phi)) - x[i] ) for i = 0 to n-1
  //                        f1[i] =
  // ( (y0 + R*sin(phi)) - y[i] ) for i = n to dim*n - 1
  //                    and f2[i] =
  // ( (z0 + b*phi     ) - z[i] )
  // That means, minimise the three functions f0[i], f1[i] and f2[i]

  float x0 = gsl_vector_get(par, 0);
  float y0 = gsl_vector_get(par, 1);
  float z0 = gsl_vector_get(par, 2);
  float R = gsl_vector_get(par, 3);
  float b = gsl_vector_get(par, 4);
  int n = ((struct data*)d)->n;
  float* x = ((struct data*)d)->x;
  float* y = ((struct data*)d)->y;
  float* z = ((struct data*)d)->z;
  float fi = 0.0;
  float phii = 0.0;

  // first dimension
  for (int i = 0; i < n; i++) {
    phii = atan2(y[i] - y0, x[i] - x0);
    fi = (x0 + R * cos(phii)) - x[i];
    gsl_vector_set(f, i, fi);
  }
  // second dimension
  for (int i = 0; i < n; i++) {
    phii = atan2(y[i] - y0, x[i] - x0);
    fi = (y0 + R * sin(phii)) - y[i];
    gsl_vector_set(f, i + n, fi);
  }
  // third dimension
  for (int i = 0; i < n; i++) {
    phii = atan2(y[i] - y0, x[i] - x0);
    fi = (z0 + b * phii) - z[i];
    gsl_vector_set(f, i + 2 * n, fi);
  }

  return GSL_SUCCESS;
}

int dfunctParametrisation2(const gsl_vector* par, void* d, gsl_matrix* J) {
  //     For helix fitting
  float x0 = gsl_vector_get(par, 0);
  float y0 = gsl_vector_get(par, 1);
  float R = gsl_vector_get(par, 3);
  float b = gsl_vector_get(par, 4);
  int n = ((struct data*)d)->n;
  float* x = ((struct data*)d)->x;
  float* y = ((struct data*)d)->y;
  float phii = 0.0;

  // calculate Jacobi's matrix J[i][j] = dfi/dparj

  // part of Jacobi's matrix corresponding to first dimension
  for (int i = 0; i < n; i++) {
    phii = atan2(y[i] - y0, x[i] - x0);

    gsl_matrix_set(J, i, 0,
                   1 - R * sin(phii) * ((y[i] - y0) / ((x[i] - x0) * (x[i] - x0) + (y[i] - y0) * (y[i] - y0))));
    gsl_matrix_set(J, i, 1, R * sin(phii) * ((x[i] - x0) / ((x[i] - x0) * (x[i] - x0) + (y[i] - y0) * (y[i] - y0))));
    gsl_matrix_set(J, i, 2, 0);
    gsl_matrix_set(J, i, 3, cos(phii));
    gsl_matrix_set(J, i, 4, 0);
  }

  // part of Jacobi's matrix corresponding to second dimension
  for (int i = 0; i < n; i++) {
    phii = atan2(y[i] - y0, x[i] - x0);

    gsl_matrix_set(J, i + n, 0,
                   R * cos(phii) * ((y[i] - y0) / ((x[i] - x0) * (x[i] - x0) + (y[i] - y0) * (y[i] - y0))));
    gsl_matrix_set(J, i + n, 1,
                   1 + R * cos(phii) * ((x[i] - x0) / ((x[i] - x0) * (x[i] - x0) + (y[i] - y0) * (y[i] - y0))));
    gsl_matrix_set(J, i + n, 2, 0);
    gsl_matrix_set(J, i + n, 3, sin(phii));
    gsl_matrix_set(J, i + n, 4, 0);
  }

  // part of Jacobi's matrix corresponding to third dimension
  for (int i = 0; i < n; i++) {
    phii = atan2(y[i] - y0, x[i] - x0);

    gsl_matrix_set(J, i + 2 * n, 0, b * ((y[i] - y0) / ((x[i] - x0) * (x[i] - x0) + (y[i] - y0) * (y[i] - y0))));
    gsl_matrix_set(J, i + 2 * n, 1, b * ((x[i] - x0) / ((x[i] - x0) * (x[i] - x0) + (y[i] - y0) * (y[i] - y0))));
    gsl_matrix_set(J, i + 2 * n, 2, 1);
    gsl_matrix_set(J, i + 2 * n, 3, 0);
    gsl_matrix_set(J, i + 2 * n, 4, phii);
  }

  return GSL_SUCCESS;
}

int fdfParametrisation2(const gsl_vector* par, void* d, gsl_vector* f, gsl_matrix* J) {
  //     For helix fitting
  functParametrisation2(par, d, f);
  dfunctParametrisation2(par, d, J);

  return GSL_SUCCESS;
}

int functParametrisation3(const gsl_vector* par, void* d, gsl_vector* f) {
  //     For helix fitting
  // calculate fit function f0[i] =
  // ( ( ( (1/omega) - d0 )*sin(Phi0) + ( 1/fabs(omega) )*cos( ( -omega/sqrt(1+tanL^2) )*s + Phi0 +(
  // (omega*pi)/(2*fabs(omega)) ) ) ) - x[i] ) for i = 0 to n-1
  //                        f1[i] =
  // ( ( (-1.0)*( (1/omega) - d0 )*cos(Phi0) + ( 1/fabs(omega) )*sin( ( -omega/sqrt(1+tanL^2) )*s + Phi0 +(
  // (omega*pi)/(2*fabs(omega)) ) ) ) - y[i] ) for i = n to dim*n - 1
  //                    and f2[i] =
  // ( ( z0 + (tanL/sqrt(1+tanL^2))*s ) - z[i] )
  // That means, minimise the three functions f0[i], f1[i] and f2[i]

  double z0 = gsl_vector_get(par, 0);
  double Phi0 = gsl_vector_get(par, 1);
  double omega = gsl_vector_get(par, 2);
  double d0 = gsl_vector_get(par, 3);
  double tanL = gsl_vector_get(par, 4);
  int n = ((struct data*)d)->n;
  float* x = ((struct data*)d)->x;
  float* y = ((struct data*)d)->y;
  float* z = ((struct data*)d)->z;
  double phii = 0.0;
  double fi = 0.0;
  double si = 0.0;

  // first dimension
  for (int i = 0; i < n; i++) {
    phii = atan2((((double)y[i]) + ((1 / omega) - d0) * cos(Phi0)), (((double)x[i]) - ((1 / omega) - d0) * sin(Phi0)));
    fi = (((1 / omega) - d0) * sin(Phi0) + (1 / fabs(omega)) * cos(phii)) - ((double)x[i]);
    gsl_vector_set(f, i, fi);
  }
  // second dimension
  for (int i = 0; i < n; i++) {
    phii = atan2((((double)y[i]) + ((1 / omega) - d0) * cos(Phi0)), (((double)x[i]) - ((1 / omega) - d0) * sin(Phi0)));
    fi = ((-1.0) * ((1 / omega) - d0) * cos(Phi0) + (1 / fabs(omega)) * sin(phii)) - ((double)y[i]);
    gsl_vector_set(f, i + n, fi);
  }
  // third dimension
  for (int i = 0; i < n; i++) {
    phii = atan2((((double)y[i]) + ((1 / omega) - d0) * cos(Phi0)), (((double)x[i]) - ((1 / omega) - d0) * sin(Phi0)));
    si = (-1.0) * ((sqrt(1 + pow(tanL, 2))) / omega) * (phii - Phi0 - (omega * M_PI) / (2 * fabs(omega)));
    fi = (z0 + (tanL / sqrt(1 + pow(tanL, 2))) * si) - ((double)z[i]);
    gsl_vector_set(f, i + 2 * n, fi);
  }

  return GSL_SUCCESS;
}

int dfunctParametrisation3(const gsl_vector* par, void* d, gsl_matrix* J) {
  //     For helix fitting
  // double z0    = gsl_vector_get(par,0); // not needed
  double Phi0 = gsl_vector_get(par, 1);
  double omega = gsl_vector_get(par, 2);
  double d0 = gsl_vector_get(par, 3);
  double tanL = gsl_vector_get(par, 4);
  int n = ((struct data*)d)->n;
  float* x = ((struct data*)d)->x;
  float* y = ((struct data*)d)->y;
  // float* z = ((struct data*)d)->z; // not needed
  double phii = 0.0;
  double si = 0.0;

  // calculate Jacobi's matrix J[i][j] = dfi/dparj

  // part of Jacobi's matrix corresponding to first dimension
  for (int i = 0; i < n; i++) {
    phii = atan2((((double)y[i]) + ((1 / omega) - d0) * cos(Phi0)), (((double)x[i]) - ((1 / omega) - d0) * sin(Phi0)));
    si = (-1.0) * ((sqrt(1 + pow(tanL, 2))) / omega) * (phii - Phi0 - (omega * M_PI) / (2 * fabs(omega)));

    gsl_matrix_set(J, i, 0, 0);
    gsl_matrix_set(J, i, 1, ((1 / omega) - d0) * cos(Phi0) - (1 / fabs(omega)) * sin(phii));
    gsl_matrix_set(J, i, 2,
                   ((-1.0) * sin(Phi0)) / pow(omega, 2) - ((signum(omega)) / (pow(fabs(omega), 2))) * cos(phii) -
                       (1 / fabs(omega)) * sin(phii) *
                           ((((-1.0) / sqrt(1 + pow(tanL, 2))) * si) + (M_PI) / (2 * fabs(omega)) -
                            (signum(omega) * omega * M_PI) / (2 * pow(fabs(omega), 2))));
    gsl_matrix_set(J, i, 3, (-1.0) * sin(Phi0));
    gsl_matrix_set(J, i, 4,
                   ((-1.0) / fabs(omega)) * sin(phii) * ((omega * tanL * si) / sqrt(pow(1 + pow(tanL, 2), 3))));
  }

  // part of Jacobi's matrix corresponding to second dimension
  for (int i = 0; i < n; i++) {
    phii = atan2((((double)y[i]) + ((1 / omega) - d0) * cos(Phi0)), (((double)x[i]) - ((1 / omega) - d0) * sin(Phi0)));
    si = (-1.0) * ((sqrt(1 + pow(tanL, 2))) / omega) * (phii - Phi0 - (omega * M_PI) / (2 * fabs(omega)));

    gsl_matrix_set(J, i + n, 0, 0);
    gsl_matrix_set(J, i + n, 1, ((1 / omega) - d0) * sin(Phi0) + (1 / fabs(omega)) * cos(phii));
    gsl_matrix_set(J, i + n, 2,
                   cos(Phi0) / pow(omega, 2) + ((signum(omega)) / (pow(fabs(omega), 2))) * sin(phii) +
                       (1 / fabs(omega)) * cos(phii) *
                           ((((-1.0) / sqrt(1 + pow(tanL, 2))) * si) + (M_PI) / (2 * fabs(omega)) -
                            (signum(omega) * omega * M_PI) / (2 * pow(fabs(omega), 2))));
    gsl_matrix_set(J, i + n, 3, cos(Phi0));
    gsl_matrix_set(J, i + n, 4, (1 / fabs(omega)) * cos(phii) * ((omega * tanL * si) / sqrt(pow(1 + pow(tanL, 2), 3))));
  }

  // part of Jacobi's matrix corresponding to third dimension
  for (int i = 0; i < n; i++) {
    phii = atan2((((double)y[i]) + ((1 / omega) - d0) * cos(Phi0)), (((double)x[i]) - ((1 / omega) - d0) * sin(Phi0)));
    si = (-1.0) * ((sqrt(1 + pow(tanL, 2))) / omega) * (phii - Phi0 - (omega * M_PI) / (2 * fabs(omega)));

    gsl_matrix_set(J, i + 2 * n, 0, 1.0);
    gsl_matrix_set(J, i + 2 * n, 1, 0);
    gsl_matrix_set(J, i + 2 * n, 2, 0);
    gsl_matrix_set(J, i + 2 * n, 3, 0);
    gsl_matrix_set(J, i + 2 * n, 4, si / sqrt(1 + pow(tanL, 2)) - (pow(tanL, 2) * si) / sqrt(pow(1 + pow(tanL, 2), 3)));
  }

  return GSL_SUCCESS;
}

int fdfParametrisation3(const gsl_vector* par, void* d, gsl_vector* f, gsl_matrix* J) {
  //     For helix fitting
  functParametrisation3(par, d, f);
  dfunctParametrisation3(par, d, J);

  return GSL_SUCCESS;
}

ClusterShapes::ClusterShapes(int nhits, float* a, float* x, float* y, float* z)
    :

      m_nHits(nhits), m_aHit(nhits, 0.0), m_xHit(nhits, 0.0), m_yHit(nhits, 0.0), m_zHit(nhits, 0.0),
      m_exHit(nhits, 1.0), m_eyHit(nhits, 1.0), m_ezHit(nhits, 1.0),
      // all hits are assumed to be "cylindrical"
      m_ifNotGravity(1), m_ifNotInertia(1) {
  for (int i = 0; i < nhits; ++i) {
    m_aHit[i] = a[i];
    m_xHit[i] = x[i];
    m_yHit[i] = y[i];
    m_zHit[i] = z[i];
  }
}

float* ClusterShapes::getCentreOfGravity() {
  if (m_ifNotGravity == 1)
    findGravity();
  return &m_analogGravity[0];
}

float* ClusterShapes::getEigenVecInertia() {
  if (m_ifNotInertia == 1)
    findInertia();
  return &m_VecAnalogInertia[0];
}
float* ClusterShapes::getEigenVecInertiaErrors() {
  // this is a pure dummy to allow MarlinPandora development!
  if (m_ifNotInertia == 1)
    findInertia();
  return &m_VecAnalogInertia[0];
}

int ClusterShapes::FitHelix(int max_iter, int status_out, int parametrisation, float* parameter, float* dparameter,
                            float& chi2, float& distmax, int direction) {
  // Modified by Hengne Li @ LAL

  double parameterdb[5];
  double dparameterdb[5];
  double chi2db;
  double distmaxdb;
  for (int i = 0; i < 5; i++) {
    parameterdb[i] = double(parameter[i]);
    dparameterdb[i] = double(dparameter[i]);
  }
  chi2db = double(chi2);
  distmaxdb = double(distmax);

  int returnvalue =
      FitHelix(max_iter, status_out, parametrisation, parameterdb, dparameterdb, chi2db, distmaxdb, direction);

  for (int i = 0; i < 5; i++) {
    parameter[i] = float(parameterdb[i]);
    dparameter[i] = float(dparameterdb[i]);
  }
  chi2 = float(chi2db);
  distmax = float(distmaxdb);

  return returnvalue;
}

int ClusterShapes::FitHelix(int max_iter, int status_out, int parametrisation, double* parameter, double* dparameter,
                            double& chi2, double& distmax, int direction) {
  // FIXME: version with double typed parameters needed 2006/06/10 OW

  if (m_nHits < 3) {
    std::cout << "ClusterShapes : helix fit impossible, two few points";
    std::cout << std::endl;
    for (int i = 0; i < 5; ++i) {
      parameter[i] = 0.0;
      dparameter[i] = 0.0;
    }
    return 1;
  }

  // find initial parameters

  double Rmin = 1.0e+10;
  double Rmax = -1.0;
  int i1 = -1;

  // 1st loop
  for (int i = 0; i < m_nHits; ++i) {
    double Rz = sqrt(m_xHit[i] * m_xHit[i] + m_yHit[i] * m_yHit[i]);
    if (Rz < Rmin) {
      Rmin = Rz;
      i1 = i;
    }
    if (Rz > Rmax) {
      Rmax = Rz;
    }
  }

  // debug
  /*
    for (int i = 0; i < m_nHits; ++i) std::cout << i << "  " << m_xHit[i] << "  " << m_yHit[i] << "  " << m_zHit[i] <<
    std::endl; std::cout << std::endl << Rmin << "  " <<  Rmax << "  " << i1 << std::endl;
  */

  // 2nd loop
  double Upper = Rmin + 1.1 * (Rmax - Rmin);
  double Lower = Rmin + 0.9 * (Rmax - Rmin);
  double dZmin = 1.0e+20;

  int i3 = -1;

  for (int i = 0; i < m_nHits; ++i) {
    double Rz = sqrt(m_xHit[i] * m_xHit[i] + m_yHit[i] * m_yHit[i]);
    if ((Rz > Lower) && (Rz < Upper)) {
      double dZ = fabs(m_zHit[i] - m_zHit[i1]);
      if (dZ < dZmin) {
        dZmin = dZ;
        i3 = i;
      }
    }
  }

  // debug
  // std::cout << std::endl << Upper << "  " << Lower << "  " << dZmin << "  " << i3 << std::endl;

  double z1 = std::min(m_zHit[i1], m_zHit[i3]);
  double z3 = std::max(m_zHit[i1], m_zHit[i3]);

  int i2 = -1;
  double dRmin = 1.0e+20;
  double Rref = 0.5 * (Rmax + Rmin);

  // 3d loop

  for (int i = 0; i < m_nHits; ++i) {
    if (m_zHit[i] >= z1 && m_zHit[i] <= z3) {
      double Rz = sqrt(m_xHit[i] * m_xHit[i] + m_yHit[i] * m_yHit[i]);
      double dRz = fabs(Rz - Rref);
      if (dRz < dRmin) {
        i2 = i;
        dRmin = dRz;
      }
    }
  }

  // int problematic = 0;

  if (i2 < 0 || i2 == i1 || i2 == i3) {
    // problematic = 1;
    //  std::cout << "here we are " << std::endl;
    for (int i = 0; i < m_nHits; ++i) {
      if (i != i1 && i != i3) {
        i2 = i;
        if (m_zHit[i2] < z1) {
          int itemp = i1;
          i1 = i2;
          i2 = itemp;
        } else if (m_zHit[i2] > z3) {
          int itemp = i3;
          i3 = i2;
          i2 = itemp;
        }
        break;
      }
    }
    // std::cout << i1 << " " << i2 << " " << i3 << std::endl;
  }

  double x0 = 0.5 * (m_xHit[i2] + m_xHit[i1]);
  double y0 = 0.5 * (m_yHit[i2] + m_yHit[i1]);
  double x0p = 0.5 * (m_xHit[i3] + m_xHit[i2]);
  double y0p = 0.5 * (m_yHit[i3] + m_yHit[i2]);
  double ax = m_yHit[i2] - m_yHit[i1];
  double ay = m_xHit[i1] - m_xHit[i2];
  double axp = m_yHit[i3] - m_yHit[i2];
  double ayp = m_xHit[i2] - m_xHit[i3];
  double det = ax * ayp - axp * ay;
  double time;

  if (det == 0.) {
    time = 500.;
  } else {
    gsl_matrix* A = gsl_matrix_alloc(2, 2);
    gsl_vector* B = gsl_vector_alloc(2);
    gsl_vector* T = gsl_vector_alloc(2);
    gsl_matrix_set(A, 0, 0, ax);
    gsl_matrix_set(A, 0, 1, -axp);
    gsl_matrix_set(A, 1, 0, ay);
    gsl_matrix_set(A, 1, 1, -ayp);
    gsl_vector_set(B, 0, x0p - x0);
    gsl_vector_set(B, 1, y0p - y0);
    gsl_linalg_HH_solve(A, B, T);
    time = gsl_vector_get(T, 0);
    gsl_matrix_free(A);
    gsl_vector_free(B);
    gsl_vector_free(T);
  }

  double X0 = x0 + ax * time;
  double Y0 = y0 + ay * time;

  double dX = m_xHit[i1] - X0;
  double dY = m_yHit[i1] - Y0;

  double R0 = sqrt(dX * dX + dY * dY);

  /*
    if (problematic == 1) {
    std::cout << i1 << " " << i2 << " " << i3 << std::endl;
    std::cout << m_xHit[i1] << " " << m_yHit[i1] << " " << m_zHit[i1] << std::endl;
    std::cout << m_xHit[i2] << " " << m_yHit[i2] << " " << m_zHit[i2] << std::endl;
    std::cout << m_xHit[i3] << " " << m_yHit[i3] << " " << m_zHit[i3] << std::endl;
    std::cout << "R0 = " << R0 << std::endl;
    }
  */

  double phi1 = (double)atan2(m_yHit[i1] - Y0, m_xHit[i1] - X0);
  double phi2 = (double)atan2(m_yHit[i2] - Y0, m_xHit[i2] - X0);
  double phi3 = (double)atan2(m_yHit[i3] - Y0, m_xHit[i3] - X0);

  // testing bz > 0 hypothesis

  if (phi1 > phi2)
    phi2 = phi2 + 2.0 * M_PI;
  if (phi1 > phi3)
    phi3 = phi3 + 2.0 * M_PI;
  if (phi2 > phi3)
    phi3 = phi3 + 2.0 * M_PI;

  double bz_plus = (phi3 - phi1) / (m_zHit[i3] - m_zHit[i1]);
  double phi0_plus = phi1 - bz_plus * m_zHit[i1];
  double dphi_plus = fabs(bz_plus * m_zHit[i2] + phi0_plus - phi2);

  // testing bz < 0 hypothesis

  phi1 = (double)atan2(m_yHit[i1] - Y0, m_xHit[i1] - X0);
  phi2 = (double)atan2(m_yHit[i2] - Y0, m_xHit[i2] - X0);
  phi3 = (double)atan2(m_yHit[i3] - Y0, m_xHit[i3] - X0);

  if (phi1 < phi2)
    phi2 = phi2 - 2.0 * M_PI;
  if (phi1 < phi3)
    phi3 = phi3 - 2.0 * M_PI;
  if (phi2 < phi3)
    phi3 = phi3 - 2.0 * M_PI;

  double bz_minus = (phi3 - phi1) / (m_zHit[i3] - m_zHit[i1]);
  double phi0_minus = phi1 - bz_minus * m_zHit[i1];
  double dphi_minus = fabs(bz_minus * m_zHit[i2] + phi0_minus - phi2);

  double bz;
  double phi0;

  if (dphi_plus < dphi_minus) {
    bz = bz_plus;
    phi0 = phi0_plus;
  } else {
    bz = bz_minus;
    phi0 = phi0_minus;
  }

  double par_init[5];

  if (parametrisation == 1) {
    par_init[0] = (double)X0;
    par_init[1] = (double)Y0;
    par_init[2] = (double)R0;
    par_init[3] = (double)bz;
    par_init[4] = (double)phi0;
  } else if (parametrisation == 2) {
    par_init[0] = (double)X0;
    par_init[1] = (double)Y0;
    par_init[2] = (double)(-phi0 / bz);
    par_init[3] = (double)R0;
    par_init[4] = (double)(1 / bz);
  }

  else if (parametrisation == 3) { // parameter vector: (z0,Phi0,omega,d0,tanL)

    // debug
    // std::cout << std::setprecision(6) << "InitFit (X0,Y0,R0,bz,phi0) = " << "(" << X0 << "," << Y0 << "," << R0 <<
    // "," << bz << "," << phi0 << ")" << std::endl;

    // debug
    /*
          X0 = -1205.28;
          Y0 = 175.317;
          R0 = 1217.97;
          bz = 0.00326074;
          phi0 = -0.144444;
    */

    // debug
    // std::cout << std::setprecision(6) << "InitUsed (X0,Y0,R0,bz,phi0) = " << "(" << X0 << "," << Y0 << "," << R0 <<
    // "," << bz << "," << phi0 << ")" << std::endl;

    double omega = 1 / R0 * direction;
    double tanL = (-1.0) * omega / bz;

    double Phi0 = (-1.0) * atan2(X0, Y0);

    if (direction == 1) {
      if (tanL >= 0.0)
        Phi0 += M_PI; // add pi (see LC-DET-2006-004) //  >= or > ?
      else
        Phi0 -= M_PI; // < or <= ?
    }

    // double d0 = R0 - X0/sin(Phi0);
    // double d0 = R0 + Y0/cos(Phi0);

    double d0 = 0.0;
    if (true /*direction != 1*/)
      d0 = R0 - sqrt(X0 * X0 + Y0 * Y0);
    // else d0 = R0 + sqrt(X0*X0 + Y0*Y0);

    // double d0 = R0 - ( (X0-Y0)/(sqrt(2.0)*cos(pi/4 - Phi0)) );
    // double d0 = R0 - ((X0-Y0)/(sin(Phi0)+cos(Phi0)));

    // double Phi0 = asin(X0/(R0-d0));

    double z0 = (1 / bz) * ((-1.0) * phi0 + Phi0 + (omega * M_PI) / (2.0 * fabs(omega)));

    // debug
    /*
      std::cout << std::setprecision(6) << "InitFitCalculated (d0,z0,phi0,omega,tanL) = " << "(" << d0 << "," << z0 <<
      "," << Phi0 << "," << omega << "," << tanL << ")"
      << "  " << "sign(omega) = " << direction << std::endl;
    */

    // debug
    /*
      d0    = 0.00016512;
      z0    = 0.000853511;
      Phi0  = 1.11974;
      omega = -4.22171e-05;
      tanL  = -0.33436;
    */

    // debug
    // std::cout << std::setprecision(6) << "InitFitUsed (d0,z0,phi0,omega,tanL) = " << "(" << d0 << "," << z0 << "," <<
    // Phi0 << "," << omega << "," << tanL << ")" << std::endl;

    par_init[0] = z0;
    par_init[1] = Phi0;
    par_init[2] = omega;
    par_init[3] = d0;
    par_init[4] = tanL;
  } else
    return 1;

  // local variables
  int status = 0;
  int iter = 0;

  int npar = 5; // five parameters to fit
  int ndim = 0;
  if (parametrisation == 1)
    ndim = 2; // two dependent dimensions
  else if (parametrisation == 2)
    ndim = 3; // three dependent dimensions
  else if (parametrisation == 3)
    ndim = 3; // three dependent dimensions

  else
    return 1;

  double chi2_nofit = 0.0;
  int iFirst = 1;
  for (int ipoint = 0; ipoint < m_nHits; ipoint++) {
    std::array<double, 2> distRPZ;
    double Dist = DistanceHelix(m_xHit[ipoint], m_yHit[ipoint], m_zHit[ipoint], X0, Y0, R0, bz, phi0, distRPZ);
    double chi2rphi = distRPZ[0] / m_exHit[ipoint];
    chi2rphi = chi2rphi * chi2rphi;
    double chi2z = distRPZ[1] / m_ezHit[ipoint];
    chi2z = chi2z * chi2z;
    chi2_nofit = chi2_nofit + chi2rphi + chi2z;
    if (Dist > distmax || iFirst == 1) {
      distmax = Dist;
      iFirst = 0;
    }
  }
  chi2_nofit = chi2_nofit / double(m_nHits);

  if (status_out == 1) {
    for (int i = 0; i < 5; ++i) {
      parameter[i] = (double)par_init[i];
      dparameter[i] = 0.0;
    }
    chi2 = chi2_nofit;
    return 0;
  }

  // converging criteria
  const double abs_error = 1e-4;
  const double rel_error = 1e-4;

  gsl_multifit_function_fdf fitfunct;

  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;

  gsl_multifit_fdfsolver* s = gsl_multifit_fdfsolver_alloc(T, ndim * m_nHits, npar);

  gsl_matrix* covar = gsl_matrix_alloc(npar, npar); // covariance matrix

  data d;
  d.n = m_nHits;
  d.x = &m_xHit[0];
  d.y = &m_yHit[0];
  d.z = &m_zHit[0];
  d.ex = &m_exHit[0];
  d.ey = &m_eyHit[0];
  d.ez = &m_ezHit[0];

  if (parametrisation == 1) {
    fitfunct.f = &functParametrisation1;
    fitfunct.df = &dfunctParametrisation1;
    fitfunct.fdf = &fdfParametrisation1;
    fitfunct.n = ndim * m_nHits;
    fitfunct.p = npar;
    fitfunct.params = &d;
  } else if (parametrisation == 2) {
    fitfunct.f = &functParametrisation2;
    fitfunct.df = &dfunctParametrisation2;
    fitfunct.fdf = &fdfParametrisation2;
    fitfunct.n = ndim * m_nHits;
    fitfunct.p = npar;
    fitfunct.params = &d;
  } else if (parametrisation == 3) {
    fitfunct.f = &functParametrisation3;
    fitfunct.df = &dfunctParametrisation3;
    fitfunct.fdf = &fdfParametrisation3;
    fitfunct.n = ndim * m_nHits;
    fitfunct.p = npar;
    fitfunct.params = &d;
  } else
    return 1;

  gsl_vector_view pinit = gsl_vector_view_array(par_init, npar);
  gsl_multifit_fdfsolver_set(s, &fitfunct, &pinit.vector);

  // perform fit
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(s);

    if (status)
      break;
    status = gsl_multifit_test_delta(s->dx, s->x, abs_error, rel_error);

  } while (status == GSL_CONTINUE && iter < max_iter);

  // fg: jacobian has been dropped from gsl_multifit_fdfsolver in gsl 2:
  gsl_matrix* J = gsl_matrix_alloc(s->fdf->n, s->fdf->p);
  gsl_multifit_fdfsolver_jac(s, J);
  gsl_multifit_covar(J, rel_error, covar);
  //  gsl_multifit_covar (s->J, rel_error, covar);

  chi2 = 0.0;

  if (parametrisation == 1) {
    X0 = (double)gsl_vector_get(s->x, 0);
    Y0 = (double)gsl_vector_get(s->x, 1);
    R0 = (double)gsl_vector_get(s->x, 2);
    bz = (double)gsl_vector_get(s->x, 3);
    phi0 = (double)gsl_vector_get(s->x, 4);
  } else if (parametrisation == 2) {
    X0 = (double)gsl_vector_get(s->x, 0);
    Y0 = (double)gsl_vector_get(s->x, 1);
    R0 = (double)gsl_vector_get(s->x, 3);
    bz = (double)(1 / gsl_vector_get(s->x, 4));
    phi0 = (double)(-gsl_vector_get(s->x, 2) / gsl_vector_get(s->x, 4));
  } else if (parametrisation == 3) { // (parameter vector: (z0,phi0,omega,d0,tanL)

    double z0 = gsl_vector_get(s->x, 0);
    double Phi0 = gsl_vector_get(s->x, 1);
    double omega = gsl_vector_get(s->x, 2);
    double d0 = gsl_vector_get(s->x, 3);
    double tanL = gsl_vector_get(s->x, 4);

    X0 = (double)(((1 / omega) - d0) * sin(Phi0));
    Y0 = (double)((-1) * ((1 / omega) - d0) * cos(Phi0));
    R0 = (double)(1 / fabs(omega));
    bz = (double)((-1) * (omega / tanL));
    phi0 = (double)(((z0 * omega) / tanL) + Phi0 + ((omega * M_PI) / (2 * fabs(omega))));
  } else
    return 1;

  iFirst = 1;
  double ddmax = 0.0;
  for (int ipoint = 0; ipoint < m_nHits; ipoint++) {
    std::array<double, 2> distRPZ;
    double Dist = DistanceHelix(m_xHit[ipoint], m_yHit[ipoint], m_zHit[ipoint], X0, Y0, R0, bz, phi0, distRPZ);
    double chi2rphi = distRPZ[0] / m_exHit[ipoint];
    chi2rphi = chi2rphi * chi2rphi;
    double chi2z = distRPZ[1] / m_ezHit[ipoint];
    chi2z = chi2z * chi2z;
    chi2 = chi2 + chi2rphi + chi2z;
    if (Dist > ddmax || iFirst == 1) {
      iFirst = 0;
      ddmax = Dist;
    }
  }

  chi2 = chi2 / double(m_nHits);
  if (chi2 < chi2_nofit) {
    for (int i = 0; i < npar; i++) {
      parameter[i] = gsl_vector_get(s->x, i);
      dparameter[i] = sqrt(gsl_matrix_get(covar, i, i));
    }
    distmax = ddmax;
  } else {
    chi2 = chi2_nofit;
    for (int i = 0; i < npar; i++) {
      parameter[i] = (double)par_init[i];
      dparameter[i] = 0.0;
    }
  }

  //  if (problematic == 1)
  //    std::cout << "chi2 = " << chi2 << std::endl;

  gsl_multifit_fdfsolver_free(s);
  gsl_matrix_free(covar);
  return 0;
}

// ##########################################
// #####                                #####
// #####        private methods         #####
// #####                                #####
// ##########################################

void ClusterShapes::findElipsoid() {
  /**   Elipsoid parameter calculations see cluster_proper.f  */
  float cx, cy, cz;
  float dx, dy, dz;
  float r_hit_max, d_begn, d_last, r_max, proj;
  if (m_ifNotInertia == 1)
    findInertia();

  // Find Minumal and Maximal Lenght for Principal axis
  r_hit_max = -100000.;
  d_begn = 100000.;
  d_last = -100000.;
  cx = m_VecAnalogInertia[0];
  cy = m_VecAnalogInertia[1];
  cz = m_VecAnalogInertia[2];
  for (int i = 0; i < m_nHits; ++i) {
    dx = m_xHit[i] - m_xgr;
    dy = m_yHit[i] - m_ygr;
    dz = m_zHit[i] - m_zgr;
    r_max = sqrt(dx * dx + dy * dy + dz * dz);
    if (r_max > r_hit_max)
      r_hit_max = r_max;
    proj = dx * cx + dy * cy + dz * cz;
    if (proj < d_begn)
      d_begn = proj;
    //            lad_begn = ladc(L)
    if (proj > d_last)
      d_last = proj;
    //            lad_last = ladc(L)
  }
  //        if (r_hit_max > 0.0)
  //	  m_r1 = 1.05*r_hit_max; // + 5% of length
}

void ClusterShapes::findGravity() {
  m_totAmpl = 0.;
  for (int i = 0; i < 3; ++i) {
    m_analogGravity[i] = 0.0;
  }
  for (int i = 0; i < m_nHits; ++i) {
    m_totAmpl += m_aHit[i];
    m_analogGravity[0] += m_aHit[i] * m_xHit[i];
    m_analogGravity[1] += m_aHit[i] * m_yHit[i];
    m_analogGravity[2] += m_aHit[i] * m_zHit[i];
  }
  for (int i = 0; i < 3; ++i) {
    m_analogGravity[i] /= m_totAmpl;
  }
  m_xgr = m_analogGravity[0];
  m_ygr = m_analogGravity[1];
  m_zgr = m_analogGravity[2];
  m_ifNotGravity = 0;
}

void ClusterShapes::findInertia() {
  double aIne[3][3];
  //  float radius1;
  float radius2 = 0.0;

  findGravity();

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      aIne[i][j] = 0.0;
    }
  }

  for (int i = 0; i < m_nHits; ++i) {
    float dX = m_xHit[i] - m_analogGravity[0];
    float dY = m_yHit[i] - m_analogGravity[1];
    float dZ = m_zHit[i] - m_analogGravity[2];
    aIne[0][0] += m_aHit[i] * (dY * dY + dZ * dZ);
    aIne[1][1] += m_aHit[i] * (dX * dX + dZ * dZ);
    aIne[2][2] += m_aHit[i] * (dX * dX + dY * dY);
    aIne[0][1] -= m_aHit[i] * dX * dY;
    aIne[0][2] -= m_aHit[i] * dX * dZ;
    aIne[1][2] -= m_aHit[i] * dY * dZ;
  }

  for (int i = 0; i < 2; ++i) {
    for (int j = i + 1; j < 3; ++j) {
      aIne[j][i] = aIne[i][j];
    }
  }
  //****************************************
  // analog Inertia
  //****************************************

  gsl_matrix_view aMatrix = gsl_matrix_view_array((double*)aIne, 3, 3);
  gsl_vector* aVector = gsl_vector_alloc(3);
  gsl_matrix* aEigenVec = gsl_matrix_alloc(3, 3);
  gsl_eigen_symmv_workspace* wa = gsl_eigen_symmv_alloc(3);
  gsl_eigen_symmv(&aMatrix.matrix, aVector, aEigenVec, wa);
  gsl_eigen_symmv_free(wa);
  gsl_eigen_symmv_sort(aVector, aEigenVec, GSL_EIGEN_SORT_ABS_ASC);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      m_VecAnalogInertia[i + 3 * j] = gsl_matrix_get(aEigenVec, i, j);
    }
  }

  // Main principal points away from IP

  m_radius = 0.;
  radius2 = 0.;

  for (int i = 0; i < 3; ++i) {
    m_radius += m_analogGravity[i] * m_analogGravity[i];
    radius2 += (m_analogGravity[i] + m_VecAnalogInertia[i]) * (m_analogGravity[i] + m_VecAnalogInertia[i]);
  }

  if (radius2 < m_radius) {
    for (int i = 0; i < 3; ++i)
      m_VecAnalogInertia[i] = -m_VecAnalogInertia[i];
  }

  m_radius = sqrt(m_radius);
  m_ifNotInertia = 0;

  // The final job
  findWidth();
  findElipsoid();

  gsl_vector_free(aVector);
  gsl_matrix_free(aEigenVec);
}

void ClusterShapes::findWidth() {
  if (m_ifNotInertia == 1)
    findInertia();
}

/**
   Function sdist(xp,yp,zp,cx,cy,cz,xv,yv,zv)
   c----------------------------------------------------------------------
   c        Distance from line to point
   c       xp, yp, zp -- point is at the line
   c       xv, yv, zv -- point is out of line
   ********************************************************************
   *     Last update     V.L.Morgunov     08-Apr-2002                 *
   ********************************************************************
   real xp,yp,zp,cx,cy,cz,xv,yv,zv,t1,t2,t3,tt,sdist

   t1 = cy*(zp-zv)-cz*(yp-yv)
   t2 = cz*(xp-xv)-cx*(zp-zv)
   t3 = cx*(yp-yv)-cy*(xp-xv)
   tt = sqrt(cx**2+cy**2+cz**2)
   sdist = sqrt(t1**2+t2**2+t3**2)/tt

   return
   end
*/

double ClusterShapes::DistanceHelix(double x, double y, double z, double X0, double Y0, double R0, double bz,
                                    double phi0, std::array<double, 2>& distRPZ) {
  double phi = atan2(y - Y0, x - X0);
  double R = sqrt((y - Y0) * (y - Y0) + (x - X0) * (x - X0));
  double dXY2 = (R - R0) * (R - R0);
  double _const_2pi = 2.0 * M_PI;
  double xN = (bz * z + phi0 - phi) / _const_2pi;

  int n1 = 0;
  int n2 = 0;
  int nSpirals = 0;

  if (xN > 0) {
    n1 = (int)xN;
    n2 = n1 + 1;
  } else {
    n1 = (int)xN - 1;
    n2 = n1 + 1;
  }

  if (fabs(n1 - xN) < fabs(n2 - xN)) {
    nSpirals = n1;
  } else {
    nSpirals = n2;
  }

  double dZ = (phi + _const_2pi * nSpirals - phi0) / bz - z;
  double dZ2 = dZ * dZ;

  distRPZ[0] = sqrt(dXY2);
  distRPZ[1] = sqrt(dZ2);

  return sqrt(dXY2 + dZ2);
}
