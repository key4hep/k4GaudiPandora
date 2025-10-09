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
