#ifndef K4GAUDIPANDORA_CLUSTERSHAPES_H
#define K4GAUDIPANDORA_CLUSTERSHAPES_H

#include <array>
#include <vector>

/**
 *    Utility class to derive properties of clusters, such as centre of gravity,
 *    axes of inertia, fits of the cluster shape and so on. All the details are
 *    explained in the documentation of the methods. Several classes of the GSL
 *    (GNU Scientific Library) are needed in this class.
 *
 *    @authors V. Morgunov (ITEP/DESY), A. Raspereza (DESY), O. Wendt (DESY)
 *    @version $Id: ClusterShapes.h,v 1.14 2007-04-27 13:56:53 owendt Exp $
 *
 */
class ClusterShapes {
public:
  /**
   *    Constructor
   *    @param nhits : number of hits in the cluster
   *    @param a     : amplitudes of elements ('cells') of the cluster. Stored in
   *                   an array, with one entry for each element ('cell'). Each entry
   *                   is depending on coordinates x,y,z (Cartesian), which are stored
   *                   in the arrays x,y,z.
   *    @param x,y,z : array of coordinates corresponding to the array of amplitudes a.
   *
   *
   */
  ClusterShapes(int nhits, float* a, float* x, float* y, float* z);

  /**
   * returns an array, which represents a vector from the origin of the
   * coordinate system, i.\ e.\ IP, to the centre of gravity of the cluster. The centre
   * of gravity is calculated with the energy of the entries of the cluster.
   */
  float* getCentreOfGravity();

  /**
   * array of the three main axes of inertia (9 entries) starting
   * with the axis corresponding to the smallest inertia of mass
   * (main principal axis). All axes are normalised to a length
   * of 1.
   */
  float* getEigenVecInertia();
  // this is (for now) a pure dummy to allow MarlinPandora development!
  float* getEigenVecInertiaErrors();

  /**
   * performs a least square fit on a helix path in space, which
   * which is defined as (Cartesian coordiantes):
   *
   * 1. parametrisation:
   * x[i] = x0 + R*cos(b*z[i] + phi0)
   * y[i] = y0 + R*sin(b*z[i] + phi0)
   * z[i] = z[i],
   * where x0,y0,R,b and phi0 are the parameters to be fitted and
   * x[i],y[i],z[i] are the (Cartesian) coordiantes of the space
   * points.
   *
   * 2. parametrisation:
   * x[i] = x0 + R*cos(phi)
   * y[i] = y0 + R*sin(phi)
   * z[i] = z0 + b*phi
   * and phi = atan2( y[i]-y0 , x[i]-x0 ),
   * where x0,y0,z0,R and b are the parameters to be fitted and
   * x[i],y[i],z[i] are the (Cartesian) coordiantes of the space
   * points.
   *
   * The method returns 1 if an error occured and 0 if not.
   *
   * The following output/input parameters are returned/needed:
   *
   * OUTPUTS:
   * @param parameter     : array of parameters to be fitted.
   *                        For parametrisation 1: parameter[5] = {x0,y0,R,b,phi0}
   *                        For parametrisation 2: parameter[5] = {x0,y0,z0,R,b}
   * @param dparameter    : error on the parameters, that means:
   *                        dparameter[i] = sqrt( CovarMatrix[i][i] )
   * @param chi2          : chi2 of the fit
   * @param distmax       : maximal distance between the points x[i],y[i]
   *                        z[i] and the fitted function
   *
   * INPUTS:
   * @param parametrisation : 1 for first and 2 for second parametrisation (see above)
   * @param max_iter        : maximal number of iterations, before fit cancels
   * @param status_out      : if set to 1, only the initial parameters of
   *                          the fit are calculated and are stored in
   *                          parameter. The entries of dparameter are
   *                          set to 0.0
   */
  int FitHelix(int max_iter, int status_out, int parametrisation, double* parameter, double* dparameter, double& chi2,
               double& distmax, int direction = 1);

  int FitHelix(int max_iter, int status_out, int parametrisation, float* parameter, float* dparameter, float& chi2,
               float& distmax, int direction = 1);

private:
  int m_nHits;

  std::vector<float> m_aHit;
  std::vector<float> m_xHit;
  std::vector<float> m_yHit;
  std::vector<float> m_zHit;
  std::vector<float> m_exHit;
  std::vector<float> m_eyHit;
  std::vector<float> m_ezHit;
  std::vector<float> m_xl;
  std::vector<float> m_xt;
  std::vector<float> m_t;
  std::vector<float> m_s;
  std::vector<int> m_types;

  int m_ifNotGravity = 1;
  float m_totAmpl = 0.0;
  float m_radius = 0.0;
  float m_xgr = 0.0;
  float m_ygr = 0.0;
  float m_zgr = 0.0;
  std::array<float, 3> m_analogGravity{0.0, 0.0, 0.0};

  int m_ifNotInertia = 1;
  float m_VecAnalogInertia[9];

  void findElipsoid();
  void findGravity();
  void findInertia();
  void findWidth();
  double DistanceHelix(double x, double y, double z, double X0, double Y0, double R0, double bz, double phi0,
                       double* distRPhiZ);

  // private methods for non-linear, multidim. fitting (helix)
  // static int functParametrisation1(const gsl_vector* par, void* data, gsl_vector* f);
  // static int dfunctParametrisation1(const gsl_vector* par, void* d, gsl_matrix* J);
  // static int fdfParametrisation1(const gsl_vector* par, void* d, gsl_vector* f, gsl_matrix* J);
};

#endif
