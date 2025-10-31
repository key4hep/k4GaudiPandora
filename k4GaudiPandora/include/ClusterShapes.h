/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
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

private:
  int m_nHits;

  std::vector<float> m_aHit;
  std::vector<float> m_xHit;
  std::vector<float> m_yHit;
  std::vector<float> m_zHit;
  std::vector<float> m_exHit;
  std::vector<float> m_eyHit;
  std::vector<float> m_ezHit;

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
};

#endif
