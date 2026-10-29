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

#include "PfoCreatorIdea.h"

#include <edm4hep/CalorimeterHit.h>
#include <edm4hep/MutableCluster.h>
#include <edm4hep/MutableReconstructedParticle.h>

#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Track.h"

#include <cmath>
#include <limits>

PfoCreatorIdea::PfoCreatorIdea(const Settings& settings, pandora::Pandora& pandora, const Gaudi::Algorithm* algorithm)
    : PfoCreatorBase(pandora, algorithm), m_settings(settings) {}

pandora::StatusCode
PfoCreatorIdea::CreateParticleFlowObjects(edm4hep::ClusterCollection& clusterColl,
                                          edm4hep::ReconstructedParticleCollection& aPfoColl) const {
  const pandora::PfoList* pandoraPfoList = nullptr;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(m_pandora, pandoraPfoList))

  for (const auto* aPfo : *pandoraPfoList) {
    auto reconstructedParticle = aPfoColl.create();

    // No reference point and no start vertex: IDEA uses its own vertex finder and fitter, not
    // pandora's track-based reference point.
    this->AddTracksToRecoParticle(aPfo, reconstructedParticle);
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                             this->AddClustersToRecoParticle(aPfo, clusterColl, reconstructedParticle))
    this->SetRecoParticlePropertiesFromPFO(aPfo, reconstructedParticle);
  } // loop pandora PFO list

  return pandora::STATUS_CODE_SUCCESS;
}

pandora::StatusCode
PfoCreatorIdea::AddClustersToRecoParticle(const pandora::ParticleFlowObject* const pPandoraPfo,
                                          edm4hep::ClusterCollection& clusterColl,
                                          edm4hep::MutableReconstructedParticle& reconstructedParticle) const {
  const pandora::ClusterList& clusterList(pPandoraPfo->GetClusterList());

  for (const auto* pPandoraCluster : clusterList) {
    pandora::CaloHitList pandoraCaloHitList;
    pPandoraCluster->GetOrderedCaloHitList().FillCaloHitList(pandoraCaloHitList);

    // isolated hits will only contribute to the energy
    // but not the position and covariance matrix
    const auto isolatedCaloHitList = pPandoraCluster->GetIsolatedCaloHitList();

    double hitE = 0., hitX = 0., hitY = 0., hitZ = 0.;
    auto cluster = clusterColl.create();

    for (const auto* pCaloHit : pandoraCaloHitList) {
      const auto* hit = static_cast<const edm4hep::CalorimeterHit*>(pCaloHit->GetParentAddress());
      cluster.addToHits(*hit);
      const auto pos = hit->getPosition();
      const double e = hit->getEnergy();
      hitE += e;
      hitX += pos[0] * e;
      hitY += pos[1] * e;
      hitZ += pos[2] * e;
    } // loop calo hits in cluster

    for (const auto* pIsoHit : isolatedCaloHitList) {
      cluster.addToHits(*static_cast<const edm4hep::CalorimeterHit*>(pIsoHit->GetParentAddress()));
    } // loop isolated hits in cluster

    // The cluster energy is the calorimeter's own measurement of this cluster: the dual-readout
    // combination E_DR = (S - chi C) / (1 - chi), summed over ECAL and HCAL.  That is what the
    // registered DualReadoutCorrection plugin returns, and it sits in the HADRONIC slot, so it is
    // taken unconditionally -- the electromagnetic slot has no plugin registered and would hand
    // back the bare scintillation-plus-Cherenkov sum, which is not an energy.  The pfo is free to
    // use a different estimator (photons the scintillation sum, charged pfos the track momentum);
    // that is a pfo-level choice and does not belong here.
    cluster.setEnergy(pPandoraCluster->GetCorrectedHadronicEnergy(m_pandora));

    // Propagate the per-hit energy errors through the same combination, ECAL and HCAL separately
    // because their chi differ.  Same hit accounting as DualReadoutCorrection: both readout
    // channels, isolated hits included.
    double varEcalS = 0., varEcalC = 0., varHcalS = 0., varHcalC = 0.;
    for (const pandora::CaloHitList& hitList : {pandoraCaloHitList, isolatedCaloHitList}) {
      for (const auto* pCaloHit : hitList) {
        const auto* hit = static_cast<const edm4hep::CalorimeterHit*>(pCaloHit->GetParentAddress());
        const double eErr = hit->getEnergyError();
        const double var = eErr * eErr;
        const bool isCherenkov = pCaloHit->GetHitType() == pandora::DRC_CHEREN;

        if (pCaloHit->GetElectromagneticEnergy() > 0.f) {
          if (isCherenkov)
            varEcalC += var;
          else
            varEcalS += var;
        } else {
          if (isCherenkov)
            varHcalC += var;
          else
            varHcalS += var;
        }
      }
    }
    const double ecalScale = 1. / (1. - m_settings.m_chiEcal);
    const double hcalScale = 1. / (1. - m_settings.m_chiHcal);
    const double varEcal = (varEcalS + m_settings.m_chiEcal * m_settings.m_chiEcal * varEcalC) * ecalScale * ecalScale;
    const double varHcal = (varHcalS + m_settings.m_chiHcal * m_settings.m_chiHcal * varHcalC) * hcalScale * hcalScale;
    cluster.setEnergyError(float(std::sqrt(varEcal + varHcal)));

    if (hitE <= std::numeric_limits<float>::epsilon()) {
      m_algorithm.warning() << "PfoCreatorIdea::AddClustersToRecoParticle: invalid cluster energy " << hitE << endmsg;
      return pandora::STATUS_CODE_FAILURE;
    }

    const double xBar = hitX / hitE, yBar = hitY / hitE, zBar = hitZ / hitE;
    cluster.setPosition({float(xBar), float(yBar), float(zBar)});

    // Covariance of the barycentre propagated from the per-hit energy errors.  With
    // xBar = sum(E_i x_i) / E_tot we have d(xBar)/dE_i = (x_i - xBar) / E_tot, hence
    // Cov_ab = sum_i sigma_Ei^2 (a_i - aBar)(b_i - bBar) / E_tot^2.  Subtracting the barycentre is
    // what makes it translation invariant; without it the covariance grows with the distance to
    // the origin.  This needs a second pass, since the barycentre is only known after the first.
    double cxx = 0., cxy = 0., cyy = 0., cxz = 0., cyz = 0., czz = 0.;
    for (const auto* pCaloHit : pandoraCaloHitList) {
      const auto* hit = static_cast<const edm4hep::CalorimeterHit*>(pCaloHit->GetParentAddress());
      const double eErr = hit->getEnergyError();
      const double var = eErr * eErr;
      const auto pos = hit->getPosition();
      const double dx = pos[0] - xBar, dy = pos[1] - yBar, dz = pos[2] - zBar;
      cxx += var * dx * dx;
      cxy += var * dx * dy;
      cyy += var * dy * dy;
      cxz += var * dx * dz;
      cyz += var * dy * dz;
      czz += var * dz * dz;
    }
    const double invE2 = 1. / (hitE * hitE);
    // edm4hep packs the covariance as the lower triangle: {xx, xy, yy, xz, yz, zz}
    cluster.setPositionError({float(cxx * invE2), float(cxy * invE2), float(cyy * invE2), float(cxz * invE2),
                              float(cyz * invE2), float(czz * invE2)});

    reconstructedParticle.addToClusters(cluster);
  } // loop clusters in PFO

  return pandora::STATUS_CODE_SUCCESS;
}
