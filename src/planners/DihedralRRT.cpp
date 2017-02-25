/*KGSX: Biomolecular Kino-geometric Sampling and Fitting of Experimental Data*/

/* DihedralRRT.cpp
 *
 *  Created on: 19.08.2015
 *      Author: StDoBudd
 */

#include "DihedralRRT.h"

#include <iomanip>
#include <stack>
#include <IO.h>

#include "core/Molecule.h"
#include "core/Chain.h"
#include "Logger.h"
#include "core/Transformation.h"
#include "CTKTimer.h"
#include "metrics/RMSD.h"
#include "metrics/Dihedral.h"


using namespace std;

//const double RMSD_DIFF_THRESHOLD = 0.01;
const double MOV_DIH_THRESHOLD = 0.001;

extern double selectNodeTime;


DihedralRRT::DihedralRRT(
    Molecule *protein,
    Direction *direction,
    int numSamples,
    double maxDistance,
    bool sampleRandom
) :
    SamplingPlanner(),
    m_protein(protein),
    m_direction(direction),
    m_numSamples(numSamples),
    m_maxDistance(maxDistance),
    m_sampleRandom(sampleRandom){
  m_protein = protein;
  m_numDOFs = m_protein->m_spanningTree->getNumDOFs();

  Configuration *pSmp = new Configuration(m_protein);
  m_protein->m_conf = pSmp;
  m_target = nullptr;
  m_samples.push_back(pSmp);
  pSmp->m_vdwEnergy = 99999;
  pSmp->m_id = 0; // m_root

  m_deform_mag = 0.25;
  m_rand_radius = 2;

  m_top_min_rmsd = 99999;
  m_top_min_rmsd_id = -1;
  m_minMovDihDistance = 99999;
  m_minMovDihDistance_id = -1;
}

DihedralRRT::~DihedralRRT() {
  for (list<Configuration *>::iterator iter = m_samples.begin(); iter != m_samples.end(); iter++) {
    delete *iter;
  }
}

void DihedralRRT::generateSamples() {

  int sample_id = 0, max_depth = 0, failed_trials = 0, total_trials = 0;
  Configuration *pTarget = nullptr, *pClosestSmp, *pNewSmp = nullptr;
  gsl_vector *gradient = gsl_vector_alloc(m_protein->totalDofNum());

  //Save initial file (for movie)
  IO::writeNewSample(m_samples.front(), m_samples.front(), sample_id, m_workingDir, m_saveData);

  bool createNewTarget = false;

  log("debug")<<"Dihedral RRT, generating "<<m_numSamples<<" samples."<<endl;
  while (sample_id < m_numSamples) {
    ++total_trials;

    CTKTimer timer;
    double start_time = timer.getTimeNow();

    ///Sample Randomly or follow a path (if possible and target has not been reached)
    if (m_sampleRandom || pTarget == nullptr || createNewTarget) {
      log("debug") << "Generating new target, getting new seed" << endl;
      pTarget = GenerateRandConf(); //generates random conformation within sampling hyper-sphere
      createNewTarget = false;
      pClosestSmp = SelectClosestNode(pTarget); //pick closest existing node as seed
      double end_time = timer.getTimeNow();
      selectNodeTime += end_time - start_time;
    } else {
      log("debug") << "Using latest sample as seed" << endl;
      pClosestSmp = m_samples.back();
    }
    log("debug") << "Direction" << endl;
    log("debug") << "Target " << pTarget<< endl;
    m_direction->gradient(pClosestSmp, pTarget, gradient); ///desired direction depending on specified gradient
    log("debug") << "Move" << endl;
    pNewSmp = m_move->move(pClosestSmp, gradient); ///obtaining new conformation, depending on move and direction

    if (!pNewSmp->updatedMolecule()->inCollision()) { //non-colliding conformation
      log("debug") << "Valid conf" << endl;
      ++sample_id;

      pNewSmp->m_distanceToIni = m_metric->distance(pNewSmp, m_samples.front());
      pNewSmp->m_distanceToParent = m_metric->distance(pNewSmp, pClosestSmp);
      pNewSmp->m_id = sample_id;
      m_samples.push_back(pNewSmp);

      IO::writeNewSample(pNewSmp, m_samples.front(), sample_id, m_workingDir, m_saveData);

      if (pNewSmp->m_treeDepth > max_depth)
        max_depth = pNewSmp->m_treeDepth; ///keep track of tree-depth

      double distToRandGoal = m_metric->distance(pNewSmp, pTarget);

      log("samplingStatus") << "> New structure: " << m_protein->getName() << "_new_" << sample_id << ".pdb";
      log("samplingStatus") << " .. Distance to initial: " << setprecision(6) << pNewSmp->m_distanceToIni;
      log("samplingStatus") << " .. Distance to current target: " << setprecision(3) << distToRandGoal;
      log("samplingStatus") << " .. accessible dofs: " << pNewSmp->m_clashFreeDofs << endl << endl;

      if (distToRandGoal <= MOV_DIH_THRESHOLD) {//current target reached
        delete pTarget;
        createNewTarget = true;
      }
    } else {
      log("debug")<<"Dihedral RRT: Failed due to collision."<<endl;
      ++failed_trials;
      delete pTarget;
      createNewTarget = true;
    }
    log("debug") << "After sample" << endl;
  }

  gsl_vector_free(gradient);
}

Configuration *DihedralRRT::GenerateRandConf() {
  Configuration *pNewSmp = new Configuration(m_protein);

  double length = 0.0;
  for (int i = 0; i < m_numDOFs; ++i) {
    pNewSmp->m_dofs[i] = Math3D::dPi * RandomN1P1();
    length += pNewSmp->m_dofs[i] * pNewSmp->m_dofs[i];
  }
  length = sqrt(length);
//  if (ExploreOptions::getOptions()->scaleToRadius) {
    double factor = pow(Random01(), 1.0 / m_numDOFs) * m_maxDistance / length; ///ToDo: Do we still want this scaling?
    for (int i = 0; i < m_numDOFs; ++i) {
      pNewSmp->m_dofs[i] = factor * pNewSmp->m_dofs[i];
    }
//  }

  pNewSmp->m_id = -1;
  return pNewSmp;
}

Configuration *DihedralRRT::SelectClosestNode(Configuration *pTarget) {
  Configuration *pSmp, *pMinSmp;
  double min_distance = 1000000.0;
  double distance;
  pMinSmp = nullptr;

  for (list<Configuration *>::iterator iter = m_samples.begin(); iter != m_samples.end(); ++iter) {
    pSmp = *iter;
    distance = m_metric->distance(pSmp, pTarget);
    if (distance < min_distance) {
      min_distance = distance;
      pMinSmp = pSmp;
    }
  }

  return pMinSmp;
}



