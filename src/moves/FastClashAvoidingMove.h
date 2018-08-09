/*

Excited States software: KGS
Contributors: See CONTRIBUTORS.txt
Contact: kgs-contact@simtk.org

Copyright (C) 2009-2017 Stanford University

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

This entire text, including the above copyright notice and this permission notice
shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.

*/
//
// Created by Dominik Budday on 15.03.16.
//

#ifndef KGS_FASTCLASHAVOIDINGMOVE_H
#define KGS_FASTCLASHAVOIDINGMOVE_H

#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "metrics/Metric.h"
#include "metrics/RMSD.h"
#include "moves/Move.h"
#include "core/Configuration.h"

class FastClashAvoidingMove : public Move
{
 public:
  FastClashAvoidingMove( double maxRotation,
                         int trialSteps,
                         const std::string& atomTypes,
                         bool projectConstraints
  );

 protected:
  Configuration* performMove(Configuration* current, gsl_vector* gradient);

 private:
  Configuration* projectOnClashNullspace(Configuration *conf,
                                      gsl_vector *gradient,
                                      std::set<std::pair<Atom *, Atom *> > &collisions);

  gsl_matrix* computeClashAvoidingJacobian( Configuration* conf,
                                            std::map<int,int>& dofMap,
                                            std::set< std::pair<Atom*,Atom*> >& collisions);

  /** Return a map that associates cycle-dofs and constrained dofs with a general dofs. */
  std::map<int,int> collectConstrainedDofMap(Configuration* conf, std::set< std::pair<Atom*,Atom*> >& allCollisions);
  const int m_trialSteps;
  const bool m_projectConstraints;
  const std::string m_collisionCheckAtomTypes;
};


#endif //KGS_FASTCLASHAVOIDINGMOVE_H
