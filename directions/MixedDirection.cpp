/*
    KGSX: Biomolecular Kino-geometric Sampling and Fitting of Experimental Data
    Yao et al, Proteins. 2012 Jan;80(1):25-43
    e-mail: latombe@cs.stanford.edu, vdbedem@slac.stanford.edu, julie.bernauer@inria.fr

        Copyright (C) 2011-2013 Stanford University

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


#include <gsl/gsl_vector_double.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include "MixedDirection.h"

MixedDirection::MixedDirection():
    m_totalWeight(0)
{ }

void MixedDirection::addDirection(Direction* dir, double weight)
{
  m_directions.push_back(dir);
  m_weights.push_back(weight);
  m_totalWeight+=weight;
}

void MixedDirection::computeGradient(Configuration* conf, Configuration* target, gsl_vector* ret)
{
  assert(!m_directions.empty());

  m_directions[0]->gradient(conf,target,ret);
  //std::cout<<"MixedDirection::computeGradient - rand norm: "<<gsl_blas_dnrm2(ret)<<std::endl;

  gsl_vector* tmp = gsl_vector_calloc(ret->size);
  for(size_t i=1;i<m_directions.size();i++){
    m_directions[i]->gradient(conf,target,tmp);
    gsl_vector_add(ret,tmp);
  }

}
