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

#include <cassert>
#include "DOF.h"

DOF::DOF(const KinEdge* edge):
    m_edge(edge),
    m_value(0),
    m_index(-1),
    m_cycleIndex(-1)
{ }

DOF::~DOF(){};


void DOF::setValue(double val)
{
  m_value = val;
}

double DOF::getValue() const
{
  return m_value;
}

unsigned int DOF::getIndex() const
{
  return (unsigned int)m_index;
}

unsigned int DOF::getCycleIndex() const
{
  return (unsigned int)m_cycleIndex;
}

void DOF::setIndex(unsigned int idx)
{
  m_index = idx;
}

void DOF::setCycleIndex(unsigned int idx)
{
  m_cycleIndex = idx;
}

bool DOF::isDOFligand() const{
    if(m_edge->StartVertex->isligand() or m_edge->EndVertex->isligand()){
        return true;
    }
    else{
        return false;
    }
}

bool DOF::isDOFnull() const {
    if (m_edge->StartVertex->isnullligand() and m_edge->EndVertex->isnullligand()) {
        return true;
    }
    else{
        return false;
    }
}

bool DOF::isDOFbinding() const{
    int t1=0;
    int t2=0;
    if(m_edge->StartVertex->isligand() && m_edge->StartVertex->m_rigidbody!=nullptr){
        t1=1;
    }
    if(m_edge->EndVertex->isligand() && m_edge->EndVertex->m_rigidbody!=nullptr){
        t2=1;
    }

    if(t1+t2==1){
        return true;
    }
    else{
        return false;
    }
}

std::pair<int,int> DOF::dofatomid() const{
    std::pair<int,int> atom_id;
    atom_id.first=m_edge->getBond()->m_atom1->getId();
    atom_id.second=m_edge->getBond()->m_atom2->getId();
    return atom_id;
}