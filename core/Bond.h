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
#ifndef MOLECULEBOND_H
#define MOLECULEBOND_H

#include <string>

class Atom;

class Bond {
  public:
	Atom* Atom1; // Atom1 Id is always smaller than Atom2 Id (is this correct?! @D)
	Atom* Atom2;
	std::string BondType;
	int Bars;
	bool constrained;

	Bond(Atom* atom1, Atom* atom2, std::string bond_type);
	Bond(Bond & bond);
	Bond();
	~Bond();
	void print();

	bool isLocked ();
	bool isPeptideBond ();
	bool isHbond ();
	
  double getTorsion();
};

std::ostream& operator<<(std::ostream& os, const Bond & b);
std::ostream& operator<<(std::ostream& os, const Bond * b);

#endif
