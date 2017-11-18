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




#include <string>
#include <iostream>
#include <metrics/Dihedral.h>

#include "core/Molecule.h"
#include "core/Chain.h"
#include "core/Grid.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "Logger.h"
#include "metrics/RMSD.h"
#include "CTKTimer.h"

using namespace std;

extern double prevRMSDTime;
extern double currentRMSDTime;

Molecule * myReadFile(string pdbFile){
  char* tmp = realpath(pdbFile.c_str(), nullptr);
  if(tmp==nullptr){ cerr<<pdbFile<<" is not a valid PDB-file"<<endl; exit(-1); }

  Molecule* protein = IO::readPdb(pdbFile);
  return protein;
}


int main( int argc, char* argv[] ){
  enableLogger("rmsd");
  if(argc<3){ cerr<<"Too few arguments. Please specify PDB-file in arguments"<<endl; exit(-1);}
  
  Selection sel("all");
//  Configuration* reference = new Configuration(myReadFile(argv[1]));
  Molecule* reference = myReadFile(argv[1]);
  for(int i=2;i<argc;i++){
    Molecule * p = myReadFile(argv[i]);
//    Configuration* c = new Configuration(p);

    CTKTimer timer;
    timer.Reset();
    double start_time = timer.LastElapsedTime();

    double rmsd = p->alignReferencePositionsTo(reference, sel);
//        metric->align(p, reference->getMolecule());
    log("rmsd")<<argv[i]<<" : "<<rmsd<<endl;

    double end_time = timer.ElapsedTime();

    string outfile = p->getName()+"_aligned.pdb";
    IO::writePdb(p,outfile);

    cout<<"Alignment took "<<end_time-start_time<<endl;
//    cout<<"RMSD only took "<<end_time_2-end_time<<endl;
//    delete c;
    delete p;
  }

}


