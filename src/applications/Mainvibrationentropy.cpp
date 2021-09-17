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
#include <list>
#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_vector.h>
#include <math/gsl_helpers.h>
#include <math/NullspaceSVD.h>
#include <math/Eigenvalue.h>

#include "core/Molecule.h"
#include "core/Grid.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "Logger.h"
#include "applications/options/VibrationentropyOptions.h"
#include "../core/Configuration.h"
#include "CTKTimer.h"

extern double jacobianAndNullspaceTime;
extern double rigidityTime;

using namespace std;

int main( int argc, char* argv[] ){

  CTKTimer timer;
  timer.Reset();
  double start_time = timer.LastElapsedTime();

  if(argc<2){ cerr<<"Too few arguments. Please specify PDB-file in arguments"<<endl; exit(-1);}

  VibrationentropyOptions::createOptions(argc,argv);
  VibrationentropyOptions& options = *(VibrationentropyOptions::getOptions());

  enableLogger("rigidity");
  enableLogger("default");

  enableLogger("so"); //print options
  options.print();

  string out_path = options.workingDirectory;
  Selection movingResidues(options.residueNetwork);
  Molecule* protein = IO::readPdb(
      options.initialStructureFile,
      options.extraCovBonds,
      options.hydrogenbondMethod,
      options.hydrogenbondFile
  );

    Molecule* equilibrium = IO::readPdb(
            options.equilibriumStructureFile,
            options.extraCovBonds,
            options.hydrogenbondMethod,
            options.hydrogenbondFile
    );
    
    if (options.setligand!=""){
        protein->setatomligand(options.setligand);
    }
    
  protein->initializeTree(movingResidues,1.0,options.roots);
  string name = protein->getName();

  log("rigidity")<<"Molecule with ligand has:"<<endl;
  log("rigidity") << "> " << protein->getAtoms().size() << " atoms" << endl;
  log("rigidity")<<"> "<<protein->getInitialCollisions().size()<<" initial collisions"<<endl;
  log("rigidity")<<"> "<<protein->m_spanningTree->m_cycleAnchorEdges.size()<<" total bond constraints"<<endl;
  log("rigidity")<<"> "<<protein->getHBonds().size()<<" hydrogen bonds"<<endl;
  log("rigidity")<<"> "<<protein->getHydrophobicBonds().size()<<" hydrophobic bonds"<<endl;
  log("rigidity")<<"> "<<protein->getDBonds().size()<<" distance bonds"<<endl;
  log("rigidity") << "> " << protein->m_spanningTree->getNumDOFs() << " DOFs of which " << protein->m_spanningTree->getNumCycleDOFs() << " are cycle-DOFs" << endl;
  log("rigidity") << "> " << protein->m_spanningTree->ligand_dof_id.size() << " DOFs of which " << protein->m_spanningTree->ligand_cycledof_id.size() << " are cycle-DOFs from ligand\n" << endl;

  Configuration* conf = protein->m_conf;
  NullspaceSVD ns = *(dynamic_cast<NullspaceSVD*>(conf->getNullspace()));
  int numRows = ns.getMatrix()->size1;
  int numCols = ns.getMatrix()->size2;
  int nullspaceCols = ns.getNullspaceSize();
  int rankJacobian = numCols - nullspaceCols;
  int numRedundantCons = numRows-rankJacobian;

    //Site transfer DOF analysis
    ///Only if source and sink are provided, compute mutual information
    bool mutualInformation = options.source !="";
    if(mutualInformation) {
        Selection source(options.source);
        Selection sink(options.sink);
        /*double mutInfo_noligand = protein_noligand->m_conf->siteDOFTransfer(source, sink,
                                                          ns_noligand.getBasis()); /// change this to V-matrix for whole sliding mechanism*/
        double mutInfo = protein->m_conf->siteDOFTransfer(source, sink,
                                                          ns.getBasis()); /// change this to V-matrix for whole sliding mechanism
    }

  /// Create larger rigid substructures for rigid cluster decomposition
    //Molecule* rigidified_noligand = protein_noligand->collapseRigidBonds(options.collapseRigid);
    Molecule* rigidified = protein->collapseRigidBonds(options.collapseRigid);

  ///Write PDB File for pyMol usage
  int sample_id = 1;
  string out_file = out_path + "output/" + name + "_new_" +
                    std::to_string((long long)sample_id)
                    //static_cast<ostringstream*>( &(ostringstream() << sample_id) )->str()
                    + ".pdb";

  rigidified->writeRigidbodyIDToBFactor();
  rigidified->m_conf->m_vdwEnergy = protein->vdwEnergy();
  IO::writePdb(rigidified, out_file);

  //if(conf->checknocoupling(conf->getNullspace(),conf->getNullspacenocoupling()) && options.nocoupling == "true"){std::cout<<"this protein has no binding bond between ligand and protein"<<std::endl;exit(-1);}

  /*rigidified_noligand->writeRigidbodyIDToBFactor();
  rigidified_noligand->m_conf->m_vdwEnergy = protein_noligand->vdwEnergy();
  IO::writePdb(rigidified_noligand, out_file_noligand);*/

    string proteinonlyname="noprotonly";
    if(options.proteinonly){proteinonlyname="protonly";}

for(int i=0;i<options.entropycutoff.size();i++) {
    for(int j=0;j<options.vdwenergycutoff.size();j++) {
        conf->Hessianmatrixentropy(options.entropycutoff[i],options.coefficient,options.vdwenergycutoff[j],conf->getNullspace(), equilibrium,options.proteinonly,"false");
        if(options.getHessian) {
            ///save Jacobian and Nullspace to file
            string outJac1 =
                    out_path + "output/" + name + "_Hessian_" + std::to_string((int) options.vdwenergycutoff[j]) + "_" +
                    std::to_string((int) options.entropycutoff[i]) + "_" + proteinonlyname + "_test_" +
                    std::to_string((long long) sample_id)
                    + ".txt";
            gsl_matrix_outtofile(conf->geteigenvalue()->getHessiantorsionangle(), outJac1);
        }
        string outNull1 =
                out_path + "output/" + name + "_eigen_" + std::to_string((int) options.vdwenergycutoff[j]) + "_" +
                std::to_string((int) options.entropycutoff[i]) + "_" + proteinonlyname + "_test_" +
                std::to_string((long long) sample_id)
                + ".txt";
        gsl_vector_outtofile(conf->geteigenvalue()->getSingularvalue(), outNull1);

        if (!conf->checknocoupling() && options.nocoupling == "true") {
            conf->Hessianmatrixentropy(options.entropycutoff[i], options.coefficient,options.vdwenergycutoff[j],conf->getNullspacenocoupling(),equilibrium, options.proteinonly,"true");
            ///save Jacobian and Nullspace to file
            if(options.getHessian){
                string outJac2 =
                        out_path + "output/" + name + "_nocoupling_Hessian_" +
                        std::to_string((int) options.vdwenergycutoff[j]) +
                        "_" + std::to_string((int) options.entropycutoff[i]) + "_" + proteinonlyname + "_test_" +
                        std::to_string((long long) sample_id)
                        + ".txt";
                gsl_matrix_outtofile(conf->geteigenvalue()->getHessiantorsionangle(), outJac2);
            }
            string outNull2 =
                    out_path + "output/" + name + "_nocoupling_eigen_" + std::to_string((int) options.vdwenergycutoff[j]) +
                    "_" + std::to_string((int) options.entropycutoff[i]) + "_" + proteinonlyname + "_test_" +
                    std::to_string((long long) sample_id)
                    + ".txt";
            /*FILE * f2 = fopen (outNull2.c_str(), "w");
            gsl_vector_fwrite (f2, Evnocoupling->getSingularvalue());
            fclose (f2);
            FILE * f3 = fopen (outJac2.c_str(), "w");
            gsl_matrix_fwrite (f3, Evnocoupling->getHessiantorsionangle());
            fclose (f3);*/

            gsl_vector_outtofile(conf->geteigenvalue()->getSingularvalue(), outNull2);
            //gsl_vector_outtofile(conf->getNullspacenocoupling()->buildDofRigid(), rigiddof2);
        }
    }
}

    //Print final status
    double end_time = timer.ElapsedTime();
    log("rigidity")<< "Took "<<(end_time-start_time)<<" seconds to perform rigidity analysis\n";

    if(options.saveData <= 0) return 0;

    ///save pyMol coloring script
    string pyMol = out_path + "output/" + name + "_pyMol_" +
                   std::to_string((long long) sample_id)
                   + ".pml";
    string statFile = out_path + "output/" + name + "_stats_" +
                      std::to_string((long long) sample_id)
                      + ".txt";
    ///Write statistics
    IO::writeStats(protein, statFile, rigidified); //original protein with all bonds etc, rigidified one for cluster info

    ///Write pyMol script
    IO::writePyMolScript(rigidified, out_file, pyMol, protein);



    return 0;
}


