/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/


#include "LoadCaseData.h"

#include <iomanip>
#include <iostream>

using namespace std;


CLoadCaseData :: ~CLoadCaseData()
{
	delete [] node_load;
	delete [] dof;
	delete [] load;
	delete [] node_nbc;
	delete [] dof_nbc;
	delete [] nbc;
}

// allocate storage
void CLoadCaseData :: Allocate(unsigned int num_load, unsigned int num_nbc)
{
	nloads = num_load; 
	node_load = new unsigned int[nloads]; 
	dof = new unsigned int[nloads];
	load = new double[nloads];
	nnbc = num_nbc;
	node_nbc = new unsigned int[2*nnbc];
	dof_nbc = new unsigned int[2*nnbc];
	nbc = new double[2*nnbc];
}; 

//	Read load case data from stream Input
bool CLoadCaseData :: Read(ifstream& Input)
{
//	Load case number (LL) and number of concentrated loads in this load case(NL)
	
	// number of load, number of natural bounce condition
	unsigned int NL, NNBC;

	Input >> NL >> NNBC;
	
	Allocate(NL, NNBC); //load number in this one
	for (unsigned int i = 0; i < NL; i++)
		Input >> node_load[i] >> dof[i] >> load[i];

	for (unsigned int i = 0; i < NNBC; i++){
		Input >> node_nbc[2*i] >> node_nbc[2*i+1] >> dof_nbc[i] >> nbc[2*i] >> nbc[2*i+1] >> EleGrp_num >> Ele_num; //node1,node2,dof,nbc load1,load2,element
		
		//require 2 dof are the same
		//if(dof_nbc[j]!=dof_nbc[j+1]){
			//cerr << "*** Error *** 2 dof should be the same for Natural Bounce Condition!" << endl
			//<< "    Wrong dof input : " << i+1 << endl;
			//return false;
		//}
	}

	return true;
}

//	Write load case data to stream
void CLoadCaseData::Write(COutputter& output) //output as a COutputter
{
	for (unsigned int i = 0; i < nloads; i++){
		output << setw(17) << node_load[i] << setw(13) << dof[i]  << setw(19) << load[i] << endl; 
	}
	for (unsigned int i = 0; i < nnbc; i++){
		output << setw(7) << node_nbc[2*i] << setw(3) << node_nbc[2*i+1] << setw(7) << dof_nbc[i] << setw(19) 
		<< nbc[2*i] << setw(13) << nbc[2*i+1] << setw(9) << Ele_num << endl;
	}
	
}
