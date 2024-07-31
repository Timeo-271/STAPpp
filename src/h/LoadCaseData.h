/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Outputter.h"

using namespace std;

//! Class LoadData is used to store load data
class CLoadCaseData
{
public:

	unsigned int nloads;	//!< Number of concentrated loads in this load case
	unsigned int* node_load;		//!< Node number to which this load is applied
	unsigned int* dof;		//!< Degree of freedom number for this load component
	double* load;			//!< Magnitude of load
	unsigned int nnbc;
	unsigned int* node_nbc;
	unsigned int* dof_nbc;
	double* nbc;
	unsigned int Ele_num;
	unsigned int EleGrp_num;

public:

	CLoadCaseData():nloads(0),nnbc(0),Ele_num(0),EleGrp_num(0),node_load(nullptr),dof(nullptr),load(nullptr),
	node_nbc(nullptr),dof_nbc(nullptr),nbc(nullptr){};

	~CLoadCaseData();

//!	Set nloads, and new array node, dof and load
	void Allocate(unsigned int num_load, unsigned int num_nbc);

//!	Read load case data from stream Input
	bool Read(ifstream& Input);

//!	Write load case data to stream
	void Write(COutputter& output);
};
