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

class CElementGroup;

using namespace std;

//!	Node class
class CNode
{
public:

//!	Maximum number of degrees of freedom per node
/*!	For 3D bar and solid elements, NDF = 3. For 3D beam or shell elements, NDF = 5 or 6 */
	unsigned int NDF_;

//!	Node numer
	unsigned int NodeNumber;

//!	x, y and z coordinates of the node
	double* XYZ;

	double* BODY;

//!	Boundary code of each degree of freedom of the node
/*!		0: The corresponding degree of freedom is active (defined in the global system) */
/*!		1: The corresponding degree of freedom in nonactive (not defined) */
/*!	After call Domain::CalculateEquationNumber(), bcode stores the global equation number */
/*!	corresponding to each degree of freedom of the node */
	unsigned int* bcode_;

//!	Constructor
	CNode(): NodeNumber(0),NDF_(3), XYZ(new double[3]), bcode_(new unsigned int[3]), BODY(new double[3]) 
	{for (int i = 0; i < 3; ++i) {
        XYZ[i] = 0.0;
        bcode_[i] = 0;
        BODY[i] = 0;}}

// !Desconstructor
	~CNode();

//!	Read nodal point data from stream Input
	bool Read(ifstream& Input);

//!	Output nodal point data to stream
	void Write(COutputter& output);

//!	Output equation numbers of nodal point to stream OutputFile
	void WriteEquationNo(COutputter& OutputFile);

//!	Write nodal displacement
	void WriteNodalDisplacement(COutputter& OutputFile, double* Displacement);
};
