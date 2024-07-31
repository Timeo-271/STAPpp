/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/


#include <iostream>
#include <iomanip>

#include "Node.h"
#include "ElementGroup.h"
#include "LoadCaseData.h"


//Desconstructor
CNode::~CNode() {
        if (XYZ)
            delete [] XYZ;
        
        if (bcode_)
            delete [] bcode_;

		if (BODY)
            delete [] BODY;
    }

//	Read element data from stream Input
bool CNode::Read(ifstream& Input)
{
	Input >> NodeNumber;	// node number
	Input >> bcode_[0] >> bcode_[1] >> bcode_[2]
		  >> XYZ[0] >> XYZ[1] >> XYZ[2] 
		  >> BODY[0] >> BODY[1] >>BODY[2];
	return true;
}

//	Output nodal point data to stream
void CNode::Write(COutputter& output)
{
	output << setw(9) << NodeNumber << setw(5) << bcode_[0] << setw(5) << bcode_[1] << setw(5) << bcode_[2]
		   << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
	
}

//	Output equation numbers of nodal point to stream
void CNode::WriteEquationNo(COutputter& output)
{
	output << setw(9) << NodeNumber << "       ";

	for (unsigned int dof = 0; dof < NDF_; dof++)	// Loop over for DOFs of node np
	{
		output << setw(5) << bcode_[dof];
	}

	output << endl;
}

//	Write nodal displacement
void CNode::WriteNodalDisplacement(COutputter& output, double* Displacement)
{
	output << setw(5) << NodeNumber << "        ";

	for (unsigned int j = 0; j < NDF_; j++)
	{
		if (bcode_[j] == 0)
		{
			output << setw(18) << 0.0;
		}
		else
		{
			output << setw(18) << Displacement[bcode_[j] - 1];
		}
	}

	output << endl;
}
