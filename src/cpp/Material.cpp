/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

bool C2DMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> mu >> thickness;	// Young's modulus, Poisson ratio and thickness of the plate

	Input >> plane_stress;	// Plane stress indicator, True if plane stress, False if plane strain

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}

void C2DMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << mu << setw(16) << thickness << setw(12) << plane_stress << endl;
}
