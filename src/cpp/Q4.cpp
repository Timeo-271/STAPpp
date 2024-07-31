#include "Q4.h"
#include "Gauss.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//  Constructor
CQ4::CQ4()
{
	NEN_ = 4;
	nodes_ = new CNode * [NEN_];

	ND_ = 12;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
	
	Nmat = new double[8];

	detJ = new double[9];
}

//  Read element data from stream Input
bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int Mset;
	unsigned int N1, N2, N3, N4;

	Input >> N1 >> N2 >> N3 >> N4 >> Mset;
	ElementMaterial_ = dynamic_cast<C2DMaterial*>(MaterialSets) + Mset - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//  Write element data to stream
void CQ4::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(9) << nodes_[2]->NodeNumber
		   << setw(9) << nodes_[3]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

//  Calculate element stiffness matrix
//  Upper triangular matrix stored as an array column by colum starting from the diagonal element
void CQ4::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	C2DMaterial* material_ = dynamic_cast<C2DMaterial*>(ElementMaterial_);

	double E = material_->E;
	double v = material_->mu;
	double t = material_->thickness;
	bool plane_stress = material_->plane_stress;

	if (plane_stress)	// plane strain condition
	{
		E = E / (1 - v * v);
		v = v / (1 - v);
	}
	(*this).Nmat = new double[8];

	double w[3] = {0};
	double gp[3] = {0};
	GaussianQuadrature Q4gauss;
	unsigned int ngp = Q4gauss.ngp;
	Q4gauss.Gauss_Calculate();
	for(unsigned int i = 0; i<ngp; i++){
		w[i] = Q4gauss.g_w[i];
		gp[i] = Q4gauss.g_p[i];
	}

	for (unsigned int i = 0; i < ngp; i++) {
		for (unsigned int j = 0; j < ngp; j++) {

			// Gaussian quadrature weights and points
			double W_i = w[i];
			double W_j = w[j];
			double P_i = gp[i];
			double P_j = gp[j];

			// Jacobian matrix and det
			double J_11 = 1.0 / 4 * ((P_j - 1) * (nodes_[0]->XYZ[0]) + (1 - P_j) * (nodes_[1]->XYZ[0]) + (1 + P_j) * (nodes_[2]->XYZ[0]) + (-P_j - 1) * (nodes_[3]->XYZ[0]));
			double J_12 = 1.0 / 4 * ((P_j - 1) * (nodes_[0]->XYZ[1]) + (1 - P_j) * (nodes_[1]->XYZ[1]) + (1 + P_j) * (nodes_[2]->XYZ[1]) + (-P_j - 1) * (nodes_[3]->XYZ[1]));
			double J_21 = 1.0 / 4 * ((P_i - 1) * (nodes_[0]->XYZ[0]) + (-P_i - 1) * (nodes_[1]->XYZ[0]) + (1 + P_i) * (nodes_[2]->XYZ[0]) + (1 - P_i) * (nodes_[3]->XYZ[0]));
			double J_22 = 1.0 / 4 * ((P_i - 1) * (nodes_[0]->XYZ[1]) + (-P_i - 1) * (nodes_[1]->XYZ[1]) + (1 + P_i) * (nodes_[2]->XYZ[1]) + (1 - P_i) * (nodes_[3]->XYZ[1]));
			double det_J = J_11 * J_22 - J_12 * J_21;

			// Inv of J
			double J_inv_11 = J_22 / det_J;
			double J_inv_12 = -J_12 / det_J;
			double J_inv_21 = -J_21 / det_J;
			double J_inv_22 = J_11 / det_J;

			// Grad of shape function
			double N1_x = 1.0 / 4 * (J_inv_11 * (P_j - 1) + J_inv_12 * (P_i - 1));
			double N2_x = 1.0 / 4 * (J_inv_11 * (1 - P_j) + J_inv_12 * (-P_i - 1));
			double N3_x = 1.0 / 4 * (J_inv_11 * (1 + P_j) + J_inv_12 * (1 + P_i));
			double N4_x = 1.0 / 4 * (J_inv_11 * (-P_j - 1) + J_inv_12 * (1 - P_i));
			double N1_y = 1.0 / 4 * (J_inv_21 * (P_j - 1) + J_inv_22 * (P_i - 1));
			double N2_y = 1.0 / 4 * (J_inv_21 * (1 - P_j) + J_inv_22 * (-P_i - 1));
			double N3_y = 1.0 / 4 * (J_inv_21 * (1 + P_j) + J_inv_22 * (1 + P_i));
			double N4_y = 1.0 / 4 * (J_inv_21 * (-P_j - 1) + J_inv_22 * (1 - P_i));

			Matrix[0] += W_i * W_j * (N1_x * N1_x * E / (1 - v * v) + N1_y * N1_y * E / 2 / (1 + v)) * det_J * t;

			Matrix[1] += W_i * W_j * (N1_y * N1_y * E / (1 - v * v) + N1_x * N1_x * E / 2 / (1 + v)) * det_J * t;
			Matrix[2] += W_i * W_j * (N1_x * N1_y * E / 2 / (1 - v)) * det_J * t;

			Matrix[6] += W_i * W_j * (N2_x * N2_x * E / (1 - v * v) + N2_y * N2_y * E / 2 / (1 + v)) * det_J * t;
			Matrix[8] += W_i * W_j * (N1_y * N2_x * E * v / (1 - v * v) + N1_x * N2_y * E / 2 / (1 + v)) * det_J * t;
			Matrix[9] += W_i * W_j * (N1_x * N2_x * E / (1 - v * v) + N1_y * N2_y * E / 2 / (1 + v)) * det_J * t;

			Matrix[10] += W_i * W_j * (N2_y * N2_y * E / (1 - v * v) + N2_x * N2_x * E / 2 / (1 + v)) * det_J * t;
			Matrix[11] += W_i * W_j * (N2_x * N2_y * E / 2 / (1 - v)) * det_J * t;
			Matrix[13] += W_i * W_j * (N1_y * N2_y * E / (1 - v * v) + N1_x * N2_x * E / 2 / (1 + v)) * det_J * t;
			Matrix[14] += W_i * W_j * (N1_x * N2_y * E * v / (1 - v * v) + N1_y * N2_x * E / 2 / (1 + v)) *det_J * t;

			Matrix[21] += W_i * W_j * (N3_x * N3_x * E / (1 - v * v) + N3_y * N3_y * E / 2 / (1 + v)) * det_J * t;
			Matrix[23] += W_i * W_j * (N2_y * N3_x * E * v / (1 - v * v) + N2_x * N3_y * E / 2 / (1 + v)) * det_J * t;
			Matrix[24] += W_i * W_j * (N2_x * N3_x * E / (1 - v * v) + N2_y * N3_y * E / 2 / (1 + v)) * det_J * t;
			Matrix[26] += W_i * W_j * (N1_y * N3_x * E * v / (1 - v * v) + N1_x * N3_y * E / 2 / (1 + v)) * det_J * t;
			Matrix[27] += W_i * W_j * (N1_x * N3_x * E / (1 - v * v) + N1_y * N3_y * E / 2 / (1 + v)) * det_J * t;
			
			Matrix[28] += W_i * W_j * (N3_y * N3_y * E / (1 - v * v) + N3_x * N3_x * E / 2 / (1 + v)) * det_J * t;
			Matrix[29] += W_i * W_j * (N3_x * N3_y * E / 2 / (1 - v)) * det_J * t;
			Matrix[31] += W_i * W_j * (N2_y * N3_y * E / (1 - v * v) + N2_x * N3_x * E / 2 / (1 + v)) * det_J * t;
			Matrix[32] += W_i * W_j * (N2_x * N3_y * E * v / (1 - v * v) + N2_y * N3_x * E / 2 / (1 + v)) * det_J * t;
			Matrix[34] += W_i * W_j * (N1_y * N3_y * E / (1 - v * v) + N1_x * N3_x * E / 2 / (1 + v)) * det_J * t;
			Matrix[35] += W_i * W_j * (N1_x * N3_y * E * v / (1 - v * v) + N1_y * N3_x * E / 2 / (1 + v)) * det_J * t;
			
			Matrix[45] += W_i * W_j * (N4_x * N4_x * E / (1 - v * v) + N4_y * N4_y * E / 2 / (1 + v)) * det_J * t;
			Matrix[47] += W_i * W_j * (N3_y * N4_x * E * v / (1 - v * v) + N3_x * N4_y * E / 2 / (1 + v)) * det_J * t;
			Matrix[48] += W_i * W_j * (N3_x * N4_x * E / (1 - v * v) + N3_y * N4_y * E / 2 / (1 + v)) * det_J * t;
			Matrix[50] += W_i * W_j * (N2_y * N4_x * E * v / (1 - v * v) + N2_x * N4_y * E / 2 / (1 + v)) * det_J * t;
			Matrix[51] += W_i * W_j * (N2_x * N4_x * E / (1 - v * v) + N2_y * N4_y * E / 2 / (1 + v)) * det_J * t;
			Matrix[53] += W_i * W_j * (N1_y * N4_x * E * v / (1 - v * v) + N1_x * N4_y * E / 2 / (1 + v)) * det_J * t;
			Matrix[54] += W_i * W_j * (N1_x * N4_x * E / (1 - v * v) + N1_y * N4_y * E / 2 / (1 + v)) * det_J * t;

			Matrix[55] += W_i * W_j * (N4_y * N4_y * E / (1 - v * v) + N4_x * N4_x * E / 2 / (1 + v)) * det_J * t;
			Matrix[56] += W_i * W_j * (N4_x * N4_y * E / 2 / (1 - v)) * det_J * t;
			Matrix[58] += W_i * W_j * (N3_y * N4_y * E / (1 - v * v) + N3_x * N4_x * E / 2 / (1 + v)) * det_J * t;
			Matrix[59] += W_i * W_j * (N3_x * N4_y * E * v / (1 - v * v) + N3_y * N4_x * E / 2 / (1 + v)) * det_J * t;
			Matrix[61] += W_i * W_j * (N2_y * N4_y * E / (1 - v * v) + N2_x * N4_x * E / 2 / (1 + v)) * det_J * t;
			Matrix[62] += W_i * W_j * (N2_x * N4_y * E * v / (1 - v * v) + N2_y * N4_x * E / 2 / (1 + v)) * det_J * t;
			Matrix[64] += W_i * W_j * (N1_y * N4_y * E / (1 - v * v) + N1_x * N4_x * E / 2 / (1 + v)) * det_J * t;
			Matrix[65] += W_i * W_j * (N1_x * N4_y * E * v / (1 - v * v) + N1_y * N4_x * E / 2 / (1 + v)) * det_J * t;
		}
	}
	
}

//  Calculate element stress
void CQ4::ElementStress(double* stress, double* Displacement)
{
	C2DMaterial* material_ = dynamic_cast<C2DMaterial*>(ElementMaterial_);

	double E = material_->E;
	double v = material_->mu;
	bool plane_stress = material_->plane_stress;

	if (plane_stress)	// plane strain condition
	{
		E = E / (1 - v * v);
		v = v / (1 - v);
	}

	for (unsigned int i = 0; i < 24; i++) {
		stress[i] = 0.0;
	}

	// Rebuild local displacement
	double displacement[12];

	for (unsigned int i = 0; i < 12; i++)
	{
		if (LocationMatrix_[i])
		{
			displacement[i] = Displacement[LocationMatrix_[i] - 1];
		}
		else
		{
			displacement[i] = 0;
		}
	}
	
	double w[3] = {0};
	double gp[3] = {0};
	GaussianQuadrature Q4gauss;
	unsigned int ngp = Q4gauss.ngp;
	Q4gauss.Gauss_Calculate();
	for(unsigned int i = 0; i<ngp; i++){
		w[i] = Q4gauss.g_w[i];
		gp[i] = Q4gauss.g_p[i];
	}

	for (unsigned i = 0; i < ngp; i++)
	{
		for (unsigned j = 0; j < ngp; j++)
		{
			// Coordinates of the stress superconvergence point
			double P_i = gp[i];
			double P_j = gp[j];

			// Jacobian matrix and det
			double J_11 = 1.0 / 4 * ((P_j - 1) * (nodes_[0]->XYZ[0]) + (1 - P_j) * (nodes_[1]->XYZ[0]) + (1 + P_j) * (nodes_[2]->XYZ[0]) + (-P_j - 1) * (nodes_[3]->XYZ[0]));
			double J_12 = 1.0 / 4 * ((P_j - 1) * (nodes_[0]->XYZ[1]) + (1 - P_j) * (nodes_[1]->XYZ[1]) + (1 + P_j) * (nodes_[2]->XYZ[1]) + (-P_j - 1) * (nodes_[3]->XYZ[1]));
			double J_21 = 1.0 / 4 * ((P_i - 1) * (nodes_[0]->XYZ[0]) + (-P_i - 1) * (nodes_[1]->XYZ[0]) + (1 + P_i) * (nodes_[2]->XYZ[0]) + (1 - P_i) * (nodes_[3]->XYZ[0]));
			double J_22 = 1.0 / 4 * ((P_i - 1) * (nodes_[0]->XYZ[1]) + (-P_i - 1) * (nodes_[1]->XYZ[1]) + (1 + P_i) * (nodes_[2]->XYZ[1]) + (1 - P_i) * (nodes_[3]->XYZ[1]));
			double det_J = J_11 * J_22 - J_12 * J_21;

			// Inv of J
			double J_inv_11 = J_22 / det_J;
			double J_inv_12 = -J_12 / det_J;
			double J_inv_21 = -J_21 / det_J;
			double J_inv_22 = J_11 / det_J;

			// Grad of shape function
			double N1_x = 1.0 / 4 * (J_inv_11 * (P_j - 1) + J_inv_12 * (P_i - 1));
			double N2_x = 1.0 / 4 * (J_inv_11 * (1 - P_j) + J_inv_12 * (-P_i - 1));
			double N3_x = 1.0 / 4 * (J_inv_11 * (1 + P_j) + J_inv_12 * (1 + P_i));
			double N4_x = 1.0 / 4 * (J_inv_11 * (-P_j - 1) + J_inv_12 * (1 - P_i));
			double N1_y = 1.0 / 4 * (J_inv_21 * (P_j - 1) + J_inv_22 * (P_i - 1));
			double N2_y = 1.0 / 4 * (J_inv_21 * (1 - P_j) + J_inv_22 * (-P_i - 1));
			double N3_y = 1.0 / 4 * (J_inv_21 * (1 + P_j) + J_inv_22 * (1 + P_i));
			double N4_y = 1.0 / 4 * (J_inv_21 * (-P_j - 1) + J_inv_22 * (1 - P_i));

			stress[(i * 2 + j) * 6 + 0] = E / (1 - v * v) *
				(N1_x * displacement[0] + N2_x * displacement[3] + N3_x * displacement[6] + N4_x * displacement[9])
				+ E * v / (1 - v * v) *
				(N1_y * displacement[1] + N2_y * displacement[4] + N3_y * displacement[7] + N4_y * displacement[10]);

			stress[(i * 2 + j) * 6 + 1] = E * v / (1 - v * v) *
				(N1_x * displacement[0] + N2_x * displacement[3] + N3_x * displacement[6] + N4_x * displacement[9])
				+ E / (1 - v * v) *
				(N1_y * displacement[1] + N2_y * displacement[4] + N3_y * displacement[7] + N4_y * displacement[10]);

			if (plane_stress)
			{
				stress[(i * 2 + j) * 6 + 2] = v / (1 + v) * (stress[0] + stress[1]);
			}
			else
			{
				stress[(i * 2 + j) * 6 + 2] = 0;
			}

			stress[(i * 2 + j) * 6 + 3] = E / 2 / (1 + v) *
				(N1_y * displacement[0] + N2_y * displacement[3] + N3_y * displacement[6] + N4_y * displacement[9]
					+ N1_x * displacement[1] + N2_x * displacement[4] + N3_x * displacement[7] + N4_x * displacement[10]);
		}
	}
}

// N define
void CQ4::Calculate_Nmat(double psi, double eta){
	
	Nmat[0] = 0.25*(1-psi)*(1-eta);
	Nmat[1] = 0.25*(1-psi)*(1-eta);
	Nmat[2] = 0.25*(1+psi)*(1-eta);
	Nmat[3] = 0.25*(1+psi)*(1-eta);
	Nmat[4] = 0.25*(1+psi)*(1+eta);
	Nmat[5] = 0.25*(1+psi)*(1+eta);
	Nmat[6] = 0.25*(1-psi)*(1+eta);
	Nmat[7] = 0.25*(1-psi)*(1+eta);
}

void CQ4::Calculate_detJ(){
	
	double w[3] = {0};
	double gp[3] = {0};
	GaussianQuadrature Q4gauss;
	unsigned int ngp = Q4gauss.ngp;
	Q4gauss.Gauss_Calculate();
	for(unsigned int i = 0; i<ngp; i++){
		w[i] = Q4gauss.g_w[i];
		gp[i] = Q4gauss.g_p[i];
	}

	for(unsigned int i = 0;i<ngp;i++){
		for(unsigned int j=0; j<ngp; j++){
			double P_i = gp[i];
			double P_j = gp[j];
			double J_11 = 1.0 / 4 * ((P_j - 1) * (nodes_[0]->XYZ[0]) + (1 - P_j) * (nodes_[1]->XYZ[0]) + (1 + P_j) * (nodes_[2]->XYZ[0]) + (-P_j - 1) * (nodes_[3]->XYZ[0]));
			double J_12 = 1.0 / 4 * ((P_j - 1) * (nodes_[0]->XYZ[1]) + (1 - P_j) * (nodes_[1]->XYZ[1]) + (1 + P_j) * (nodes_[2]->XYZ[1]) + (-P_j - 1) * (nodes_[3]->XYZ[1]));
			double J_21 = 1.0 / 4 * ((P_i - 1) * (nodes_[0]->XYZ[0]) + (-P_i - 1) * (nodes_[1]->XYZ[0]) + (1 + P_i) * (nodes_[2]->XYZ[0]) + (1 - P_i) * (nodes_[3]->XYZ[0]));
			double J_22 = 1.0 / 4 * ((P_i - 1) * (nodes_[0]->XYZ[1]) + (-P_i - 1) * (nodes_[1]->XYZ[1]) + (1 + P_i) * (nodes_[2]->XYZ[1]) + (1 - P_i) * (nodes_[3]->XYZ[1]));
			detJ[ngp*i+j] = J_11 * J_22 - J_12 * J_21;
		}
	}
}

void CQ4::Calculate_NBC(double* Force, double t1, double t2){
	clear(Force,2);
	double w[3] = {0};
	double gp[3] = {0};
	GaussianQuadrature Q4gauss;
	Q4gauss.Gauss_Calculate();
	unsigned int ngp = Q4gauss.ngp;
	//std::cout<<ngp<<endl;
	for(unsigned int i = 0; i<ngp; i++){
		w[i] = Q4gauss.g_w[i];
		gp[i] = Q4gauss.g_p[i];
	}
	for(unsigned int i = 0;i<ngp;i++){
		Calculate_Nmat(gp[i], -1);
		double t = 0.5*(1-gp[i])*t1 + 0.5*(1+gp[i])*t2;
		Force[0] += w[i]*Nmat[0]*t;
		Force[1] += w[i]*Nmat[2]*t;
	}
}

void CQ4::Calculate_BODY(double* Force){
	clear(Force,12);

	double b[8]={0};
	b[0] = nodes_[0]->BODY[0];
	b[1] = nodes_[0]->BODY[1];
	b[2] = nodes_[1]->BODY[0];
	b[3] = nodes_[1]->BODY[1];
	b[4] = nodes_[2]->BODY[0];
	b[5] = nodes_[2]->BODY[1];
	b[6] = nodes_[3]->BODY[0];
	b[7] = nodes_[3]->BODY[1];

	double w[3] = {0};
	double gp[3] = {0};
	GaussianQuadrature Q4gauss;
	Q4gauss.Gauss_Calculate();
	unsigned int ngp = Q4gauss.ngp;
	double be[2] = {0};
	double F[8] = {0};

	Calculate_detJ();

	for(unsigned int i = 0;i<ngp;i++){
		for(unsigned int j=0; j<ngp; j++){
			Calculate_Nmat(gp[i],gp[j]);
			for(unsigned int k=0;k<4;k=k+2){
				be[0] += Nmat[k]*b[k];
				be[1] += Nmat[k+1]*b[k+1];
			}
			for(unsigned int k=0;k<4;k=k+2){
			F[k] += w[i]*w[j]*detJ[ngp*i+j]*(be[0]*Nmat[k]);
			F[k+1] += w[i]*w[j]*detJ[ngp*i+j]*(be[1]*Nmat[k+1]);
			}
		}
	}

	unsigned int row[] = {0,1,3,4,6,7,9,10};
	for(unsigned int i = 0;i<8;i++){
		Force[row[i]] = F[i];
	}
}