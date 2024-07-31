# pragma once

# include "Element.h"

using namespace std;

//! Q4 element class
class CQ4 : public CElement
{
public:

	double* Nmat;

	double* detJ;

public:

//! Constructor
	CQ4();

//! Deconstructor
	~CQ4(){
		if(Nmat)
			delete[] Nmat;
		if(detJ)
			delete[] detJ;
	}

//! Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//! Write element data to stream
	virtual void Write(COutputter& output);

//! Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//! Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

	virtual void Calculate_Nmat(double csi, double eta);

	virtual void Calculate_detJ();

	virtual void Calculate_NBC(double* Force, double t1, double t2);

	virtual void Calculate_BODY(double* Force);
};
