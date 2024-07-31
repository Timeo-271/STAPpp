#pragma once

class GaussianQuadrature {

public:

	static unsigned int ngp;

	static double weights_1p[1];

	static double points_1p[1];

	static double weights_2p[2];

	static double points_2p[2];

	static double weights_3p[3];

	static double points_3p[3];

	double* g_w; //weight for gauss

	double* g_p; // gauss point

public:

	GaussianQuadrature();

	~GaussianQuadrature(){
		if(g_w)
			delete[] g_w;
		if(g_p)
			delete[] g_p;
	}

	void Gauss_Calculate();

};