#include "Gauss.h"

unsigned int GaussianQuadrature::ngp = 2;

double GaussianQuadrature::weights_1p[1] = { 2.0 };
double GaussianQuadrature::points_1p[1] = { 0.0 };

double GaussianQuadrature::weights_2p[2] = { 1.0, 1.0 };
double GaussianQuadrature::points_2p[2] = { -0.57735026919, 0.57735026919 };

double GaussianQuadrature::weights_3p[3] = { 0.555555555556, 0.888888888889, 0.555555555556 };
double GaussianQuadrature::points_3p[3] = { -0.77459666924, 0.0, 0.77459666924 };

GaussianQuadrature::GaussianQuadrature(){
    g_p = new double[ngp];
    g_w = new double[ngp];
}

void GaussianQuadrature::Gauss_Calculate(){
    if(ngp == 1){
        g_w[0] =  weights_1p[0];
        g_p[0] =  points_1p[0];
    }
    else if(ngp == 2){
        g_w[0] =  weights_2p[0];
        g_w[1] =  weights_2p[1];
        g_p[0] =  points_2p[0];
        g_p[1] =  points_2p[1];
    }
    else if(ngp == 3){
        g_w[0] =  weights_3p[0];
        g_w[1] =  weights_3p[1];
        g_w[2] =  weights_3p[2];
        g_p[0] =  points_3p[0];
        g_p[1] =  points_3p[1];
        g_p[2] =  points_3p[2];
    }
}