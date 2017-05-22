#pragma once
#include "Math/MathUtil.h"

namespace Gregory
{
using GFloat = double;
//////////////////////////////////////////////////////////////////////////
// Coner Point p0
// p0 = (n-3)/(n+5)*v + 4/(n*(n+5))*sum(mi+ci)
//    = coef1 * v + coef2 * sum(mi+ci)
inline GFloat corner_coef1(int N) {
	return GFloat(N - 3) / (N + 5);
}
inline GFloat corner_coef2(int N) {
	return 4.0 / (N * (N + 5));
}
const GFloat pce2[] = { corner_coef1(2), corner_coef2(2) };
const GFloat pce3[] = { corner_coef1(3), corner_coef2(3) };
const GFloat pce4[] = { corner_coef1(4), corner_coef2(4) };
const GFloat pce5[] = { corner_coef1(5), corner_coef2(5) };
const GFloat pce6[] = { corner_coef1(6), corner_coef2(6) };
const GFloat pce7[] = { corner_coef1(7), corner_coef2(7) };
const GFloat pce8[] = { corner_coef1(8), corner_coef2(8) };
const GFloat pce9[] = { corner_coef1(9), corner_coef2(9) };
inline const GFloat* corner_coef(int valence)
{
	switch (valence)
	{
	case 2: return pce3;
	case 3: return pce3;
	case 4: return pce4;
	case 5: return pce5;
	case 6: return pce6;
	case 7: return pce7;
	case 8: return pce8;
	case 9: return pce9;
	default: break;
	}
}
//////////////////////////////////////////////////////////////////////////
inline GFloat theta_valence(int N) {
	return 1.0 / sqrt(4 + sqr(cos(M_PI / N)));
}
inline GFloat lamda_valence(int N) {
	GFloat invN = 1.0 / GFloat(N);
	return 0.0625 * (
		5.0 + cos(M_TWOPI * invN)
		+ cos(M_PI * invN) * sqrt(18.0 + 2.0 * cos(M_TWOPI * invN))
		);
}

// Precomputed theta for limit tangent
// at different vertex valence
const GFloat theta_v3 = theta_valence(3);
const GFloat theta_v4 = theta_valence(4);
const GFloat theta_v5 = theta_valence(5);
const GFloat theta_v6 = theta_valence(6);
const GFloat theta_v7 = theta_valence(7);
const GFloat theta_v8 = theta_valence(8);
const GFloat theta_v9 = theta_valence(9);
// query function for lambda
inline GFloat edge_theta(int valence)
{
	switch (valence)
	{
	//case 2: return theta_v2;
	case 3: return theta_v3;
	case 4: return theta_v4;
	case 5: return theta_v5;
	case 6: return theta_v6;
	case 7: return theta_v7;
	case 8: return theta_v8;
	case 9: return theta_v9;
	default: return theta_valence(valence);
	}
}
// Precomputed lambda for calculating edge points
const GFloat eigen_v2 = lamda_valence(2);//TODO: not sure how to handle v2 yet
const GFloat eigen_v3 = lamda_valence(3);
const GFloat eigen_v4 = lamda_valence(4);
const GFloat eigen_v5 = lamda_valence(5);
const GFloat eigen_v6 = lamda_valence(6);
const GFloat eigen_v7 = lamda_valence(7);
const GFloat eigen_v8 = lamda_valence(8);
const GFloat eigen_v9 = lamda_valence(9);
// query function for lambda
inline GFloat edge_eigen_val(int valence)
{
	switch (valence)
	{
	case 2: return eigen_v2;
	case 3: return eigen_v3;
	case 4: return eigen_v4;
	case 5: return eigen_v5;
	case 6: return eigen_v6;
	case 7: return eigen_v7;
	case 8: return eigen_v8;
	case 9: return eigen_v9;
	default: return lamda_valence(valence);
	}
}
// cosPi_N[i] = cos(Pi * i / N), 0<=i<2N 
const GFloat cosPi_3[] = { 1, 0.5, -0.5, -1, -0.5, 0.5 };
const GFloat cosPi_4[] = {
	 1,  sqrt(2) * 0.5, 0, -sqrt(2) * 0.5,
	-1, -sqrt(2) * 0.5, 0,  sqrt(2) * 0.5
};
const GFloat cosPi_5[] = {
	 1, cos(M_PI*0.2), cos(M_PI*0.4), cos(M_PI*0.6), cos(M_PI*0.8),
	-1, cos(M_PI*1.2), cos(M_PI*1.4), cos(M_PI*1.6), cos(M_PI*1.8)
};
const GFloat cosPi_6[] = {
	 1,  sqrt(3) * 0.5,  0.5, 0, -0.5, -sqrt(3) * 0.5,
	-1, -sqrt(3) * 0.5, -0.5, 0,  0.5,  sqrt(3) * 0.5
};
const GFloat cosPi_7[] = {
	 1, cos(M_PI / 7.0),     cos(M_PI * 2 / 7.0), cos(M_PI * 3 / 7.0),
	    cos(M_PI * 4 / 7.0), cos(M_PI * 5 / 7.0), cos(M_PI * 6 / 7.0),
	-1, cos(M_PI * 6 / 7.0), cos(M_PI * 5 / 7.0), cos(M_PI * 4 / 7.0),
	    cos(M_PI * 3 / 7.0), cos(M_PI * 2 / 7.0), cos(M_PI / 7.0)
};
const GFloat cosPi_8[] = {
	 1,  cos(M_PI * 0.125),  sqrt(2) * 0.5,  cos(M_PI * 0.375), // 0-3
	 0, -cos(M_PI * 0.375), -sqrt(2) * 0.5, -cos(M_PI * 0.125), //4-7
	-1, -cos(M_PI * 0.125), -sqrt(2) * 0.5, -cos(M_PI * 0.375),//8-11
	 0,  cos(M_PI * 0.375),  sqrt(2) * 0.5,  cos(M_PI * 0.125),//12-15
};
const GFloat cosPi_9[] = {
	   1, cos(M_PI / 9.0),     cos(M_PI * 2 / 9.0),
	 0.5, cos(M_PI * 4 / 9.0), cos(M_PI * 5 / 9.0),
	-0.5, cos(M_PI * 7 / 9.0), cos(M_PI * 8 / 9.0),
	  -1, cos(M_PI * 8 / 90),  cos(M_PI * 7 / 9.0),
	-0.5, cos(M_PI * 5 / 9.0), cos(M_PI * 4 / 9.0),
	 0.5, cos(M_PI * 2 / 9.0), cos(M_PI / 9.0)
};
// read in numerator, denominator to return value
inline GFloat cosPi(int num, int den)
{
	while (num >= 2 * den) num -= 2 * den;

	switch (den)
	{
	//case 2: return cosPi_2[num];
	case 3: return cosPi_3[num];
	case 4: return cosPi_4[num];
	case 5: return cosPi_5[num];
	case 6: return cosPi_6[num];
	case 7: return cosPi_7[num];
	case 8: return cosPi_8[num];
	case 9: return cosPi_9[num];
	default: return cos(M_PI * num / den);
	}
}
inline GFloat edgeP_coefMi(int i, int N)
{
	return (1.0 - edge_theta(N)*cosPi(1, N)) * cosPi(i * 2, N);
}

inline GFloat edgeP_coefCi(int i, int N)
{
	return 2 * edge_theta(N) * cosPi(i * 2 + 1, N);
}
}