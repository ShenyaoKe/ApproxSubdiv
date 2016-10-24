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
	return 0.0625 * (
		5 + cos(M_TWOPI / N)
		+ cos(M_PI / N) * sqrt(18 + 2 * cos(M_TWOPI / N))
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
// Precomputed lambda for calculating edge points
const GFloat lambda_v3 = lamda_valence(3);
const GFloat lambda_v4 = lamda_valence(4);
const GFloat lambda_v5 = lamda_valence(5);
const GFloat lambda_v6 = lamda_valence(6);
const GFloat lambda_v7 = lamda_valence(7);
const GFloat lambda_v8 = lamda_valence(8);
const GFloat lambda_v9 = lamda_valence(9);

const GFloat cosPi_3[] = { 1, 0.5, -0.5, -1, -0.5, 0.5 };
const GFloat cosPi_4[] = {
	1, sqr(2) * 0.5, 0, -sqr(2) * 0.5,
	-1, -sqr(2) * 0.5, 0, sqr(2) * 0.5
};
const GFloat cosPi_5[] = {
	1, cos(M_PI*0.2), cos(M_PI*0.4), cos(M_PI*0.6), cos(M_PI*0.8),
	-1, cos(M_PI*1.2), cos(M_PI*1.4), cos(M_PI*1.6), cos(M_PI*1.8)
};
const GFloat cosPi_6[] = {
	1, sqrt(3) * 0.5, 0.5, 0, -0.5, -sqrt(3) * 0.5,
	-1, -sqrt(3) * 0.5, -0.5, 0, 0.5, sqrt(3) * 0.5
};
const GFloat cosPi_7[] = {
	1, cos(M_PI / 7.0), cos(M_PI * 2 / 7.0), cos(M_PI * 3 / 7.0),
	cos(M_PI * 4 / 7.0), cos(M_PI * 5 / 7.0), cos(M_PI * 6 / 7.0),
	-1, cos(M_PI * 6 / 7.0), cos(M_PI * 5 / 7.0), cos(M_PI * 4 / 7.0),
	cos(M_PI * 3 / 7.0), cos(M_PI * 2 / 7.0), cos(M_PI / 7.0),
};
const GFloat cosPi_8[] = {
	1, cos(M_PI * 0.125), sqr(2) * 0.5, cos(M_PI * 0.375), 0,//0-4
	-cos(M_PI * 0.125), -sqr(2) * 0.5, -cos(M_PI * 0.375), //5-7
	-1, -cos(M_PI * 0.375), -sqr(2) * 0.5, -cos(M_PI * 0.125),//8-11
	0, cos(M_PI * 0.375), sqr(2) * 0.5, cos(M_PI * 0.125),//12-15
};
const GFloat cosPi_9[] = {
	1, cos(M_PI / 9.0), cos(M_PI * 2 / 9.0),
	0.5, cos(M_PI * 4 / 9.0), cos(M_PI * 5 / 9.0),
	-0.5, cos(M_PI * 7 / 9.0), cos(M_PI * 8 / 9.0),
	-1, cos(M_PI * 8 / 90), cos(M_PI * 7 / 9.0),
	-0.5, cos(M_PI * 5 / 9.0), cos(M_PI * 4 / 9.0),
	0.5, cos(M_PI * 2 / 9.0), cos(M_PI / 9.0)
};

}