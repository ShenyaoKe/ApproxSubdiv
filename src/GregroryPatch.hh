#pragma once
#include "Math/MathUtil.h"

namespace Gregory
{

	inline double theta_valence(int N) { return sqrt(4 + sqr(cos(M_PI / N))); }
	inline double lamda_valence(int N) {
		return 0.0625 * (
			5 + cos(M_TWOPI / N)
			+ cos(M_PI / N) * sqrt(18 + 2 * cos(M_TWOPI / N))
			);
	}
const double theta_valence3 = theta_valence(3);
const double theta_valence5 = theta_valence(5);
const double theta_valence6 = theta_valence(6);
const double theta_valence7 = theta_valence(7);
const double theta_valence8 = theta_valence(8);
const double theta_valence9 = theta_valence(9);
const double lambda_valence3 = lamda_valence(3);
const double lambda_valence5 = lamda_valence(5);
const double lambda_valence6 = lamda_valence(6);
const double lambda_valence7 = lamda_valence(7);
const double lambda_valence8 = lamda_valence(8);
const double lambda_valence9 = lamda_valence(9);

}