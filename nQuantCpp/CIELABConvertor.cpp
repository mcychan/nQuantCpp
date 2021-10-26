#include "stdafx.h"
#include "CIELABConvertor.h"

#define _USE_MATH_DEFINES
#include <algorithm>
#include <math.h>

using namespace std;

#ifdef _DEBUG
#define new DEBUG_NEW
#endif	
	
void CIELABConvertor::RGB2LAB(const Color& c1, Lab& lab)
{
	auto r = c1.GetR() / 255.0, g = c1.GetG() / 255.0, b = c1.GetB() / 255.0;
	double x, y, z;

	r = (r > 0.04045) ? pow((r + 0.055) / 1.055, 2.4) : r / 12.92;
	g = (g > 0.04045) ? pow((g + 0.055) / 1.055, 2.4) : g / 12.92;
	b = (b > 0.04045) ? pow((b + 0.055) / 1.055, 2.4) : b / 12.92;

	x = (r * 0.4124 + g * 0.3576 + b * 0.1805) / 0.95047;
	y = (r * 0.2126 + g * 0.7152 + b * 0.0722) / 1.00000;
	z = (r * 0.0193 + g * 0.1192 + b * 0.9505) / 1.08883;

	x = (x > 0.008856) ? cbrt(x) : (7.787 * x) + 16.0 / 116.0;
	y = (y > 0.008856) ? cbrt(y) : (7.787 * y) + 16.0 / 116.0;
	z = (z > 0.008856) ? cbrt(z) : (7.787 * z) + 16.0 / 116.0;

	lab.alpha = c1.GetA();
	lab.L = (116 * y) - 16;
	lab.A = 500 * (x - y);
	lab.B = 200 * (y - z);
}
	
ARGB CIELABConvertor::LAB2RGB(const Lab& lab){
	auto y = (lab.L + 16) / 116.0;
	auto x = lab.A / 500 + y;
	auto z = y - lab.B / 200.0;
	double r, g, b;

	x = 0.95047 * ((x * x * x > 0.008856) ? x * x * x : (x - 16.0 / 116.0) / 7.787);
	y = 1.00000 * ((y * y * y > 0.008856) ? y * y * y : (y - 16.0 / 116.0) / 7.787);
	z = 1.08883 * ((z * z * z > 0.008856) ? z * z * z : (z - 16.0 / 116.0) / 7.787);

	r = x *  3.2406 + y * -1.5372 + z * -0.4986;
	g = x * -0.9689 + y *  1.8758 + z *  0.0415;
	b = x *  0.0557 + y * -0.2040 + z *  1.0570;

	r = (r > 0.0031308) ? (1.055 * pow(r, 1.0 / 2.4) - 0.055) : 12.92 * r;
	g = (g > 0.0031308) ? (1.055 * pow(g, 1.0 / 2.4) - 0.055) : 12.92 * g;
	b = (b > 0.0031308) ? (1.055 * pow(b, 1.0 / 2.4) - 0.055) : 12.92 * b;

	return Color::MakeARGB(clamp((int)lab.alpha, 0, BYTE_MAX), clamp((int)rint(r * BYTE_MAX), 0, BYTE_MAX),
		clamp((int)rint(g * BYTE_MAX), 0, BYTE_MAX), clamp((int)rint(b * BYTE_MAX), 0, BYTE_MAX));
}

/*******************************************************************************
* Conversions.
******************************************************************************/

inline const double deg2Rad(const double deg)
{
	return deg * M_PI / 180.0;
}

double CIELABConvertor::L_prime_div_k_L_S_L(const Lab& lab1, const Lab& lab2)
{
	const auto k_L = 1.0; // kL
	auto deltaLPrime = lab2.L - lab1.L;
	auto barLPrime = (lab1.L + lab2.L) / 2.0;
	auto S_L = 1 + ((0.015 * sqr(barLPrime - 50.0)) / _sqrt(20 + sqr(barLPrime - 50.0)));
	return deltaLPrime / (k_L * S_L);
}

double CIELABConvertor::C_prime_div_k_L_S_L(const Lab& lab1, const Lab& lab2, double& a1Prime, double& a2Prime, double& CPrime1, double& CPrime2)
{
	const auto k_C = 1.0, K1 = 0.045; // K1
	const auto pow25To7 = 6103515625.0; /* pow(25, 7) */
	auto C1 = _sqrt((lab1.A * lab1.A) + (lab1.B * lab1.B));
	auto C2 = _sqrt((lab2.A * lab2.A) + (lab2.B * lab2.B));
	auto barC = (C1 + C2) / 2.0;
	auto G = 0.5 * (1 - _sqrt(pow(barC, 7) / (pow(barC, 7) + pow25To7)));
	a1Prime = (1.0 + G) * lab1.A;
	a2Prime = (1.0 + G) * lab2.A;

	CPrime1 = _sqrt((a1Prime * a1Prime) + (lab1.B * lab1.B));
	CPrime2 = _sqrt((a2Prime * a2Prime) + (lab2.B * lab2.B));
	auto deltaCPrime = CPrime2 - CPrime1;
	auto barCPrime = (CPrime1 + CPrime2) / 2.0;
	
	auto S_C = 1 + (K1 * barCPrime);
	return deltaCPrime / (k_C * S_C);
}

double CIELABConvertor::H_prime_div_k_L_S_L(const Lab& lab1, const Lab& lab2, const double a1Prime, const double a2Prime, const double CPrime1, const double CPrime2, double& barCPrime, double& barhPrime)
{
	const auto k_H = 1.0, K2 = 0.015; // K2
	const auto deg360InRad = deg2Rad(360.0);
	const auto deg180InRad = deg2Rad(180.0);
	auto CPrimeProduct = CPrime1 * CPrime2;
	double hPrime1;
	if (lab1.B == 0 && a1Prime == 0)
		hPrime1 = 0.0;
	else {
		hPrime1 = atan2(lab1.B, a1Prime);
		/*
		* This must be converted to a hue angle in degrees between 0
		* and 360 by addition of 2π to negative hue angles.
		*/
		if (hPrime1 < 0)
			hPrime1 += deg360InRad;
	}
	double hPrime2;
	if (lab2.B == 0 && a2Prime == 0)
		hPrime2 = 0.0;
	else {
		hPrime2 = atan2(lab2.B, a2Prime);
		/*
		* This must be converted to a hue angle in degrees between 0
		* and 360 by addition of 2π to negative hue angles.
		*/
		if (hPrime2 < 0)
			hPrime2 += deg360InRad;
	}
	double deltahPrime;
	if (CPrimeProduct == 0)
		deltahPrime = 0;
	else {
		/* Avoid the fabs() call */
		deltahPrime = hPrime2 - hPrime1;
		if (deltahPrime < -deg180InRad)
			deltahPrime += deg360InRad;
		else if (deltahPrime > deg180InRad)
			deltahPrime -= deg360InRad;
	}

	auto deltaHPrime = 2.0 * _sqrt(CPrimeProduct) * sin(deltahPrime / 2.0);
	auto hPrimeSum = hPrime1 + hPrime2;
	if ((CPrime1 * CPrime2) == 0) {
		barhPrime = hPrimeSum;
	}
	else {
		if (fabs(hPrime1 - hPrime2) <= deg180InRad)
			barhPrime = hPrimeSum / 2.0;
		else {
			if (hPrimeSum < deg360InRad)
				barhPrime = (hPrimeSum + deg360InRad) / 2.0;
			else
				barhPrime = (hPrimeSum - deg360InRad) / 2.0;
		}
	}

	barCPrime = (CPrime1 + CPrime2) / 2.0;
	auto T = 1.0 - (0.17 * cos(barhPrime - deg2Rad(30.0))) +
		(0.24 * cos(2.0 * barhPrime)) +
		(0.32 * cos((3.0 * barhPrime) + deg2Rad(6.0))) -
		(0.20 * cos((4.0 * barhPrime) - deg2Rad(63.0)));
	auto S_H = 1 + (K2 * barCPrime * T);
	return deltaHPrime / (k_H * S_H);
}

double CIELABConvertor::R_T(const double barCPrime, const double barhPrime, const double C_prime_div_k_L_S_L, const double H_prime_div_k_L_S_L)
{
	const double pow25To7 = 6103515625.0; /* pow(25, 7) */
	auto deltaTheta = deg2Rad(30.0) * exp(-pow((barhPrime - deg2Rad(275.0)) / deg2Rad(25.0), 2.0));
	auto R_C = 2.0 * _sqrt(pow(barCPrime, 7.0) / (pow(barCPrime, 7.0) + pow25To7));
	auto R_T = -sin(2.0 * deltaTheta) * R_C;
	return R_T * C_prime_div_k_L_S_L * H_prime_div_k_L_S_L;
}

double CIELABConvertor::CIEDE2000(const Lab& lab1, const Lab& lab2)
{
	auto deltaL_prime_div_k_L_S_L = L_prime_div_k_L_S_L(lab1, lab2);
	double a1Prime, a2Prime, CPrime1, CPrime2;
	auto deltaC_prime_div_k_L_S_L = C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
	double barCPrime, barhPrime;
	auto deltaH_prime_div_k_L_S_L = H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
	auto deltaR_T = R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
	return
		pow(deltaL_prime_div_k_L_S_L, 2.0) +
		pow(deltaC_prime_div_k_L_S_L, 2.0) +
		pow(deltaH_prime_div_k_L_S_L, 2.0) +
		deltaR_T;
}
