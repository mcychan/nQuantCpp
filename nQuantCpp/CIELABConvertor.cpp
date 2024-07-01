#include "stdafx.h"
#include "CIELABConvertor.h"

#define _USE_MATH_DEFINES
#include <algorithm>
#include <math.h>
#include <iostream>

using namespace std;

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

static const double XYZ_WHITE_REFERENCE_X = 95.047;
static const double XYZ_WHITE_REFERENCE_Y = 100;
static const double XYZ_WHITE_REFERENCE_Z = 108.883;
static const double XYZ_EPSILON = 0.008856;
static const double XYZ_KAPPA = 903.3;

static float pivotXyzComponent(double component) {
	return component > XYZ_EPSILON
		? (float) cbrt(component)
		: (float)((XYZ_KAPPA * component + 16) / 116.0);
}

static double gammaToLinear(int channel)
{
	auto c = channel / 255.0;
	return c < 0.04045 ? c / 12.92 : pow((c + 0.055) / 1.055, 2.4);
}
	
void CIELABConvertor::RGB2LAB(const Color& c1, Lab& lab)
{
	auto sr = gammaToLinear(c1.GetR());
	auto sg = gammaToLinear(c1.GetG());
	auto sb = gammaToLinear(c1.GetB());
	const auto x = pivotXyzComponent(100 * (sr * 0.4124 + sg * 0.3576 + sb * 0.1805) / XYZ_WHITE_REFERENCE_X);
	const auto y = pivotXyzComponent(100 * (sr * 0.2126 + sg * 0.7152 + sb * 0.0722) / XYZ_WHITE_REFERENCE_Y);
	const auto z = pivotXyzComponent(100 * (sr * 0.0193 + sg * 0.1192 + sb * 0.9505) / XYZ_WHITE_REFERENCE_Z);

	lab.alpha = c1.GetA();
	lab.L = max(0, 116 * y - 16);
	lab.A = 500 * (x - y);
	lab.B = 200 * (y - z);
}
	
ARGB CIELABConvertor::LAB2RGB(const Lab& lab){
	const auto fy = (lab.L + 16.0) / 116.0;
	const auto fx = lab.A / 500 + fy;
	const auto fz = fy - lab.B / 200.0;
	double tmp = fx * fx * fx;
	const auto xr = tmp > XYZ_EPSILON ? tmp : (116.0 * fx - 16) / XYZ_KAPPA;
	const auto yr = lab.L > XYZ_KAPPA * XYZ_EPSILON ? fy * fy * fy : lab.L / XYZ_KAPPA;
	tmp = fz * fz * fz;
	const auto zr = tmp > XYZ_EPSILON ? tmp : (116.0 * fz - 16) / XYZ_KAPPA;
	const auto x = xr * XYZ_WHITE_REFERENCE_X;
	const auto y = yr * XYZ_WHITE_REFERENCE_Y;
	const auto z = zr * XYZ_WHITE_REFERENCE_Z;

	double r = (x * 3.2406 + y * -1.5372 + z * -0.4986) / 100.0;
	double g = (x * -0.9689 + y * 1.8758 + z * 0.0415) / 100.0;
	double b = (x * 0.0557 + y * -0.2040 + z * 1.0570) / 100.0;
	r = r > 0.0031308 ? 1.055 * pow(r, 1 / 2.4) - 0.055 : 12.92 * r;
	g = g > 0.0031308 ? 1.055 * pow(g, 1 / 2.4) - 0.055 : 12.92 * g;
	b = b > 0.0031308 ? 1.055 * pow(b, 1 / 2.4) - 0.055 : 12.92 * b;

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

double CIELABConvertor::Y_Diff(const Color& c1, const Color& c2)
{
	auto color2Y = [](const Color& c) -> double {
		auto sr = gammaToLinear(c.GetR());
		auto sg = gammaToLinear(c.GetG());
		auto sb = gammaToLinear(c.GetB());
		return sr * 0.2126 + sg * 0.7152 + sb * 0.0722;
	};
		
	auto y = color2Y(c1);
	auto y2 = color2Y(c2);
	return abs(y2 - y) * XYZ_WHITE_REFERENCE_Y;
}

double CIELABConvertor::U_Diff(const Color& c1, const Color& c2)
{
	auto color2U = [](const Color& c) -> double {
		return -0.09991 * c.GetR() - 0.33609 * c.GetG() + 0.436 * c.GetB();
	};

	auto u = color2U(c1);
	auto u2 = color2U(c2);
	return abs(u2 - u);
}
