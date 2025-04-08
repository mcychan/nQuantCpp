#pragma once
#ifndef _WIN32
	#include <guiddef.h>
	#include <gdiplus.h>
	using namespace Gdiplus;

	#ifndef BYTE
		typedef unsigned char BYTE;
	#endif
	#ifndef BYTE_MAX
		#include <climits>
		#define BYTE_MAX UCHAR_MAX
	#endif
#endif
#define _USE_MATH_DEFINES
#include <math.h>

class CIELABConvertor
{

public:
	struct Lab {
		BYTE alpha = BYTE_MAX;
		float A = 0.0;
		float B = 0.0;
		float L = 0.0;
	};
	
	static ARGB LAB2RGB(const Lab& lab);
	static void RGB2LAB(const Color& c1, Lab& lab);
	static float L_prime_div_k_L_S_L(const Lab& lab1, const Lab& lab2);
	static float C_prime_div_k_L_S_L(const Lab& lab1, const Lab& lab2, double& a1Prime, double& a2Prime, double& CPrime1, double& CPrime2);
	static float H_prime_div_k_L_S_L(const Lab& lab1, const Lab& lab2, const double a1Prime, const double a2Prime, const double CPrime1, const double CPrime2, double& barCPrime, double& barhPrime);
	static float R_T(const double barCPrime, const double barhPrime, const double C_prime_div_k_L_S_L, const double H_prime_div_k_L_S_L);
	
	/* From the paper "The CIEDE2000 Color-Difference Formula: Implementation Notes, */
	/* Supplementary Test Data, and Mathematical Observations", by */
	/* Gaurav Sharma, Wencheng Wu and Edul N. Dalal, */
	/* Color Res. Appl., vol. 30, no. 1, pp. 21-30, Feb. 2005. */
	/* Return the CIEDE2000 Delta E color difference measure squared, for two Lab values */
	static float CIEDE2000(const Lab& lab1, const Lab& lab2);
	
	static double Y_Diff(const Color& c1, const Color& c2);
	static double U_Diff(const Color& c1, const Color& c2);
};