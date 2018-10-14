#pragma once

class CIELABConvertor
{

public:
	struct Lab {
		byte alpha = BYTE_MAX;
		double A = 0.0;
		double B = 0.0;
		double L = 0.0;
	};
	
	static ARGB LAB2RGB(const Lab& lab);
	static void RGB2LAB(const Color& c1, Lab& lab);
	static double L_prime_div_k_L_S_L(const Lab& lab1, const Lab& lab2);
	static double C_prime_div_k_L_S_L(const Lab& lab1, const Lab& lab2, double& a1Prime, double& a2Prime, double& CPrime1, double& CPrime2);
	static double H_prime_div_k_L_S_L(const Lab& lab1, const Lab& lab2, const double a1Prime, const double a2Prime, const double CPrime1, const double CPrime2, double& barCPrime, double& barhPrime);
	static double R_T(const double barCPrime, const double barhPrime, const double C_prime_div_k_L_S_L, const double H_prime_div_k_L_S_L);
	
	/* From the paper "The CIEDE2000 Color-Difference Formula: Implementation Notes, */
	/* Supplementary Test Data, and Mathematical Observations", by */
	/* Gaurav Sharma, Wencheng Wu and Edul N. Dalal, */
	/* Color Res. Appl., vol. 30, no. 1, pp. 21-30, Feb. 2005. */
	/* Return the CIEDE2000 Delta E color difference measure squared, for two Lab values */
	static double CIEDE2000(const Lab& lab1, const Lab& lab2);
};