/*************************************************
* 人工晶状体屈光度计算
* Calculation of IOL Power
* 
* @By RyotaKirkland
* 
* 目前已支持的公式：SRK-I、SRK-II、SRK-T、Binkhorst-I、HofferQ、Holladay-I、Haigis、Haigis-L、Shammas.
* Formulae Available: SRK-I、SRK-II、SRK-T、Binkhorst-I、HofferQ、Holladay-I、Haigis、Haigis-L、Shammas.
* 
* 支持正常眼IOL屈光度计算、给定IOL屈光度时的术后预留屈光度计算。
* Support the calculation of IOL power for normal eyes and the calculation of refractive error for a given IOL power.
* 
* 关于A、SF（HofferQ）、a[i]（Haigis，i=0，1，2）值：请参阅相关IOL规格书获取
* About A, SF (HofferQ), a[i] (Haigis, i=0, 1, 2) values: please refer to the relevant IOL specifications.
*************************************************/

#include <iostream>
#include "math.h"
#define pi 3.141592653589793238462643383279
//第一部分：正常眼IOL屈光度计算
//Part 1. The Calculation of IOL Power for Normal Eyes

double SRK_I(double AL, double K) {
	//SRK-I公式，输入变量为眼轴长度AL以及角膜曲率K
	//SRK-I Formula, variations input are Axial Length (AL), and Corneal Curvature (K).
	
	double A = 118.3;
	//A常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//A const provided here is a sample value.
	double P;
	//P为IOL屈光度
	//P stands for IOL power.
	P = A - 2.5 * AL - 0.9 * K;
	return P;
}

double SRK_II(double AL, double K) {
	//SRK-II公式，输入变量为眼轴长度AL以及角膜曲率K
	//SRK-II Formula, variations input are Axial Length (AL), and Corneal Curvature (K).
	
	double A = 118.3;
	//A常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//A const provided here is a sample value.
	double Ac;
	//Ac为修正后的A常数
	//Ac stands for corrected A const.
	if (AL < 20) {
		Ac = A + 3;
	}
	else if (AL >= 20 && AL < 21) {
		Ac = A + 2;
	}
	else if (AL >= 21 && AL < 22) {
		Ac = A + 1;
	}
	else if (AL >= 22 && AL < 24.5) {
		Ac = A;
	}
	else {
		Ac = A - 0.5;
	}
	double P;
	//P为IOL屈光度
	//P stands for IOL power.
	P = Ac - 2.5 * AL - 0.9 * K;
	return P;
}

double SRK_T(double AL, double K)
{
	//SRK-T公式，输入变量为眼轴长度AL以及角膜曲率K
	//SRK-T Formula, variations input are Axial Length (AL), and Corneal Curvature (K).
	
	double A = 118.3;
	//A常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//A const provided here is a sample value.
	double r = 337.5 / K;
	//计算角膜曲率半径
	//Calculation of Corneal Radius of Curvature.
	double L_COR;
	//计算眼轴长修正值
	//Calculation of Corrected Axial Length.
	if (AL > 24.2) {
		L_COR = -3.446 + 1.716 * AL - 0.02378 * AL * AL;
	}
	else {
		L_COR = AL;
	}

	double cornea_widith = -5.41 + 0.58412 * L_COR + 0.098 * K;
	//计算角膜宽度
	//Calculation of Corneal Widith.
	double cornea_height = r - sqrt(r * r - (cornea_widith * cornea_widith) / 4);
	//计算角膜高度
	//Calculation of Corneal Height.
	double ACD_const = 0.62467 * A - 68.74709;
	//计算前房深度
	//Calculation of Anterior Chamber Depth (ACD).
	double offset = ACD_const - 3.336;
	//计算偏移值
	//Calculation of ACD offset.
	double ACD_est = cornea_height + offset;
	//计算术后有效人工晶体位置估计值
	//Calculation of estimated postoperative effective IOL position.
	double retina_thick = 0.65696 - 0.02029 * AL;
	//计算视网膜厚度
	//Calculation of Retina Thickness.
	double L_optical = AL + retina_thick;
	//计算光轴长度
	//Calculation of Optical Axial Length.
	double V, na, nc, ncm1;
	//V为眼睛到镜片的距离，na为前房折射率，nc为角膜折射率，ncm1为nc-1的值
	//V stands for the distance of eye and spectacle lens.
	//na stands for the refractive index of Anterior Chamber.
	//nc stands for the refractive index of Cornea.
	//ncm1 equals to nc-1.
	V = 12;
	na = 1.336;
	nc = 1.333;
	ncm1 = 0.333;

	double P;
	P = 1000 * na * (na * r - ncm1 * L_optical) / ((L_optical - ACD_est) * (na * r - ncm1 * ACD_est));
	//P为IOL屈光度
	//P stands for IOL power.
	/*double IOL_Refr_is_Zero = (1000 * na * (na * r - ncm1 * L_optical)) / ((L_optical - ACD_est) * (na * r - ncm1 * ACD_est));
	return IOL_Refr_is_Zero;*/
	return P;
}

double Binkhorst_I(double AL, double K, double ACD) {
	//Binkhorst-I公式，输入变量为眼轴长度AL、角膜曲率K以及前房深度ACD
	//Binkhorst-I Formula, variations input are Axial Length (AL), Corneal Curvature (K), and Anterior Chamber Depth (ACD).

	double r = 337.5 / K;
	//计算角膜曲率半径
	//Calculation of Corneal Radius of Curvature.
	double P;
	P = 1336 * (4 * r - AL) / ((AL - ACD) * (4 * r - ACD));
	//P为IOL屈光度
	//P stands for IOL power.
	return P;
}


double HofferQ(double AL, double K) {
	//HofferQ公式，输入变量为眼轴长度AL以及角膜曲率K
	//HofferQ Formula, variations input are Axial Length (AL) and Corneal Curvature (K).
	
	double M, G, pACD, ACD;
	//判断ACD计算参量M和G的取值
	if (AL <= 23) {
		M = 1;
		G = 28;
	}
	else if (AL > 23) {
		M = -1;
		G = 23.5;
	}
	if (AL > 31) {
		AL = 31;
	}
	else if (AL < 18.5) {
		AL = 18.5;
	}
	double Krad = K * pi / 180;
	pACD = 5.15;
	//ACD常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//ACD const provided here is a sample value.
	ACD = pACD + 0.3 * (AL - 23.5) + (tan(Krad)) * (tan(Krad)) + 0.1 * M * (23.5 - AL) * (23.5 - AL) * tan((0.1 * (G - AL) * (G - AL)) * pi / 180) - 0.99166;
	//计算ACD预测值
	//Calculaton of predicted ACD.
	double P;
	P = (1336 / (AL - ACD - 0.05)) - (1.336 / ((1.336 / K ) - ((ACD + 0.05) / 1000)));
	//P为IOL屈光度
	//P stands for IOL power.
	return P;
}

double Holladay_I(double AL, double K) {
	//Holladay-I公式，输入变量为眼轴长度AL、角膜曲率K以及前房深度ACD
	//Holladay-I Formula, variations input are Axial Length (AL), Corneal Curvature (K), and Anterior Chamber Depth (ACD).
	
	double R;
	R = 337.5 / K;
	//计算角膜曲率半径
	//Calculation of Corneal Radius of Curvature.
	double na, nc;
	na = 1.336;
	nc = 1.333;
	//na为前房折射率，nc为角膜折射率
	//na stands for the refractive index of Anterior Chamber.
	//nc stands for the refractive index of Cornea.
	double AL_modified = AL + 0.2;
	//修正眼轴长测量值
	//Modified Axial Length
	double Rag, AG, ACD;
	//计算ACD预测值
	//Calculaton of predicted ACD.
	if (R >= 7) {
		Rag = R;
	}
	else {
		Rag = 7;
	}
	AG = 12.5 * AL / 23.45;
	if (AG > 13.5) {
		AG = 13.5;
	}
	ACD = 0.56 + Rag - sqrt(Rag * Rag - (AG * AG / 4));
	double SF = 1.37;
	//手术因子常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//Surgery Factor (SF) provided here is a sample value.
	double P;
	P = 1000 * na * (na * R - (nc - 1) * AL_modified) / ((AL_modified - ACD - SF) * (na * R - (nc - 1) * (ACD + SF)));
	//P为IOL屈光度
	//P stands for IOL power.
	return P;
}

double Haigis(double AL, double K, double ACD) {
	//Haigis公式，输入变量为眼轴长度AL、角膜曲率K以及前房深度ACD
	//Haigis Formula, variations input are Axial Length (AL), Corneal Curvature (K), and Anterior Chamber Depth (ACD).

	double a0, a1, a2;
	a0 = 0.325;
	a1 = 0.255;
	a2 = 0.141;
	//Haigis常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//Haigis constants provided here are sample values.
	double R;
	R = 337.5 / K;
	//计算角膜曲率半径
	//Calculation of Corneal Radius of Curvature.
	double na, nc;
	na = 1.336;
	nc = 1.3315;
	//na为前房折射率，nc为角膜折射率
	//na stands for the refractive index of Anterior Chamber.
	//nc stands for the refractive index of Cornea.
	double ACD_optical = a0 + a1 * ACD + a2 * AL;
	//光学ACD值
	//Optical ACD
	double Dc = (nc - 1) / R;
	//角膜屈光度
	//Refraction Power of cornea.
	double P;
	P = 1000 * (na / (AL - ACD_optical) - na / (na / Dc - ACD_optical));
	//P为IOL屈光度
	//P stands for IOL power.
	return P;
}

double Haigis_L(double AL, double K, double ACD) {
	//Haigis-L公式，输入变量为眼轴长度AL、角膜曲率K以及前房深度ACD
	//Haigis-L Formula, variations input are Axial Length (AL), Corneal Curvature (K), and Anterior Chamber Depth (ACD).

	double a0, a1, a2;
	a0 = 0.325;
	a1 = 0.255;
	a2 = 0.141;
	//Haigis常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//Haigis constants provided here are sample values.
	double R;
	R = 337.5 / K;
	//计算角膜曲率半径
	//Calculation of Corneal Radius of Curvature.
	double Rm;
	Rm = 331.5 / (-5.1625 * R + 82.2603 - 0.35);
	//计算修正角膜曲率半径
	//Calculation of modified Corneal Radius of Curvature.
	double na, nc;
	na = 1.336;
	nc = 1.3315;
	//na为前房折射率，nc为角膜折射率
	//na stands for the refractive index of Anterior Chamber.
	//nc stands for the refractive index of Cornea.
	double ACD_optical = a0 + a1 * ACD + a2 * AL;
	//光学ACD值
	//Optical ACD
	double Dc = (nc - 1) / Rm;
	//角膜屈光度
	//Refraction Power of cornea.
	double P;
	P = 1000 * (na / (AL - ACD_optical) - na / (na / Dc - ACD_optical));
	//P为IOL屈光度
	//P stands for IOL power.
	return P;
}

double Shammas(double AL, double K) {
	//Shammas公式，输入变量为眼轴长度AL以及角膜曲率K
	//Shammas Formula, variations input are Axial Length (AL), and Corneal Curvature (K).

	double A = 118.3;
	//A常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//A const provided here is a sample value.
	double Km = 1.14 * K - 6.8;
	//计算修正角膜曲率
	//Calculation of modified Corneal Curvature.
	double ACD_predict = 0.5835 * A - 64.4;
	//计算ACD预测值
	//Calculation of predicted ACD.
	double P;
	P = 1336 / (AL - 0.1 * (AL - 23) - (ACD_predict + 0.05)) - 1 / (1.0125 / Km - (ACD_predict + 0.05) / 1336);
	//P为IOL屈光度
	//P stands for IOL power.
	return P;
}

//第二部分：给定IOL屈光度时的术后预留屈光度、以及其他参数计算
//Part 2. The Calculation of refractive error for a given IOL power, and other parameters.

double SRK_I_Refr(double P0, double P) {
	//SRK-I公式，输入变量为正视IOL屈光度P0以及屈光不正IOL屈光度P
	//SRK-I Formula, variations input are Emmetropia Diopter (P0), and IOL diopter to be implanted (P).

	double Rf;
	if (P0 > 14) {
		Rf = 1.25;
	}
	else {
		Rf = 1;
	}
	//屈光度修正系数。
	//Diopter Correction Factor.
	double Refr;
	//Refr为术后预留屈光度
	//Refr stands for refractive error.
	Refr = (P0 - P) / Rf;
	return Refr;
}

double SRK_II_Refr(double P0, double P) {
	//SRK-II公式，输入变量为正视IOL屈光度P0以及屈光不正IOL屈光度P
	//SRK-II Formula, variations input are Emmetropia Diopter (P0), and IOL diopter to be implanted (P).

	double Rf;
	if (P0 > 14) {
		Rf = 1.25;
	}
	else {
		Rf = 1;
	}
	//屈光度修正系数。
	//Diopter Correction Factor.
	double Refr;
	//Refr为术后预留屈光度
	//Refr stands for refractive error.
	Refr = (P0 - P) / Rf;
	return Refr;
}

double SRK_T_Refr(double AL, double K, double P)
{
	//SRK-T公式，输入变量为眼轴长度AL、角膜曲率K以及屈光不正IOL度数P
	//SRK-T Formula, variations input are Axial Length (AL), Corneal Curvature (K), and IOL diopter to be implanted (P).

	double A = 118.3;
	//A常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//A const provided here is a sample value.
	double r = 337.5 / K;
	//计算角膜曲率半径
	//Calculation of Corneal Radius of Curvature.
	double L_COR;
	//计算眼轴长修正值
	//Calculation of Corrected Axial Length.
	if (AL > 24.2) {
		L_COR = -3.446 + 1.716 * AL - 0.02378 * AL * AL;
	}
	else {
		L_COR = AL;
	}

	double cornea_widith = -5.41 + 0.58412 * L_COR + 0.098 * K;
	//计算角膜宽度
	//Calculation of Corneal Widith.
	double cornea_height = r - sqrt(r * r - (cornea_widith * cornea_widith) / 4);
	//计算角膜高度
	//Calculation of Corneal Height.
	double ACD_const = 0.62467 * A - 68.74709;
	//计算前房深度
	//Calculation of Anterior Chamber Depth (ACD).
	double offset = ACD_const - 3.336;
	//计算偏移值
	//Calculation of ACD offset.
	double ACD_est = cornea_height + offset;
	//计算术后有效人工晶体位置估计值
	//Calculation of estimated postoperative effective IOL position.
	double retina_thick = 0.65696 - 0.02029 * AL;
	//计算视网膜厚度
	//Calculation of Retina Thickness.
	double L_optical = AL + retina_thick;
	//计算光轴长度
	//Calculation of Optical Axial Length.
	double V, na, nc, ncm1;
	//V为眼睛到镜片的距离，na为前房折射率，nc为角膜折射率，ncm1为nc-1的值
	//V stands for the distance of eye and spectacle lens.
	//na stands for the refractive index of Anterior Chamber.
	//nc stands for the refractive index of Cornea.
	//ncm1 equals to nc-1.
	V = 12;
	na = 1.336;
	nc = 1.333;
	ncm1 = 0.333;

	double Refr;
	//Refr为术后预留屈光度
	//Refr stands for refractive error.
	Refr = (1000 * na * (na * r - ncm1 * L_optical) - P * (L_optical - ACD_est) * (na * r - ncm1 * ACD_est)) / (na * (V * (na * r - ncm1 * L_optical) + L_optical * r) - 0.001 * P * (L_optical - ACD_est) * (V * (na * r - ncm1 * ACD_est) + ACD_est * r));
	return Refr;
}

double Binkhorst_I_Refr(double P0, double P) {
	//Binkhorst-I公式，输入变量为正视IOL屈光度P0以及屈光不正IOL屈光度P
	//Binkhorst-I Formula, variations input are Emmetropia Diopter (P0), and IOL diopter to be implanted (P).

	double Rf;
	if (P0 > 14) {
		Rf = 1.25;
	}
	else {
		Rf = 1;
	}
	//屈光度修正系数。
	//Diopter Correction Factor.
	double Refr;
	//Refr为术后预留屈光度
	//Refr stands for refractive error.
	Refr = (P0 - P) / Rf;
	return Refr;
}

double HofferQ_Refr(double AL, double K, double P) {
	//HofferQ公式，输入变量为眼轴长度AL、角膜曲率K以及前房深度ACD
	//HofferQ Formula, variations input are Axial Length (AL), Corneal Curvature (K), and Anterior Chamber Depth (ACD).
	double M, G, pACD, ACD;
	//判断ACD计算参量M和G的取值
	if (AL <= 23) {
		M = 1;
		G = 28;
	}
	else if (AL > 23) {
		M = -1;
		G = 23.5;
	}
	if (AL > 31) {
		AL = 31;
	}
	else if (AL < 18.5) {
		AL = 18.5;
	}
	double Krad = K * pi / 180;
	pACD = 5.15;
	//ACD常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//ACD const provided here is a sample value.
	ACD = pACD + 0.3 * (AL - 23.5) + (tan(Krad)) * (tan(Krad)) + 0.1 * M * (23.5 - AL) * (23.5 - AL) * tan((0.1 * (G - AL) * (G - AL)) * pi / 180) - 0.99166;
	//计算ACD预测值
	//Calculaton of predicted ACD.
	double R, Refr;
	R = 1.336 / (1.336 / (1336 / (AL - ACD - 0.05) - P) + (ACD + 0.05) / 1000) - K;
	Refr = R / (1 + 0.012 * R);
	//Refr为术后预留屈光度
	//Refr stands for refractive error.
	return Refr;
}

double Holladay_I_Refr(double AL, double K, double P) {
	//Holladay-I公式，输入变量为眼轴长度AL、角膜曲率K以及前房深度ACD
	//Holladay-I Formula, variations input are Axial Length (AL), Corneal Curvature (K), and Anterior Chamber Depth (ACD).

	double R;
	R = 337.5 / K;
	//计算角膜曲率半径
	//Calculation of Corneal Radius of Curvature.
	double na, nc;
	na = 1.336;
	nc = 1.333;
	//na为前房折射率，nc为角膜折射率
	//na stands for the refractive index of Anterior Chamber.
	//nc stands for the refractive index of Cornea.
	double AL_modified = AL + 0.2;
	//修正眼轴长测量值
	//Modified Axial Length
	double Rag, AG, ACD;
	//计算ACD预测值
	//Calculaton of predicted ACD.
	if (R >= 7) {
		Rag = R;
	}
	else {
		Rag = 7;
	}
	AG = 12.5 * AL / 23.45;
	if (AG > 13.5) {
		AG = 13.5;
	}
	ACD = 0.56 + Rag - sqrt(Rag * Rag - (AG * AG / 4));
	double SF = 1.37;
	//手术因子常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//Surgery Factor (SF) provided here is a sample value.
	double Refr, Refr1, Refr2;
	Refr1 = 1000 * na * (na * R - (nc - 1) * AL_modified) - P * (AL_modified - ACD - SF) * (na * R - (nc - 1) * (ACD + SF));
	Refr2 = na * (12 * (na * R - (nc - 1) * AL_modified) + AL_modified * R) - 0.001 * P * (AL_modified - ACD - SF) * (12 * (na * R - (nc - 1) * (ACD + SF)) + (ACD + SF) * R);
	Refr = Refr1 / Refr2;
	//Refr为术后预留屈光度
	//Refr stands for refractive error.
	return Refr;
}

double Haigis_Refr(double AL, double K, double ACD, double P) {
	//Haigis公式，输入变量为眼轴长度AL、角膜曲率K以及前房深度ACD
	//Haigis Formula, variations input are Axial Length (AL), Corneal Curvature (K), and Anterior Chamber Depth (ACD).

	double a0, a1, a2;
	a0 = 0.95;
	a1 = 0.4;
	a2 = 0.1;
	//Haigis常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//Haigis constants provided here are sample values.
	double na = 1.336;
	double nc = 1.3315;
	//na为前房折射率，nc为角膜折射率
	//na stands for the refractive index of Anterior Chamber.
	//nc stands for the refractive index of Cornea.
	double R = 337.5 / K;
	//计算角膜曲率半径
	//Calculation of Corneal Radius of Curvature.
	double ACD_optical = a0 + a1 * ACD + a2 * AL;
	//光学ACD值
	//Optical ACD
	double Dc = (nc - 1) / R;
	//角膜屈光度
	//Refraction Power of cornea.
	double V = 12;
	double Refr;
	Refr = 1000 / (1 / (na / (na / (na / (AL - ACD_optical) - P / 1000) + ACD_optical) - Dc) + V);
	//Refr为术后预留屈光度
	//Refr stands for refractive error.
	return Refr;
}

double Haigis_L_Refr(double AL, double K, double ACD, double P) {
	//Haigis-L公式，输入变量为眼轴长度AL、角膜曲率K以及前房深度ACD
	//Haigis-L Formula, variations input are Axial Length (AL), Corneal Curvature (K), and Anterior Chamber Depth (ACD).

	double a0, a1, a2;
	a0 = 0.325;
	a1 = 0.255;
	a2 = 0.141;
	//Haigis常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//Haigis constants provided here are sample values.
	double R;
	R = 337.5 / K;
	//计算角膜曲率半径
	//Calculation of Corneal Radius of Curvature.
	double Rm;
	Rm = 331.5 / (-5.1625 * R + 82.2603 - 0.35);
	//计算修正角膜曲率半径
	//Calculation of modified Corneal Radius of Curvature.
	double na, nc;
	na = 1.336;
	nc = 1.3315;
	//na为前房折射率，nc为角膜折射率
	//na stands for the refractive index of Anterior Chamber.
	//nc stands for the refractive index of Cornea.
	double ACD_optical = a0 + a1 * ACD + a2 * AL;
	//光学ACD值
	//Optical ACD
	double Dc = (nc - 1) / Rm;
	//角膜屈光度
	//Refraction Power of cornea.
	double V = 12;
	double Refr;
	Refr = 1000 / (1 / (na / (na / (na / (AL - ACD_optical) - P / 1000) + ACD_optical) - Dc) + V);
	//Refr为术后预留屈光度
	//Refr stands for refractive error.
	return Refr;
}

double Shammas_Refr(double AL, double K, double P) {
	//Shammas公式，输入变量为眼轴长度AL以及角膜曲率K
	//Shammas Formula, variations input are Axial Length (AL), and Corneal Curvature (K).

	double A = 118.3;
	//A常数，本文件提供示例值，实际请参阅相关IOL规格书等资料获取。
	//A const provided here is a sample value.
	double Km = 1.14 * K - 6.8;
	//计算修正角膜曲率
	//Calculation of modified Corneal Curvature.
	double ACD_predict = 0.5835 * A - 64.4;
	//计算ACD预测值
	//Calculation of predicted ACD.
	double Refr;
	Refr = 1.0125 / (1 / (1336 / (AL - 0.1 * (AL - 23) - (ACD_predict + 0.05)) - P) + (ACD_predict + 0.05) / 1336) - Km;
	//P为IOL屈光度
	//P stands for IOL power.
	return Refr;
}
