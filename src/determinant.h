
#ifndef _DETERMINANT_H_
#define _DETERMINANT_H_

inline double det2(
	double a11, double a12, 
	double a21, double a22)
{
	return a11*a22-a12*a21;
}

inline double det3(
	double a11, double a12, double a13, 
	double a21, double a22, double a23, 
	double a31, double a32, double a33)
{
	return a11*det2(a22, a23, a32, a33)
	      -a12*det2(a21, a23, a31, a33)
	      +a13*det2(a21, a22, a31, a32);
}

inline double det4(
	double a11, double a12, double a13, double a14, 
	double a21, double a22, double a23, double a24,
	double a31, double a32, double a33, double a34,
	double a41, double a42, double a43, double a44)
{
	return a11*det3(a22, a23, a24, a32, a33, a34, a42, a43, a44)
	      -a12*det3(a21, a23, a24, a31, a33, a34, a41, a43, a44)
	      +a13*det3(a21, a22, a24, a31, a32, a34, a41, a42, a44)
	      -a14*det3(a21, a22, a23, a31, a32, a33, a41, a42, a43);
}

#endif // _DETERMINANT_H_
