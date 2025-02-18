#if !defined( MESH_GEOMETRY_H )
#define MESH_GEOMETRY_H
#include <math.h>
#include <float.h>
#include <algorithm>

#include "../../config.h"
#if defined(_NO_STD_MINMAX) || defined(WIN32)
	#include "minmaxpatch.h"
#endif

inline
int inCircle(double x[4], double y[4])
{
	double x10 = x[1] - x[0];
	double x20 = x[2] - x[0];
	double x30 = x[3] - x[0];
	double x21 = x[2] - x[1];
	double x31 = x[3] - x[1];
	double y10 = y[1] - y[0];
	double y20 = y[2] - y[0];
	double y30 = y[3] - y[0];
	double y21 = y[2] - y[1];
	double y31 = y[3] - y[1];

	double a1 = x10 * y20, b1 = x20 * y10;
	double a2 = x10 * y30, b2 = x30 * y10;
	double a3 = x20 * x21, b3 = y20 * y21;
	double a4 = x30 * x31, b4 = y30 * y31;

	double t1 = a1 - b1;
	double t2 = a2 - b2;
	double t3 = a3 + b3;
	double t4 = a4 + b4;
	double t14 = t1 * t4, t23 = t2 * t3;
	double det = t14 - t23;

	if (det >= 0.0)
		return 1;

	// rounding error evaluation:
	//		sum:     a = b + c -> a_error <= b_error + c_error + max(|a|, |b|, |c|) * DBL_EPSILON
	//		product: a = b * c -> a_error <= b_error * |c| + c_error * |b| + b_error * c_error + |a| * DBL_EPSILON
	//             the term (b_error * c_error) is very small and can therefore be ignored

	// maximum rounding error of coordinate differences
	double maxcoord = std::max(std::max(std::max(fabs(x[0]), fabs(x[1])),
	                                    std::max(fabs(x[2]), fabs(x[3]))),
	                           std::max(std::max(fabs(y[0]), fabs(y[1])),
	                                    std::max(fabs(y[2]), fabs(y[3]))));

	// maximum rounding error of t1..t4
	x10 = fabs(x10);
	y10 = fabs(y10);

	x20 = fabs(x20); y20 = fabs(y20);
	a1 = fabs(a1); b1 = fabs(b1); t1 = fabs(t1);
	
	double e1 = maxcoord * (x10 + x20 + y10 + y20) +
	            a1 + b1 + std::max(std::max(a1, b1), t1);

	x30 = fabs(x30); y30 = fabs(y30);
	a2 = fabs(a2); b2 = fabs(b2); t2 = fabs(t2);

	double e2 = maxcoord * (x10 + x30 + y10 + y30) +
	            a2 + b2 + std::max(std::max(a2, b2), t2);

	x21 = fabs(x21); y21 = fabs(y21);
	a3 = fabs(a3);	b3 = fabs(b3);	t3 = fabs(t3);

	double e3 = maxcoord * (x21 + x20 + y21 + y20) +
	            a3 + b3 + std::max(std::max(a3, b3), t3);

	x31 = fabs(x31); y31 = fabs(y31);
	a4 = fabs(a4); b4 = fabs(b4);	t4 = fabs(t4);

	double e4 = maxcoord * (x31 + x30 + y31 + y30) +
	            a4 + b4 + std::max(std::max(a4, b4), t4);

	t14 = fabs(t14); t23 = fabs(t23);
	// maximum rounding error of det
	double eps = e1 * t4 + e4 * t1 + e2 * t3 + e3 * t2 +
	             t14 + t23 + std::max(std::max(t14, t23), fabs(det));

	eps *= 2.0 * DBL_EPSILON;

	if(det < -eps)
		return 0;

	return 1;
}

inline
double sideDet(double x1,double y1,double x2,double y2,double x3,double y3, double *err = NULL)
{
	double d23 = x2 * y3, d32 = x3 * y2;
	double d13 = x1 * y3, d31 = x3 * y1;
	double d12 = x1 * y2, d21 = x2 * y1;
	double d1 = d23 - d32, d2 = d13 - d31, d3 = d12 - d21;

	if (err)
	{
		// maximum rounding errors of d1, d2 and d3
		//          product 1   product 2   difference
		double e1 = fabs(d23) + fabs(d32) + std::max(std::max(fabs(d23), fabs(d32)), fabs(d1));
		double e2 = fabs(d13) + fabs(d31) + std::max(std::max(fabs(d13), fabs(d31)), fabs(d2));
		double e3 = fabs(d12) + fabs(d21) + std::max(std::max(fabs(d12), fabs(d21)), fabs(d3));

		// maximum rounding error of sideDet
		*err = (e1 + e2 + e3 + std::max(std::max(fabs(d1), fabs(d2)), fabs(d1 - d2)) +
		        std::max(std::max(fabs(d1 - d2), fabs(d3)), d1 - d2 + d3)) * 2.0 * DBL_EPSILON;
	}

	return d1 - d2 + d3;
}

inline
int orientation(double x1,double y1,double x2,double y2,double x3,double y3)
{
	double val;
	double eps;
	val = sideDet(x1,y1,x2,y2,x3,y3, &eps);

	if(val <= eps)
		return 0;
	return 1;
}

inline
double distance(double x1, double y1, double x2, double y2)
{
	double dx = x2 - x1;
	double dy = y2 - y1;
	return sqrt(dx*dx + dy*dy);
}

inline
void vertexLocation(double p1x, double p1y, double p2x, double p2y,
						  double p3x, double p3y, double *vx, double *vy)
{
	double a,b,c,d;
	double r1,r2;
	double det;
	
	a = 2 * (p2x - p1x);
	b = 2 * (p2y - p1y);
	c = 2 * (p3x - p2x);
	d = 2 * (p3y - p2y);
	r1 = -(-p2x * p2x + p1x * p1x - p2y * p2y + p1y * p1y);
	r2 = -(-p3x * p3x + p2x * p2x - p3y * p3y + p2y * p2y);
	det = a * d - b * c;
	*vx = (1 / det) * (d * r1 - b * r2);
	*vy = (1 / det) * (-c * r1 + a * r2);
}

//#define MAP(a) (a)
//#define UNMAP(a) (a)
#define MAP(a) log(a)
#define UNMAP(a) exp(a)

#endif /* MESH_GEOMETRY_H */
