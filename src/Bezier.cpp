#include <stdio.h>
#include <math.h>

#include "mathFunc.h"
#include "initial.hpp"


double Rx(double* b, double* q, int j) {
	int i;
	double ret = 0.0;
	double temp = 0.0;

	for (i = 0; i <= BEZIER_ORDER; i++) {
		ret += b[i] * (factorial(BEZIER_ORDER) / (factorial(BEZIER_ORDER - i) * factorial(i))) * pow(q[j], i) * pow((1 - q[j]), BEZIER_ORDER - i);
		temp += (factorial(BEZIER_ORDER) / (factorial(BEZIER_ORDER - i) * factorial(i)));
		//printf("%.7lf\n", temp);
	}
	return(ret);
}
double Ry(double* b, double* q, int j) {
	int i;
	double ret = 0.0;

	for (i = 0; i <= BEZIER_ORDER; i++) {
		ret += b[i] * (factorial(BEZIER_ORDER) / (factorial(BEZIER_ORDER - i) * factorial(i))) * pow(q[j], i) * pow((1 - q[j]), BEZIER_ORDER - i);

	}
	return(ret);
}



double d1Rxdq1(double* b, double* q, int j) {
	int i;
	double ret = 0.0;

	for (i = 1; i <= BEZIER_ORDER; i++) {
		ret += (-b[i - 1] + b[i]) * (factorial(BEZIER_ORDER) / (factorial(BEZIER_ORDER - i) * factorial(i))) * i * pow(q[j], i - 1) * pow((1 - q[j]), BEZIER_ORDER - i);

	}
	return(ret);
}
double d1Rydq1(double* b, double* q, int j) {
	int i;
	double ret = 0.0;

	for (i = 1; i <= BEZIER_ORDER; i++) {
		ret += (-b[i - 1] + b[i]) * (factorial(BEZIER_ORDER) / (factorial(BEZIER_ORDER - i) * factorial(i))) * i * pow(q[j], i - 1) * pow((1 - q[j]), BEZIER_ORDER - i);

	}
	return(ret);
}



double d2Rxdq2(double* b, double* q, int j) {
	int i;
	double ret = 0.0;

	for (i = 2; i <= BEZIER_ORDER; i++) {
		ret += (b[i] - 2 * b[i - 1] + b[i - 2])
			* (factorial(BEZIER_ORDER) / (factorial(BEZIER_ORDER - i) * factorial(i))) * i * (i - 1) * pow(q[j], i - 2) * pow((1 - q[j]), BEZIER_ORDER - i);

	}
	return(ret);
}

double d2Rydq2(double* b, double* q, int j) {
	int i;
	double ret = 0.0;

	for (i = 2; i <= BEZIER_ORDER; i++) {
		ret += (b[i] - 2 * b[i - 1] + b[i - 2])
			* (factorial(BEZIER_ORDER) / (factorial(BEZIER_ORDER - i) * factorial(i))) * i * (i - 1) * pow(q[j], i - 2) * pow((1 - q[j]), BEZIER_ORDER - i);

	}
	return(ret);
}


double d3Rxdq3(double* b, double* q, int j) {
	int i;
	double ret = 0.0;

	for (i = 3; i <= BEZIER_ORDER; i++) {
		ret += (b[i] - 3 * b[i - 1] + 3 * b[i - 2] - b[i - 3])
			* (factorial(BEZIER_ORDER) / (factorial(BEZIER_ORDER - i) * factorial(i))) * i * (i - 1) * (i - 2) * pow(q[j], i - 3) * pow((1 - q[j]), BEZIER_ORDER - i);

	}
	return(ret);
}
double d3Rydq3(double* b, double* q, int j) {
	int i;
	double ret = 0.0;

	for (i = 3; i <= BEZIER_ORDER; i++) {
		ret += (b[i] - 3 * b[i - 1] + 3 * b[i - 2] - b[i - 3])
			* (factorial(BEZIER_ORDER) / (factorial(BEZIER_ORDER - i) * factorial(i))) * i * (i - 1) * (i - 2) * pow(q[j], i - 3) * pow((1 - q[j]), BEZIER_ORDER - i);

	}
	return(ret);
}



double d4Rxdq4(double* b, double* q, int j) {
	int i;
	double ret = 0.0;

	for (i = 4; i <= BEZIER_ORDER; i++) {
		ret += (b[i] - 4 * b[i - 1] + 6 * b[i - 2] - 4 * b[i - 3] + b[i - 4])
			* (factorial(BEZIER_ORDER) / (factorial(BEZIER_ORDER - i) * factorial(i))) * i * (i - 1) * (i - 2) * (i - 3) * pow(q[j], i - 4) * pow((1 - q[j]), BEZIER_ORDER - i);

	}
	return(ret);
}
double d4Rydq4(double* b, double* q, int j) {
	int i;
	double ret = 0.0;

	for (i = 4; i <= BEZIER_ORDER; i++) {
		ret += (b[i] - 4 * b[i - 1] + 6 * b[i - 2] - 4 * b[i - 3] + b[i - 4])
			* (factorial(BEZIER_ORDER) / (factorial(BEZIER_ORDER - i) * factorial(i))) * i * (i - 1) * (i - 2) * (i - 3) * pow(q[j], i - 4) * pow((1 - q[j]), BEZIER_ORDER - i);

	}
	return(ret);
}

// 5階微分：Rx 成分
double d5Rxdq5(double* b, double* q, int j) {
    int i;
    double ret = 0.0;

    for (i = 5; i <= BEZIER_ORDER; i++) {
        // binomial 差分: b[i] - 5 b[i-1] + 10 b[i-2] - 10 b[i-3] + 5 b[i-4] - b[i-5]
        double diff =
            b[i]
          - 5.0 * b[i - 1]
          + 10.0 * b[i - 2]
          - 10.0 * b[i - 3]
          + 5.0 * b[i - 4]
          - 1.0 * b[i - 5];

        // binomial(C(BEZIER_ORDER,i))
        double coef = factorial(BEZIER_ORDER)
                    / (factorial(BEZIER_ORDER - i) * factorial(i));

        // i*(i-1)*(i-2)*(i-3)*(i-4)
        double perm = (double)i * (i - 1) * (i - 2) * (i - 3) * (i - 4);

        ret += diff
             * coef
             * perm
             * pow(q[j], i - 5)
             * pow((1.0 - q[j]), BEZIER_ORDER - i);
    }
    return ret;
}

// 5階微分：Ry 成分
double d5Rydq5(double* b, double* q, int j) {
    int i;
    double ret = 0.0;

    for (i = 5; i <= BEZIER_ORDER; i++) {
        double diff =
            b[i]
          - 5.0 * b[i - 1]
          + 10.0 * b[i - 2]
          - 10.0 * b[i - 3]
          + 5.0 * b[i - 4]
          - 1.0 * b[i - 5];

        double coef = factorial(BEZIER_ORDER)
                    / (factorial(BEZIER_ORDER - i) * factorial(i));

        double perm = (double)i * (i - 1) * (i - 2) * (i - 3) * (i - 4);

        ret += diff
             * coef
             * perm
             * pow(q[j], i - 5)
             * pow((1.0 - q[j]), BEZIER_ORDER - i);
    }
    return ret;
}


double norm(double* bx, double* by, double* q, int j) {
	double ret = 0.0;
	ret = sqrt(pow(d1Rxdq1(bx, q, j), 2) + pow(d1Rydq1(by, q, j), 2));

	return(ret);
}