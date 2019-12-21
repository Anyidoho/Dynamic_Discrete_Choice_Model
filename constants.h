#ifndef constants_h
#define constants_h

#define DEBUG 1 //comment to turn off debug mode
#define CONST_NUM_HOUSEHOLDS_FLORENCE 1027
//#define CONST_NUM_HOUSEHOLDS_MICHAEL  330
#define CONST_NUM_TIME_STEPS 24
#define CONST_NUM_COVARIATES 19 //Includes constant
#define CONST_NUM_SD_COVARIATES 11 //(Should contain attributes that are time-invariant)

#define CONST_SPACING 10
#define CONST_WIDE_SPACING 15
#define CONST_EXTRA_WIDE_SPACING 25

#define CONST_MAX_ITERATIONS 100
#define CONST_CONVERGENCE 1E-5

#define CONST_EXP 2.71828
#define CONST_INFTY 1E20
#define CONST_EPSILON 10E-15
#define CONST_LARGE_EPSILON 10E-12
#define CONST_EULER 0.577215664901532

#define CONST_NO_TRANSFORM 1
//#define CONST_LOGIT_TRANSFORM 1
//#define CONST_SQUARE_TRANSFORM 1

//#define CONST_OPTIMIZE_DFACT 1;//Comment to fix discount factor (Then make sure that only CONST_NO_TRANSFORM is turned on)

#endif 