#ifndef _EMD_H
#define _EMD_H
/*
    emd.h

    Last update: 3/24/98
    Modified by Simon Urbanek: 2011/02/28
    - adapted to R interface
    - add "extrapolate" parameter allowing asymmetric extrapolation

    An implementation of the Earth Movers Distance.
    Based of the solution for the Transportation problem as described in
    "Introduction to Mathematical Programming" by F. S. Hillier and 
    G. J. Lieberman, McGraw-Hill, 1990.

    Copyright (C) 1998 Yossi Rubner
    Computer Science Department, Stanford University
    E-Mail: rubner@cs.stanford.edu   URL: http://vision.stanford.edu/~rubner

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in
	all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
	THE SOFTWARE.

*/


/* DEFINITIONS */
#ifndef MAX_SIG_SIZE
#define MAX_SIG_SIZE   511
#endif
#define MAX_ITERATIONS 3000
#define EMD_INFINITY   1e20
#define EPSILON        1e-6

#define FDIM  4

typedef struct
{
  float loc[FDIM];
} feature_t;

typedef struct
{
  int n;                /* Number of features in the signature */
  feature_t *Features;  /* Pointer to the features vector */
  float *Weights;       /* Pointer to the weights of the features */
} signature_t;


typedef struct
{
  int from;             /* Feature number in signature 1 */
  int to;               /* Feature number in signature 2 */
  float amount;         /* Amount of flow from "from" to "to" */
} flow_t;

/* NEW.SU - define a function that calculates all base distances between signatures */
/*          calculates cost matrix _C and returns the maximum */
typedef float (dist_fn_t)(signature_t *, signature_t *);
/* a less-efficient method is to set default_dist and use calc_dist_default */
typedef float (dist_t)(feature_t *, feature_t *);

/* END.SU */

float emd_rubner(signature_t *Signature1, signature_t *Signature2,
		 flow_t *Flow, int *FlowSize, int extrapolate, dist_fn_t *dfn);

float calc_dist_L2(signature_t *Signature1, signature_t *Signature2);
float calc_dist_L1(signature_t *Signature1, signature_t *Signature2);
float calc_dist_default(signature_t *Signature1, signature_t *Signature2);

void set_default_dist(dist_t *fn);

#endif
