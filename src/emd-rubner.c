/*
 emd.c
Last update: 3/14/98
Modified by Simon Urbanek: 2011/02/28
- added extrapolation support
- add pluggable cost computing functions
Modified by Simon Urbanek: 2011/02/28
- changed code to use dynamically allocated structures
in order to remove the limit on signature sizes    
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EMD_RUBNER_MAIN 1
#include "emd-rubner.h"

#define R_NO_REMAP 1
#include <R_ext/Error.h>
#include <R_ext/Print.h>

#ifndef  DEBUG_LEVEL
#define DEBUG_LEVEL 0
#endif
/*
DEBUG_LEVEL:
0 = NO MESSAGES
1 = PRINT THE NUMBER OF ITERATIONS AND THE FINAL RESULT
2 = PRINT THE RESULT AFTER EVERY ITERATION
3 = PRINT ALSO THE FLOW AFTER EVERY ITERATION
4 = PRINT A LOT OF INFORMATION (PROBABLY USEFUL ONLY FOR THE AUTHOR)
*/


#define MAX_SIG_SIZE1  (MAX_SIG_SIZE+1)  /* FOR THE POSIBLE DUMMY FEATURE */

/* NEW TYPES DEFINITION */

/* node1_t IS USED FOR SINGLE-LINKED LISTS */
typedef struct node1_t {
  int i;
  double val;
  struct node1_t *Next;
} node1_t;

/* node1_t IS USED FOR DOUBLE-LINKED LISTS */
typedef struct node2_t {
  int i, j;
  double val;
  struct node2_t *NextC;               /* NEXT COLUMN */
struct node2_t *NextR;               /* NEXT ROW */
} node2_t;

static void *mem_alloc(long size) {
  void *m = malloc(size);
  if (!m) Rf_error("Unable to memory (%.1f MB) in emdist", ((double)size) / (1024.0 * 1024.0));
  return m;
}

static void mem_free(void *m) {
  free(m);
}

#define SMALL_SIG 512  /* signatures smaller than this are handled in
static memory, for larger ones we'll need
to perform memory allocation (slower). It
doesn't really make sense to make this too
big, since longer signatures also take
longer to compute so memory allocation is
the lesser problem. */

static node1_t local_U[SMALL_SIG], local_V[SMALL_SIG];

#define MAT(A, M, N) A[(M) + ((N) * max_sig)]

/* GLOBAL VARIABLE DECLARATION */
static int max_sig; /* size of all structures */
static int _n1, _n2;                          /* SIGNATURES SIZES */
static float *_C;                             /* THE COST MATRIX */
static float local_C[SMALL_SIG * SMALL_SIG];

static node2_t *_X;                           /* THE BASIC VARIABLES VECTOR */
static node2_t local_X[SMALL_SIG * 2];

/* VARIABLES TO HANDLE _X EFFICIENTLY */
static node2_t *_EndX, *_EnterX;
static char *_IsX;
static char local_IsX[SMALL_SIG * SMALL_SIG];
static node2_t **_RowsX, **_ColsX;
static node2_t *local_RowsX[SMALL_SIG], *local_ColsX[SMALL_SIG];
static node2_t **Loop;
static node2_t *local_Loop[SMALL_SIG * 2];

static char *IsUsed;
static char local_IsUsed[SMALL_SIG * 2];
static node1_t local_Ur[SMALL_SIG], local_Vr[SMALL_SIG];
static double local_Delta[SMALL_SIG * SMALL_SIG];
static double *Delta = local_Delta;
static node1_t *Ur = local_Ur, *Vr = local_Vr;


static double _maxW;
static float _maxC;

/* DECLARATION OF FUNCTIONS */
static float init(signature_t *Signature1, signature_t *Signature2, int extrapolate, dist_fn_t *dfn);
static void findBasicVariables(node1_t *U, node1_t *V);
static int isOptimal(node1_t *U, node1_t *V);
static int findLoop(node2_t **Loop);
static void newSol();
static void russel(double *S, double *D);
static void addBasicVariable(int minI, int minJ, double *S, double *D, 
                             node1_t *PrevUMinI, node1_t *PrevVMinJ,
                             node1_t *UHead);
#if DEBUG_LEVEL > 0
static void printSolution();
#endif


static void free_globals() {
  if (_C != local_C) mem_free(_C);
  if (_X != local_X) mem_free(_X);
  if (_IsX != local_IsX) mem_free(_IsX);
  if (_RowsX != local_RowsX) mem_free(_RowsX);
  if (Loop != local_Loop) mem_free(Loop);
  if (IsUsed != local_IsUsed) mem_free(IsUsed);
  if (Delta != local_Delta) mem_free(Delta);
  if (Ur != local_Ur) mem_free(Ur);
}

/******************************************************************************
float emd(signature_t *Signature1, signature_t *Signature2,
          flow_t *Flow, int *FlowSize, int extrapolate, dist_fn_t *dfn)

where
Signature1, Signature2  Pointers to signatures that their distance we want
to compute.
Dist       Pointer to the ground distance. i.e. the function that computes
the distance between two features.
Flow       (Optional) Pointer to a vector of flow_t (defined in emd.h) 
where the resulting flow will be stored. Flow must have n1+n2-1
elements, where n1 and n2 are the sizes of the two signatures
respectively.
If NULL, the flow is not returned.
FlowSize   (Optional) Pointer to an integer where the number of elements in
Flow will be stored
extrapolate  if set to 1 or 2 extrapolates the distance to the full mass
of 1 or 2 assuming truncated signature in the other sig. This has
any effect only if the mass of the other signature is larger

******************************************************************************/

float emd_rubner(signature_t *Signature1, signature_t *Signature2,
                 flow_t *Flow, int *FlowSize, int extrapolate, dist_fn_t *dfn)
{
  int itr;
  double totalCost;
  float w;
  node2_t *XP;
  flow_t *FlowP;
  node1_t *U = local_U, *V = local_V;
  max_sig = ((Signature1->n > Signature2->n) ? Signature1->n : Signature2->n) + 1;
  if (max_sig > SMALL_SIG) {
    U = (node1_t*) mem_alloc(sizeof(node1_t) * (max_sig * 2));
    V = U + max_sig;
  }
  
  w = init(Signature1, Signature2, extrapolate, dfn);
  
#if DEBUG_LEVEL > 1
  Rprintf("\nINITIAL SOLUTION:\n");
  printSolution();
#endif
  
  if (_n1 > 1 && _n2 > 1)  /* IF _n1 = 1 OR _n2 = 1 THEN WE ARE DONE */
  {
    for (itr = 1; itr < MAX_ITERATIONS; itr++)
    {
      /* FIND BASIC VARIABLES */
      findBasicVariables(U, V);
      
      /* CHECK FOR OPTIMALITY */
      if (isOptimal(U, V))
        break;
      
      /* IMPROVE SOLUTION */
      newSol();
      
#if DEBUG_LEVEL > 1
      Rprintf("\nITERATION # %d \n", itr);
      printSolution();
#endif
    }
    
    if (itr == MAX_ITERATIONS)
      Rf_warning("emd: Maximum number of iterations has been reached (%d)",
                 MAX_ITERATIONS);
  }
  
  /* COMPUTE THE TOTAL FLOW */
  totalCost = 0;
  if (Flow != NULL)
    FlowP = Flow;
  for(XP=_X; XP < _EndX; XP++)
  {
    if (XP == _EnterX)  /* _EnterX IS THE EMPTY SLOT */
  continue;
    if (XP->i == Signature1->n || XP->j == Signature2->n)  /* DUMMY FEATURE */
  continue;
    
    if (XP->val == 0)  /* ZERO FLOW */
  continue;
    
    totalCost += (double)XP->val * MAT(_C, XP->i, XP->j);
    if (Flow != NULL)
    {
      FlowP->from = XP->i;
      FlowP->to = XP->j;
      FlowP->amount = XP->val;
      FlowP++;
    }
  }
  if (Flow != NULL)
    *FlowSize = FlowP-Flow;
  
#if DEBUG_LEVEL > 0
  Rprintf("\n*** OPTIMAL SOLUTION (%d ITERATIONS): %f ***\n", itr, totalCost);
#endif
  
  if (U != local_U) mem_free(U);
  free_globals();
  
  /* RETURN THE NORMALIZED COST == EMD */
  /*return (float)(totalCost / w);*/
  printf("%f", totalCost)
}

/* BEGIN NEW.SU */
static float dist_L2(feature_t *a, feature_t *b) {
  float d = 0.0;
  int i = 0;
  while (i < FDIM) {
    float s = a->loc[i] - b->loc[i];
    d += s * s;
    i++;
  }
  return sqrtf(d);
}

static float dist_L1(feature_t *a, feature_t *b) {
  float d = 0.0;
  int i = 0;
  while (i < FDIM) {
    float s = a->loc[i] - b->loc[i];
    if (s > 0) 
      d += s;
    else
      d -= s;
    i++;
  }
  return d;
}

float calc_dist_L2(signature_t *Signature1, signature_t *Signature2) {
  feature_t *P1, *P2;
  float x = 0;
  int i, j;
  int n1 = Signature1->n, n2 = Signature2->n;
  for(i = 0, P1 = Signature1->Features; i < n1; i++, P1++)
    for(j = 0, P2 = Signature2->Features; j < n2; j++, P2++) 
    {
      MAT(_C, i, j) = dist_L2(P1, P2);
      if (MAT(_C, i, j) > x)
        x = MAT(_C, i, j);
    }
    return x;
}

float calc_dist_L1(signature_t *Signature1, signature_t *Signature2) {
  feature_t *P1, *P2;
  float x = 0;
  int i, j;
  int n1 = Signature1->n, n2 = Signature2->n;
  for(i = 0, P1 = Signature1->Features; i < n1; i++, P1++)
    for(j = 0, P2 = Signature2->Features; j < n2; j++, P2++) 
    {
      MAT(_C, i, j) = dist_L1(P1, P2);
      if (MAT(_C, i, j) > x)
        x = MAT(_C, i, j);
    }
    return x;
}

/* this is much slower since it involves function calls, but then it's more flexible ... */
static dist_t *default_dist;

void set_default_dist(dist_t *fn) {
  default_dist = fn;
}

float calc_dist_default(signature_t *Signature1, signature_t *Signature2) {
  feature_t *P1, *P2;
  float x = 0;
  int i, j;
  int n1 = Signature1->n, n2 = Signature2->n;
  for(i = 0, P1 = Signature1->Features; i < n1; i++, P1++)
    for(j = 0, P2 = Signature2->Features; j < n2; j++, P2++) 
    {
      MAT(_C, i, j) = default_dist(P1, P2);
      if (MAT(_C, i, j) > x)
        x = MAT(_C, i, j);
    }
    return x;
}


/* END NEW.SU */



/**********************
init
**********************/

static double local_S[SMALL_SIG], local_D[SMALL_SIG];

static float init(signature_t *Signature1, signature_t *Signature2, int extrapolate, dist_fn_t *dfn)
{
  int i, j;
  double sSum, dSum, diff;
  double *S, *D;
  
  _n1 = Signature1->n;
  _n2 = Signature2->n;
  
  max_sig = ((_n1 > _n2) ? _n1 : _n2) + 1;
  /* if the signatures are too big, use allocated memory, otherwise use local static memory */
  if (max_sig > SMALL_SIG) {
    S  = (double*) mem_alloc(sizeof(double) * max_sig * 2);
    D  = S + max_sig;
    _C = (float*) mem_alloc(sizeof(float) * max_sig * max_sig);
    _X = (node2_t*) mem_alloc(sizeof(node2_t) * max_sig * 2);
    Loop = (node2_t**) mem_alloc(sizeof(node2_t*) * max_sig * 2);
    _IsX = (char*) mem_alloc(max_sig * max_sig);
    _RowsX = (node2_t**) mem_alloc(sizeof(node2_t*) * max_sig * 2);
    _ColsX = _RowsX + max_sig;
    IsUsed = (char*) mem_alloc(max_sig * 2);
    Delta = (double*) mem_alloc(sizeof(double) * max_sig * max_sig);
    Ur = (node1_t*) mem_alloc(sizeof(node1_t) * max_sig * 2);
    Vr = Ur + max_sig;
  } else {
    S  = local_S;
    D  = local_D;
    _C = local_C;
    _X = local_X;
    Loop = local_Loop;
    _IsX = local_IsX;
    _RowsX = local_RowsX;
    _ColsX = local_ColsX;
    IsUsed = local_IsUsed;
    Delta = local_Delta;
    Ur = local_Ur;
    Vr = local_Vr;
  }
  
  /* COMPUTE THE DISTANCE MATRIX */
  _maxC = dfn(Signature1, Signature2);
  
  /* SUM UP THE SUPPLY AND DEMAND */
  sSum = 0.0;
  for(i=0; i < _n1; i++)
  {
    S[i] = Signature1->Weights[i];
    sSum += Signature1->Weights[i];
    _RowsX[i] = NULL;
  }
  dSum = 0.0;
  for(j=0; j < _n2; j++)
  {
    D[j] = Signature2->Weights[j];
    dSum += Signature2->Weights[j];
    _ColsX[j] = NULL;
  }
  
  /* IF SUPPLY DIFFERENT THAN THE DEMAND, ADD A ZERO-COST DUMMY CLUSTER */
  diff = sSum - dSum;
  if (fabs(diff) >= EPSILON * sSum)
  {
    if (diff < 0.0)
    {
      for (j=0; j < _n2; j++)
        MAT(_C, _n1, j) = 0;
      S[_n1] = -diff;
      _RowsX[_n1] = NULL;
      _n1++;
    }
    else
    {
      for (i=0; i < _n1; i++)
        MAT(_C, i, _n2) = 0;
      D[_n2] = diff;
      _ColsX[_n2] = NULL;
      _n2++;
    }
  }
  
  /* INITIALIZE THE BASIC VARIABLE STRUCTURES */
  memset(_IsX, 0, max_sig * max_sig);
  
  _EndX = _X;
  
  _maxW = sSum > dSum ? sSum : dSum;
  
  /* FIND INITIAL SOLUTION */
  russel(S, D);
  
  _EnterX = _EndX++;  /* AN EMPTY SLOT (ONLY _n1+_n2-1 BASIC VARIABLES) */
  
  if (S != local_S) mem_free(S);
  
  if (extrapolate == 1 && sSum > dSum) return dSum * dSum / sSum;
  if (extrapolate == 2 && dSum > sSum) return sSum * sSum / dSum;
  return sSum > dSum ? dSum : sSum;
}


/**********************
findBasicVariables
**********************/
static void findBasicVariables(node1_t *U, node1_t *V)
{
  int i, j, found;
  int UfoundNum, VfoundNum;
  node1_t u0Head, u1Head, *CurU, *PrevU;
  node1_t v0Head, v1Head, *CurV, *PrevV;
  
  /* INITIALIZE THE ROWS LIST (U) AND THE COLUMNS LIST (V) */
  u0Head.Next = CurU = U;
  for (i=0; i < _n1; i++)
  {
    CurU->i = i;
    CurU->Next = CurU+1;
    CurU++;
  }
  (--CurU)->Next = NULL;
  u1Head.Next = NULL;
  
  CurV = V+1;
  v0Head.Next = _n2 > 1 ? V+1 : NULL;
  for (j=1; j < _n2; j++)
  {
    CurV->i = j;
    CurV->Next = CurV+1;
    CurV++;
  }
  (--CurV)->Next = NULL;
  v1Head.Next = NULL;
  
  /* THERE ARE _n1+_n2 VARIABLES BUT ONLY _n1+_n2-1 INDEPENDENT EQUATIONS,
  SO SET V[0]=0 */
  V[0].i = 0;
  V[0].val = 0;
  v1Head.Next = V;
  v1Head.Next->Next = NULL;
  
  /* LOOP UNTIL ALL VARIABLES ARE FOUND */
  UfoundNum=VfoundNum=0;
  while (UfoundNum < _n1 || VfoundNum < _n2)
  {
    
#if DEBUG_LEVEL > 3
    Rprintf("UfoundNum=%d/%d,VfoundNum=%d/%d\n",UfoundNum,_n1,VfoundNum,_n2);
    Rprintf("U0=");
    for(CurU = u0Head.Next; CurU != NULL; CurU = CurU->Next)
      Rprintf("[%d]",CurU-U);
    Rprintf("\n");
    Rprintf("U1=");
    for(CurU = u1Head.Next; CurU != NULL; CurU = CurU->Next)
      Rprintf("[%d]",CurU-U);
    Rprintf("\n");
    Rprintf("V0=");
    for(CurV = v0Head.Next; CurV != NULL; CurV = CurV->Next)
      Rprintf("[%d]",CurV-V);
    Rprintf("\n");
    Rprintf("V1=");
    for(CurV = v1Head.Next; CurV != NULL; CurV = CurV->Next)
      Rprintf("[%d]",CurV-V);
    Rprintf("\n\n");
#endif
    
    found = 0;
    if (VfoundNum < _n2)
    {
      /* LOOP OVER ALL MARKED COLUMNS */
      PrevV = &v1Head;
      for (CurV=v1Head.Next; CurV != NULL; CurV=CurV->Next)
      {
        j = CurV->i;
        /* FIND THE VARIABLES IN COLUMN j */
        PrevU = &u0Head;
        for (CurU=u0Head.Next; CurU != NULL; CurU=CurU->Next)
        {
          i = CurU->i;
          if (MAT(_IsX, i, j))
          {
            /* COMPUTE U[i] */
            CurU->val = MAT(_C, i, j) - CurV->val;
            /* ...AND ADD IT TO THE MARKED LIST */
            PrevU->Next = CurU->Next;
            CurU->Next = u1Head.Next != NULL ? u1Head.Next : NULL;
            u1Head.Next = CurU;
            CurU = PrevU;
          }
          else
            PrevU = CurU;
        }
        PrevV->Next = CurV->Next;
        VfoundNum++;
        found = 1;
      }
    }
    if (UfoundNum < _n1)
    {
      /* LOOP OVER ALL MARKED ROWS */
      PrevU = &u1Head;
      for (CurU=u1Head.Next; CurU != NULL; CurU=CurU->Next)
      {
        i = CurU->i;
        /* FIND THE VARIABLES IN ROWS i */
        PrevV = &v0Head;
        for (CurV=v0Head.Next; CurV != NULL; CurV=CurV->Next)
        {
          j = CurV->i;
          if (MAT(_IsX, i, j))
          {
            /* COMPUTE V[j] */
            CurV->val = MAT(_C, i, j) - CurU->val;
            /* ...AND ADD IT TO THE MARKED LIST */
            PrevV->Next = CurV->Next;
            CurV->Next = v1Head.Next != NULL ? v1Head.Next: NULL;
            v1Head.Next = CurV;
            CurV = PrevV;
          }
          else
            PrevV = CurV;
        }
        PrevU->Next = CurU->Next;
        UfoundNum++;
        found = 1;
      }
    }
    if (! found)
      Rf_error("emd: Unexpected error in findBasicVariables!\nThis typically happens when the EPSILON defined in emd-rubner.h is not right for the scale of the problem.");
  }
}




/**********************
isOptimal
**********************/
static int isOptimal(node1_t *U, node1_t *V)
{    
  double delta, deltaMin;
  int i, j, minI, minJ;
  
  /* FIND THE MINIMAL Cij-Ui-Vj OVER ALL i,j */
  deltaMin = EMD_INFINITY;
  for(i=0; i < _n1; i++)
    for(j=0; j < _n2; j++)
      if (! MAT(_IsX, i, j))
      {
        delta = MAT(_C, i, j) - U[i].val - V[j].val;
        if (deltaMin > delta)
        {
          deltaMin = delta;
          minI = i;
          minJ = j;
        }
      }
      
#if DEBUG_LEVEL > 3
      Rprintf("deltaMin=%f\n", deltaMin);
#endif
      
      if (deltaMin == EMD_INFINITY)
        Rf_error("emd: Unexpected error in isOptimal.");
      
      _EnterX->i = minI;
      _EnterX->j = minJ;
      
      /* IF NO NEGATIVE deltaMin, WE FOUND THE OPTIMAL SOLUTION */
      return deltaMin >= -EPSILON * _maxC;
      
      /*
      return deltaMin >= -EPSILON;
      */
}

/**********************
newSol
**********************/
static void newSol()
{
  int i, j, k;
  double xMin;
  int steps;
  node2_t *CurX, *LeaveX;
  
#if DEBUG_LEVEL > 3
  Rprintf("EnterX = (%d,%d)\n", _EnterX->i, _EnterX->j);
#endif
  
  /* ENTER THE NEW BASIC VARIABLE */
  i = _EnterX->i;
  j = _EnterX->j;
  MAT(_IsX, i, j) = 1;
  _EnterX->NextC = _RowsX[i];
  _EnterX->NextR = _ColsX[j];
  _EnterX->val = 0;
  _RowsX[i] = _EnterX;
  _ColsX[j] = _EnterX;
  
  /* FIND A CHAIN REACTION */
  steps = findLoop(Loop);
  
  /* FIND THE LARGEST VALUE IN THE LOOP */
  xMin = EMD_INFINITY;
  for (k=1; k < steps; k+=2)
  {
    if (Loop[k]->val < xMin)
    {
      LeaveX = Loop[k];
      xMin = Loop[k]->val;
    }
  }
  
  /* UPDATE THE LOOP */
  for (k=0; k < steps; k+=2)
  {
    Loop[k]->val += xMin;
    Loop[k+1]->val -= xMin;
  }
  
#if DEBUG_LEVEL > 3
  Rprintf("LeaveX = (%d,%d)\n", LeaveX->i, LeaveX->j);
#endif
  
  /* REMOVE THE LEAVING BASIC VARIABLE */
  i = LeaveX->i;
  j = LeaveX->j;
  MAT(_IsX, i, j) = 0;
  if (_RowsX[i] == LeaveX)
    _RowsX[i] = LeaveX->NextC;
  else
    for (CurX=_RowsX[i]; CurX != NULL; CurX = CurX->NextC)
      if (CurX->NextC == LeaveX)
      {
        CurX->NextC = CurX->NextC->NextC;
        break;
      }
      if (_ColsX[j] == LeaveX)
        _ColsX[j] = LeaveX->NextR;
      else
        for (CurX=_ColsX[j]; CurX != NULL; CurX = CurX->NextR)
          if (CurX->NextR == LeaveX)
          {
            CurX->NextR = CurX->NextR->NextR;
            break;
          }
          
          /* SET _EnterX TO BE THE NEW EMPTY SLOT */
          _EnterX = LeaveX;
}



/**********************
findLoop
**********************/
static int findLoop(node2_t **Loop)
{
  int steps;
  node2_t **CurX, *NewX;
  
  memset(IsUsed, 0, _n1 + _n2);
  
  CurX = Loop;
  NewX = *CurX = _EnterX;
  IsUsed[_EnterX-_X] = 1;
  steps = 1;
  
  do
  {
    if (steps%2 == 1)
    {
      /* FIND AN UNUSED X IN THE ROW */
      NewX = _RowsX[NewX->i];
      while (NewX != NULL && IsUsed[NewX-_X])
        NewX = NewX->NextC;
    }
    else
    {
      /* FIND AN UNUSED X IN THE COLUMN, OR THE ENTERING X */
      NewX = _ColsX[NewX->j];
      while (NewX != NULL && IsUsed[NewX-_X] && NewX != _EnterX)
        NewX = NewX->NextR;
      if (NewX == _EnterX)
        break;
    }
    
    if (NewX != NULL)  /* FOUND THE NEXT X */
    {
      /* ADD X TO THE LOOP */
      *++CurX = NewX;
      IsUsed[NewX-_X] = 1;
      steps++;
#if DEBUG_LEVEL > 3
      Rprintf("steps=%d, NewX=(%d,%d)\n", steps, NewX->i, NewX->j);    
#endif
    }
    else  /* DIDN'T FIND THE NEXT X */
    {
      /* BACKTRACK */
      do
      {
        NewX = *CurX;
        do 
        {
          if (steps%2 == 1)
            NewX = NewX->NextR;
          else
            NewX = NewX->NextC;
        } while (NewX != NULL && IsUsed[NewX-_X]);
        
        if (NewX == NULL)
        {
          IsUsed[*CurX-_X] = 0;
          CurX--;
          steps--;
        }
      } while (NewX == NULL && CurX >= Loop);
      
#if DEBUG_LEVEL > 3
      Rprintf("BACKTRACKING TO: steps=%d, NewX=(%d,%d)\n",
              steps, NewX->i, NewX->j);    
#endif
      IsUsed[*CurX-_X] = 0;
      *CurX = NewX;
      IsUsed[NewX-_X] = 1;
    }     
  } while(CurX >= Loop);
  
  if (CurX == Loop)
    Rf_error("emd: Unexpected error in findLoop!");
#if DEBUG_LEVEL > 3
{
  int i;
  Rprintf("FOUND LOOP:\n");
  for (i=0; i < steps; i++)
    Rprintf("%d: (%d,%d)\n", i, Loop[i]->i, Loop[i]->j);
}
#endif

return steps;
}

/**********************
russel
**********************/
static void russel(double *S, double *D)
{
  int i, j, found, minI, minJ;
  double deltaMin, oldVal, diff;
  node1_t uHead, *CurU, *PrevU;
  node1_t vHead, *CurV, *PrevV;
  node1_t *PrevUMinI = PrevUMinI, *PrevVMinJ, *Remember;
  
  
  
  /* INITIALIZE THE ROWS LIST (Ur), AND THE COLUMNS LIST (Vr) */
  uHead.Next = CurU = Ur;
  for (i=0; i < _n1; i++)
  {
    CurU->i = i;
    CurU->val = -EMD_INFINITY;
    CurU->Next = CurU+1;
    CurU++;
  }
  (--CurU)->Next = NULL;
  
  vHead.Next = CurV = Vr;
  for (j=0; j < _n2; j++)
  {
    CurV->i = j;
    CurV->val = -EMD_INFINITY;
    CurV->Next = CurV+1;
    CurV++;
  }
  (--CurV)->Next = NULL;
  
  /* FIND THE MAXIMUM ROW AND COLUMN VALUES (Ur[i] AND Vr[j]) */
  for(i=0; i < _n1 ; i++)
    for(j=0; j < _n2 ; j++)
    {
      float v;
      v = MAT(_C, i, j);
      if (Ur[i].val <= v)
        Ur[i].val = v;
      if (Vr[j].val <= v)
        Vr[j].val = v;
    }
    
    /* COMPUTE THE Delta MATRIX */
    for(i=0; i < _n1 ; i++)
      for(j=0; j < _n2 ; j++)
        MAT(Delta, i, j) = MAT(_C, i, j) - Ur[i].val - Vr[j].val;
  
  /* FIND THE BASIC VARIABLES */
  do
  {
#if DEBUG_LEVEL > 3
    Rprintf("Ur=");
    for(CurU = uHead.Next; CurU != NULL; CurU = CurU->Next)
      Rprintf("[%d]",CurU-Ur);
    Rprintf("\n");
    Rprintf("Vr=");
    for(CurV = vHead.Next; CurV != NULL; CurV = CurV->Next)
      Rprintf("[%d]",CurV-Vr);
    Rprintf("\n");
    Rprintf("\n\n");
#endif
    
    /* FIND THE SMALLEST Delta[i][j] */
    found = 0; 
    deltaMin = EMD_INFINITY;      
    PrevU = &uHead;
    for (CurU=uHead.Next; CurU != NULL; CurU=CurU->Next)
    {
      int i;
      i = CurU->i;
      PrevV = &vHead;
      for (CurV=vHead.Next; CurV != NULL; CurV=CurV->Next)
      {
        int j;
        j = CurV->i;
        if (deltaMin > MAT(Delta, i, j))
        {
          deltaMin = MAT(Delta, i, j);
          minI = i;
          minJ = j;
          PrevUMinI = PrevU;
          PrevVMinJ = PrevV;
          found = 1;
        }
        PrevV = CurV;
      }
      PrevU = CurU;
    }
    
    if (! found)
      break;
    
    /* ADD X[minI][minJ] TO THE BASIS, AND ADJUST SUPPLIES AND COST */
    Remember = PrevUMinI->Next;
    addBasicVariable(minI, minJ, S, D, PrevUMinI, PrevVMinJ, &uHead);
    
    /* UPDATE THE NECESSARY Delta[][] */
    if (Remember == PrevUMinI->Next)  /* LINE minI WAS DELETED */
    {
      for (CurV=vHead.Next; CurV != NULL; CurV=CurV->Next)
      {
        int j;
        j = CurV->i;
        if (CurV->val == MAT(_C, minI, j))  /* COLUMN j NEEDS UPDATING */
        {
          /* FIND THE NEW MAXIMUM VALUE IN THE COLUMN */
          oldVal = CurV->val;
          CurV->val = -EMD_INFINITY;
          for (CurU=uHead.Next; CurU != NULL; CurU=CurU->Next)
          {
            int i;
            i = CurU->i;
            if (CurV->val <= MAT(_C, i, j))
              CurV->val = MAT(_C, i, j);
          }
          
          /* IF NEEDED, ADJUST THE RELEVANT Delta[*][j] */
          diff = oldVal - CurV->val;
          if (fabs(diff) < EPSILON * _maxC)
            for (CurU=uHead.Next; CurU != NULL; CurU=CurU->Next)
              MAT(Delta, CurU->i, j) += diff;
        }
      }
    }
    else  /* COLUMN minJ WAS DELETED */
    {
      for (CurU=uHead.Next; CurU != NULL; CurU=CurU->Next)
      {
        int i;
        i = CurU->i;
        if (CurU->val == MAT(_C, i, minJ))  /* ROW i NEEDS UPDATING */
        {
          /* FIND THE NEW MAXIMUM VALUE IN THE ROW */
          oldVal = CurU->val;
          CurU->val = -EMD_INFINITY;
          for (CurV=vHead.Next; CurV != NULL; CurV=CurV->Next)
          {
            int j;
            j = CurV->i;
            if(CurU->val <= MAT(_C, i, j))
              CurU->val = MAT(_C, i, j);
          }
          
          /* If NEEDED, ADJUST THE RELEVANT Delta[i][*] */
          diff = oldVal - CurU->val;
          if (fabs(diff) < EPSILON * _maxC)
            for (CurV=vHead.Next; CurV != NULL; CurV=CurV->Next)
              MAT(Delta, i, CurV->i) += diff;
        }
      }
    }
  } while (uHead.Next != NULL || vHead.Next != NULL);
}




/**********************
addBasicVariable
**********************/
static void addBasicVariable(int minI, int minJ, double *S, double *D, 
                             node1_t *PrevUMinI, node1_t *PrevVMinJ,
                             node1_t *UHead)
{
  double T;
  
  if (fabs(S[minI]-D[minJ]) <= EPSILON * _maxW)  /* DEGENERATE CASE */
  {
    T = S[minI];
    S[minI] = 0;
    D[minJ] -= T; 
  }
  else if (S[minI] < D[minJ])  /* SUPPLY EXHAUSTED */
  {
    T = S[minI];
    S[minI] = 0;
    D[minJ] -= T; 
  }
  else  /* DEMAND EXHAUSTED */
  {
    T = D[minJ];
    D[minJ] = 0; 
    S[minI] -= T; 
  }
  
  /* X(minI,minJ) IS A BASIC VARIABLE */
  MAT(_IsX, minI, minJ) = 1; 
  
  _EndX->val = T;
  _EndX->i = minI;
  _EndX->j = minJ;
  _EndX->NextC = _RowsX[minI];
  _EndX->NextR = _ColsX[minJ];
  _RowsX[minI] = _EndX;
  _ColsX[minJ] = _EndX;
  _EndX++;
  
  /* DELETE SUPPLY ROW ONLY IF THE EMPTY, AND IF NOT LAST ROW */
  if (S[minI] == 0 && UHead->Next->Next != NULL)
    PrevUMinI->Next = PrevUMinI->Next->Next;  /* REMOVE ROW FROM LIST */
  else
    PrevVMinJ->Next = PrevVMinJ->Next->Next;  /* REMOVE COLUMN FROM LIST */
}


#if DEBUG_LEVEL > 0

/**********************
 printSolution
 **********************/
static void printSolution()
{
  node2_t *P;
  double totalCost;
  
  totalCost = 0;
  
#if DEBUG_LEVEL > 2
  Rprintf("SIG1\tSIG2\tFLOW\tCOST\n");
#endif
  for(P=_X; P < _EndX; P++)
    if (P != _EnterX && MAT(_IsX, P->i, P->j))
    {
#if DEBUG_LEVEL > 2
      Rprintf("%d\t%d\t%f\t%f\n", P->i, P->j, P->val, MAT(_C, P->i, P->j));
#endif
      totalCost += (double)P->val * MAT(_C, P->i, P->j);
    }
    
    Rprintf("COST = %f\n", totalCost);
}
#endif
