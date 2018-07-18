/* this version of the EMD uses Rubner's code */
#include "emd-rubner.h"

#include <stdlib.h>
#include <string.h>

#define R_NO_REMAP 1
#define USE_RINTERNALS 1
#include <Rinternals.h>

static SEXP dist_clos, cf1, cf2; /* dist_clos is the closure, cf1/2 are cached vectors that we allocate only once and re-use */

static float eval_dist(feature_t *f1, feature_t *f2) {
    double *x = REAL(cf1), *y = REAL(cf2);
    int i;
    for (i = 0; i < FDIM; i++) {
	x[i] = f1->loc[i];
	y[i] = f2->loc[i];
    }
    SEXP res = Rf_eval(Rf_lang3(dist_clos, cf1, cf2), R_GlobalEnv);
    if (TYPEOF(res) == INTSXP && LENGTH(res) == 1)
	return (float) (INTEGER(res)[0]);
    if (TYPEOF(res) != REALSXP || LENGTH(res) != 1)
	Rf_error("invalid distance result - must be a numeric vector of length one");
    return (float)(REAL(res)[0]);
}


SEXP emd_r(SEXP sBase, SEXP sCur, SEXP sExtra, SEXP sFlows, SEXP sDist) {
  SEXP sBaseDim = Rf_getAttrib(sBase, R_DimSymbol);
  SEXP sCurDim = Rf_getAttrib(sCur, R_DimSymbol);
  if (sBaseDim == R_NilValue || LENGTH(sBaseDim) != 2) Rf_error("base must be a matrix");
  if (sCurDim  == R_NilValue || LENGTH(sCurDim)  != 2) Rf_error("cur must be a matrix");
  int *baseDim = INTEGER(sBaseDim);
  int *curDim = INTEGER(sCurDim);
  int baseRows = baseDim[0], baseCol = baseDim[1];
  int curRows = curDim[0], curCol = curDim[1];
  if (TYPEOF(sDist) != CLOSXP && (TYPEOF(sDist) != STRSXP || LENGTH(sDist) != 1)) Rf_error("invalid distance specification");
  const char *distName = (TYPEOF(sDist) == STRSXP) ? CHAR(STRING_ELT(sDist, 0)) : 0;
  dist_fn_t *dist_fn = 0;
  if (!distName) {
      dist_fn = calc_dist_default;
      set_default_dist(eval_dist);
      dist_clos = sDist;
      cf1 = PROTECT(Rf_allocVector(REALSXP, FDIM));
      cf2 = PROTECT(Rf_allocVector(REALSXP, FDIM));
  } else {
      if (!strcmp(distName, "euclidean")) dist_fn = calc_dist_L2;
      if (!strcmp(distName, "manhattan")) dist_fn = calc_dist_L1;
  }
  if (!dist_fn)
      Rf_error("invalid distance specification");
  sBase = Rf_coerceVector(sBase, REALSXP);
  sCur = Rf_coerceVector(sCur, REALSXP);
  double *baseVal = REAL(sBase);
  double *curVal = REAL(sCur);
  flow_t *flows = NULL;
  int n_flows = 0;

  if (baseCol != curCol) Rf_error("base and current sets must have the same dimensionality");
  if (baseCol < 2) Rf_error("at least two columns (weight and location) are required");
  if (baseCol > FDIM + 1) Rf_warning("more than %d dimensions are used, those will be ignored", FDIM);

  signature_t baseSig, curSig;
  
  baseSig.n = baseRows;
  baseSig.Features = (feature_t*) R_alloc(baseRows, sizeof(feature_t));
  baseSig.Weights  = (float*) R_alloc(baseRows, sizeof(float));
  curSig.n = curRows;
  curSig.Features = (feature_t*) R_alloc(curRows, sizeof(feature_t));
  curSig.Weights  = (float*) R_alloc(curRows, sizeof(float));

  int i, j;
  for (i = 0; i < baseRows; i++) {
    for (j = 0; j < FDIM; j++)
      baseSig.Features[i].loc[j] = (j + 1 < baseCol) ? baseVal[i + (j + 1) * baseRows] : 0.0;
    baseSig.Weights[i] = baseVal[i];
  }
  for (i = 0; i < curRows; i++) {
    for (j = 0; j < FDIM; j++)
      curSig.Features[i].loc[j] = (j + 1 < curCol) ? curVal[i + (j + 1) * curRows] : 0.0;
    curSig.Weights[i] = curVal[i];
  }
  
  if (Rf_asLogical(sFlows) == TRUE) {
      flows = malloc(sizeof(flow_t) * (baseRows + curRows - 1));
      if (!flows)
	  Rf_error("unable to allocate memory for flows");
  }

  double d = emd_rubner(&baseSig, &curSig, flows, flows ? &n_flows : NULL, Rf_asInteger(sExtra), dist_fn);

  if (!distName) /* cf1, cf2 */
      UNPROTECT(2);
  
  if (!flows)
      return Rf_ScalarReal(d);

  SEXP res = PROTECT(Rf_ScalarReal(d));
  SEXP fl = PROTECT(Rf_allocVector(VECSXP, 3)); /* must protect due to install() */
  Rf_setAttrib(res, Rf_install("flows"), fl);
  UNPROTECT(1);
  SEXP f_from = Rf_allocVector(INTSXP, n_flows);  SET_VECTOR_ELT(fl, 0, f_from);
  SEXP f_to   = Rf_allocVector(INTSXP, n_flows);  SET_VECTOR_ELT(fl, 1, f_to);
  SEXP f_amt  = Rf_allocVector(REALSXP, n_flows); SET_VECTOR_ELT(fl, 2, f_amt);
  int * i_from = INTEGER(f_from), * i_to = INTEGER(f_to);
  double * r_amt = REAL(f_amt);
  
  for (i = 0; i < n_flows; i++) {
      i_from[i] = flows[i].from;
      i_to[i]   = flows[i].to;
      r_amt[i]  = flows[i].amount;
  }
  free(flows);
  
  UNPROTECT(1);
  return res;
}
