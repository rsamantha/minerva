#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  The following name(s) appear with different usages
  e.g., with different numbers of arguments:

    mineRonevar

  This needs to be resolved in the tables and any declarations.
*/

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP mineRall(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mineRonevar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"mineRall",    (DL_FUNC) &mineRall,    7},
    {"mineRonevar", (DL_FUNC) &mineRonevar, 6},
    {NULL, NULL, 0}
};

void R_init_minerva(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

