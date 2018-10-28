/*
 * To register the functions implemented in C++, see
 * http://cran.r-project.org/doc/manuals/R-exts.html#Registering-native-routines
 *
 * But it seems not to work as it is supposed to be in that
 * they are still working if not registered.
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rversion.h>

static const R_CallMethodDef CallEntries[] = {
  {NULL, NULL, 0}
};

void attribute_visible R_init_rstanarm(DllInfo *dll) {
  // next line is necessary to avoid a NOTE from R CMD check
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE); // necessary for .onLoad() to work
}
