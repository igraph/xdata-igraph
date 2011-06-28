/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2011  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "rinterface.h"

igraph_error_handler_t *R_igraph_oldhandler;

void R_igraph_myhandler (const char *reason, const char *file,
			 int line, int igraph_errno) {
  IGRAPH_FINALLY_FREE();
  error("At %s:%i : %s, %s", file, line, reason, 
	igraph_strerror(igraph_errno));
}

void R_igraph_warning_handler(const char *reason, const char *file,
			      int line, int igraph_errno) {
  warning("At %s:%i :%s", file, line, reason);
}

extern int R_interrupts_pending;

int R_igraph_interrupt_handler(void *data) {
#if  ( defined(HAVE_AQUA) || defined(Win32) )
  R_CheckUserInterrupt();
#else
  if (R_interrupts_pending) {
    IGRAPH_FINALLY_FREE();
    R_CheckUserInterrupt();
  }
#endif
  return 0;
}

int R_igraph_progress_handler(const char *message, igraph_real_t percent,
			      void * data) {
  SEXP ec;
  int ecint;
  PROTECT(ec=EVAL(lang3(install(".igraph.progress"), ScalarReal(percent), 
			ScalarString(mkChar(message)))));
  ecint=INTEGER(ec)[0];
  UNPROTECT(1);
  return ecint;
}

int R_igraph_status_handler(const char *message, void *data) {
  SEXP ec;
  int ecint;
  PROTECT(ec=EVAL(lang2(install(".igraph.status"), 
			ScalarString(mkChar(message)))));
  ecint=INTEGER(ec)[0];
  UNPROTECT(1);
  return 0;
}

SEXP R_igraph_init(SEXP progress, SEXP status) {
  igraph_set_error_handler(R_igraph_myhandler);
  igraph_set_warning_handler(R_igraph_warning_handler);
  igraph_set_interruption_handler(R_igraph_interrupt_handler);
  igraph_i_set_attribute_table(&R_igraph_attribute_table);
  if (LOGICAL(status)[0]) { 
    igraph_set_status_handler(R_igraph_status_handler);
  }
  if (LOGICAL(progress)[0]) {
    igraph_set_progress_handler(R_igraph_progress_handler);
  }
  return R_NilValue;
}

SEXP R_igraph_set_verbose(SEXP verbose) {
  if (LOGICAL(verbose)[0]) {
    igraph_set_status_handler(R_igraph_status_handler);    
    igraph_set_progress_handler(R_igraph_progress_handler);
  } else {
    igraph_set_status_handler(0);
    igraph_set_progress_handler(0);
  }
  return R_NilValue;
}

SEXP R_igraph_finalizer() {
  EVAL(lang4(install(".igraph.progress"), ScalarReal(0.0), 
	     ScalarString(mkChar("")), ScalarLogical(1)));
  IGRAPH_FINALLY_FREE();
  return R_NilValue;
}
