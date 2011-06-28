/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
   Rue de l'Industrie 5, Lausanne 1005, Switzerland
   
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

#ifndef R_IGRAPH_RINTERFACE_H
#define R_IGRAPH_RINTERFACE_H

#include "igraph.h"
#include "igraph_error.h"

#include "config.h"

#define USE_RINTERNALS
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <stdio.h>

/* R global tables */

extern igraph_attribute_table_t R_igraph_attribute_table;
extern igraph_attribute_table_t *R_igraph_attribute_oldtable;

/* Conversions */

SEXP R_igraph_vector_to_SEXP(const igraph_vector_t *v);
SEXP R_igraph_vector_to_SEXPp1(const igraph_vector_t *v);
SEXP R_igraph_vector_bool_to_SEXP(igraph_vector_bool_t *v);
SEXP R_igraph_vector_long_to_SEXP(igraph_vector_long_t *v);
SEXP R_igraph_0orvector_to_SEXP(igraph_vector_t *v);
SEXP R_igraph_0orvector_bool_to_SEXP(igraph_vector_bool_t *v);
SEXP R_igraph_0orvector_long_to_SEXP(igraph_vector_long_t *v);
SEXP R_igraph_0orvector_to_SEXPp1(igraph_vector_t *v);
SEXP R_igraph_matrix_to_SEXP(igraph_matrix_t *m);
SEXP R_igraph_0ormatrix_to_SEXP(igraph_matrix_t *m);
SEXP R_igraph_array3_to_SEXP(igraph_array3_t *a);
SEXP R_igraph_0orarray3_to_SEXP(igraph_array3_t *a);
SEXP R_igraph_strvector_to_SEXP(const igraph_strvector_t *m);
SEXP R_igraph_to_SEXP(igraph_t *graph);
SEXP R_igraph_vectorlist_to_SEXP(const igraph_vector_ptr_t *ptr);
SEXP R_igraph_vectorlist_to_SEXP_p1(const igraph_vector_ptr_t *ptr);
void R_igraph_vectorlist_destroy(igraph_vector_ptr_t *ptr);
SEXP R_igraph_hrg_to_SEXP(const igraph_hrg_t *hrg);
SEXP R_igraph_arpack_options_to_SEXP(igraph_arpack_options_t *opt);
SEXP R_igraph_bliss_info_to_SEXP(igraph_bliss_info_t *info);

int R_igraph_SEXP_to_strvector(SEXP rval, igraph_strvector_t *sv);
int R_igraph_SEXP_to_strvector_copy(SEXP rval, igraph_strvector_t *sv);
int R_SEXP_to_vector(SEXP sv, igraph_vector_t *v);
int R_SEXP_to_vector_copy(SEXP sv, igraph_vector_t *v);
int R_SEXP_to_matrix(SEXP pakl, igraph_matrix_t *akl);
int R_igraph_SEXP_to_array3(SEXP rval, igraph_array3_t *a);
int R_igraph_SEXP_to_array3_copy(SEXP rval, igraph_array3_t *a);
int R_SEXP_to_igraph_matrix_copy(SEXP pakl, igraph_matrix_t *akl);
int R_SEXP_to_igraph(SEXP graph, igraph_t *res);
int R_SEXP_to_igraph_copy(SEXP graph, igraph_t *res);
int R_SEXP_to_igraph_vs(SEXP rit, igraph_t *graph, igraph_vs_t *it);
int R_SEXP_to_igraph_es(SEXP rit, igraph_t *graph, igraph_es_t *it);
int R_SEXP_to_igraph_adjlist(SEXP vectorlist, igraph_adjlist_t *ptr);
int R_SEXP_to_vector_bool(SEXP sv, igraph_vector_bool_t *v);
int R_SEXP_to_vector_int(SEXP sv, igraph_vector_int_t *v);
int R_SEXP_to_hrg(SEXP shrg, igraph_hrg_t *hrg);
int R_SEXP_to_hrg_copy(SEXP shrg, igraph_hrg_t *hrg);
int R_SEXP_to_igraph_arpack_options(SEXP in, igraph_arpack_options_t *opt);
int R_SEXP_to_igraph_layout_drl_options(SEXP in, 
					igraph_layout_drl_options_t *opt);

/* Utility functions */

SEXP R_igraph_i_list5(SEXP i1, SEXP i2, SEXP i3, SEXP i4, SEXP i5);
SEXP R_igraph_i_lang7(SEXP s, SEXP t, SEXP u, SEXP v, SEXP w, SEXP x, SEXP y);
SEXP R_igraph_getListElement(SEXP list, const char *str);

/* Hashing vertex names */

#define R_IGRAPH_NAME "name"

SEXP igraph_Rhashtable_init(int minsize);
SEXP igraph_Rhashtable_insert(SEXP hashtable, const char *key, int value);
int igraph_Rhashtable_search(const SEXP hashtable, const char *key, 
			     SEXP names);
SEXP igraph_Rhashtable_remove(SEXP hashtable, const char *key, SEXP names);
int igraph_Rhashtable_count(const SEXP hashtable);
int igraph_Rhashtable_update(SEXP hashtable, const char *key, int newvalue,
			     SEXP names);

unsigned int igraph_Rhashtable_string_hash(const char *key);

SEXP R_igraph_hash_create(igraph_attribute_record_t *attrrec, int first_id);
SEXP R_igraph_hash_create2(SEXP names);
SEXP R_igraph_hash_add(SEXP hash, SEXP names, SEXP pfirst_id);
SEXP R_igraph_hash_match(SEXP graph, SEXP query);
SEXP R_igraph_hash_permute(SEXP hash, SEXP oldnames, SEXP newnames, SEXP idx);
SEXP R_igraph_hash_update(SEXP hash, SEXP index, SEXP oldval, SEXP newval);

/* Attributes */

int R_SEXP_to_attr_comb(SEXP input, igraph_attribute_combination_t *comb);
int R_igraph_attribute_init(igraph_t *graph, igraph_vector_ptr_t *attr);
void R_igraph_attribute_destroy(igraph_t *graph);
int R_igraph_attribute_copy(igraph_t *to, const igraph_t *from,
			    igraph_bool_t ga, igraph_bool_t va, 
			    igraph_bool_t ea);
int R_igraph_attribute_add_vertices(igraph_t *graph, long int nv, 
				    igraph_vector_ptr_t *nattr);
int R_igraph_attribute_permute_vertices(const igraph_t *graph,
					igraph_t *newgraph,
					const igraph_vector_t *idx);
int R_igraph_attribute_add_edges(igraph_t *graph, 
				 const igraph_vector_t *edges,
				 igraph_vector_ptr_t *nattr);
int R_igraph_attribute_permute_edges(const igraph_t *graph,
				     igraph_t *newgraph,
				     const igraph_vector_t *idx);
igraph_bool_t R_igraph_attribute_has_attr(const igraph_t *graph,
					  igraph_attribute_elemtype_t type,
					  const char *name);
int R_igraph_attribute_gettype(const igraph_t *graph,
			       igraph_attribute_type_t *type,
			       igraph_attribute_elemtype_t elemtype,
			       const char *name);
int R_igraph_attribute_get_numeric_graph_attr(const igraph_t *graph,
					      const char *name, 
					      igraph_vector_t *value);
int R_igraph_attribute_get_string_graph_attr(const igraph_t *graph,
					     const char *name,
					     igraph_strvector_t *value);
int R_igraph_attribute_get_numeric_vertex_attr(const igraph_t *graph, 
					       const char *name,
					       igraph_vs_t vs,
					       igraph_vector_t *value);
int R_igraph_attribute_get_string_vertex_attr(const igraph_t *graph, 
					      const char *name,
					      igraph_vs_t vs,
					      igraph_strvector_t *value);
int R_igraph_attribute_get_numeric_edge_attr(const igraph_t *graph,
					     const char *name,
					     igraph_es_t es,
					     igraph_vector_t *value);
int R_igraph_attribute_get_string_edge_attr(const igraph_t *graph,
					    const char *name,
					    igraph_es_t es,
					    igraph_strvector_t *value);

#endif
