/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge MA, 02139 USA
   
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

#ifndef IGRAPH_HASHTABLE_H
#define IGRAPH_HASHTABLE_H

#include "igraph_vector.h"
#include "igraph_vector_ptr.h"

#include <stdio.h>

typedef unsigned int igraph_hash_function_t (const void *key);
typedef igraph_bool_t igraph_hash_eq_function_t (const void *key1, 
						 const void *key2);
typedef int igraph_hash_printkey_function_t (const void *key, FILE *file);

typedef struct {
  igraph_vector_int_t first;
  igraph_vector_int_t next;
  igraph_vector_int_t value;
  igraph_vector_uint_t hashvalue;
  igraph_vector_ptr_t key;
  int free;
  int tablelength;
  int entrycount;
  int loadlimit;
  int primeindex;
  igraph_hash_function_t *hashfn;
  igraph_hash_eq_function_t *eqfn;
} igraph_hashtable_t;

int igraph_hashtable_init(igraph_hashtable_t *hashtable,
			  int minsize, igraph_hash_function_t *hashfunction,
			  igraph_hash_eq_function_t *eqfunction);
void igraph_hashtable_destroy(igraph_hashtable_t *hashtable);
int igraph_hashtable_insert(igraph_hashtable_t *hashtable,
			    const void *key, int value);
int igraph_hashtable_search(const igraph_hashtable_t *hashtable, 
			    const void *key);
int igraph_hashtable_remove(igraph_hashtable_t *hashtable, const void *key);
int igraph_hashtable_count(const igraph_hashtable_t *hashtable);
int igraph_hashtable_update(igraph_hashtable_t *hashtable, const void *key, 
			    int newvalue);
int igraph_hashtable_fprint(const igraph_hashtable_t *hashtable, 
			    igraph_hash_printkey_function_t *printfunction, 
			    FILE *file);
int igraph_hashtable_print(const igraph_hashtable_t *hashtable, 
			   igraph_hash_printkey_function_t *printfunction);

unsigned int igraph_hashtable_string_hash(const void *key);
igraph_bool_t igraph_hashtable_string_eq(const void *key1, const void *key2);

int igraph_hashtable_string_init(igraph_hashtable_t *hashtable, int minsize);
int igraph_hashtable_string_fprint(const igraph_hashtable_t *hashtable,
				   FILE *file);
int igraph_hashtable_string_print(const igraph_hashtable_t *hashtable);

#endif

