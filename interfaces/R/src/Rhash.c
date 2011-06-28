/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2011  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph.h"
#include "igraph_error.h"
#include "rinterface.h"

#include "config.h"

#define USE_RINTERNALS
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

static const unsigned int igraph_Rhashtable_primes[] = {
  5, 11, 27, 53, 97, 193, 389,
  769, 1543, 3079, 6151,
  12289, 24593, 49157, 98317,
  196613, 393241, 786433, 1572869,
  3145739, 6291469, 12582917, 25165843,
  50331653, 100663319, 201326611, 402653189,
  805306457, 1610612741
};
const unsigned int igraph_Rhashtable_prime_table_length = 
  sizeof(igraph_Rhashtable_primes)/sizeof(igraph_Rhashtable_primes[0]);
const double igraph_Rhashtable_max_load_factor = 0.65;

#define FIRST(x)       (INTEGER(VECTOR_ELT((x), 0)))
#define NEXT(x)        (INTEGER(VECTOR_ELT((x), 1)))
#define VALUE(x)       (INTEGER(VECTOR_ELT((x), 2)))
#define HASHVALUE(x)   ((unsigned int*)(INTEGER(VECTOR_ELT((x), 3))))
#define FREE(x)        (INTEGER(VECTOR_ELT((x), 4))[0])
#define TABLELENGTH(x) (INTEGER(VECTOR_ELT((x), 4))[1])
#define PRIMEINDEX(x)  (INTEGER(VECTOR_ELT((x), 4))[2])
#define ENTRYCOUNT(x)  (INTEGER(VECTOR_ELT((x), 4))[3])
#define LOADLIMIT(x)   (INTEGER(VECTOR_ELT((x), 4))[4])
#define SIZE(x)        (GET_LENGTH(VECTOR_ELT((x), 1)))

SEXP igraph_Rhashtable_init(int minsize) {

  SEXP result, meta;
  int pindex, size=igraph_Rhashtable_primes[0];
  if (minsize > (1u << 30)) {
    error("igraph hash table too big");
  }
  for (pindex=0; pindex < igraph_Rhashtable_prime_table_length; pindex++) {
    if (igraph_Rhashtable_primes[pindex] > minsize) {
      size = igraph_Rhashtable_primes[pindex]; 
      break;
    }
  }
  
  PROTECT(result = NEW_LIST(5));
  SET_VECTOR_ELT(result, 0, NEW_INTEGER(size)); /* first */
  memset(INTEGER(VECTOR_ELT(result, 0)), 0, size * sizeof(int));
  SET_VECTOR_ELT(result, 1, NEW_INTEGER(0));    /* next */
  SET_VECTOR_ELT(result, 2, NEW_INTEGER(0));    /* value */
  SET_VECTOR_ELT(result, 3, NEW_INTEGER(0));    /* hashvalue */
  SET_VECTOR_ELT(result, 4, NEW_INTEGER(5));	/* meta */
  meta = VECTOR_ELT(result, 4);

  INTEGER(meta)[0] = 0;		/* free */
  INTEGER(meta)[1] = size; 	/* tablelength */
  INTEGER(meta)[2] = pindex;	/* primeindex */
  INTEGER(meta)[3] = 0;		/* entrycount */
  INTEGER(meta)[4] = 		/* loadlimit */
    (unsigned int) ceil(size * igraph_Rhashtable_max_load_factor);

  UNPROTECT(1);
  return result;
}

unsigned int igraph_Rhashtable_hash(const char *key) {
  /* Aim to protect against poor hash functions by adding logic here
   * - logic taken from java 1.4 hashtable source */
  unsigned int i = igraph_Rhashtable_string_hash(key);
  i += ~(i << 9);
  i ^=  ((i >> 14) | (i << 18)); /* >>> */
  i +=  (i << 4);
  i ^=  ((i >> 10) | (i << 22)); /* >>> */
  return i;
}

SEXP igraph_Rhashtable_expand(SEXP hashtable) {

  SEXP oldfirstsexp;
  int *oldfirst, *newfirst;
  int oldsize, newsize, e, mynext, i;
  int primeindex = PRIMEINDEX(hashtable);
  int tablelength = TABLELENGTH(hashtable);
  unsigned int *hashvalue = HASHVALUE(hashtable);
  int *next = NEXT(hashtable);

  PROTECT(oldfirstsexp = VECTOR_ELT(hashtable, 0));
  oldfirst=INTEGER(oldfirstsexp);

  if (primeindex == igraph_Rhashtable_prime_table_length - 1) {
    error("igraph hash table too big");
  }
  oldsize = igraph_Rhashtable_primes[ primeindex ];
  newsize = igraph_Rhashtable_primes[ primeindex + 1 ];
  SET_VECTOR_ELT(hashtable, 0, NEW_INTEGER(newsize));
  newfirst=FIRST(hashtable);
  memset(newfirst, 0, newsize * sizeof(int));
  PRIMEINDEX(hashtable) += 1;

  /* We might need to rehash some elements, as we are hashing into a
     larger interval now. */
  for (i=0; i<tablelength; i++) {
    for (e=oldfirst[i]; e != 0; e = mynext) {
      int index = hashvalue[e-1] % newsize;
      mynext = next[e-1];
      next[e-1] = newfirst[index];
      newfirst[index] = e;
    }
  }
  
  TABLELENGTH(hashtable) = newsize;
  LOADLIMIT(hashtable) =
    (unsigned int) ceil(newsize * igraph_Rhashtable_max_load_factor);

  UNPROTECT(1);
  return hashtable;
}

SEXP igraph_Rhashtable_resize(SEXP hashtable) {
  
  int *oldvec, *next;
  unsigned int *olduvec;
  int size=SIZE(hashtable);
  int newsize= size ? size * 2 : 1;  
  int i;
  
  oldvec = NEXT(hashtable);
  SET_VECTOR_ELT(hashtable, 1, NEW_INTEGER(newsize));
  memcpy(NEXT(hashtable), oldvec, size * sizeof(int));
  
  oldvec = VALUE(hashtable);
  SET_VECTOR_ELT(hashtable, 2, NEW_INTEGER(newsize));
  memcpy(VALUE(hashtable), oldvec, size * sizeof(int));
  
  olduvec = HASHVALUE(hashtable);
  SET_VECTOR_ELT(hashtable, 3, NEW_INTEGER(newsize));
  memcpy(HASHVALUE(hashtable), olduvec, size * sizeof(int));

  /* Need to update the list of free cells */
  next=NEXT(hashtable);
  for (i=size; i<newsize; i++) {
    next[i] = i+2;
  }
  next[newsize-1] = FREE(hashtable);
  FREE(hashtable) = size+1;

  return hashtable;
}

/* TODO: PROTECT() properly. */

SEXP igraph_Rhashtable_insert(SEXP hashtable,
			     const char *key, int value) {

  int pos, index;
  int *next;

  /* First vector full? */
  if (ENTRYCOUNT(hashtable) + 1 > LOADLIMIT(hashtable)) {
    hashtable = igraph_Rhashtable_expand(hashtable);
  }
  /* Others full? */
  if (!FREE(hashtable)) {
    hashtable = igraph_Rhashtable_resize(hashtable);
  }
  next = NEXT(hashtable);
  pos = FREE(hashtable) - 1;
  FREE(hashtable) = next[pos];
  
  ENTRYCOUNT(hashtable) += 1;
  HASHVALUE(hashtable)[pos] = igraph_Rhashtable_hash(key);
  index = HASHVALUE(hashtable)[pos] % TABLELENGTH(hashtable);
  VALUE(hashtable)[pos] = value;
  NEXT(hashtable)[pos] = FIRST(hashtable)[index];
  FIRST(hashtable)[index] = pos+1;
  
  return hashtable;
}

int igraph_Rhashtable_search(SEXP hashtable, const char *key, SEXP names) {
  int *value = VALUE(hashtable);
  unsigned int hashvalue = igraph_Rhashtable_hash(key);
  unsigned int index = hashvalue % TABLELENGTH(hashtable);
  unsigned int e = FIRST(hashtable)[index];
  while (e) {
    --e;
    if (hashvalue == HASHVALUE(hashtable)[e]) {
      const char *key2 = CHAR(STRING_ELT(names, value[e]));
      if (!strcmp(key, key2)) {
	return VALUE(hashtable)[e];
      }
    }
    e = NEXT(hashtable)[e];
  }
  
  return NA_INTEGER;		/* TODO: why doesn't this work? */
}

int igraph_Rhashtable_update(SEXP hashtable, const char *key, int newvalue, 
			     SEXP names) {
  int *value = VALUE(hashtable);
  unsigned int hashvalue = igraph_Rhashtable_hash(key);
  unsigned int index = hashvalue % TABLELENGTH(hashtable);
  unsigned int e = FIRST(hashtable)[index];
  while (e) {
    const char *key2 = CHAR(STRING_ELT(names, value[--e]));
    if (hashvalue == HASHVALUE(hashtable)[e] && !strcmp(key, key2)) {
      VALUE(hashtable)[e] = newvalue;
    }
    e = NEXT(hashtable)[e];
  }
  
  return 0;
}

SEXP igraph_Rhashtable_remove(SEXP hashtable, const char *key, SEXP names) {
  int *value = VALUE(hashtable);
  unsigned int hashvalue = igraph_Rhashtable_hash(key);
  unsigned int index = hashvalue % TABLELENGTH(hashtable);
  unsigned int e = FIRST(hashtable)[index];
  unsigned int prev = 0, next;
  
  while (e) {
    const char *key2 = CHAR(STRING_ELT(names, value[--e]));
    next = NEXT(hashtable)[e];
    if (hashvalue == HASHVALUE(hashtable)[e] && !strcmp(key, key2)) {
      ENTRYCOUNT(hashtable) -= 1;
      
      if (prev) {
	NEXT(hashtable)[prev] = next;
      } else {
	FIRST(hashtable)[index] = next;
      }
      
      NEXT(hashtable)[e] = FREE(hashtable);
      FREE(hashtable) = e+1;

      return hashtable;
    }
    prev = e;
    e = next;
  }
  
  return hashtable;
}

int igraph_Rhashtable_count(const SEXP hashtable) {
  return ENTRYCOUNT(hashtable);
}

/* --------------- */
/* Strings as keys */ 
/* --------------- */

#define PER_LOOP (8)

// macro HASH_SDBM is the equivalent of
//   hash(i) = hash(i - 1) * 65599 + str[i];
#define HASH_INIT (0)
#define HASH_SDBM(h) (h << 6) + (h << 16) - h

// macro HASH_DJB2M is the equivalent of
//   hash(i) = hash(i - 1) * 33 + str[i];
//#define HASH_INIT (5381)
//#define HASH_DJB2(h) ((h << 5) + h)

#define HASH_STEP(h, s) h = *s++ + HASH_SDBM(h);

// loop unwound version using a duff device
unsigned int igraph_Rhashtable_string_hash(const char *key) {
  register char *s = (char*) key;
  unsigned int l=strlen(s);
  register unsigned int h = HASH_INIT;
  register unsigned int n = (l + PER_LOOP - 1) / PER_LOOP;

  switch (l % PER_LOOP) {
  case 0: do { 
      HASH_STEP(h, s)
    case 7:      HASH_STEP(h, s)
    case 6:      HASH_STEP(h, s)
    case 5:      HASH_STEP(h, s)
    case 4:      HASH_STEP(h, s)
    case 3:      HASH_STEP(h, s)
    case 2:      HASH_STEP(h, s)
    case 1:      HASH_STEP(h, s)
	} while (--n);
  }

  return h;
}

/* TODO: PROTECT() properly. */

SEXP R_igraph_hash_create(igraph_attribute_record_t *attrrec, int first_id) {
  igraph_strvector_t *strvec = (igraph_strvector_t *) attrrec->value;
  int i, n=igraph_strvector_size(strvec);
  SEXP hash;

  PROTECT(hash = igraph_Rhashtable_init(n));
  for (i=0; i<n; i++, first_id++) {
    hash = igraph_Rhashtable_insert(hash, STR(*strvec, i), first_id);
  }
  
  UNPROTECT(1);
  return hash;
}

/* TODO: PROTECT() properly. */

SEXP R_igraph_hash_update(SEXP hash, SEXP index, SEXP oldval, SEXP newval) {
  int i, n=GET_LENGTH(index);
  int *rindex=INTEGER(index);

  for (i=0; i<n; i++) {
    hash = igraph_Rhashtable_remove(hash, CHAR(STRING_ELT(oldval, i)),
				    oldval);
  } 
  for (i=0; i<n; i++) {
    hash = igraph_Rhashtable_insert(hash, CHAR(STRING_ELT(newval, i)), 
				    rindex[i]-1);
  }
  
  return hash;
}

/* TODO: PROTECT() properly. */

SEXP R_igraph_hash_create2(SEXP names) {
  int i, n=GET_LENGTH(names);
  SEXP hash;
  
  PROTECT(hash = igraph_Rhashtable_init(n));
  for (i=0; i<n; i++) {
    hash = igraph_Rhashtable_insert(hash, CHAR(STRING_ELT(names, i)), i);
  }
  
  UNPROTECT(1);
  return hash;
}

SEXP R_igraph_hash_match(SEXP graph, SEXP query) {
  SEXP hashtable=VECTOR_ELT(VECTOR_ELT(graph, 8), 4);
  SEXP names=R_igraph_getListElement(VECTOR_ELT(VECTOR_ELT(graph, 8), 2), 
				     "name");
  SEXP result;
  int i, len=GET_LENGTH(query);
  int *iresult;
  int *value = VALUE(hashtable);
  unsigned int tablelength = TABLELENGTH(hashtable);
  int *first = FIRST(hashtable);
  unsigned int *hashvalue = HASHVALUE(hashtable);
  int *next = NEXT(hashtable);
  
  PROTECT(result = NEW_INTEGER(len));
  iresult=INTEGER(result);
  for (i=0; i<len; i++) {
    const char *key = CHAR(STRING_ELT(query, i));
    unsigned int myhashvalue = igraph_Rhashtable_hash(key);
    unsigned int index = myhashvalue % tablelength;
    unsigned int e = first[index];
    iresult[i] = NA_INTEGER;	/* TODO */
    while (e) {
      --e; 
      if (myhashvalue == hashvalue[e]) { 
	const char *key2 = CHAR(STRING_ELT(names, value[e]));
	if (!strcmp(key, key2)) { iresult[i] = value[e]+1; break; }
      }
      e = next[e];
    }
  }
  
  UNPROTECT(1);
  return result;
}

/* The new vertex names were already added to 'names' */
/* TODO: PROTECT() properly. */

SEXP R_igraph_hash_add(SEXP hash, SEXP names, SEXP pfirst_id) {
  int first_id=INTEGER(pfirst_id)[0];
  int i, len=GET_LENGTH(names) - first_id;
  
  for (i=0; i<len; i++, first_id++) {
    hash = igraph_Rhashtable_insert(hash, CHAR(STRING_ELT(names, first_id)), 
				    first_id);
  }
  
  return hash;
}

/* TODO: PROTECT() properly. */

SEXP R_igraph_hash_permute(SEXP hash, SEXP oldnames, SEXP newnames, 
			   SEXP idx) {

  /* We rather craete a new hash, TODO: optimize this. */

  SEXP result;
  int i, n=GET_LENGTH(newnames);
  
  PROTECT(result = igraph_Rhashtable_init(n));
  for (i=0; i<n; i++) {
    result = igraph_Rhashtable_insert(result, CHAR(STRING_ELT(newnames, i)), 
				      i);
  }
  
  UNPROTECT(1);
  return result;
}
