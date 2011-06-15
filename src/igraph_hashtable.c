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

#include "igraph_hashtable.h"

#include <string.h>

static const unsigned int igraph_hashtable_primes[] = {
  5, 11, 27, 53, 97, 193, 389,
  769, 1543, 3079, 6151,
  12289, 24593, 49157, 98317,
  196613, 393241, 786433, 1572869,
  3145739, 6291469, 12582917, 25165843,
  50331653, 100663319, 201326611, 402653189,
  805306457, 1610612741
};
const unsigned int igraph_hashtable_prime_table_length = 
  sizeof(igraph_hashtable_primes)/sizeof(igraph_hashtable_primes[0]);
const double igraph_hashtable_max_load_factor = 0.65;

int igraph_hashtable_init(igraph_hashtable_t *hashtable,
			  int minsize, igraph_hash_function_t *hashfunction,
			  igraph_hash_eq_function_t *eqfunction) {

  int pindex, size=igraph_hashtable_primes[0];
  if (minsize > (1u << 30)) {
    IGRAPH_ERROR("Hash table too big", IGRAPH_EINVAL);
  }
  for (pindex=0; pindex < igraph_hashtable_prime_table_length; pindex++) {
    if (igraph_hashtable_primes[pindex] > minsize) {
      size = igraph_hashtable_primes[pindex]; 
      break;
    }
  }
  
  IGRAPH_CHECK(igraph_vector_int_init(&hashtable->first, size));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &hashtable->first);
  IGRAPH_CHECK(igraph_vector_int_init(&hashtable->next, 0));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &hashtable->next);
  IGRAPH_CHECK(igraph_vector_int_init(&hashtable->value, 0));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &hashtable->value);
  IGRAPH_CHECK(igraph_vector_uint_init(&hashtable->hashvalue, 0));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &hashtable->hashvalue);
  IGRAPH_CHECK(igraph_vector_ptr_init(&hashtable->key, 0));
  IGRAPH_FINALLY(igraph_vector_ptr_destroy, &hashtable->key);

  hashtable->free = 0;
  hashtable->tablelength = size;
  hashtable->primeindex = pindex;
  hashtable->entrycount = 0;
  hashtable->loadlimit = 
    (unsigned int) ceil(size * igraph_hashtable_max_load_factor);
  hashtable->hashfn = hashfunction;
  hashtable->eqfn = eqfunction;

  IGRAPH_FINALLY_CLEAN(5);
  return 0;
}

void igraph_hashtable_destroy(igraph_hashtable_t *hashtable) {
  igraph_vector_int_destroy(&hashtable->first);
  igraph_vector_int_destroy(&hashtable->next);
  igraph_vector_int_destroy(&hashtable->value);
  igraph_vector_uint_destroy(&hashtable->hashvalue);  
}

unsigned int igraph_hashtable_hash(const igraph_hashtable_t *hashtable, 
				   const void *key) {
  /* Aim to protect against poor hash functions by adding logic here
   * - logic taken from java 1.4 hashtable source */
  unsigned int i = hashtable->hashfn(key);
  i += ~(i << 9);
  i ^=  ((i >> 14) | (i << 18)); /* >>> */
  i +=  (i << 4);
  i ^=  ((i >> 10) | (i << 22)); /* >>> */
  return i;
}

int igraph_hashtable_expand(igraph_hashtable_t *hashtable) {
  int newsize, e, next, i;

  if (hashtable->primeindex == igraph_hashtable_prime_table_length - 1) {
    IGRAPH_ERROR("Hash table too big", IGRAPH_ENOMEM);
  }
  newsize = igraph_hashtable_primes[ hashtable->primeindex + 1 ];
  IGRAPH_CHECK(igraph_vector_int_resize(&hashtable->first, newsize));
  hashtable->primeindex += 1;
  for (i=hashtable->tablelength; i < newsize; i++) {
    VECTOR(hashtable->first)[i] = 0;
  }

  /* We might need to rehash some elements, as we are hashing into a
     larger interval now. */
  for (i=0; i<hashtable->tablelength; i++) {
    for (e=VECTOR(hashtable->first)[i]; e != 0; e = next) {
      int index = VECTOR(hashtable->hashvalue)[e-1] % newsize;
      next = VECTOR(hashtable->next)[e-1];
      if (index != i) { 
	VECTOR(hashtable->next)[e-1] = VECTOR(hashtable->first)[index];
	VECTOR(hashtable->first)[index] = e;
      }
    }
  }
  
  hashtable->tablelength = newsize;
  hashtable->loadlimit = 
    (unsigned int) ceil(newsize * igraph_hashtable_max_load_factor);

  return 0;
}

/* TODO: proper error handling */

int igraph_hashtable_resize(igraph_hashtable_t *hashtable) {
  
  int size=igraph_vector_int_size(&hashtable->next);
  int newsize= size ? size * 2 : 1;  
  int i;
  
  IGRAPH_CHECK(igraph_vector_int_resize(&hashtable->next, newsize));
  IGRAPH_CHECK(igraph_vector_int_resize(&hashtable->value, newsize));
  IGRAPH_CHECK(igraph_vector_uint_resize(&hashtable->hashvalue, newsize));
  IGRAPH_CHECK(igraph_vector_ptr_resize(&hashtable->key, newsize));

  /* Need to update the list of free cells */
  for (i=size; i<newsize; i++) {
    VECTOR(hashtable->next)[i] = i+2;
  }
  VECTOR(hashtable->next)[newsize-1] = hashtable->free;
  hashtable->free = size+1;

  return 0;
}

int igraph_hashtable_insert(igraph_hashtable_t *hashtable,
			    const void *key, int value) {

  int pos, index;

  /* First vector full? */
  if (hashtable->entrycount + 1 > hashtable->loadlimit) {
    IGRAPH_CHECK(igraph_hashtable_expand(hashtable));
  }
  /* Others full? */
  if (!hashtable->free) {
    IGRAPH_CHECK(igraph_hashtable_resize(hashtable));
  }
  pos = hashtable->free - 1;
  hashtable->free = VECTOR(hashtable->next)[pos];
  
  hashtable->entrycount += 1;
  VECTOR(hashtable->hashvalue)[pos] = igraph_hashtable_hash(hashtable, key);
  index = VECTOR(hashtable->hashvalue)[pos] % hashtable->tablelength;
  VECTOR(hashtable->key)[pos] = (void*) key; /* !!!!! */
  VECTOR(hashtable->value)[pos] = value;
  VECTOR(hashtable->next)[pos] = VECTOR(hashtable->first)[index];
  VECTOR(hashtable->first)[index] = pos+1;
  
  return 0;
}

int igraph_hashtable_search(const igraph_hashtable_t *hashtable,
			    const void *key) {
  unsigned int hashvalue = igraph_hashtable_hash(hashtable, key);
  unsigned int index = hashvalue % hashtable->tablelength;
  unsigned int e = VECTOR(hashtable->first)[index];
  while (e) {
    e--;
    if (hashvalue == VECTOR(hashtable->hashvalue)[e] && 
	hashtable->eqfn(key, VECTOR(hashtable->key)[e])) {
      return VECTOR(hashtable->value)[e];
    }
    e = VECTOR(hashtable->next)[e];
  }
  
  return 0;
}

int igraph_hashtable_update(igraph_hashtable_t *hashtable, const void *key, 
			    int newvalue) {
  unsigned int hashvalue = igraph_hashtable_hash(hashtable, key);
  unsigned int index = hashvalue % hashtable->tablelength;
  unsigned int e = VECTOR(hashtable->first)[index];
  while (e) {
    e--;
    if (hashvalue == VECTOR(hashtable->hashvalue)[e] && 
	hashtable->eqfn(key, VECTOR(hashtable->key)[e])) {
      VECTOR(hashtable->value)[e] = newvalue;
    }
    e = VECTOR(hashtable->next)[e];
  }
  
  return 0;
}

int igraph_hashtable_remove(igraph_hashtable_t *hashtable, const void *key) {

  unsigned int hashvalue = igraph_hashtable_hash(hashtable, key);
  unsigned int index = hashvalue % hashtable->tablelength;
  unsigned int e = VECTOR(hashtable->first)[index];
  unsigned int prev = 0, next;
  
  while (e) {
    e--;
    next = VECTOR(hashtable->next)[e];
    if (hashvalue == VECTOR(hashtable->hashvalue)[e] &&
	hashtable->eqfn(key, VECTOR(hashtable->key)[e])) {
      hashtable->entrycount -= 1;
      
      if (prev) {
	VECTOR(hashtable->next)[prev] = next;
      } else {
	VECTOR(hashtable->first)[index] = 0;
      }

      return VECTOR(hashtable->value)[e];
    }
    prev = e;
    e = next;
  }
  
  return 0;
}

int igraph_hashtable_count(const igraph_hashtable_t *hashtable) {
  return hashtable->entrycount;
}

int igraph_hashtable_fprint(const igraph_hashtable_t *hashtable, 
			   igraph_hash_printkey_function_t *printfunction, 
			   FILE *file) {
  int i, e, ret;

  /*
  fprintf(file, "Free: %i, TableLength: %i, EntryCount: %i, Loadlimit: %i, "
	  "PrimeIndex: %i\n", hashtable->free, hashtable->tablelength,
	  hashtable->entrycount, hashtable->loadlimit, hashtable->primeindex);

  fprintf(file, "FIRST: "); igraph_vector_int_fprint(&hashtable->first, file);
  fprintf(file, "NEXT: ");  igraph_vector_int_fprint(&hashtable->next, file);
  fprintf(file, "VALUE: "); igraph_vector_int_fprint(&hashtable->value, file);
  fprintf(file, "HASH: ");  igraph_vector_uint_fprint(&hashtable->hashvalue,
						      file);
  */
  
  for (i=0; i<hashtable->tablelength; i++) {
    for (e=VECTOR(hashtable->first)[i]; e != 0; 
	 e = VECTOR(hashtable->next)[e-1]) {
      IGRAPH_CHECK(printfunction(VECTOR(hashtable->key)[e-1], file));
      ret=fprintf(file, ": %i\n", VECTOR(hashtable->value)[e-1]);
      if (ret < 0) { 
	IGRAPH_ERROR("Cannot print hashtable", IGRAPH_EFILE);
      }
    }
  }
  return 0;
}

int igraph_hashtable_print(const igraph_hashtable_t *hashtable, 
			   igraph_hash_printkey_function_t *printfunction) {
  return igraph_hashtable_fprint(hashtable, printfunction, stdout);
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
unsigned int igraph_hashtable_string_hash(const void *key) {
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

igraph_bool_t igraph_hashtable_string_eq(const void *key1, const void *key2) {
  return ! strcmp( (const char*) key1, (const char*) key2);
}

int igraph_hashtable_string_init(igraph_hashtable_t *hashtable, int minsize) {
  return igraph_hashtable_init(hashtable, minsize, 
			       igraph_hashtable_string_hash,
			       igraph_hashtable_string_eq);
}

int igraph_hashtable_string_printkey(const void *key, FILE *file) {
  const char *str = (const char *) key;
  int ret = fprintf(file, "%s", str);
  if (ret < 0) { IGRAPH_ERROR("Cannot print hash table", IGRAPH_EFILE); }
  return 0;
}

int igraph_hashtable_string_fprint(const igraph_hashtable_t *hashtable,
				   FILE *file) {
  return igraph_hashtable_fprint(hashtable, 
				 igraph_hashtable_string_printkey, file);
}

int igraph_hashtable_string_print(const igraph_hashtable_t *hashtable) {
  return igraph_hashtable_fprint(hashtable, 
				 igraph_hashtable_string_printkey, stdout);
}
