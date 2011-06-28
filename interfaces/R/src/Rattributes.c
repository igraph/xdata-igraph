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

/******************************************************
 * Attributes                                         *
 *****************************************************/

int R_SEXP_to_attr_comb(SEXP input, igraph_attribute_combination_t *comb) {
  long int i, n=GET_LENGTH(input);
  SEXP names=GET_NAMES(input);
  
  igraph_attribute_combination_init(comb);
  
  for (i=0; i<n; i++) {
    const char *n;
    igraph_attribute_combination_type_t type;
    void *func;
    
    /* Name */
    if (!isNull(names)) {
      n=CHAR(STRING_ELT(names, i));
    }
    if (isNull(names) || strlen(n)==0) {
      n=0;
    }

    /* Type and function, if any */
    if (isFunction(VECTOR_ELT(input, i))) {
      type=IGRAPH_ATTRIBUTE_COMBINE_FUNCTION;
      func=VECTOR_ELT(input, i);      
    } else {
      type=REAL(AS_NUMERIC(VECTOR_ELT(input, i)))[0];
      func=0;
    }
    igraph_attribute_combination_add(comb, n, type, func);
  }
  
  return 0;
}

int R_igraph_attribute_init(igraph_t *graph, igraph_vector_ptr_t *attr) {
  SEXP result, names, gal;
  long int i;
  long int attrno;
  PROTECT(result=NEW_LIST(5));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(3));
  REAL(VECTOR_ELT(result, 0))[0]=0; /* R objects */
  REAL(VECTOR_ELT(result, 0))[1]=1; /* igraph_t objects */
  REAL(VECTOR_ELT(result, 0))[2]=1; /* whether the graph is safe */
  for (i=2; i<=3; i++) {
    SET_VECTOR_ELT(result, i, NEW_LIST(0)); /* val, eal */
  }
  SET_VECTOR_ELT(result, i, R_NilValue); /* hash */
  graph->attr=result;

  /* Add graph attributes */
  attrno= attr==NULL ? 0 : igraph_vector_ptr_size(attr);
  SET_VECTOR_ELT(result, 1, NEW_LIST(attrno));
  gal=VECTOR_ELT(result, 1);
  PROTECT(names=NEW_CHARACTER(attrno));
  for (i=0; i<attrno; i++) {
    igraph_attribute_record_t *rec=VECTOR(*attr)[i];
    igraph_vector_t *vec;
    igraph_strvector_t *strvec;
    SET_STRING_ELT(names, i, mkChar(rec->name));
    SET_VECTOR_ELT(gal, i, R_NilValue);
    switch (rec->type) {
    case IGRAPH_ATTRIBUTE_NUMERIC:
      vec=(igraph_vector_t*) rec->value;
      if (igraph_vector_size(vec) > 0) {
	SET_VECTOR_ELT(gal, i, NEW_NUMERIC(1));
	REAL(VECTOR_ELT(gal, i))[0]=VECTOR(*vec)[0];
      }
      break;
    case IGRAPH_ATTRIBUTE_STRING:
      strvec=(igraph_strvector_t*) rec->value;
      if (igraph_strvector_size(strvec) > 0) {
	SET_VECTOR_ELT(gal, i, NEW_CHARACTER(1));
	SET_STRING_ELT(VECTOR_ELT(gal,i), 0, mkChar(STR(*strvec, 0)));
      }
      break;
    case IGRAPH_ATTRIBUTE_R_OBJECT:
      IGRAPH_ERROR("R_objects not implemented yet", IGRAPH_UNIMPLEMENTED);
      break;
    }
  }
  SET_NAMES(gal, names);
  
  UNPROTECT(1);
  return 0;
}

void R_igraph_attribute_destroy(igraph_t *graph) {
  SEXP attr=graph->attr;
  REAL(VECTOR_ELT(attr, 0))[1] -= 1; /* refcount for igraph_t */
  if (REAL(VECTOR_ELT(attr, 0))[1]==0 && 
      REAL(VECTOR_ELT(attr, 0))[2]==1) {
    UNPROTECT_PTR(attr);
  }
  graph->attr=0;
}

/* If not copying all three attribute kinds are requested, then 
   we don't refcount, but really copy the requested ones, because 
   1) we can only refcount all three at the same time, and 
   2) the not-copied attributes will be set up by subsequent calls 
      to permute_vertices and/or permute/edges anyway. */

int R_igraph_attribute_copy(igraph_t *to, const igraph_t *from,
			    igraph_bool_t ga, igraph_bool_t va, igraph_bool_t ea) {
  SEXP fromattr=from->attr;
  if (ga && va && ea) {
    to->attr=from->attr;
    REAL(VECTOR_ELT(fromattr, 0))[1] += 1; /* refcount only */
    if (REAL(VECTOR_ELT(fromattr, 0))[1] == 1) {
      PROTECT(to->attr);
    }
  } else {
    R_igraph_attribute_init(to,0); /* Sets up many things */
    SEXP toattr=to->attr;
    if (ga) {
      SET_VECTOR_ELT(toattr, 1, duplicate(VECTOR_ELT(fromattr, 1)));
    } 
    if (va) {
      SET_VECTOR_ELT(toattr, 2, duplicate(VECTOR_ELT(fromattr, 2)));      
      SET_VECTOR_ELT(toattr, 4, duplicate(VECTOR_ELT(fromattr, 4)));
    } 
    if (ea) {
      SET_VECTOR_ELT(toattr, 3, duplicate(VECTOR_ELT(fromattr, 3)));      
    }
  }
  return 0;
}

int R_igraph_attribute_add_vertices(igraph_t *graph, long int nv, 
				    igraph_vector_ptr_t *nattr) {
  SEXP attr=graph->attr;
  SEXP val, rep=0, names, newnames;
  igraph_vector_t news;
  long int valno, i, origlen, nattrno, newattrs;
  igraph_bool_t had_names=0, added_names=0;

  if (REAL(VECTOR_ELT(attr, 0))[0]+REAL(VECTOR_ELT(attr, 0))[1] > 1) {
    SEXP newattr;
    PROTECT(newattr=duplicate(attr));
    REAL(VECTOR_ELT(attr, 0))[1] -= 1;
    if (REAL(VECTOR_ELT(attr, 0))[1] == 0) {
      UNPROTECT_PTR(attr);
    }
    REAL(VECTOR_ELT(newattr, 0))[0] = 0;
    REAL(VECTOR_ELT(newattr, 0))[1] = 1;
    attr=graph->attr=newattr;
  }

  val=VECTOR_ELT(attr, 2);
  valno=GET_LENGTH(val);  
  names=GET_NAMES(val);
  if (nattr==NULL) { 
    nattrno=0;
  } else {
    nattrno=igraph_vector_ptr_size(nattr); 
  }
  origlen=igraph_vcount(graph)-nv;

  /* First add the new attributes, if any */
  newattrs=0;
  IGRAPH_VECTOR_INIT_FINALLY(&news, 0);
  for (i=0; i<valno; i++) {
    if (!strcmp(CHAR(STRING_ELT(names, i)), R_IGRAPH_NAME)) {
      had_names = i+1;
    }
  }
  for (i=0; i<nattrno; i++) {
    igraph_attribute_record_t *nattr_entry=VECTOR(*nattr)[i];
    const char *nname=nattr_entry->name;
    long int j; 
    igraph_bool_t l=0;
    if (!strcmp(nname, R_IGRAPH_NAME)) {
      added_names=i+1;
    }
    for (j=0; !l && j<valno; j++) {
      l=!strcmp(nname, CHAR(STRING_ELT(names, j)));
    }
    if (!l) {
      newattrs++;
      IGRAPH_CHECK(igraph_vector_push_back(&news, i));
    }
  }
  if (newattrs != 0) {
    SEXP app, newval;
    PROTECT(app=NEW_LIST(newattrs));
    PROTECT(newnames=NEW_CHARACTER(newattrs));
    PROTECT(rep=EVAL(lang3(install("rep"), ScalarLogical(NA_LOGICAL), 
			   ScalarInteger(origlen))));
    for (i=0; i<newattrs; i++) {
      igraph_attribute_record_t *tmp=
	VECTOR(*nattr)[(long int)VECTOR(news)[i]];
      SET_VECTOR_ELT(app, i, rep);
      SET_STRING_ELT(newnames, i, CREATE_STRING_VECTOR(tmp->name));
    }
    UNPROTECT(1); 		/* rep */
    PROTECT(newval=EVAL(lang3(install("c"), val, app)));
    PROTECT(newnames=EVAL(lang3(install("c"), names, newnames)));
    SET_NAMES(newval, newnames);
    SET_VECTOR_ELT(attr, 2, newval);
    val=VECTOR_ELT(attr, 2);    
    valno=GET_LENGTH(val);  
    names=GET_NAMES(val);
    UNPROTECT(4);
    rep=0;
  }
  igraph_vector_destroy(&news);
  IGRAPH_FINALLY_CLEAN(1);	/* news */

  /* Now append the new values */
  for (i=0; i<valno; i++) {
    SEXP oldva=VECTOR_ELT(val, i), newva;
    const char *sexpname=CHAR(STRING_ELT(names,i));
    igraph_bool_t l=0;
    long int j;
    for (j=0; !l && j<nattrno; j++) {
      igraph_attribute_record_t *tmp=VECTOR(*nattr)[j];
      l=!strcmp(sexpname, tmp->name);
    }
    if (l) {
      /* This attribute is present in nattr */
      SEXP app=0;
      igraph_attribute_record_t *tmprec=VECTOR(*nattr)[j-1];
      switch (tmprec->type) {
      case IGRAPH_ATTRIBUTE_NUMERIC:
	if (nv != igraph_vector_size(tmprec->value)) {
	  IGRAPH_ERROR("Invalid attribute length", IGRAPH_EINVAL);
	}
	PROTECT(app=NEW_NUMERIC(nv));
	igraph_vector_copy_to(tmprec->value, REAL(app));
	break;
      case IGRAPH_ATTRIBUTE_STRING:
	if (nv != igraph_strvector_size(tmprec->value)) {
	  IGRAPH_ERROR("Invalid attribute length", IGRAPH_EINVAL);
	}
	PROTECT(app=R_igraph_strvector_to_SEXP(tmprec->value));
	break;
      case IGRAPH_ATTRIBUTE_R_OBJECT:
	/* TODO */
	IGRAPH_ERROR("R_objects not implemented yet", IGRAPH_UNIMPLEMENTED);
	break;
      default:
	warning("Ignoring unknown attribute type");
	break;
      }
      if (app!=0) {
	PROTECT(newva=EVAL(lang3(install("c"), oldva, app)));
	SET_VECTOR_ELT(val, i, newva);
	UNPROTECT(2);		/* app & newva */
      }
    } else {
      /* No such attribute, append NA's */
      if (rep==0) {
	PROTECT(rep=EVAL(lang3(install("rep"), ScalarLogical(NA_LOGICAL), 
			       ScalarInteger(nv))));
      }
      PROTECT(newva=EVAL(lang3(install("c"), oldva, rep)));
      SET_VECTOR_ELT(val, i, newva);
      UNPROTECT(1); 		/* newva */
    }
  }

  if (had_names || added_names) {
    if (had_names && ! added_names) {
      /* The new vertices have no names, nothing to do */
    } else if (!had_names && added_names) {
      SET_VECTOR_ELT(attr, 4,
		     R_igraph_hash_create(VECTOR(*nattr)[added_names-1], 
					  origlen));
    } else if (had_names && added_names) {
      SET_VECTOR_ELT(attr, 4, 
		     R_igraph_hash_add(VECTOR_ELT(attr, 4), 
				       VECTOR_ELT(val, had_names-1),
				       ScalarInteger(origlen)));
    }
  }

  if (rep != 0) {
    UNPROTECT(1);
  } 
  
  return 0;
}

/* void R_igraph_attribute_delete_vertices(igraph_t *graph,  */
/* 					const igraph_vector_t *eidx, */
/* 					const igraph_vector_t *vidx) { */
/*   SEXP attr=graph->attr; */
/*   SEXP eal, val; */
/*   long int valno, ealno, i; */
/*   if (REAL(VECTOR_ELT(attr, 0))[0]+REAL(VECTOR_ELT(attr, 0))[1] > 1) { */
/*     SEXP newattr; */
/*     PROTECT(newattr=duplicate(attr)); */
/*     REAL(VECTOR_ELT(attr, 0))[1] -= 1; */
/*     if (REAL(VECTOR_ELT(attr, 0))[1] == 0) { */
/*       UNPROTECT_PTR(attr); */
/*     } */
/*     REAL(VECTOR_ELT(newattr, 0))[0] = 0; */
/*     REAL(VECTOR_ELT(newattr, 0))[1] = 1; */
/*     attr=graph->attr=newattr; */
/*   } */

/*   /\* Vertices *\/ */
/*   val=VECTOR_ELT(attr, 2); */
/*   valno=GET_LENGTH(val); */
/*   for (i=0; i<valno; i++) { */
/*     SEXP oldva=VECTOR_ELT(val, i), newva, ss; */
/*     long int origlen=GET_LENGTH(oldva); */
/*     long int newlen=0, j; */
/*     for (j=0; j<igraph_vector_size(vidx); j++) { */
/*       if (VECTOR(*vidx)[j] > 0) { */
/* 	newlen++; */
/*       } */
/*     } */
/*     PROTECT(ss=NEW_NUMERIC(newlen)); */
/*     for (j=0; j<origlen; j++) { */
/*       if (VECTOR(*vidx)[j]>0) { */
/* 	REAL(ss)[(long int)VECTOR(*vidx)[j]-1]=j+1; */
/*       } */
/*     } */
/*     PROTECT(newva=EVAL(lang3(install("["), oldva, ss))); */
/*     SET_VECTOR_ELT(val, i, newva); */
/*     UNPROTECT(2); */
/*   }     */

/*   /\* Edges *\/ */
/*   eal=VECTOR_ELT(attr, 3); */
/*   ealno=GET_LENGTH(eal); */
/*   for (i=0; i<ealno; i++) { */
/*     SEXP oldea=VECTOR_ELT(eal, i), newea, ss; */
/*     long int origlen=GET_LENGTH(oldea); */
/*     long int newlen=0, j; */
/*     /\* calculate new length *\/ */
/*     for (j=0; j<origlen; j++) { */
/*       if (VECTOR(*eidx)[j] > 0) { */
/* 	newlen++; */
/*       } */
/*     }     */
/*     PROTECT(ss=NEW_NUMERIC(newlen)); */
/*     for (j=0; j<origlen; j++) { */
/*       if (VECTOR(*eidx)[j]>0) { */
/* 	REAL(ss)[(long int)VECTOR(*eidx)[j]-1]=j+1; */
/*       } */
/*     } */
/*     PROTECT(newea=EVAL(lang3(install("["), oldea, ss))); */
/*     SET_VECTOR_ELT(eal, i, newea); */
/*     UNPROTECT(2); */
/*   } */
/* } */

int R_igraph_attribute_permute_vertices(const igraph_t *graph,
					igraph_t *newgraph,
					const igraph_vector_t *idx) {

  igraph_bool_t has_names;

  if (graph == newgraph) {

    SEXP attr=newgraph->attr;
    SEXP val, toval;
    SEXP names;
    long int i, valno;
    long int idxlen=igraph_vector_size(idx);
    SEXP ss;

    /* We copy if we need to */
    if (REAL(VECTOR_ELT(attr, 0))[0]+REAL(VECTOR_ELT(attr, 0))[1] > 1) {
      SEXP newattr;
      PROTECT(newattr=duplicate(attr));
      REAL(VECTOR_ELT(attr, 0))[1] -= 1;
      if (REAL(VECTOR_ELT(attr, 0))[1] == 0) {
	UNPROTECT_PTR(attr);
      }
      REAL(VECTOR_ELT(newattr, 0))[0] = 0;
      REAL(VECTOR_ELT(newattr, 0))[1] = 1;
      attr=newgraph->attr=newattr;
    }

    val=VECTOR_ELT(attr,2);
    names=GET_NAMES(val);
    valno=GET_LENGTH(val);
    
    /* If we have no vertex attributes, then we don't need to do anything */
    if (valno==0) { return 0; }
    
    /* Convert idx to an R object, we will use this for indexing */
    PROTECT(ss=NEW_INTEGER(idxlen));
    for (i=0; i<idxlen; i++) {
      INTEGER(ss)[i] = VECTOR(*idx)[i]+1;
    }
    
    for (i=0; i<valno; i++) {
      SEXP oldva=VECTOR_ELT(val, i);
      SEXP newva;

      /* We do NOT do any copying, the attributes were already copied, 
	 we're doing this in place. */
      PROTECT(newva=EVAL(lang3(install("["), oldva, ss)));

      if (!strcmp(CHAR(VECTOR_ELT(names, i)), R_IGRAPH_NAME)) {
	SET_VECTOR_ELT(attr, 4, R_igraph_hash_permute(VECTOR_ELT(attr, 4),
						      oldva, newva, ss));
      }

      SET_VECTOR_ELT(val, i, newva);
      UNPROTECT(1);
    }
    
    UNPROTECT(1);    

  } else {

    SEXP attr=graph->attr;
    SEXP toattr=newgraph->attr;
    SEXP val, toval;
    SEXP names;
    long int i, valno;
    long int idxlen=igraph_vector_size(idx);
    SEXP ss;

    val=VECTOR_ELT(attr,2);
    valno=GET_LENGTH(val);
    
    /* If we have no vertex attributes, then we don't need to do anything */
    if (valno==0) { return 0; }
    
    /* Convert idx to an R object, we will use this for indexing */
    PROTECT(ss=NEW_INTEGER(idxlen));
    for (i=0; i<idxlen; i++) {
      INTEGER(ss)[i] = VECTOR(*idx)[i]+1;
    }
    
    /* Resize the vertex attribute list in 'newgraph' */
    PROTECT(toval=NEW_LIST(valno));
    PROTECT(names=GET_NAMES(val));
    SET_NAMES(toval, names);
    UNPROTECT(1);
    
    for (i=0; i<valno; i++) {
      SEXP oldva=VECTOR_ELT(val, i);
      SEXP newva;

      PROTECT(newva=EVAL(lang3(install("["), oldva, ss)));

      if (!strcmp(CHAR(VECTOR_ELT(names, i)), R_IGRAPH_NAME)) {
	SET_VECTOR_ELT(toattr, 4,
		       R_igraph_hash_permute(VECTOR_ELT(attr, 4),
					     oldva, newva, ss));
      }
      
      SET_VECTOR_ELT(toval, i, newva);
      UNPROTECT(1);
    }
    
    SET_VECTOR_ELT(toattr, 2, toval);
    UNPROTECT(2);

  }
  
  return 0;
}

int R_igraph_attribute_add_edges(igraph_t *graph, 
				 const igraph_vector_t *edges,
				 igraph_vector_ptr_t *nattr) {
  SEXP attr=graph->attr;
  SEXP eal, rep=0, names, newnames;
  igraph_vector_t news;
  long int ealno, i, origlen, nattrno, newattrs;  
  long int ne=igraph_vector_size(edges)/2;
  if (REAL(VECTOR_ELT(attr, 0))[0]+REAL(VECTOR_ELT(attr, 0))[1] > 1) {
    SEXP newattr;
    PROTECT(newattr=duplicate(attr));
    REAL(VECTOR_ELT(attr, 0))[1] -= 1;
    if (REAL(VECTOR_ELT(attr, 0))[1] == 0) {
      UNPROTECT_PTR(attr);
    }
    REAL(VECTOR_ELT(newattr, 0))[0] = 0;
    REAL(VECTOR_ELT(newattr, 0))[1] = 1;
    attr=graph->attr=newattr;
  }

  eal=VECTOR_ELT(attr, 3);
  ealno=GET_LENGTH(eal);
  names=GET_NAMES(eal);
  if (nattr==NULL) {
    nattrno=0; 
  } else {
    nattrno=igraph_vector_ptr_size(nattr);
  }
  origlen=igraph_ecount(graph)-ne;

  /* First add the new attributes, if any */
  newattrs=0;
  IGRAPH_VECTOR_INIT_FINALLY(&news, 0);
  for (i=0; i<nattrno; i++) {
    igraph_attribute_record_t *nattr_entry=VECTOR(*nattr)[i];
    const char *nname=nattr_entry->name;
    long int j;
    igraph_bool_t l=0;
    for (j=0; !l && j<ealno; j++) {
      l=!strcmp(nname, CHAR(STRING_ELT(names, j)));
    }
    if (!l) {
      newattrs++;
      IGRAPH_CHECK(igraph_vector_push_back(&news, i));
    }
  }
  if (newattrs != 0) {
    SEXP app, neweal;
    PROTECT(app=NEW_LIST(newattrs));
    PROTECT(newnames=NEW_CHARACTER(newattrs));
    PROTECT(rep=EVAL(lang3(install("rep"), ScalarLogical(NA_LOGICAL),
			   ScalarInteger(origlen))));
    for (i=0; i<newattrs; i++) {
      igraph_attribute_record_t *tmp=
	VECTOR(*nattr)[ (long int) VECTOR(news)[i]];
      SET_VECTOR_ELT(app, i, rep);
      SET_STRING_ELT(newnames, i, CREATE_STRING_VECTOR(tmp->name));
    }
    UNPROTECT(1);		/* rep */
    PROTECT(neweal=EVAL(lang3(install("c"), eal, app)));
    PROTECT(newnames=EVAL(lang3(install("c"), names, newnames)));
    SET_NAMES(neweal, newnames);
    SET_VECTOR_ELT(attr, 3, neweal);
    eal=VECTOR_ELT(attr, 3);
    ealno=GET_LENGTH(eal);
    names=GET_NAMES(eal);
    UNPROTECT(4);
    rep=0;
  }
  igraph_vector_destroy(&news);
  IGRAPH_FINALLY_CLEAN(1);

  /* Now append the new values */
  for (i=0; i<ealno; i++) {
    SEXP oldea=VECTOR_ELT(eal, i), newea;
    const char *sexpname=CHAR(STRING_ELT(names, i));
    igraph_bool_t l=0;
    long int j;
    for (j=0; !l && j<nattrno; j++) {
      igraph_attribute_record_t *tmp=VECTOR(*nattr)[j];
      l=!strcmp(sexpname, tmp->name);
    }
    if (l) {
      /* This attribute is present in nattr */
      SEXP app=0;
      igraph_attribute_record_t *tmprec=VECTOR(*nattr)[j-1];
      switch (tmprec->type) {
      case IGRAPH_ATTRIBUTE_NUMERIC:
	if (ne != igraph_vector_size(tmprec->value)) {
	  IGRAPH_ERROR("Invalid attribute length", IGRAPH_EINVAL);
	}
	PROTECT(app=NEW_NUMERIC(ne));
	igraph_vector_copy_to(tmprec->value, REAL(app));
	break;
      case IGRAPH_ATTRIBUTE_STRING:
	if (ne != igraph_strvector_size(tmprec->value)) {
	  IGRAPH_ERROR("Invalid attribute length", IGRAPH_EINVAL);
	}
	PROTECT(app=R_igraph_strvector_to_SEXP(tmprec->value));
	break;
      case IGRAPH_ATTRIBUTE_R_OBJECT:
	/* TODO */
	IGRAPH_ERROR("R objects not implemented yet", IGRAPH_UNIMPLEMENTED);
	break;
      default:
	warning("Ignoring unknown attribute type");
	break;
      }
      if (app!=0) {
	PROTECT(newea=EVAL(lang3(install("c"), oldea, app)));
	SET_VECTOR_ELT(eal, i, newea);
	UNPROTECT(2);		/* app & newea */
      }
    } else {
      /* No such attribute, append NA's */
      if (rep==0) {
	PROTECT(rep=EVAL(lang3(install("rep"), ScalarLogical(NA_LOGICAL),
			       ScalarInteger(ne))));
      }
      PROTECT(newea=EVAL(lang3(install("c"), oldea, rep)));
      SET_VECTOR_ELT(eal, i, newea);
      UNPROTECT(1);		/* newea */
    }
  }
  if (rep != 0) {
    UNPROTECT(1);
  }

  return 0;
}

/* void R_igraph_attribute_delete_edges(igraph_t *graph,  */
/* 				     const igraph_vector_t *idx) { */
/*   SEXP attr=graph->attr; */
/*   SEXP eal; */
/*   long int ealno, i; */
/*   if (REAL(VECTOR_ELT(attr, 0))[0]+REAL(VECTOR_ELT(attr, 0))[1] > 1) { */
/*     SEXP newattr; */
/*     PROTECT(newattr=duplicate(attr)); */
/*     REAL(VECTOR_ELT(attr, 0))[1] -= 1; */
/*     if (REAL(VECTOR_ELT(attr, 0))[1] == 0) { */
/*       UNPROTECT_PTR(attr); */
/*     } */
/*     REAL(VECTOR_ELT(newattr, 0))[0] = 0; */
/*     REAL(VECTOR_ELT(newattr, 0))[1] = 1; */
/*     attr=graph->attr=newattr; */
/*   } */

/*   eal=VECTOR_ELT(attr, 3); */
/*   ealno=GET_LENGTH(eal); */
/*   for (i=0; i<ealno; i++) { */
/*     SEXP oldea=VECTOR_ELT(eal, i), newea, ss; */
/*     long int origlen=GET_LENGTH(oldea); */
/*     long int newlen=0, j; */
/*     /\* create subscript vector *\/ */
/*     for (j=0; j<origlen; j++) { */
/*       if (VECTOR(*idx)[j] > 0) { */
/* 	newlen++; */
/*       } */
/*     } */
/*     PROTECT(ss=NEW_NUMERIC(newlen)); */
/*     for (j=0; j<origlen; j++) { */
/*       if (VECTOR(*idx)[j] > 0) { */
/* 	REAL(ss)[(long int)VECTOR(*idx)[j]-1] = j+1; */
/*       } */
/*     } */
/*     PROTECT(newea=EVAL(lang3(install("["), oldea, ss))); */
/*     SET_VECTOR_ELT(eal, i, newea); */
/*     UNPROTECT(2); */
/*   } */
/* } */

int R_igraph_attribute_permute_edges(const igraph_t *graph,
				     igraph_t *newgraph,
				     const igraph_vector_t *idx) {

  if (graph==newgraph) {
    
    SEXP attr=newgraph->attr;
    SEXP eal;
    long int i, ealno;
    long int idxlen=igraph_vector_size(idx);
    SEXP ss;
    
    /* We copy if we need to */
    if (REAL(VECTOR_ELT(attr, 0))[0]+REAL(VECTOR_ELT(attr, 0))[1] > 1) {
      SEXP newattr;
      PROTECT(newattr=duplicate(attr));
      REAL(VECTOR_ELT(attr, 0))[1] -= 1;
      if (REAL(VECTOR_ELT(attr, 0))[1] == 0) {
      UNPROTECT_PTR(attr);
      }
      REAL(VECTOR_ELT(newattr, 0))[0] = 0;
      REAL(VECTOR_ELT(newattr, 0))[1] = 1;
      attr=newgraph->attr=newattr;
    }
    
    eal=VECTOR_ELT(attr,3);
    ealno=GET_LENGTH(eal);
    
    /* If we have no edge attributes, then we don't need to do anything */
    if (ealno==0) { return 0; }
    
    /* Convert idx to an R object, we will use this for indexing */
    PROTECT(ss=NEW_INTEGER(idxlen));
    for (i=0; i<idxlen; i++) {
      INTEGER(ss)[i] = VECTOR(*idx)[i]+1;
    }
    
    for (i=0; i<ealno; i++) {
      SEXP oldea=VECTOR_ELT(eal, i);
      SEXP newea;
      
      /* We do NOT do any copying, the attributes were already copied, 
	 we're doing this in place. */
      PROTECT(newea=EVAL(lang3(install("["), oldea, ss)));
      SET_VECTOR_ELT(eal, i, newea);
      UNPROTECT(1);
    }
    
    UNPROTECT(1);

  } else { 

    SEXP attr=graph->attr;
    SEXP toattr=newgraph->attr;
    SEXP eal, toeal;
    SEXP names;
    long int i, ealno;
    long int idxlen=igraph_vector_size(idx);
    SEXP ss;

    eal=VECTOR_ELT(attr,3);
    ealno=GET_LENGTH(eal);
    
    /* If we have no vertex attributes, then we don't need to do anything */
    if (ealno==0) { return 0; }
    
    /* Convert idx to an R object, we will use this for indexing */
    PROTECT(ss=NEW_INTEGER(idxlen));
    for (i=0; i<idxlen; i++) {
      INTEGER(ss)[i] = VECTOR(*idx)[i]+1;
    }
    
    /* Resize the vertex attribute list in 'newgraph' */
    PROTECT(toeal=NEW_LIST(ealno));
    PROTECT(names=GET_NAMES(eal));
    SET_NAMES(toeal, names);
    UNPROTECT(1);
    
    for (i=0; i<ealno; i++) {
      SEXP oldea=VECTOR_ELT(eal, i);
      SEXP newea;
      
      PROTECT(newea=EVAL(lang3(install("["), oldea, ss)));
      SET_VECTOR_ELT(toeal, i, newea);
      UNPROTECT(1);
    }
    
    SET_VECTOR_ELT(toattr, 3, toeal);
    UNPROTECT(2);
  }
  
  return 0;
}

int R_igraph_attribute_get_info(const igraph_t *graph,
				igraph_strvector_t *gnames,
				igraph_vector_t *gtypes,
				igraph_strvector_t *vnames,
				igraph_vector_t *vtypes,
				igraph_strvector_t *enames,
				igraph_vector_t *etypes) {
  igraph_strvector_t *names[3] = { gnames, vnames, enames };
  igraph_vector_t *types[3] = { gtypes, vtypes, etypes };
  long int i, j;

  SEXP attr=graph->attr;

  for (i=0; i<3; i++) {
    igraph_strvector_t *n=names[i];
    igraph_vector_t *t=types[i];
    SEXP al=VECTOR_ELT(attr, i+1);

    if (n) {			/* return names */
      SEXP names=GET_NAMES(al);
      R_igraph_SEXP_to_strvector_copy(names, n);
    }

    if (t) {			/* return types */
      igraph_vector_resize(t, GET_LENGTH(al));
      for (j=0; j<GET_LENGTH(al); j++) {
	SEXP a=VECTOR_ELT(al, j);
	if (TYPEOF(a)==REALSXP || TYPEOF(a)==INTSXP) {
	  VECTOR(*t)[j]=IGRAPH_ATTRIBUTE_NUMERIC;
	} else if (IS_CHARACTER(a)) {
	  VECTOR(*t)[j]=IGRAPH_ATTRIBUTE_STRING;
	} else {
	  VECTOR(*t)[j]=IGRAPH_ATTRIBUTE_R_OBJECT;
	}
      }
    }
  }

  return 0;
}

igraph_bool_t R_igraph_attribute_has_attr(const igraph_t *graph,
					  igraph_attribute_elemtype_t type,
					  const char *name) {
  long int attrnum;
  SEXP res;
  
  switch (type) {
  case IGRAPH_ATTRIBUTE_GRAPH:
    attrnum=1;
    break;
  case IGRAPH_ATTRIBUTE_VERTEX:
    attrnum=2;
    break;
  case IGRAPH_ATTRIBUTE_EDGE:
    attrnum=3;
    break;
  default:
    IGRAPH_ERROR("Unkwown attribute element type", IGRAPH_EINVAL);
    break;
  }
  
  res=R_igraph_getListElement(VECTOR_ELT(graph->attr, attrnum), name);
  return res != R_NilValue;
}

int R_igraph_attribute_gettype(const igraph_t *graph,
			       igraph_attribute_type_t *type,
			       igraph_attribute_elemtype_t elemtype,
			       const char *name) {
  long int attrnum;
  SEXP res;
  
  switch (elemtype) {
  case IGRAPH_ATTRIBUTE_GRAPH:
    attrnum=1;
    break;
  case IGRAPH_ATTRIBUTE_VERTEX:
    attrnum=2;
    break;
  case IGRAPH_ATTRIBUTE_EDGE:
    attrnum=3;
    break;
  default:
    IGRAPH_ERROR("Unkwown attribute element type", IGRAPH_EINVAL);
    break;
  }
  
  res=R_igraph_getListElement(VECTOR_ELT(graph->attr, attrnum), name);
  if (IS_NUMERIC(res) || IS_INTEGER(res)) {
    *type=IGRAPH_ATTRIBUTE_NUMERIC;
  } else if (IS_CHARACTER(res)) {
    *type=IGRAPH_ATTRIBUTE_STRING;
  } else {
    *type=IGRAPH_ATTRIBUTE_R_OBJECT;
  }
  return 0;
}

int R_igraph_attribute_get_numeric_graph_attr(const igraph_t *graph,
					      const char *name, 
					      igraph_vector_t *value) {
  SEXP gal=VECTOR_ELT(graph->attr, 1);
  SEXP ga=R_igraph_getListElement(gal, name);
  
  if (ga == R_NilValue) {
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  }

  PROTECT(ga=AS_NUMERIC(ga));
  
  IGRAPH_CHECK(igraph_vector_resize(value, 1));
  VECTOR(*value)[0]=REAL(ga)[0];

  UNPROTECT(1);

  return 0;
}

int R_igraph_attribute_get_string_graph_attr(const igraph_t *graph,
					     const char *name,
					     igraph_strvector_t *value) {
  /* TODO: serialization */
  SEXP gal=VECTOR_ELT(graph->attr, 1);
  SEXP ga=R_igraph_getListElement(gal, name);
  
  if (ga == R_NilValue) {
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  }

  PROTECT(ga=AS_CHARACTER(ga));

  IGRAPH_CHECK(igraph_strvector_resize(value, 1));
  IGRAPH_CHECK(igraph_strvector_set(value, 0, CHAR(STRING_ELT(ga, 0))));

  UNPROTECT(1);

  return 0;
}

int R_igraph_attribute_get_numeric_vertex_attr(const igraph_t *graph, 
					       const char *name,
					       igraph_vs_t vs,
					       igraph_vector_t *value) {
  /* TODO: serialization */
  SEXP val=VECTOR_ELT(graph->attr, 2);
  SEXP va=R_igraph_getListElement(val, name);
  igraph_vector_t newvalue;

  if (va == R_NilValue) {
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  }
  PROTECT(va=AS_NUMERIC(va));

  if (igraph_vs_is_all(&vs)) {
    R_SEXP_to_vector_copy(va, &newvalue);
    igraph_vector_destroy(value);
    *value=newvalue;
  } else {
    igraph_vit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
    IGRAPH_FINALLY(igraph_vit_destroy, &it);
    IGRAPH_CHECK(igraph_vector_resize(value, IGRAPH_VIT_SIZE(it)));
    while (!IGRAPH_VIT_END(it)) {
      long int v=IGRAPH_VIT_GET(it);
      VECTOR(*value)[i]=REAL(va)[v];
      IGRAPH_VIT_NEXT(it);
      i++;
    }
    igraph_vit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }

  UNPROTECT(1);
  return 0;
}

int R_igraph_attribute_get_string_vertex_attr(const igraph_t *graph, 
					      const char *name,
					      igraph_vs_t vs,
					      igraph_strvector_t *value) {
  /* TODO: serialization */
  SEXP val, va;

  val=VECTOR_ELT(graph->attr, 2);  
  va=R_igraph_getListElement(val, name);
  if (va == R_NilValue) {
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  }

  PROTECT(va=AS_CHARACTER(va));
  
  if (igraph_vs_is_all(&vs)) {
    R_igraph_SEXP_to_strvector_copy(va, value);
  } else {
    igraph_vit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
    IGRAPH_FINALLY(igraph_vit_destroy, &it);
    IGRAPH_CHECK(igraph_strvector_resize(value, IGRAPH_VIT_SIZE(it)));
    while (!IGRAPH_VIT_END(it)) {
      long int v=IGRAPH_VIT_GET(it);
      const char *str=CHAR(STRING_ELT(va, v));
      IGRAPH_CHECK(igraph_strvector_set(value, i, str));
      IGRAPH_VIT_NEXT(it);
      i++;
    }
    igraph_vit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }

  UNPROTECT(1);

  return 0;
}

int R_igraph_attribute_get_numeric_edge_attr(const igraph_t *graph,
					     const char *name,
					     igraph_es_t es,
					     igraph_vector_t *value) {
  /* TODO: serialization */
  SEXP eal=VECTOR_ELT(graph->attr, 3);
  SEXP ea=R_igraph_getListElement(eal, name);
  igraph_vector_t newvalue;

  if (ea == R_NilValue) {
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  }
  PROTECT(ea=AS_NUMERIC(ea));
  
  if (igraph_es_is_all(&es)) {    
    R_SEXP_to_vector_copy(AS_NUMERIC(ea), &newvalue);
    igraph_vector_destroy(value);
    *value=newvalue;
  } else {
    igraph_eit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);
    IGRAPH_CHECK(igraph_vector_resize(value, IGRAPH_EIT_SIZE(it)));
    while (!IGRAPH_EIT_END(it)) {
      long int e=IGRAPH_EIT_GET(it);
      VECTOR(*value)[i]=REAL(ea)[e];
      IGRAPH_EIT_NEXT(it);
      i++;
    }
    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }

  UNPROTECT(1);
  return 0;
}

int R_igraph_attribute_get_string_edge_attr(const igraph_t *graph,
					    const char *name,
					    igraph_es_t es,
					    igraph_strvector_t *value) {
  /* TODO: serialization */
  SEXP eal=VECTOR_ELT(graph->attr, 3);
  SEXP ea=R_igraph_getListElement(eal, name);
  
  if (ea == R_NilValue) {
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  }

  PROTECT(ea=AS_CHARACTER(ea));
  
  if (igraph_es_is_all(&es)) {
    R_igraph_SEXP_to_strvector_copy(ea, value);
  } else {
    igraph_eit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);
    IGRAPH_CHECK(igraph_strvector_resize(value, IGRAPH_EIT_SIZE(it)));
    while (!IGRAPH_EIT_END(it)) {
      long int e=IGRAPH_EIT_GET(it);
      const char *str=CHAR(STRING_ELT(ea, e));
      IGRAPH_CHECK(igraph_strvector_set(value, i, str));
      IGRAPH_EIT_NEXT(it);
      i++;
    }
    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  UNPROTECT(1);
  
  return 0;
}

SEXP R_igraph_ac_sum_numeric(SEXP attr, 
			     const igraph_vector_ptr_t *merges) {
  SEXP res;
  SEXP attr2;
  long int i, len=igraph_vector_ptr_size(merges);
  
  PROTECT(attr2=AS_NUMERIC(attr));
  PROTECT(res=NEW_NUMERIC(len));
  
  for (i=0; i<len; i++) {
    igraph_vector_t *v=VECTOR(*merges)[i];
    long int j, n=igraph_vector_size(v);
    igraph_real_t s=0.0;
    for (j=0; j<n; j++) {
      long int src=VECTOR(*v)[j];
      s += REAL(attr2)[src];
    }
    REAL(res)[i] = s;
  }
  
  UNPROTECT(2);
  return res;
}

SEXP R_igraph_ac_prod_numeric(SEXP attr, 
			      const igraph_vector_ptr_t *merges) {
  SEXP res;
  SEXP attr2;
  long int i, len=igraph_vector_ptr_size(merges);
  
  PROTECT(attr2=AS_NUMERIC(attr));
  PROTECT(res=NEW_NUMERIC(len));
  
  for (i=0; i<len; i++) {
    igraph_vector_t *v=VECTOR(*merges)[i];
    long int j, n=igraph_vector_size(v);
    igraph_real_t s=1.0;
    for (j=0; j<n; j++) {
      long int src=VECTOR(*v)[j];
      s *= REAL(attr2)[src];
    }
    REAL(res)[i] = s;
  }
  
  UNPROTECT(2);
  return res;
}

SEXP R_igraph_ac_min_numeric(SEXP attr, 
			     const igraph_vector_ptr_t *merges) {
  SEXP res;
  SEXP attr2;
  long int i, len=igraph_vector_ptr_size(merges);
  
  PROTECT(attr2=AS_NUMERIC(attr));
  PROTECT(res=NEW_NUMERIC(len));
  
  for (i=0; i<len; i++) {
    igraph_vector_t *v=VECTOR(*merges)[i];
    long int j, n=igraph_vector_size(v);
    igraph_real_t m= n > 0 ? REAL(attr2)[(long) VECTOR(*v)[0] ] : NA_REAL;
    for (j=1; j<n; j++) {
      long int src=VECTOR(*v)[j];
      igraph_real_t val= REAL(attr2)[src];
      if (val < m) {
	m=val;
      }
    }
    REAL(res)[i] = m;
  }
  
  UNPROTECT(2);
  return res;
}

SEXP R_igraph_ac_max_numeric(SEXP attr, 
			     const igraph_vector_ptr_t *merges) {
  SEXP res;
  SEXP attr2;
  long int i, len=igraph_vector_ptr_size(merges);
  
  PROTECT(attr2=AS_NUMERIC(attr));
  PROTECT(res=NEW_NUMERIC(len));
  
  for (i=0; i<len; i++) {
    igraph_vector_t *v=VECTOR(*merges)[i];
    long int j, n=igraph_vector_size(v);
    igraph_real_t m= n > 0 ? REAL(attr2)[(long) VECTOR(*v)[0] ] : NA_REAL;
    for (j=1; j<n; j++) {
      long int src=VECTOR(*v)[j];
      igraph_real_t val= REAL(attr2)[src];
      if (val > m) {
	m=val;
      }
    }
    REAL(res)[i] = m;
  }
  
  UNPROTECT(2);
  return res;
}

SEXP R_igraph_ac_random_numeric(SEXP attr, 
				const igraph_vector_ptr_t *merges) {
  SEXP res;
  SEXP attr2;
  long int i, len=igraph_vector_ptr_size(merges);
  
  PROTECT(attr2=AS_NUMERIC(attr));
  PROTECT(res=NEW_NUMERIC(len));

  RNG_BEGIN();

  for (i=0; i<len; i++) {
    igraph_vector_t *v=VECTOR(*merges)[i];
    long int j, n=igraph_vector_size(v);
    if (n==0) {
      REAL(res)[i]=NA_REAL;
    } else if (n==1) {
      REAL(res)[i]=REAL(attr2)[(long) VECTOR(*v)[0] ];
    } else {
      long int r=RNG_INTEGER(0,n-1);
      REAL(res)[i]=REAL(attr2)[(long) VECTOR(*v)[r] ];
    }
  }

  RNG_END();

  UNPROTECT(2);
  return res;
}

SEXP R_igraph_ac_first_numeric(SEXP attr, 
			       const igraph_vector_ptr_t *merges) {
  SEXP res;
  SEXP attr2;
  long int i, len=igraph_vector_ptr_size(merges);
  
  PROTECT(attr2=AS_NUMERIC(attr));
  PROTECT(res=NEW_NUMERIC(len));
  
  for (i=0; i<len; i++) {
    igraph_vector_t *v=VECTOR(*merges)[i];
    long int j, n=igraph_vector_size(v);
    if (n==0) {
      REAL(res)[i]=NA_REAL;
    } else {
      REAL(res)[i]=REAL(attr2)[(long) VECTOR(*v)[0] ];
    }
  }
  
  UNPROTECT(2);
  return res;
}

SEXP R_igraph_ac_last_numeric(SEXP attr, 
			      const igraph_vector_ptr_t *merges) {
  SEXP res;
  SEXP attr2;
  long int i, len=igraph_vector_ptr_size(merges);
  
  PROTECT(attr2=AS_NUMERIC(attr));
  PROTECT(res=NEW_NUMERIC(len));
  
  for (i=0; i<len; i++) {
    igraph_vector_t *v=VECTOR(*merges)[i];
    long int j, n=igraph_vector_size(v);
    if (n==0) {
      REAL(res)[i]=NA_REAL;
    } else {
      REAL(res)[i]=REAL(attr2)[(long) VECTOR(*v)[n-1] ];
    }
  }
  
  UNPROTECT(2);
  return res;
}

SEXP R_igraph_ac_mean_numeric(SEXP attr, 
			      const igraph_vector_ptr_t *merges) {
  SEXP res;
  SEXP attr2;
  long int i, len=igraph_vector_ptr_size(merges);
  
  PROTECT(attr2=AS_NUMERIC(attr));
  PROTECT(res=NEW_NUMERIC(len));
  
  for (i=0; i<len; i++) {
    igraph_vector_t *v=VECTOR(*merges)[i];
    long int j, n=igraph_vector_size(v);
    igraph_real_t s= n>0 ? 0.0 : NA_REAL;
    for (j=0; j<n; j++) {
      long int src=VECTOR(*v)[j];
      s += REAL(attr2)[src];
    }
    if (n>0) { s=s/n; }
    REAL(res)[i] = s;
  }
  
  UNPROTECT(2);
  return res;
}

SEXP R_igraph_ac_median_numeric(SEXP attr, 
				const igraph_vector_ptr_t *merges) {
  SEXP res;
  SEXP attr2;
  long int i, len=igraph_vector_ptr_size(merges);
  
  PROTECT(attr2=AS_NUMERIC(attr));
  PROTECT(res=NEW_NUMERIC(len));

  for (i=0; i<len; i++) {
    igraph_vector_t *v=VECTOR(*merges)[i];
    long int j, n=igraph_vector_size(v);
    SEXP tmp, call, tmp2;
    if (n==0) {
      REAL(res)[i] = NA_REAL;
    } else if (n==1) {
      REAL(res)[i] = REAL(attr2)[ (long) VECTOR(*v)[0] ];
    } else {
      PROTECT(tmp=NEW_NUMERIC(n));
      for (j=0; j<n; j++) {
	long int src=VECTOR(*v)[j];
	REAL(tmp)[j] = REAL(attr2)[src];
      }
      PROTECT(call=lang2(install("median"), tmp));
      PROTECT(tmp2=EVAL(call));
      REAL(res)[i] = REAL(tmp2)[0];
      UNPROTECT(3);
    }
  }
  
  UNPROTECT(2);
  return res;
}

SEXP R_igraph_ac_all_other(SEXP attr,
			   const igraph_vector_ptr_t *merges,
			   const char *function_name,
			   SEXP arg) {
  SEXP res, res2;
  long int i, len=igraph_vector_ptr_size(merges);

  PROTECT(res=NEW_LIST(len));
  
  for (i=0; i<len; i++) {
    igraph_vector_t *v=VECTOR(*merges)[i];
    long int j, n=igraph_vector_size(v);
    SEXP tmp;
    PROTECT(tmp=NEW_NUMERIC(n));
    for (j=0; j<n; j++) {
      long int src=VECTOR(*v)[j];
      REAL(tmp)[j] = src+1;
    }
    if (! arg) {
      SET_VECTOR_ELT(res, i, 
		     EVAL(lang2(install(function_name), 
				EVAL(lang3(install("["), attr, tmp)))));
    } else {
      SET_VECTOR_ELT(res, i, 
		     EVAL(lang3(install(function_name), 
				EVAL(lang3(install("["), attr, tmp)),
				arg)));
    }
    UNPROTECT(1);
  }
    
  if (IS_VECTOR(attr)) {
    /* try to simplify it */
    igraph_bool_t good=1;
    for (i=0; i<len; i++) {
      if (GET_LENGTH(VECTOR_ELT(res, i)) != 1) {
	good=0;
	break;
      }
    }
    if (good) {
      res2=EVAL(lang3(install("unlist"), res, 
		      ScalarLogical(0)));
      UNPROTECT(1);
      return res2;
    }
  }
  
  UNPROTECT(1);
  return res;
}

SEXP R_igraph_ac_func(SEXP attr, 
		      const igraph_vector_ptr_t *merges,
		      SEXP func) {
  SEXP res, res2;
  long int i, len=igraph_vector_ptr_size(merges);
  
  PROTECT(res=NEW_LIST(len));
  
  for (i=0; i<len; i++) {
    igraph_vector_t *v=VECTOR(*merges)[i];
    long int j, n=igraph_vector_size(v);
    SEXP tmp;
    PROTECT(tmp=NEW_NUMERIC(n));
    for (j=0; j<n; j++) {
      long int src=VECTOR(*v)[j];
      REAL(tmp)[j] = src+1;
    }
    SET_VECTOR_ELT(res, i, 
		   EVAL(lang2(func, 
			      EVAL(lang3(install("["), attr, tmp)))));
    UNPROTECT(1);
  }
  
  if (IS_VECTOR(attr)) {
    /* try to simplify it */
    igraph_bool_t good=1;
    for (i=0; i<len; i++) {
      if (GET_LENGTH(VECTOR_ELT(res, i)) != 1) {
	good=0;
	break;
      }
    }
    if (good) {
      res2=EVAL(lang3(install("unlist"), res, 
		      ScalarLogical(0)));
      UNPROTECT(1);
      return res2;
    }
  }
  
  UNPROTECT(1);
  return res;
}

/* TODO: update the hash of vertex names */

int R_igraph_attribute_combine_vertices(const igraph_t *graph,
			 igraph_t *newgraph,
			 const igraph_vector_ptr_t *merges,
			 const igraph_attribute_combination_t *comb) {
  
  SEXP attr=graph->attr;
  SEXP toattr=newgraph->attr;
  SEXP val=VECTOR_ELT(attr, 2);
  SEXP toval;
  long int i, j, valno=GET_LENGTH(val);
  SEXP names, newnames;
  SEXP res;
  int keepno=0;
  int *TODO;
  void **funcs;
  
  /* Create the TODO list first */
  PROTECT(names=GET_NAMES(val));
  TODO=igraph_Calloc(valno, int);
  if (!TODO) {
    IGRAPH_ERROR("Cannot combine vertex attributes", 
		 IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, TODO);
  funcs=igraph_Calloc(valno, void*);
  if (!funcs) {
    IGRAPH_ERROR("Cannot combine vertex attributes",
		 IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, funcs);
  for (i=0; i<valno; i++) {
    const char *name=CHAR(STRING_ELT(names, i));
    igraph_attribute_combination_type_t todo;
    void *voidfunc;
    igraph_attribute_combination_query(comb, name, &todo, &voidfunc);
    TODO[i]=todo;
    funcs[i]=voidfunc;
    if (todo != IGRAPH_ATTRIBUTE_COMBINE_IGNORE) {
      keepno++;
    }
  }

  /* Not safe to UNPROTECT attributes */
  REAL(VECTOR_ELT(attr, 0))[2]=0;
  REAL(VECTOR_ELT(toattr, 0))[2]=0;
  
  PROTECT(res=NEW_LIST(keepno));
  PROTECT(newnames=NEW_CHARACTER(keepno));
  for (i=0, j=0; i<valno; i++) {
    SEXP va=VECTOR_ELT(val, i);
    const char *name=CHAR(STRING_ELT(names, i));
    igraph_attribute_combination_type_t todo=TODO[i];
    igraph_attribute_type_t type;
    void *voidfunc=funcs[i];
    SEXP func;

    /* What kind of attribute */
    R_igraph_attribute_gettype(graph, &type, 
			       IGRAPH_ATTRIBUTE_VERTEX, name);
    
    switch (todo) {
    case IGRAPH_ATTRIBUTE_COMBINE_DEFAULT:
      /* Never happens from R */
    case IGRAPH_ATTRIBUTE_COMBINE_IGNORE:
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
      func=(SEXP)voidfunc;
      SET_VECTOR_ELT(res, j, R_igraph_ac_func(va, merges, func));
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_SUM:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_sum_numeric(va, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(va, merges, "sum", NULL));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_PROD:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_prod_numeric(va, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(va, merges, "prod", NULL));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_MIN:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_min_numeric(va, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(va, merges, "min", NULL));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_MAX:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_max_numeric(va, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(va, merges, "max", NULL));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_RANDOM:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_random_numeric(va, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(va, merges, "sample", 
						     ScalarInteger(1)));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_FIRST:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_first_numeric(va, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(va, merges, "head", 
						     ScalarInteger(1)));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_LAST:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_last_numeric(va, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(va, merges, "tail", 
						     ScalarInteger(1)));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_MEAN:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_mean_numeric(va, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(va, merges, "mean", NULL));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_MEDIAN:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_median_numeric(va, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(va, merges, "median", NULL));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_CONCAT:
      SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(va, merges, "c", NULL));
      break;
    default:
      IGRAPH_ERROR("Unknown attribute combination", 
		   IGRAPH_UNIMPLEMENTED);
      break;
    }

    if (todo != IGRAPH_ATTRIBUTE_COMBINE_IGNORE) {
      SET_STRING_ELT(newnames, j, STRING_ELT(names, i));
      j++;
    }
  } 

  /* It is now safe to UNPROTECT attributes */
  REAL(VECTOR_ELT(attr, 0))[2]=1;
  REAL(VECTOR_ELT(toattr, 0))[2]=1;

  igraph_free(funcs);
  igraph_free(TODO);
  IGRAPH_FINALLY_CLEAN(2);

  SET_NAMES(res, newnames);
  SET_VECTOR_ELT(toattr, 2, res);
  UNPROTECT(3);
  
  return 0;
}

int R_igraph_attribute_combine_edges(const igraph_t *graph,
			 igraph_t *newgraph,
			 const igraph_vector_ptr_t *merges,
			 const igraph_attribute_combination_t *comb) {

  SEXP attr=graph->attr;
  SEXP toattr=newgraph->attr;
  SEXP eal=VECTOR_ELT(attr, 3);
  SEXP toeal;
  long int i, j, ealno=GET_LENGTH(eal);
  SEXP names, newnames;
  SEXP res;
  int keepno=0;
  int *TODO;
  void **funcs;

  /* Create the TODO list first */
  PROTECT(names=GET_NAMES(eal));
  TODO=igraph_Calloc(ealno, int);
  if (!TODO) {
    IGRAPH_ERROR("Cannot combine edge attributes", 
		 IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, TODO);
  funcs=igraph_Calloc(ealno, void*);
  if (!funcs) {
    IGRAPH_ERROR("Cannot combine edge attributes",
		 IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, funcs);
  for (i=0; i<ealno; i++) {
    const char *name=CHAR(STRING_ELT(names, i));
    igraph_attribute_combination_type_t todo;
    void *voidfunc;
    igraph_attribute_combination_query(comb, name, &todo, &voidfunc);
    TODO[i]=todo;
    funcs[i]=voidfunc;
    if (todo != IGRAPH_ATTRIBUTE_COMBINE_IGNORE) {
      keepno++;
    }
  }

  /* Not safe to UNPROTECT attributes */
  REAL(VECTOR_ELT(attr, 0))[2]=0;
  REAL(VECTOR_ELT(toattr, 0))[2]=0;

  PROTECT(res=NEW_LIST(keepno));
  PROTECT(newnames=NEW_CHARACTER(keepno));
  for (i=0, j=0; i<ealno; i++) {
    SEXP ea=VECTOR_ELT(eal, i);
    const char *name=CHAR(STRING_ELT(names, i));
    igraph_attribute_combination_type_t todo=TODO[i];
    igraph_attribute_type_t type;
    void *voidfunc=funcs[i];
    SEXP func;

    /* What kind of attribute */
    R_igraph_attribute_gettype(graph, &type, IGRAPH_ATTRIBUTE_EDGE, name);
    
    switch (todo) {
    case IGRAPH_ATTRIBUTE_COMBINE_DEFAULT:
      /* Never happens from R */
    case IGRAPH_ATTRIBUTE_COMBINE_IGNORE:
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_FUNCTION:
      func=(SEXP)voidfunc;
      SET_VECTOR_ELT(res, j, R_igraph_ac_func(ea, merges, func));
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_SUM:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_sum_numeric(ea, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(ea, merges, "sum", NULL));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_PROD:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_prod_numeric(ea, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(ea, merges, "prod", NULL));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_MIN:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_min_numeric(ea, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(ea, merges, "min", NULL));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_MAX:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_max_numeric(ea, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(ea, merges, "max", NULL));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_RANDOM:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_random_numeric(ea, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(ea, merges, "sample", 
						     ScalarInteger(1)));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_FIRST:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_first_numeric(ea, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(ea, merges, "head", 
						     ScalarInteger(1)));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_LAST:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_last_numeric(ea, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(ea, merges, "tail", 
						     ScalarInteger(1)));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_MEAN:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_mean_numeric(ea, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(ea, merges, "mean", NULL));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_MEDIAN:
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	SET_VECTOR_ELT(res, j, R_igraph_ac_median_numeric(ea, merges));
      } else {
	SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(ea, merges, "median", NULL));
      }
      break;
    case IGRAPH_ATTRIBUTE_COMBINE_CONCAT:
      SET_VECTOR_ELT(res, j, R_igraph_ac_all_other(ea, merges, "c", NULL));
      break;
    default:
      IGRAPH_ERROR("Unknown attribute combination", 
		   IGRAPH_UNIMPLEMENTED);
      break;
    }

    if (todo != IGRAPH_ATTRIBUTE_COMBINE_IGNORE) {
      SET_STRING_ELT(newnames, j, STRING_ELT(names, i));
      j++;
    }
  }   

  /* It is now safe to UNPROTECT attributes */
  REAL(VECTOR_ELT(attr, 0))[2]=1;
  REAL(VECTOR_ELT(toattr, 0))[2]=1;

  igraph_free(funcs);
  igraph_free(TODO);
  IGRAPH_FINALLY_CLEAN(2);

  SET_NAMES(res, newnames);
  SET_VECTOR_ELT(toattr, 3, res);
  UNPROTECT(3);

  return 0;
}

igraph_attribute_table_t R_igraph_attribute_table={
  &R_igraph_attribute_init, &R_igraph_attribute_destroy,
  &R_igraph_attribute_copy, &R_igraph_attribute_add_vertices,
  &R_igraph_attribute_permute_vertices, 
  &R_igraph_attribute_combine_vertices,
  &R_igraph_attribute_add_edges,
  &R_igraph_attribute_permute_edges,
  &R_igraph_attribute_combine_edges,
  &R_igraph_attribute_get_info,
  &R_igraph_attribute_has_attr, &R_igraph_attribute_gettype,
  &R_igraph_attribute_get_numeric_graph_attr,
  &R_igraph_attribute_get_string_graph_attr,
  &R_igraph_attribute_get_numeric_vertex_attr,
  &R_igraph_attribute_get_string_vertex_attr,
  &R_igraph_attribute_get_numeric_edge_attr,
  &R_igraph_attribute_get_string_edge_attr
};

igraph_attribute_table_t *R_igraph_attribute_oldtable;

