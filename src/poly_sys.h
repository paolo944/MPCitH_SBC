#ifndef POLY_SYS_H
#define POLY_SYS_H

#include "flint/fq_nmod_mpoly.h"

void create_poly_sys(fq_nmod_mpoly_t g, fq_nmod_mpoly_t **system, const fq_nmod_mpoly_ctx_t mpoly_ring, slong k, const char **x);

void clear_sys(fq_nmod_mpoly_t *system, fq_nmod_mpoly_ctx_t mpoly_ring, slong k);

void fprint_sys(fq_nmod_mpoly_t *system, const char **x, const fq_nmod_mpoly_ctx_t mpoly_ring, slong k, const char *fn);

#endif