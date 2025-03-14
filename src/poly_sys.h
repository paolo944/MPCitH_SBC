#ifndef POLY_SYS_H
#define POLY_SYS_H

#include "flint/fq_nmod_mpoly.h"
#include "flint/nmod_mpoly.h"

void init_system(nmod_mpoly_t **system, const nmod_mpoly_ctx_t mpoly_ring, slong k);

void create_poly_system(fq_nmod_mpoly_t g, nmod_mpoly_t **system, const fq_nmod_mpoly_ctx_t mpoly_ring,
    const nmod_mpoly_ctx_t system_mpoly_ring);

void clear_system(nmod_mpoly_t **system, const nmod_mpoly_ctx_t mpoly_ring, slong k);

void fprint_system(nmod_mpoly_t *system, const char **x, const nmod_mpoly_ctx_t mpoly_ring, const char *fn, slong nvars, slong k);

#endif