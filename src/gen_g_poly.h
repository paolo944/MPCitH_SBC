#ifndef GEN_G_POLY_H
#define GEN_G_POLY_H

#include "flint/fq_nmod_mpoly.h"

void gen_g_poly(fq_nmod_mpoly_t *g, fq_nmod_struct *u, fq_nmod_struct *v, const fq_nmod_mpoly_ctx_t mpoly_ring);

#endif