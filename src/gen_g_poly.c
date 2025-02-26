#include "gen_g_poly.h"
#include <stdlib.h>

void gen_g_poly(fq_nmod_mpoly_t g, fq_nmod_struct *u, fq_nmod_struct *v, 
                const fq_nmod_mpoly_ctx_t mpoly_ring, int vec_size)
{
    fq_nmod_mpoly_t ux, uy, vx, vy, tmp1, tmp2;
    fq_nmod_mpoly_init(ux, mpoly_ring);
    fq_nmod_mpoly_init(uy, mpoly_ring);
    fq_nmod_mpoly_init(vx, mpoly_ring);
    fq_nmod_mpoly_init(vy, mpoly_ring);
    fq_nmod_mpoly_init(tmp1, mpoly_ring);
    fq_nmod_mpoly_init(tmp2, mpoly_ring);

    fq_nmod_mpoly_add_fq_nmod(ux, ux, (const fq_nmod_t){u[vec_size-2]}, mpoly_ring);
    fq_nmod_mpoly_add_fq_nmod(uy, uy, (const fq_nmod_t){u[vec_size-1]}, mpoly_ring);
    fq_nmod_mpoly_add_fq_nmod(vx, vx, (const fq_nmod_t){v[vec_size-2]}, mpoly_ring);
    fq_nmod_mpoly_add_fq_nmod(vy, vy, (const fq_nmod_t){u[vec_size-1]}, mpoly_ring);

    for(int i = 0; i < vec_size-2; i++)
    {
        fq_nmod_mpoly_gen(tmp2, i, mpoly_ring);
        
        fq_nmod_mpoly_scalar_mul_fq_nmod(tmp1, tmp2, (const fq_nmod_t){u[i]}, mpoly_ring);
        fq_nmod_mpoly_add(ux, ux, tmp1, mpoly_ring);

        fq_nmod_mpoly_scalar_mul_fq_nmod(tmp1, tmp2, (const fq_nmod_t){v[i]}, mpoly_ring);
        fq_nmod_mpoly_add(vx, vx, tmp1, mpoly_ring);


        fq_nmod_mpoly_gen(tmp2, i+vec_size-2, mpoly_ring);
        
        fq_nmod_mpoly_scalar_mul_fq_nmod(tmp1, tmp2, (const fq_nmod_t){u[i]}, mpoly_ring);
        fq_nmod_mpoly_add(ux, ux, tmp1, mpoly_ring);

        fq_nmod_mpoly_scalar_mul_fq_nmod(tmp1, tmp2, (const fq_nmod_t){v[i]}, mpoly_ring);
        fq_nmod_mpoly_add(vx, vx, tmp1, mpoly_ring);
    }

    fq_nmod_mpoly_mul(tmp1, ux, vy, mpoly_ring);
    fq_nmod_mpoly_mul(tmp2, uy, vx, mpoly_ring);

    fq_nmod_mpoly_sub(g, tmp1, tmp2, mpoly_ring);

    fq_nmod_mpoly_clear(ux, mpoly_ring);
    fq_nmod_mpoly_clear(uy, mpoly_ring);
    fq_nmod_mpoly_clear(vx, mpoly_ring);
    fq_nmod_mpoly_clear(vy, mpoly_ring);
    fq_nmod_mpoly_clear(tmp1, mpoly_ring);
    fq_nmod_mpoly_clear(tmp2, mpoly_ring);

    return;
}