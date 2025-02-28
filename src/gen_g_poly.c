#include "gen_g_poly.h"
#include <stdlib.h>

static inline void gen_monomials_str(char **monomials, slong size)
{
    for(slong i = 1; i < size+1; i++)
    {
        monomials[i-1] = (char*)flint_calloc(2, sizeof(char));
        monomials[i + size-1] = (char*)flint_calloc(2, sizeof(char));
        sprintf(monomials[i-1], "x%ld", i);
        sprintf(monomials[i+size-1], "y%ld", i);
    }
    return;
}

static inline void clear_monomials_str(char **monomials, slong size)
{
    for(slong i = 0; i < size; i++)
    {
        free(monomials[i]);
        free(monomials[i + size]);
    }
    free(monomials);
    return;
}

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

    fq_nmod_mpoly_set_fq_nmod(ux, (const fq_nmod_t){u[vec_size-2]}, mpoly_ring);
    fq_nmod_mpoly_set_fq_nmod(uy, (const fq_nmod_t){u[vec_size-1]}, mpoly_ring);
    fq_nmod_mpoly_set_fq_nmod(vx, (const fq_nmod_t){v[vec_size-2]}, mpoly_ring);
    fq_nmod_mpoly_set_fq_nmod(vy, (const fq_nmod_t){v[vec_size-1]}, mpoly_ring);

    //char **x = (char**)flint_calloc(2*(vec_size-2), sizeof(char*));
    //gen_monomials_str(x, vec_size-2);

    for(int i = 0; i < vec_size-2; i++)
    {
        fq_nmod_mpoly_gen(tmp2, i, mpoly_ring);
        
        //fq_nmod_mpoly_print_pretty(tmp2, (const char**)x, mpoly_ring);
        //printf(" ");

        fq_nmod_mpoly_scalar_mul_fq_nmod(tmp1, tmp2, (const fq_nmod_t){u[i]}, mpoly_ring);
        fq_nmod_mpoly_add(ux, ux, tmp1, mpoly_ring);

        fq_nmod_mpoly_scalar_mul_fq_nmod(tmp1, tmp2, (const fq_nmod_t){v[i]}, mpoly_ring);
        fq_nmod_mpoly_add(vx, vx, tmp1, mpoly_ring);


        fq_nmod_mpoly_gen(tmp2, i+vec_size-2, mpoly_ring);

        //fq_nmod_mpoly_print_pretty(tmp2, (const char**)x, mpoly_ring);
        //printf("\n");

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

    //clear_monomials_str(x, (vec_size-2));

    return;
}