#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fq_nmod.h"
#include "flint/fq_nmod_vec.h"
#include "flint/fq_nmod_mpoly.h"
#include "read_fq.h"
#include "gen_g_poly.h"

static inline void gen_monomials_str(char **monomials, slong size)
{
    for(slong i = 0; i < size; i++)
    {
        monomials[i] = (char*)malloc(2*sizeof(char));
        monomials[i + size] = (char*)malloc(2*sizeof(char));
        sprintf(monomials[i], "x%ld", i);
        sprintf(monomials[i + size], "y%ld", i);
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

int main(void)
{
    // Random init
    flint_rand_t state;
    flint_randinit(state);

    ulong q = 2;    // Field characteristic
    slong n = 130;  // Vectors size
    slong k = 257;  // Degree of field extension

    // Field and ring init
    const char *var = "t";

    fq_nmod_ctx_t field;
    fq_nmod_ctx_init_ui(field, q, k, var);

    fq_nmod_mpoly_ctx_t mpoly_ring;
    fq_nmod_mpoly_ctx_init(mpoly_ring, 2*(n-2), ORD_DEGREVLEX, field);

    // u and v init
    fq_nmod_struct *u, *v, *x, *y;
    u = _fq_nmod_vec_init(n, field);
    v = _fq_nmod_vec_init(n, field);
    x = _fq_nmod_vec_init(n, field);
    y = _fq_nmod_vec_init(n, field);

    // Read u and v
    printf("Reading the 4 vectors/keys\n");
    read_fq_nmod_vec(u, "keys/u.pub", n, field);
    read_fq_nmod_vec(v, "keys/v.pub", n, field);
    read_fq_nmod_vec(x, "keys/x", n, field);
    read_fq_nmod_vec(y, "keys/y", n, field);

    //Compute the polynomial g
    fq_nmod_mpoly_t g;
    fq_nmod_mpoly_init(g, mpoly_ring);

    printf("generating the g polynomial\n");
    gen_g_poly(g, u, v, mpoly_ring, n);

    printf("generated g\n");

    //char **monomials = (char**)malloc(2*(n-2)*sizeof(char*));
    //gen_monomials_str(monomials, n-2);
    //fq_nmod_mpoly_print_pretty(g, (const char**)monomials, mpoly_ring);

    // Clear polys
    _fq_nmod_vec_clear(u, n, field);
    _fq_nmod_vec_clear(v, n, field);
    _fq_nmod_vec_clear(x, n, field);
    _fq_nmod_vec_clear(y, n, field);

    fq_nmod_mpoly_clear(g, mpoly_ring);

    // Field and ring clean
    fq_nmod_mpoly_ctx_clear(mpoly_ring);
    fq_nmod_ctx_clear(field);

    flint_randclear(state);

    //clear_monomials_str(monomials, n-2);

    return 0;
}
