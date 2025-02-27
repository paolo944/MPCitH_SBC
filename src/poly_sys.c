#include "poly_sys.h"
#include "flint/fq_nmod.h"
#include "flint/nmod_poly.h"
#include <stdlib.h>
#include <stdio.h>

void create_poly_sys(fq_nmod_mpoly_t g, fq_nmod_mpoly_t **system, const fq_nmod_mpoly_ctx_t mpoly_ring, slong k, const char **x)
{
    slong i;
    fq_nmod_mpoly_t q, r, divisor;
    fq_nmod_t fq_divisor;
    nmod_poly_t poly_divisor;

    fq_nmod_init(fq_divisor, mpoly_ring->fqctx);
    nmod_poly_init(poly_divisor, 2);

    fq_nmod_mpoly_init(q, mpoly_ring);
    fq_nmod_mpoly_init(r, mpoly_ring);
    fq_nmod_mpoly_init(divisor, mpoly_ring);

    fq_nmod_mpoly_set(r, g, mpoly_ring);

    system = (fq_nmod_mpoly_t **)flint_malloc(k*sizeof(fq_nmod_mpoly_t *));

    i = k-1;

    do
    {
        system[i] = (fq_nmod_mpoly_t *)flint_malloc(sizeof(fq_nmod_mpoly_t));
        fq_nmod_mpoly_init(system[i], mpoly_ring);
        // Set to t^k-1 then t^k-1... t, 1)
        nmod_poly_zero(poly_divisor);
        nmod_poly_set_coeff_ui(poly_divisor, i, 1);
        nmod_poly_print_pretty(poly_divisor, "t");
        fq_nmod_set_nmod_poly(fq_divisor, poly_divisor, mpoly_ring->fqctx);
        fq_nmod_mpoly_set_fq_nmod(divisor, fq_divisor, mpoly_ring);

        fq_nmod_mpoly_divrem(q, r, r, divisor, mpoly_ring);
        
        printf(" q: ");
        fq_nmod_mpoly_print_pretty(q, x, mpoly_ring);
        printf(" r: ");
        fq_nmod_mpoly_print_pretty(r, x, mpoly_ring);

        fq_nmod_mpoly_set(*system[i], q, mpoly_ring);
        
        i -= 1;
    }while(i > -1);

    fq_nmod_clear(fq_divisor, mpoly_ring->fqctx);
    
    nmod_poly_clear(poly_divisor);

    fq_nmod_mpoly_clear(q, mpoly_ring);
    fq_nmod_mpoly_clear(r, mpoly_ring);
    fq_nmod_mpoly_clear(divisor, mpoly_ring);

    return;
}

void clear_sys(fq_nmod_mpoly_t *system, fq_nmod_mpoly_ctx_t mpoly_ring, slong k)
{
    for(slong i = 0; i < k; i++)
    {
        printf("clear i: %ld\n", i);
        fq_nmod_mpoly_clear(system[i], mpoly_ring);
        flint_free(system[i]);
    }
    flint_free(system);
}

void fprint_sys(fq_nmod_mpoly_t *system, const char **x, const fq_nmod_mpoly_ctx_t mpoly_ring, slong k, const char *fn)
{
    FILE *f = fopen(fn, "w");
    if (!f) {
        perror("Error while opening the file\n");
        exit(EXIT_FAILURE);
    }

    fclose(f);

    for(slong i = 0; i < k; i++)
    {
        f = fopen(fn, "a");
        if (!f) {
            perror("Error while opening the file\n");
            exit(EXIT_FAILURE);
        }
        
        fq_nmod_mpoly_fprint_pretty(f, system[i], x, mpoly_ring);
        
        fclose(f);
    }
}