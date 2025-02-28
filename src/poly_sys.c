#include "poly_sys.h"
#include "flint/fq_nmod.h"
#include "flint/nmod_poly.h"
#include <stdlib.h>
#include <stdio.h>

static inline void append_system(nmod_mpoly_t **system, fq_nmod_t c, const fq_nmod_ctx_t field, 
                    const nmod_mpoly_ctx_t mpoly_ring, slong i, slong k)
{
    slong degree = fq_nmod_ctx_degree(field);
    ulong p = fq_nmod_ctx_prime(field);

    nmod_poly_t t;
    nmod_poly_init(t, p);
    fq_nmod_get_nmod_poly(t, c, field);

    ulong coeff = 0;

    unsigned long* exp = (unsigned long*)flint_calloc(k, sizeof(unsigned long));

    if (exp == NULL) {
        fprintf(stderr, "Erreur d'allocation mémoire\n");
        return;
    }

    exp[i] = 1;
    
    for(slong j = 0; j < degree; j++)
    {
        coeff = nmod_poly_get_coeff_ui(t, j);
        nmod_mpoly_set_coeff_ui_ui((nmod_mpoly_struct *)(&(*system)[j]), coeff, exp, mpoly_ring);
    }

    flint_free(exp);

    nmod_poly_clear(t);

    return;
}

void init_system(nmod_mpoly_t **system, const nmod_mpoly_ctx_t mpoly_ring, slong k)
{
    *system = (nmod_mpoly_t *)flint_calloc(k, sizeof(nmod_mpoly_t));

    if (*system == NULL) {
        fprintf(stderr, "Erreur d'allocation mémoire pour le système de polynômes\n");
        return;
    }

    for (slong i = 0; i < k; i++)
        nmod_mpoly_init((nmod_mpoly_struct *)(&(*system)[i]), mpoly_ring);

    return;
}

void create_poly_system(fq_nmod_mpoly_t g, nmod_mpoly_t **system, const fq_nmod_mpoly_ctx_t mpoly_ring,
                    const nmod_mpoly_ctx_t system_mpoly_ring)
{
    fq_nmod_t c;
    fq_nmod_init(c, mpoly_ring->fqctx);

    slong k = fq_nmod_mpoly_ctx_nvars(mpoly_ring);
    
    for(slong i = 0; i < k; i++)
    {
        fq_nmod_mpoly_get_term_coeff_fq_nmod(c, g, i, mpoly_ring);
        append_system(system, c, mpoly_ring->fqctx, system_mpoly_ring, i, k);
    }

    fq_nmod_clear(c, mpoly_ring->fqctx);

    return;
}

void clear_system(nmod_mpoly_t **system, const nmod_mpoly_ctx_t mpoly_ring, slong k)
{
    for(slong i = 0; i < k; i++)
        nmod_mpoly_clear((nmod_mpoly_struct *)(&(*system)[i]), mpoly_ring);

    flint_free(*system);
}

void fprint_system(nmod_mpoly_t *system, const char **x, const nmod_mpoly_ctx_t mpoly_ring, const char *fn, slong k)
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
        
        nmod_mpoly_fprint_pretty(f, system[i], x, mpoly_ring);

        fprintf(f, "\n");
        
        fclose(f);
    }
}