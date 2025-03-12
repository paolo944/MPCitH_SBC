#include "poly_sys.h"
#include "flint/fq_nmod.h"
#include "flint/nmod_poly.h"
#include <stdlib.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static inline void append_system(nmod_mpoly_t **system, fq_nmod_t c, fq_nmod_mpoly_t m, 
                    const fq_nmod_ctx_t field, const nmod_mpoly_ctx_t system_mpoly_ring, 
                    const fq_nmod_mpoly_ctx_t mpoly_ring)
{
    slong degree = fq_nmod_ctx_degree(field);
    ulong p = fq_nmod_ctx_prime(field);

    nmod_poly_t t;
    nmod_poly_init(t, p);
    fq_nmod_get_nmod_poly(t, c, field);

    ulong coeff = 0;

    unsigned long *exp = (unsigned long*)calloc(fq_nmod_mpoly_ctx_nvars(mpoly_ring), sizeof(unsigned long));
    fq_nmod_mpoly_get_term_exp_ui(exp, m, 0, mpoly_ring);

    #pragma omp parallel for num_threads(8) schedule(dynamic)
    for(slong j = 0; j < degree; j++)
    {
        coeff = nmod_poly_get_coeff_ui(t, j);
        if(coeff != 0)
            nmod_mpoly_set_coeff_ui_ui((nmod_mpoly_struct *)(*system + j), coeff, exp, system_mpoly_ring);
    }

    free(exp);

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
        nmod_mpoly_init((nmod_mpoly_struct *)(*system + 1), mpoly_ring);

    return;
}

void create_poly_system(fq_nmod_mpoly_t g, nmod_mpoly_t **system, const fq_nmod_mpoly_ctx_t mpoly_ring,
                    const nmod_mpoly_ctx_t system_mpoly_ring)
{
    slong nb_monomials_g;
    fq_nmod_t c;
    fq_nmod_mpoly_t m;

    fq_nmod_init(c, mpoly_ring->fqctx);
    fq_nmod_mpoly_init(m, mpoly_ring);

    nb_monomials_g = fq_nmod_mpoly_length(g, mpoly_ring);

    printf("number of monomials: %ld\n", nb_monomials_g);

    slong dixpourcent = nb_monomials_g / 10;

    for(slong i = 0; i < nb_monomials_g; i++)
    {
        if(i % dixpourcent == 0)
            printf("%ld/%ld\n", i, nb_monomials_g);
        fq_nmod_mpoly_get_term_monomial(m, g, i, mpoly_ring);
        fq_nmod_mpoly_get_term_coeff_fq_nmod(c, g, i, mpoly_ring);
        append_system(system, c, m, mpoly_ring->fqctx, system_mpoly_ring, mpoly_ring);
    }

    fq_nmod_clear(c, mpoly_ring->fqctx);
    fq_nmod_mpoly_clear(m, mpoly_ring);

    return;
}

void clear_system(nmod_mpoly_t **system, const nmod_mpoly_ctx_t mpoly_ring, slong k)
{
    for(slong i = 0; i < k; i++)
        nmod_mpoly_clear((nmod_mpoly_struct *)(*system+i), mpoly_ring);

    flint_free(*system);
}

void fprint_system(nmod_mpoly_t *system, const char **x, const nmod_mpoly_ctx_t mpoly_ring, const char *fn, slong k)
{
    FILE *f = fopen(fn, "w");
    if (!f) {
        perror("Error while opening the file\n");
        exit(EXIT_FAILURE);
    }

    for(slong i = 0; i < k-2; i++)
    {
        fprintf(f, "%s,", x[i]);
    }
    fprintf(f, "%s\n", x[k-2]);

    ulong p = nmod_mpoly_ctx_modulus(mpoly_ring);
    fprintf(f, "%ld\n", p);

    fclose(f);

    slong nnz = 0;

    for(slong i = 0; i < k-1; i++)
    {
        f = fopen(fn, "a");
        if (!f) {
            perror("Error while opening the file\n");
            exit(EXIT_FAILURE);
        }

        if(nmod_mpoly_is_canonical(system[i], mpoly_ring))
            nnz += nmod_mpoly_length(system[i], mpoly_ring);
        else
            printf("erreur, polynôme non canonique");

        nmod_mpoly_fprint_pretty(f, system[i], x, mpoly_ring);

        fprintf(f, ",\n");
        
        fclose(f);
    }
    f = fopen(fn, "a");
    if (!f) {
        perror("Error while opening the file\n");
        exit(EXIT_FAILURE);
    }
    
    if(nmod_mpoly_is_canonical(system[k-1], mpoly_ring))
            nnz += nmod_mpoly_length(system[k-1], mpoly_ring);
        else
            printf("erreur, polynôme non canonique");

    nmod_mpoly_fprint_pretty(f, system[k-1], x, mpoly_ring);

    printf("nnz: %ld\n", nnz);

    fclose(f);
}

int is_regular_seq(nmod_mpoly_t *system, const nmod_mpoly_ctx_t mpoly_ring, slong k)
{
    nmod_mpoly_t *Q;
    nmod_mpoly_t R;

    nmod_mpoly_init(R, mpoly_ring);

    init_system(&Q, mpoly_ring, k);

    for(slong i = 1; i < k; i++)
    {
        printf("test pour k = %ld\n", i);
        nmod_mpoly_divrem_ideal(&Q, R, system[i], (nmod_mpoly_struct *const*) system, i, mpoly_ring);
        if(nmod_mpoly_is_zero(R, mpoly_ring)){
            clear_system(Q, mpoly_ring, k);
            return 0;
        }
    }

    clear_system(Q, mpoly_ring, k);

    return 1;
}