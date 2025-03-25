#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <time.h>
#include <string.h>

#include "flint/flint.h"
#include "flint/fq_nmod.h"
#include "flint/fq_nmod_vec.h"
#include "flint/fq_nmod_mpoly.h"
#include "flint/nmod_mpoly.h"
#include "read_fq.h"
#include "gen_g_poly.h"
#include "poly_sys.h"

static inline void gen_monomials_str(char **monomials, slong size)
{
    for(slong i = 0; i < size; i++)
    {
        monomials[i] = (char*)calloc(5, sizeof(char));
        monomials[i + size] = (char*)calloc(5, sizeof(char));
        sprintf(monomials[i], "x%ld", i+1);
        sprintf(monomials[i + size], "y%ld", i+1);
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

int main(int argc, char **argv)
{
    if(argc != 4)
        errx(1, "Il manque des paramètres, il faut lancer comme ceci: ./model" \
        " (taille vecteur) (format) (1 si " \
        "inclure equations du corps ou 0 sinon)");

    struct timespec start, end;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &start);

    // Random init
    flint_rand_t state;
    flint_randinit(state);

    ulong q = 2;            // Field characteristic
    slong n = atoi(argv[1]);// Vectors size
    slong k = 2*(n-2)+1;    // Degree of field extension
    slong nvars = 2*(n-2);

    int format = 0;

    if(strcmp(argv[2], "hpXbred") == 0)
        format = 0;
    else if(strcmp(argv[2], "msolve") == 0)
        format = 1;
    else if(strcmp(argv[2], "magma") == 0)
        format = 2;
    else
        errx(1, "deuxième paramètre non reconnu, soit msolve soit hpXbred, soit magma");

    int field_eq = atoi(argv[3]);

    if(field_eq != 0 && field_eq != 1)
        errx(1, "le paramètre pour les équations du corps doit être 0 ou 1");
    
    char file_name[64];

    if(format == 0)
    	sprintf(file_name, "system/hpXbred/system_bilin_%ld_%ld.in", nvars, k+field_eq*nvars);
    else if(format == 1)
    	sprintf(file_name, "system/msolve/system_bilin_%ld_%ld.ms", nvars, k+field_eq*nvars);
    else if(format == 2)
    	sprintf(file_name, "system/magma/system_bilin_%ld_%ld.magma", nvars, k+field_eq*nvars);

    printf("q = %ld\nn = %ld\nk = %ld\nnvars = %ld\n", q, n, k, nvars);

    // Field and ring init
    const char *var = "t";

    fq_nmod_ctx_t field;
    fq_nmod_ctx_init_ui(field, q, k, var);

    fq_nmod_mpoly_ctx_t mpoly_ring;
    fq_nmod_mpoly_ctx_init(mpoly_ring, nvars, ORD_LEX, field);

    // u and v init
    fq_nmod_struct *u, *v;
    u = _fq_nmod_vec_init(n, field);
    v = _fq_nmod_vec_init(n, field);
    
    fq_nmod_struct *x, *y;
    x = _fq_nmod_vec_init(n, field);
    y = _fq_nmod_vec_init(n, field);

    // Read u and v
    printf("-------Reading the 4 vectors/keys\n");
    read_fq_nmod_vec(u, "keys/u.pub", n, field);
    read_fq_nmod_vec(v, "keys/v.pub", n, field);
    read_fq_nmod_vec(x, "keys/x", n, field);
    read_fq_nmod_vec(y, "keys/y", n, field);

    //Compute the polynomial g
    fq_nmod_mpoly_t g;
    fq_nmod_mpoly_init(g, mpoly_ring);

    printf("-------Generating the g polynomial\n");
    gen_g_poly(g, u, v, mpoly_ring, n);

    printf("-------Generated g\n");

    // char **monomials = (char**)calloc(2*(n-2), sizeof(char*));
    // gen_monomials_str(monomials, n-2);
    // fq_nmod_mpoly_print_pretty(g, (const char**)monomials, mpoly_ring);
    // clear_monomials_str(monomials, n-2);
// 
    // for(slong i = 0; i < g->length; i++)
    // {
        // printf("i: %ld coeff: %ld exp: %ld\n", i, g->coeffs[i], g->exps[i]);
    // }

    nmod_mpoly_t *system;
    nmod_mpoly_ctx_t system_mpoly_ring;

    nmod_mpoly_ctx_init(system_mpoly_ring, nvars, ORD_DEGREVLEX, q);

    init_system(&system, system_mpoly_ring, k+nvars);

    printf("-------Modelising the system\n");

    create_poly_system(g, &system, mpoly_ring, system_mpoly_ring);

    printf("-------Writing the system in %s\n", file_name);

    

    char **monomials = (char**)calloc(2*(n-2), sizeof(char*));
    gen_monomials_str(monomials, n-2);
    fprint_system(system, (const char**)monomials, system_mpoly_ring, file_name, nvars, k+nvars, format, field_eq);
    clear_monomials_str(monomials, n-2);
    clear_system(&system, system_mpoly_ring, k+nvars);
 
 
    // Test the keys by evaluating g on x and y
    /*
    fq_nmod_t ev;
    fq_nmod_init(ev, field);

    slong i;

    fq_nmod_struct **vals;
    vals = (fq_nmod_struct **) flint_malloc(nvars*sizeof(fq_nmod_struct *));
    for (i = 0; i < nvars/2; i++)
    {
        vals[i] = (fq_nmod_struct *) flint_malloc(sizeof(fq_nmod_struct));
        fq_nmod_init(vals[i], field);
        fq_nmod_set(vals[i], (fq_nmod_t){x[i]}, field);
    }
    for (i = nvars/2; i < nvars; i++)
    {
        vals[i] = (fq_nmod_struct *) flint_malloc(sizeof(fq_nmod_struct));
        fq_nmod_init(vals[i], field);
        fq_nmod_set(vals[i], (fq_nmod_t){y[i - nvars/2]}, field);
    }

    fq_nmod_mpoly_evaluate_all_fq_nmod(ev, g, vals, mpoly_ring);
    
    printf("-------Testing if g(x, y) = 0\n");
    if(fq_nmod_is_zero(ev, field))
        printf("\tg(x, y) = 0\n");
    else{
        printf("\tg(x, y) != 0\n\t");
        fq_nmod_print_pretty(ev, field);
    }
   
    for (slong i = 0; i < nvars; i++)
    {
        fq_nmod_clear(vals[i], field);
        flint_free(vals[i]);
    }
    flint_free(vals);

    fq_nmod_clear(ev, field);
    */

    fq_nmod_mpoly_clear(g, mpoly_ring);

    // Clear polys
    _fq_nmod_vec_clear(u, n, field);
    _fq_nmod_vec_clear(v, n, field);
    _fq_nmod_vec_clear(x, n, field);
    _fq_nmod_vec_clear(y, n, field);

    // Field and ring clean
    fq_nmod_mpoly_ctx_clear(mpoly_ring);
    nmod_mpoly_ctx_clear(system_mpoly_ring);
    fq_nmod_ctx_clear(field);

    flint_randclear(state);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("Time: %.9f secondes\n", elapsed);


    return 0;
}
