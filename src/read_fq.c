#include "read_fq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "flint/nmod_poly.h"
#include "flint/fq_nmod.h"

void read_fq_nmod_vec(fq_nmod_struct *vec, const char *fn, int size, const fq_nmod_ctx_t field)
{
    FILE *f = fopen(fn, "r");
    if (!f) {
        perror("Error while opening the file\n");
        exit(EXIT_FAILURE);
    }

    int size2 = 0;
    fscanf(f, "%d", &size2);

    if(size2 != size){
        perror("The size of the vector and the size given doesn't match\n");
        exit(EXIT_FAILURE);
    }

    printf("Reading the vector of size %d\n", size);

    fseek(f, 2, SEEK_CUR);

    nmod_poly_t tmp_poly;
    nmod_poly_init(tmp_poly, 2);

    for(int i = 0; i < size; i++)
    {
        if(nmod_poly_fread(f, tmp_poly) <= 0 && i != size-1){
            perror("Format error while reading the file\n");
            exit(EXIT_FAILURE);
        }
        fq_nmod_set_nmod_poly(&vec[i], tmp_poly, field);
    }

    nmod_poly_clear(tmp_poly);
    fclose(f);
    return;
}