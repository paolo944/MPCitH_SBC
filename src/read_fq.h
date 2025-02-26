#ifndef READ_FQ_H
#define READ_FQ_H

#include "flint/fq_nmod_vec.h"

void read_fq_nmod_vec(fq_nmod_struct *vec, const char *fn, int size, const fq_nmod_ctx_t field);

#endif