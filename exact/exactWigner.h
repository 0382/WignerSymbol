#pragma once
#ifndef EXACTWIGNER_H
#define EXACTWIGNER_H

#include <gmp.h>

#ifdef __cplusplus
extern "C"
{
using _Bool = bool;
#endif

typedef struct _sqrt_rational
{
    // sn/sd * \sqrt{rn/rd}
    mpz_t sn, sd, rn, rd;
} sqrt_rational, qsqrt_t[1];

typedef sqrt_rational *qsqrt_ptr;
typedef const sqrt_rational *qsqrt_srcptr;

void qsqrt_init(qsqrt_ptr x);
#define qsqrt_sgn(x) mpz_sgn(x->sn)
#define qsqrt_is_zero(x) (mpz_sgn(x->sn) == 0)

void qsqrt_set_ui(qsqrt_ptr x, unsigned long n);
void qsqrt_set_ui4(qsqrt_ptr x, unsigned long sn, unsigned long sd, unsigned long rn, unsigned long rd);
void qsqrt_set_si(qsqrt_ptr x, long n);
void qsqrt_set_si4(qsqrt_ptr x, long sn, unsigned long sd, unsigned long rn, unsigned long rd);
void qsqrt_set(qsqrt_ptr x, const qsqrt_srcptr y);

void qsqrt_init_set_ui(qsqrt_ptr x, unsigned long n);
void qsqrt_init_set_ui4(qsqrt_ptr x, unsigned long sn, unsigned long sd, unsigned long rn, unsigned long rd);
void qsqrt_init_set_si(qsqrt_ptr x, long n);
void qsqrt_init_set_si4(qsqrt_ptr x, long sn, unsigned long sd, unsigned long rn, unsigned long rd);
void qsqrt_init_set(qsqrt_ptr x, const qsqrt_srcptr y);

_Bool qsqrt_eq(qsqrt_srcptr x, qsqrt_srcptr y);

void qsqrt_clear(qsqrt_ptr x);
void qsqrt_print(qsqrt_srcptr x);
double qsqrt_get_d(qsqrt_srcptr x);
// the result is in the form of "sn/sd√(rn/rd)"
// the `str` should be large enough
// if `str` is NULL, the function will allocate a new string
// the caller is responsible for freeing the string if it is allocated
char *qsqrt_get_str(qsqrt_srcptr x, char *str);

void qsqrt_simplify(qsqrt_ptr x, unsigned long hint);

_Bool check_CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
_Bool check_CG0(int j1, int j2, int j3);
_Bool check_3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
_Bool check_6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
_Bool check_9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9);
_Bool check_Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda);

// assume `ans` is initialized
// if the result of the coefficient should be zero,
// else return `hint`, and `ans` is set to the coefficient
// `hint` means the maximum prime factor of `rn` and `rd` is `<= hint`
int exact_CG(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
int exact_3j(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
int exact_CG0(qsqrt_ptr ans, int j1, int j2, int j3);
int exact_3j0(qsqrt_ptr ans, int j1, int j2, int j3);
int exact_6j(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
int exact_Racah(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
int exact_9j(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9);
int exact_norm9j(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9);
int exact_Moshinsky(qsqrt_ptr ans, int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda);

double ef_CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
double ef_3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
double ef_CG0(int j1, int j2, int j3);
double ef_3j0(int j1, int j2, int j3);
double ef_6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
double ef_Racah(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
double ef_9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9);
double ef_norm9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9);
double ef_Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda);

#ifdef __cplusplus
}
#endif

#endif // EXACTWIGNER_H