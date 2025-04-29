#include "exactWigner.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <threads.h>

static inline int maxint(int a, int b) { return (a > b) ? a : b; }
static inline int minint(int a, int b) { return (a < b) ? a : b; }
static inline _Bool isodd(int x) { return (x % 2 != 0); }
static inline _Bool iseven(int x) { return (x % 2 == 0); }
static inline int iphase(int x) { return isodd(x) ? -1 : 1; }

typedef struct _prime_table
{
    unsigned long size;
    unsigned long *primes;
} prime_table_t;

// extend the prime table to include all primes up to n
static void _extend_primes_to(prime_table_t *table, unsigned long n)
{
    if (table->primes == NULL)
    {
        table->primes = (unsigned long *)malloc(sizeof(unsigned long) * 4);
        assert(table->primes != NULL);
        table->size = 4;
        table->primes[0] = 2;
        table->primes[1] = 3;
        table->primes[2] = 5;
        table->primes[3] = 7;
    }
    unsigned long size = table->size;
    unsigned long *primes = table->primes;
    if (n < primes[size - 1])
        return;
    while (primes[size - 1] < n)
    {
        unsigned long new_size = size * 2;
        unsigned long *new_primes = (unsigned long *)malloc(sizeof(unsigned long) * new_size);
        assert(new_primes != NULL);
        memcpy(new_primes, primes, sizeof(unsigned long) * size);
        free(primes);
        primes = new_primes;
        for (unsigned long i = size; i < new_size; ++i)
        {
            unsigned long next = primes[i - 1] + 2;
            while (1)
            {
                unsigned long j;
                for (j = 0; j < size; ++j)
                {
                    if (next % primes[j] == 0)
                        break;
                }
                if (j == size)
                    break;
                next += 2;
            }
            primes[i] = next;
        }
        size = new_size;
    }
    table->size = size;
    table->primes = primes;
}

// simplify x*\sqrt{t}, move square factors of `t` into `x`
// only move square factors of primes <= hint
// `q, tt` are temporary variables
static void simplify(mpz_ptr x, mpz_ptr t, unsigned long hint, mpz_ptr q, mpz_ptr tt)
{
    static thread_local prime_table_t table = {0, NULL};
    _extend_primes_to(&table, hint);
    mpz_set_ui(tt, 1);
    mp_bitcnt_t exp2 = mpz_scan1(t, 0);
    if (exp2 > 0)
    {
        unsigned long ei = exp2 / 2;
        unsigned long ti = exp2 % 2;
        if (ei > 0)
        {
            mpz_mul_2exp(x, x, ei);
        }
        if (ti == 1)
        {
            mpz_set_ui(tt, 2);
        }
        mpz_tdiv_q_2exp(t, t, exp2);
    }
    for (unsigned long i = 1; i < table.size; ++i)
    {
        unsigned long p = table.primes[i];
        if (p > hint)
            break;
        if (mpz_cmp_ui(t, p) < 0)
            break;
        unsigned long r = mpz_tdiv_q_ui(q, t, p);
        unsigned long e = 0;
        while (r == 0)
        {
            ++e;
            mpz_set(t, q);
            r = mpz_tdiv_q_ui(q, t, p);
        }
        if (e > 0)
        {
            unsigned long xi = e / 2;
            unsigned long ti = e % 2;
            if (xi > 0)
            {
                mpz_ui_pow_ui(q, p, xi);
                mpz_mul(x, x, q);
            }
            if (ti == 1)
            {
                mpz_mul_ui(tt, tt, p);
            }
        }
    }
    mpz_mul(t, t, tt);
}

static mpz_ptr bin(mpz_ptr ans, unsigned long n, unsigned long k)
{
    mpz_bin_uiui(ans, n, k);
    return ans;
}

static mpz_ptr omega(mpz_ptr ans, mpz_ptr temp, int j1, int j3, int g)
{
    bin(ans, g, j3);
    mpz_mul(ans, ans, bin(temp, j3, g - j1));
    return ans;
}

static mpz_ptr divgcd(mpz_ptr g, mpz_ptr n, mpz_ptr d)
{
    mpz_gcd(g, n, d);
    if (mpz_cmp_ui(g, 1) == 0)
        return g;
    mpz_divexact(n, n, g);
    mpz_divexact(d, d, g);
    return g;
}

// t is buffer, do not use mpq_add
static void q_add(mpz_ptr t, mpz_ptr an, mpz_ptr ad, mpz_srcptr bn, mpz_srcptr bd)
{
    mpz_mul(an, an, bd);
    mpz_addmul(an, ad, bn);
    mpz_mul(ad, ad, bd);
    divgcd(t, an, ad);
}

// simplify sn*\sqrt(1/rd) -> sn/sd*\sqrt(rn/rd)
static void simplify2(mpz_ptr sn, mpz_ptr sd, mpz_ptr rn, mpz_ptr rd, unsigned long hint)
{
    mpz_t q, tt;
    mpz_init(q);
    mpz_init(tt);
    divgcd(rn, sn, rd);
    divgcd(sd, rn, rd);
    mpz_set_ui(sd, 1);
    simplify(sn, rn, hint, q, tt);
    simplify(sd, rd, hint, q, tt);
    mpz_clear(q);
    mpz_clear(tt);
}

// simplify sn*\sqrt(rn/rd) -> sn/sd*\sqrt(rn/rd)
static void simplify3(mpz_ptr sn, mpz_ptr sd, mpz_ptr rn, mpz_ptr rd, unsigned long hint)
{
    mpz_t q, tt;
    mpz_init(q);
    mpz_init(tt);
    divgcd(sd, rn, rd);
    divgcd(sd, sn, rd);
    mpz_mul(rn, rn, sd);
    divgcd(sd, rn, rd);
    mpz_set_ui(sd, 1);
    simplify(sn, rn, hint, q, tt);
    simplify(sd, rd, hint, q, tt);
    mpz_clear(q);
    mpz_clear(tt);
}

// simplify sn/sd*\sqrt(rn/rd) -> sn/sd*\sqrt(rn/rd), t is buffer
static void simplify4(mpz_ptr t, mpz_ptr sn, mpz_ptr sd, mpz_ptr rn, mpz_ptr rd, unsigned long hint)
{
    mpz_t q;
    mpz_init(q);
    divgcd(t, sn, sd);
    divgcd(t, rn, rd);
    simplify(sn, rn, hint, q, t);
    simplify(sd, rd, hint, q, t);
    divgcd(t, sn, sd);
    divgcd(t, sn, rd);
    mpz_mul(rn, rn, t);
    divgcd(t, sd, rn);
    mpz_mul(rd, rd, t);
    divgcd(t, rn, rd);
    mpz_clear(q);
}

void qsqrt_simplify(qsqrt_ptr x, unsigned long hint)
{
    mpz_t t;
    mpz_init(t);
    simplify4(t, x->sn, x->sd, x->rn, x->rd, hint);
    mpz_clear(t);
}

void qsqrt_init(qsqrt_ptr x)
{
    mpz_init(x->sn);
    mpz_init(x->sd);
    mpz_init(x->rn);
    mpz_init(x->rd);
}

void qsqrt_set_ui(qsqrt_ptr x, unsigned long n)
{
    mpz_set_ui(x->sn, n);
    mpz_set_ui(x->sd, 1);
    mpz_set_ui(x->rn, 1);
    mpz_set_ui(x->rd, 1);
}

void qsqrt_set_ui4(qsqrt_ptr x, unsigned long sn, unsigned long sd, unsigned long rn, unsigned long rd)
{
    mpz_set_ui(x->sn, sn);
    mpz_set_ui(x->sd, sd);
    mpz_set_ui(x->rn, rn);
    mpz_set_ui(x->rd, rd);
}

void qsqrt_set_si(qsqrt_ptr x, long n)
{
    mpz_set_si(x->sn, n);
    mpz_set_ui(x->sd, 1);
    mpz_set_ui(x->rn, 1);
    mpz_set_ui(x->rd, 1);
}

void qsqrt_set_si4(qsqrt_ptr x, long sn, unsigned long sd, unsigned long rn, unsigned long rd)
{
    mpz_set_si(x->sn, sn);
    mpz_set_ui(x->sd, sd);
    mpz_set_ui(x->rn, rn);
    mpz_set_ui(x->rd, rd);
}

void qsqrt_set(qsqrt_ptr x, const qsqrt_srcptr y)
{
    mpz_set(x->sn, y->sn);
    mpz_set(x->sd, y->sd);
    mpz_set(x->rn, y->rn);
    mpz_set(x->rd, y->rd);
}

void qsqrt_init_set_ui(qsqrt_ptr x, unsigned long n)
{
    mpz_init_set_ui(x->sn, n);
    mpz_init_set_ui(x->sd, 1);
    mpz_init_set_ui(x->rn, 1);
    mpz_init_set_ui(x->rd, 1);
}

void qsqrt_init_set_ui4(qsqrt_ptr x, unsigned long sn, unsigned long sd, unsigned long rn, unsigned long rd)
{
    mpz_init_set_ui(x->sn, sn);
    mpz_init_set_ui(x->sd, sd);
    mpz_init_set_ui(x->rn, rn);
    mpz_init_set_ui(x->rd, rd);
}

void qsqrt_init_set_si(qsqrt_ptr x, long n)
{
    mpz_init_set_si(x->sn, n);
    mpz_init_set_ui(x->sd, 1);
    mpz_init_set_ui(x->rn, 1);
    mpz_init_set_ui(x->rd, 1);
}

void qsqrt_init_set_si4(qsqrt_ptr x, long sn, unsigned long sd, unsigned long rn, unsigned long rd)
{
    mpz_init_set_si(x->sn, sn);
    mpz_init_set_ui(x->sd, sd);
    mpz_init_set_ui(x->rn, rn);
    mpz_init_set_ui(x->rd, rd);
}

void qsqrt_init_set(qsqrt_ptr x, const qsqrt_srcptr y)
{
    mpz_init_set(x->sn, y->sn);
    mpz_init_set(x->sd, y->sd);
    mpz_init_set(x->rn, y->rn);
    mpz_init_set(x->rd, y->rd);
}

_Bool qsqrt_eq(qsqrt_srcptr x, qsqrt_srcptr y)
{
    if (mpz_cmp(x->sn, y->sn) != 0)
        return 0;
    if (mpz_cmp(x->sd, y->sd) != 0)
        return 0;
    if (mpz_cmp(x->rn, y->rn) != 0)
        return 0;
    if (mpz_cmp(x->rd, y->rd) != 0)
        return 0;
    return 1;
}

void qsqrt_clear(qsqrt_ptr x)
{
    mpz_clear(x->sn);
    mpz_clear(x->sd);
    mpz_clear(x->rn);
    mpz_clear(x->rd);
}

void qsqrt_print(qsqrt_srcptr x) { gmp_printf("%Zd/%Zd√(%Zd/%Zd)\n", x->sn, x->sd, x->rn, x->rd); }

double qsqrt_get_d(qsqrt_srcptr x)
{
    signed long esn = 0, esd = 0, ern = 0, erd = 0;
    double sn = mpz_get_d_2exp(&esn, x->sn);
    double sd = mpz_get_d_2exp(&esd, x->sd);
    double rn = mpz_get_d_2exp(&ern, x->rn);
    double rd = mpz_get_d_2exp(&erd, x->rd);
    double ret = sn / sd * sqrt(rn / rd);
    long exp = 2 * (esn - esd) + (ern - erd);
    long oexp = exp & 1;
    exp = (exp - oexp) / 2;
    ret = ldexp(ret, exp);
    if (oexp != 0)
        ret *= sqrt(2);
    return ret;
}

// for 9j
static double qsqrt_get_d2(qsqrt_srcptr x)
{
    signed long esn = 0, erd = 0;
    double sn = mpz_get_d_2exp(&esn, x->sn);
    double rd = mpz_get_d_2exp(&erd, x->rd);
    double ret = sn / sqrt(rd);
    long exp = 2 * esn - erd;
    long oexp = exp & 1;
    exp = (exp - oexp) / 2;
    ret = ldexp(ret, exp);
    if (oexp != 0)
        ret *= sqrt(2);
    return ret;
}

// for 3j, 6j
static double qsqrt_get_d3(qsqrt_srcptr x)
{
    signed long esn = 0, ern = 0, erd = 0;
    double sn = mpz_get_d_2exp(&esn, x->sn);
    double rn = mpz_get_d_2exp(&ern, x->rn);
    double rd = mpz_get_d_2exp(&erd, x->rd);
    double ret = sn * sqrt(rn / rd);
    long exp = 2 * esn + ern - erd;
    long oexp = exp & 1;
    exp = (exp - oexp) / 2;
    ret = ldexp(ret, exp);
    if (oexp != 0)
        ret *= sqrt(2);
    return ret;
}

char *qsqrt_get_str(qsqrt_srcptr x, char *str)
{
    int len = 0;
    if (str == NULL)
    {
        len = mpz_sizeinbase(x->sn, 10) + mpz_sizeinbase(x->sd, 10) + mpz_sizeinbase(x->rn, 10) +
              mpz_sizeinbase(x->rd, 10) + 10;
        str = (char *)malloc(len);
        assert(str != NULL);
    }
    char *p = str;
    mpz_get_str(p, 10, x->sn);
    p += strlen(p);
    *p++ = '/';
    mpz_get_str(p, 10, x->sd);
    p += strlen(p);
    static const char temp[] = "√(";
    memcpy(p, temp, sizeof(temp) - 1);
    p += sizeof(temp) - 1;
    mpz_get_str(p, 10, x->rn);
    p += strlen(p);
    *p++ = '/';
    mpz_get_str(p, 10, x->rd);
    p += strlen(p);
    *p++ = ')';
    *p = '\0';
    return str;
}

// assume dj1, dj2, dj3 are all non-negative
static inline _Bool check_couple(int dj1, int dj2, int dj3)
{
    if (isodd(dj1 + dj2 + dj3))
        return 0;
    if (dj1 > dj2 + dj3 || dj2 > dj1 + dj3 || dj3 > dj1 + dj2)
        return 0;
    return 1;
}

_Bool check_CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    if (dj1 < 0 || dj2 < 0 || dj3 < 0)
        return 0;
    if (abs(dm1) > dj1 || abs(dm2) > dj2 || abs(dm3) > dj3)
        return 0;
    if (isodd(dj1 + dm1) || isodd(dj2 + dm2) || isodd(dj3 + dm3))
        return 0;
    if (!check_couple(dj1, dj2, dj3))
        return 0;
    if (dm1 + dm2 != dm3)
        return 0;
    return 1;
}

_Bool check_CG0(int j1, int j2, int j3)
{
    if (j1 < 0 || j2 < 0 || j3 < 0)
        return 0;
    if (isodd(j1 + j2 + j3))
        return 0;
    if (j1 > j2 + j3 || j2 > j1 + j3 || j3 > j1 + j2)
        return 0;
    return 1;
}

_Bool check_3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    if (dj1 < 0 || dj2 < 0 || dj3 < 0)
        return 0;
    if (abs(dm1) > dj1 || abs(dm2) > dj2 || abs(dm3) > dj3)
        return 0;
    if (isodd(dj1 + dm1) || isodd(dj2 + dm2) || isodd(dj3 + dm3))
        return 0;
    if (!check_couple(dj1, dj2, dj3))
        return 0;
    if (dm1 + dm2 + dm3 != 0)
        return 0;
    return 1;
}

_Bool check_6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6)
{
    if (dj1 < 0 || dj2 < 0 || dj3 < 0 || dj4 < 0 || dj5 < 0 || dj6 < 0)
        return 0;
    if (!(check_couple(dj1, dj2, dj3) && check_couple(dj1, dj5, dj6) && check_couple(dj2, dj4, dj6) &&
          check_couple(dj3, dj4, dj5)))
        return 0;
    return 1;
}

_Bool check_9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9)
{
    if (dj1 < 0 || dj2 < 0 || dj3 < 0 || dj4 < 0 || dj5 < 0 || dj6 < 0 || dj7 < 0 || dj8 < 0 || dj9 < 0)
        return 0;
    if (!(check_couple(dj1, dj2, dj3) && check_couple(dj4, dj5, dj6) && check_couple(dj7, dj8, dj9) &&
          check_couple(dj1, dj4, dj7) && check_couple(dj2, dj5, dj8) && check_couple(dj3, dj6, dj9)))
        return 0;
    return 1;
}

_Bool check_Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda)
{
    if (N < 0 || L < 0 || n < 0 || l < 0 || n1 < 0 || l1 < 0 || n2 < 0 || l2 < 0 || lambda < 0)
        return 0;
    if (lambda > L + l || L > lambda + l || l > lambda + L)
        return 0;
    if (lambda > l1 + l2 || l1 > lambda + l2 || l2 > lambda + l1)
        return 0;
    const int F = 2 * N + L;
    const int f = 2 * n + l;
    const int f1 = 2 * n1 + l1;
    const int f2 = 2 * n2 + l2;
    if (F + f != f1 + f2)
        return 0;
    return 1;
}

// assume `ans` is initialized
static int impl_CG(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    const int J = (dj1 + dj2 + dj3) / 2;
    const int jm1 = J - dj1;
    const int jm2 = J - dj2;
    const int jm3 = J - dj3;
    const int j1mm1 = (dj1 - dm1) / 2;
    const int j2mm2 = (dj2 - dm2) / 2;
    const int j3mm3 = (dj3 - dm3) / 2;
    const int j2pm2 = (dj2 + dm2) / 2;
    const int low = maxint(0, maxint(j1mm1 - jm2, j2pm2 - jm1));
    const int high = minint(jm3, minint(j1mm1, j2pm2));
    mpz_ptr sum = ans->sn;
    mpz_ptr t = ans->sd;
    mpz_ptr tz = ans->rd;
    mpz_set_ui(sum, 0);
    for (int z = low; z <= high; ++z)
    {
        bin(tz, jm3, z);
        mpz_mul(tz, tz, bin(t, jm2, j1mm1 - z));
        mpz_mul(tz, tz, bin(t, jm1, j2pm2 - z));
        mpz_sub(sum, tz, sum);
    }
    if (mpz_sgn(sum) == 0)
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    if (isodd(high))
    {
        mpz_neg(sum, sum);
    }
    mpz_ptr rn = ans->rn;
    mpz_ptr rd = ans->rd;
    bin(rn, dj1, jm2);
    mpz_mul(rn, rn, bin(t, dj2, jm3));
    bin(rd, J + 1, jm3);
    mpz_mul(rd, rd, bin(t, dj1, j1mm1));
    mpz_mul(rd, rd, bin(t, dj2, j2mm2));
    mpz_mul(rd, rd, bin(t, dj3, j3mm3));
    return J + 1;
}

static int impl_CG0(qsqrt_ptr ans, int j1, int j2, int j3)
{
    const int J = j1 + j2 + j3;
    const int g = J / 2;
    mpz_ptr sn = ans->sn;
    mpz_ptr t = ans->sd;
    bin(sn, g, j3);
    mpz_mul(sn, sn, bin(t, j3, g - j1));
    if (isodd(g - j3))
    {
        mpz_neg(sn, sn);
    }
    mpz_ptr rd = ans->rd;
    bin(rd, J + 1, 2 * j3 + 1);
    mpz_mul(rd, rd, bin(t, 2 * j3, J - 2 * j1));
    return J + 1;
}

static int impl_3j0(qsqrt_ptr ans, int j1, int j2, int j3)
{
    const int J = j1 + j2 + j3;
    const int g = J / 2;
    mpz_ptr sn = ans->sn;
    mpz_ptr t = ans->sd;
    bin(sn, g, j3);
    mpz_mul(sn, sn, bin(t, j3, g - j1));
    if (isodd(g))
    {
        mpz_neg(sn, sn);
    }
    mpz_ptr rd = ans->rd;
    mpz_set_ui(rd, (unsigned long)(2 * j3 + 1));
    mpz_mul(rd, rd, bin(t, J + 1, 2 * j3 + 1));
    mpz_mul(rd, rd, bin(t, 2 * j3, J - 2 * j1));
    return J + 1;
}

// assume `ans` is initialized
static int impl_3j(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    const int J = (dj1 + dj2 + dj3) / 2;
    const int jm1 = J - dj1;
    const int jm2 = J - dj2;
    const int jm3 = J - dj3;
    const int j1mm1 = (dj1 - dm1) / 2;
    const int j2mm2 = (dj2 - dm2) / 2;
    const int j3mm3 = (dj3 - dm3) / 2;
    const int j1pm1 = (dj1 + dm1) / 2;

    mpz_ptr sum = ans->sn;
    mpz_ptr t = ans->sd;
    mpz_ptr tz = ans->rd;
    mpz_set_ui(sum, 0);
    const int low = maxint(0, maxint(j1pm1 - jm2, j2mm2 - jm1));
    const int high = minint(jm3, minint(j1pm1, j2mm2));
    for (int z = low; z <= high; ++z)
    {
        bin(tz, jm3, z);
        mpz_mul(tz, tz, bin(t, jm2, j1pm1 - z));
        mpz_mul(tz, tz, bin(t, jm1, j2mm2 - z));
        mpz_sub(sum, tz, sum);
    }
    if (mpz_sgn(sum) == 0)
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    if (isodd(dj1 + (dj3 + dm3) / 2 + high))
    {
        mpz_neg(sum, sum);
    }
    mpz_ptr rn = ans->rn;
    mpz_ptr rd = ans->rd;
    bin(rn, dj1, jm2);
    mpz_mul(rn, rn, bin(t, dj2, jm3));
    mpz_set_ui(rd, J + 1);
    mpz_mul(rd, rd, bin(t, J, jm3));
    mpz_mul(rd, rd, bin(t, dj1, j1mm1));
    mpz_mul(rd, rd, bin(t, dj2, j2mm2));
    mpz_mul(rd, rd, bin(t, dj3, j3mm3));
    return J + 1;
}

static int impl_6j(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dj4, int dj5, int dj6)
{
    const int j123 = (dj1 + dj2 + dj3) / 2;
    const int j156 = (dj1 + dj5 + dj6) / 2;
    const int j426 = (dj4 + dj2 + dj6) / 2;
    const int j453 = (dj4 + dj5 + dj3) / 2;
    const int jpm123 = (dj1 + dj2 - dj3) / 2;
    const int jpm132 = (dj1 + dj3 - dj2) / 2;
    const int jpm231 = (dj2 + dj3 - dj1) / 2;
    const int jpm156 = (dj1 + dj5 - dj6) / 2;
    const int jpm426 = (dj4 + dj2 - dj6) / 2;
    const int jpm453 = (dj4 + dj5 - dj3) / 2;

    const int low = maxint(maxint(j123, j453), maxint(j426, j156));
    const int high = minint(jpm123 + j453, minint(jpm132 + j426, jpm231 + j156));
    mpz_ptr sum = ans->sn;
    mpz_ptr t = ans->sd;
    mpz_ptr tx = ans->rd;
    mpz_set_ui(sum, 0);
    for (int x = low; x <= high; ++x)
    {
        bin(tx, x + 1, j123 + 1);
        mpz_mul(tx, tx, bin(t, jpm123, x - j453));
        mpz_mul(tx, tx, bin(t, jpm132, x - j426));
        mpz_mul(tx, tx, bin(t, jpm231, x - j156));
        mpz_sub(sum, tx, sum);
    }
    if (mpz_sgn(sum) == 0)
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    if (isodd(high))
    {
        mpz_neg(sum, sum);
    }
    mpz_ptr rn = ans->rn;
    mpz_ptr rd = ans->rd;
    bin(rn, j123 + 1, dj1 + 1);
    mpz_mul(rn, rn, bin(t, dj1, jpm123));
    bin(rd, j156 + 1, dj1 + 1);
    mpz_mul(rd, rd, bin(t, dj1, jpm156));
    mpz_mul(rd, rd, bin(t, j453 + 1, dj4 + 1));
    mpz_mul(rd, rd, bin(t, dj4, jpm453));
    mpz_mul(rd, rd, bin(t, j426 + 1, dj4 + 1));
    mpz_mul(rd, rd, bin(t, dj4, jpm426));
    mpz_mul_ui(rd, rd, (unsigned long)(dj4 + 1));
    mpz_mul_ui(rd, rd, (unsigned long)(dj4 + 1));
    return high + 1;
}

// not simplfied
static int impl_norm_9j(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9)
{
    const int j123 = (dj1 + dj2 + dj3) / 2;
    const int j456 = (dj4 + dj5 + dj6) / 2;
    const int j789 = (dj7 + dj8 + dj9) / 2;
    const int j147 = (dj1 + dj4 + dj7) / 2;
    const int j258 = (dj2 + dj5 + dj8) / 2;
    const int j369 = (dj3 + dj6 + dj9) / 2;
    const int pm123 = (dj1 + dj2 - dj3) / 2;
    const int pm132 = (dj1 + dj3 - dj2) / 2;
    const int pm231 = (dj2 + dj3 - dj1) / 2;
    const int pm456 = (dj4 + dj5 - dj6) / 2;
    const int pm465 = (dj4 + dj6 - dj5) / 2;
    const int pm564 = (dj5 + dj6 - dj4) / 2;
    const int pm789 = (dj7 + dj8 - dj9) / 2;
    const int pm798 = (dj7 + dj9 - dj8) / 2;
    const int pm897 = (dj8 + dj9 - dj7) / 2;
    const int dtl = maxint(abs(dj2 - dj6), maxint(abs(dj4 - dj8), abs(dj1 - dj9)));
    const int dth = minint(dj2 + dj6, minint(dj4 + dj8, dj1 + dj9));
    mpz_ptr sum = ans->sn;
    mpz_ptr Pt = ans->sd;
    mpz_ptr ABC = ans->rn;
    mpz_ptr tx = ans->rd;
    mpz_set_ui(sum, 0);
    mpz_t t;
    mpz_init(t);
    for (int dt = dtl; dt <= dth; dt += 2)
    {
        const int j19t = (dj1 + dj9 + dt) / 2;
        const int j26t = (dj2 + dj6 + dt) / 2;
        const int j48t = (dj4 + dj8 + dt) / 2;
        mpz_set_ui(ABC, dt + 1);
        const int xl = maxint(maxint(j123, j369), maxint(j26t, j19t));
        const int xh = minint(pm123 + j369, minint(pm132 + j26t, pm231 + j19t));
        mpz_set_ui(Pt, 0);
        for (int x = xl; x <= xh; ++x)
        {
            bin(tx, x + 1, j19t + 1);
            mpz_mul(tx, tx, bin(t, (dj1 + dj9 - dt) / 2, x - j26t));
            mpz_mul(tx, tx, bin(t, (dj1 + dt - dj9) / 2, x - j369));
            mpz_mul(tx, tx, bin(t, (dt + dj9 - dj1) / 2, x - j123));
            mpz_sub(Pt, tx, Pt);
        }
        mpz_mul(ABC, ABC, Pt);
        const int yl = maxint(maxint(j456, j258), maxint(j26t, j48t));
        const int yh = minint(pm456 + j26t, minint(pm465 + j258, pm564 + j48t));
        mpz_set_ui(Pt, 0);
        for (int y = yl; y <= yh; ++y)
        {
            bin(tx, y + 1, j26t + 1);
            mpz_mul(tx, tx, bin(t, (dj2 + dj6 - dt) / 2, y - j48t));
            mpz_mul(tx, tx, bin(t, (dt + dj6 - dj2) / 2, y - j258));
            mpz_mul(tx, tx, bin(t, (dt + dj2 - dj6) / 2, y - j456));
            mpz_sub(Pt, tx, Pt);
        }
        mpz_mul(ABC, ABC, Pt);
        const int zl = maxint(maxint(j789, j147), maxint(j48t, j19t));
        const int zh = minint(pm789 + j19t, minint(pm798 + j48t, pm897 + j147));
        mpz_set_ui(Pt, 0);
        for (int z = zl; z <= zh; ++z)
        {
            bin(tx, z + 1, j48t + 1);
            mpz_mul(tx, tx, bin(t, (dj4 + dj8 - dt) / 2, z - j19t));
            mpz_mul(tx, tx, bin(t, (dt + dj8 - dj4) / 2, z - j147));
            mpz_mul(tx, tx, bin(t, (dt + dj4 - dj8) / 2, z - j789));
            mpz_sub(Pt, tx, Pt);
        }
        mpz_mul(ABC, ABC, Pt);
        if (isodd(xh + yh + zh))
        {
            mpz_neg(ABC, ABC);
        }
        mpz_add(sum, sum, ABC);
    }
    if (mpz_sgn(sum) == 0)
    {
        qsqrt_set_ui(ans, 0);
        mpz_clear(t);
        return 0;
    }
    if (isodd(dth))
    {
        mpz_neg(sum, sum);
    }
    mpz_ptr P0 = ans->rd;
    mpz_set_ui(P0, dj9 + 1);
    mpz_pow_ui(P0, P0, 2);
    mpz_mul(P0, P0, bin(t, j123 + 1, dj3 + 1));
    mpz_mul(P0, P0, bin(t, dj3, pm231));
    mpz_mul(P0, P0, bin(t, j456 + 1, dj6 + 1));
    mpz_mul(P0, P0, bin(t, dj6, pm564));
    mpz_mul(P0, P0, bin(t, j789 + 1, dj9 + 1));
    mpz_mul(P0, P0, bin(t, dj9, pm897));
    mpz_mul(P0, P0, bin(t, j147 + 1, dj7 + 1));
    mpz_mul(P0, P0, bin(t, dj7, (dj4 + dj7 - dj1) / 2));
    mpz_mul(P0, P0, bin(t, j258 + 1, dj8 + 1));
    mpz_mul(P0, P0, bin(t, dj8, (dj5 + dj8 - dj2) / 2));
    mpz_mul(P0, P0, bin(t, j369 + 1, dj9 + 1));
    mpz_mul(P0, P0, bin(t, dj9, (dj6 + dj9 - dj3) / 2));
    const int maxJ = maxint(maxint(j123, j456), maxint(maxint(j789, j147), maxint(j258, j369)));
    mpz_clear(t);
    return maxJ + 1;
}

static void _m9j(mpz_ptr sum, mpz_ptr temp, mpz_ptr tx, mpz_ptr Pt, mpz_ptr ABC, int j1, int j2, int j3, int j4, int j5,
                 int j6, int j7, int j8, int j9)
{
    const int j123 = j1 + j2 + j3;
    const int j456 = j4 + j5 + j6;
    const int j789 = j7 + j8 + j9;
    const int j147 = j1 + j4 + j7;
    const int j258 = j2 + j5 + j8;
    const int j369 = j3 + j6 + j9;
    const int pm123 = j1 + j2 - j3;
    const int pm132 = j1 + j3 - j2;
    const int pm231 = j2 + j3 - j1;
    const int pm456 = j4 + j5 - j6;
    const int pm465 = j4 + j6 - j5;
    const int pm564 = j5 + j6 - j4;
    const int pm789 = j7 + j8 - j9;
    const int pm798 = j7 + j9 - j8;
    const int pm897 = j8 + j9 - j7;
    mpz_set_ui(sum, 0);
    const int tl = maxint(abs(j2 - j6), maxint(abs(j4 - j8), abs(j1 - j9)));
    const int th = minint(j2 + j6, minint(j4 + j8, j1 + j9));
    for (int t = tl; t <= th; ++t)
    {
        const int j19t = j1 + j9 + t;
        const int j26t = j2 + j6 + t;
        const int j48t = j4 + j8 + t;
        mpz_set_ui(ABC, (unsigned long)(2 * t + 1));
        const int xl = maxint(maxint(j123, j369), maxint(j26t, j19t));
        const int xh = minint(pm123 + j369, minint(pm132 + j26t, pm231 + j19t));
        mpz_set_ui(Pt, 0);
        for (int x = xl; x <= xh; ++x)
        {
            bin(tx, x + 1, j19t + 1);
            mpz_mul(tx, tx, bin(temp, j1 + j9 - t, x - j26t));
            mpz_mul(tx, tx, bin(temp, j1 + t - j9, x - j369));
            mpz_mul(tx, tx, bin(temp, t + j9 - j1, x - j123));
            mpz_sub(Pt, tx, Pt);
        }
        mpz_mul(ABC, ABC, Pt);
        const int yl = maxint(maxint(j456, j258), maxint(j26t, j48t));
        const int yh = minint(pm456 + j26t, minint(pm465 + j258, pm564 + j48t));
        mpz_set_ui(Pt, 0);
        for (int y = yl; y <= yh; ++y)
        {
            bin(tx, y + 1, j26t + 1);
            mpz_mul(tx, tx, bin(temp, j2 + j6 - t, y - j48t));
            mpz_mul(tx, tx, bin(temp, t + j6 - j2, y - j258));
            mpz_mul(tx, tx, bin(temp, t + j2 - j6, y - j456));
            mpz_sub(Pt, tx, Pt);
        }
        mpz_mul(ABC, ABC, Pt);
        const int zl = maxint(maxint(j789, j147), maxint(j48t, j19t));
        const int zh = minint(pm789 + j19t, minint(pm798 + j48t, pm897 + j147));
        mpz_set_ui(Pt, 0);
        for (int z = zl; z <= zh; ++z)
        {
            bin(tx, z + 1, j48t + 1);
            mpz_mul(tx, tx, bin(temp, j4 + j8 - t, z - j19t));
            mpz_mul(tx, tx, bin(temp, t + j8 - j4, z - j147));
            mpz_mul(tx, tx, bin(temp, t + j4 - j8, z - j789));
            mpz_sub(Pt, tx, Pt);
        }
        mpz_mul(ABC, ABC, Pt);
        if (isodd(xh + yh + zh))
        {
            mpz_neg(ABC, ABC);
        }
        mpz_add(sum, sum, ABC);
    }
}

static int impl_Moshinsky(qsqrt_ptr ans, int Ncom, int Lcom, int nrel, int lrel, int n1, int l1, int n2, int l2,
                          int lam)
{
    const int F = 2 * Ncom + Lcom;
    const int f = 2 * nrel + lrel;
    const int f1 = 2 * n1 + l1;
    const int f2 = 2 * n2 + l2;
    const int chi = f1 + f2;
    const int NL = Ncom + Lcom;
    const int nl = nrel + lrel;
    const int nl1 = n1 + l1;
    const int nl2 = n2 + l2;

    mpz_ptr sum_n = ans->sn;
    mpz_ptr sum_d = ans->sd;
    mpz_ptr Rn = ans->rn;
    mpz_ptr Rd = ans->rd;

    mpz_t t, tx, Pt, ABC, M9j;
    mpz_init(t);
    mpz_init(tx);
    mpz_init(Pt);
    mpz_init(ABC);
    mpz_init(M9j);
    mpz_t FA, An, Ad, Bn, Bd, Cn, Cd, Dn, Dd;
    mpz_init(FA);
    mpz_init(An);
    mpz_init(Ad);
    mpz_init(Bn);
    mpz_init(Bd);
    mpz_init(Cn);
    mpz_init(Cd);
    mpz_init(Dn);
    mpz_init(Dd);

    bin(Rn, chi + 2, f1 + 1);
    bin(Rd, chi + 2, F + 1);
    mpz_mul(Rd, Rd, bin(t, Lcom + lrel + lam + 1, 2 * lam + 1));
    mpz_mul(Rd, Rd, bin(t, 2 * lam, lam + Lcom - lrel));
    mpz_mul(Rd, Rd, bin(t, l1 + l2 + lam + 1, 2 * lam + 1));
    mpz_mul(Rd, Rd, bin(t, 2 * lam, lam + l1 - l2));
    divgcd(t, Rn, Rd);

    mpz_mul(Rn, Rn, bin(t, 2 * NL + 1, NL));
    mpz_mul_ui(Rd, Rd, 2 * Lcom + 1);
    mpz_mul(Rd, Rd, bin(t, F + 1, Ncom));

    mpz_mul(Rn, Rn, bin(t, 2 * nl + 1, nl));
    mpz_mul_ui(Rd, Rd, 2 * lrel + 1);
    mpz_mul(Rd, Rd, bin(t, f + 1, nrel));

    mpz_mul(Rn, Rn, bin(t, 2 * nl1 + 1, nl1));
    mpz_mul_ui(Rd, Rd, 2 * l1 + 1);
    mpz_mul(Rd, Rd, bin(t, f1 + 1, n1));

    mpz_mul(Rn, Rn, bin(t, 2 * nl2 + 1, nl2));
    mpz_mul_ui(Rd, Rd, 2 * l2 + 1);
    mpz_mul(Rd, Rd, bin(t, f2 + 1, n2));
    divgcd(t, Rn, Rd);

    mpz_set_ui(t, (unsigned long)(f1 + 2));
    mpz_mul_ui(t, t, (unsigned long)(f2 + 2));
    mpz_mul_ui(t, t, (unsigned long)(2 * lam + 1));
    mpz_pow_ui(t, t, 2);
    mpz_mul_2exp(t, t, chi + (Lcom + lrel + l1 + l2));
    mpz_mul(Rd, Rd, t);
    divgcd(t, Rn, Rd);

    mpz_set_ui(sum_n, 0);
    mpz_set_ui(sum_d, 1);
    for (int fa = 0; fa <= minint(f1, F); ++fa)
    {
        const int fb = f1 - fa;
        const int fc = F - fa;
        const int fd = f2 - fc;
        if (fd < 0)
            continue;
        bin(FA, f1 + 2, fa + 1);
        mpz_mul(FA, FA, bin(t, f2 + 2, fc + 1));
        for (int la = fa & 1; la <= fa; la += 2)
        {
            const int na = (fa - la) / 2;
            const int nla = na + la;
            mpz_mul_ui(An, FA, 2 * la + 1);
            mpz_mul(An, An, bin(t, fa + 1, na));
            mpz_mul_2exp(An, An, la);
            bin(Ad, 2 * nla + 1, nla);
            divgcd(t, An, Ad);
            for (int lb = abs(l1 - la); lb <= minint(l1 + la, fb); lb += 2)
            {
                const int nb = (fb - lb) / 2;
                const int nlb = nb + lb;
                mpz_mul_ui(Bn, An, 2 * lb + 1);
                mpz_mul(Bn, Bn, bin(t, fb + 1, nb));
                mpz_mul_2exp(Bn, Bn, lb);
                mpz_mul(Bd, Ad, bin(t, 2 * nlb + 1, nlb));
                mpz_mul(Bn, Bn, omega(tx, t, la, l1, (la + lb + l1) / 2));
                // Δ(lalbl1)
                mpz_mul(Bd, Bd, bin(t, la + lb + l1 + 1, 2 * l1 + 1));
                mpz_mul(Bd, Bd, bin(t, 2 * l1, l1 + la - lb));
                divgcd(t, Bn, Bd);
                for (int lc = abs(Lcom - la); lc <= minint(Lcom + la, fc); lc += 2)
                {
                    const int nc = (fc - lc) / 2;
                    const int nlc = nc + lc;
                    mpz_mul_ui(Cn, Bn, 2 * lc + 1);
                    mpz_mul(Cn, Cn, bin(t, fc + 1, nc));
                    mpz_mul_2exp(Cn, Cn, lc);
                    mpz_mul(Cd, Bd, bin(t, 2 * nlc + 1, nlc));
                    mpz_mul(Cn, Cn, omega(tx, t, la, Lcom, (la + lc + Lcom) / 2));
                    // Δ(lalcl3)
                    mpz_mul(Cd, Cd, bin(t, la + lc + Lcom + 1, 2 * Lcom + 1));
                    mpz_mul(Cd, Cd, bin(t, 2 * Lcom, Lcom + la - lc));
                    divgcd(t, Cn, Cd);
                    const int ldmin = maxint(abs(l2 - lc), abs(lrel - lb));
                    const int ldmax = minint(fd, minint(l2 + lc, lrel + lb));
                    for (int ld = ldmin; ld <= ldmax; ld += 2)
                    {
                        const int nd = (fd - ld) / 2;
                        const int nld = nd + ld;
                        mpz_mul_ui(Dn, Cn, 2 * ld + 1);
                        mpz_mul(Dn, Dn, bin(t, fd + 1, nd));
                        mpz_mul_2exp(Dn, Dn, ld);
                        mpz_mul(Dd, Cd, bin(t, 2 * nld + 1, nld));
                        mpz_mul(Dn, Dn, omega(tx, t, lc, l2, (lc + ld + l2) / 2));
                        mpz_mul(Dn, Dn, omega(tx, t, lb, lrel, (lb + ld + lrel) / 2));
                        // Δ(lbldl4)
                        mpz_mul(Dd, Dd, bin(t, lb + ld + lrel + 1, 2 * lrel + 1));
                        mpz_mul(Dd, Dd, bin(t, 2 * lrel, lrel + lb - ld));
                        // Δ(lcldl2)
                        mpz_mul(Dd, Dd, bin(t, lc + ld + l2 + 1, 2 * l2 + 1));
                        mpz_mul(Dd, Dd, bin(t, 2 * l2, l2 + lc - ld));
                        divgcd(t, Dn, Dd);
                        _m9j(M9j, t, tx, Pt, ABC, la, lb, l1, lc, ld, l2, Lcom, lrel, lam);
                        mpz_mul(Dn, Dn, M9j);
                        divgcd(t, Dn, Dd);
                        if (isodd(ld))
                        {
                            mpz_neg(Dn, Dn);
                        }
                        q_add(t, sum_n, sum_d, Dn, Dd);
                    }
                }
            }
        }
    }
    if (mpz_sgn(sum_n) == 0)
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    int hint = maxint(chi + 2, maxint(Lcom + lrel, l1 + l2) + lam + 1);
    simplify4(t, sum_n, sum_d, Rn, Rd, hint);
    mpz_clear(t);
    mpz_clear(tx);
    mpz_clear(Pt);
    mpz_clear(ABC);
    mpz_clear(M9j);
    mpz_clear(FA);
    mpz_clear(An);
    mpz_clear(Ad);
    mpz_clear(Bn);
    mpz_clear(Bd);
    mpz_clear(Cn);
    mpz_clear(Cd);
    mpz_clear(Dn);
    mpz_clear(Dd);
    return hint;
}

static int gcdint(int a, int b)
{
    if (b == 0)
        return a;
    return gcdint(b, a % b);
}

static int divgcdint(int *a, int *b)
{
    const int g = gcdint(*a, *b);
    *a /= g;
    *b /= g;
    return g;
}

static int impl_Moshinsky_d(qsqrt_ptr ans, int Ncom, int Lcom, int nrel, int lrel, int n1, int l1, int n2, int l2,
                            int lam, int m1w1, int m2w2)
{
    const int F = 2 * Ncom + Lcom;
    const int f = 2 * nrel + lrel;
    const int f1 = 2 * n1 + l1;
    const int f2 = 2 * n2 + l2;
    const int chi = f1 + f2;
    const int NL = Ncom + Lcom;
    const int nl = nrel + lrel;
    const int nl1 = n1 + l1;
    const int nl2 = n2 + l2;
    divgcdint(&m1w1, &m2w2);

    mpz_ptr sum_n = ans->sn;
    mpz_ptr sum_d = ans->sd;
    mpz_ptr Rn = ans->rn;
    mpz_ptr Rd = ans->rd;

    mpz_t t, tx, Pt, ABC, M9j;
    mpz_init(t);
    mpz_init(tx);
    mpz_init(Pt);
    mpz_init(ABC);
    mpz_init(M9j);
    mpz_t FAn, FAd, An, Ad, Bn, Bd, Cn, Cd, Dn, Dd;
    mpz_init(FAn);
    mpz_init(FAd);
    mpz_init(An);
    mpz_init(Ad);
    mpz_init(Bn);
    mpz_init(Bd);
    mpz_init(Cn);
    mpz_init(Cd);
    mpz_init(Dn);
    mpz_init(Dd);

    bin(Rn, chi + 2, f1 + 1);
    bin(Rd, chi + 2, F + 1);
    mpz_mul(Rd, Rd, bin(t, Lcom + lrel + lam + 1, 2 * lam + 1));
    mpz_mul(Rd, Rd, bin(t, 2 * lam, lam + Lcom - lrel));
    mpz_mul(Rd, Rd, bin(t, l1 + l2 + lam + 1, 2 * lam + 1));
    mpz_mul(Rd, Rd, bin(t, 2 * lam, lam + l1 - l2));
    divgcd(t, Rn, Rd);

    mpz_mul(Rn, Rn, bin(t, 2 * NL + 1, NL));
    mpz_mul_ui(Rd, Rd, 2 * Lcom + 1);
    mpz_mul(Rd, Rd, bin(t, F + 1, Ncom));

    mpz_mul(Rn, Rn, bin(t, 2 * nl + 1, nl));
    mpz_mul_ui(Rd, Rd, 2 * lrel + 1);
    mpz_mul(Rd, Rd, bin(t, f + 1, nrel));

    mpz_mul(Rn, Rn, bin(t, 2 * nl1 + 1, nl1));
    mpz_mul_ui(Rd, Rd, 2 * l1 + 1);
    mpz_mul(Rd, Rd, bin(t, f1 + 1, n1));

    mpz_mul(Rn, Rn, bin(t, 2 * nl2 + 1, nl2));
    mpz_mul_ui(Rd, Rd, 2 * l2 + 1);
    mpz_mul(Rd, Rd, bin(t, f2 + 1, n2));
    divgcd(t, Rn, Rd);

    mpz_set_ui(t, (unsigned long)(f1 + 2));
    mpz_mul_ui(t, t, (unsigned long)(f2 + 2));
    mpz_mul_ui(t, t, (unsigned long)(2 * lam + 1));
    mpz_pow_ui(t, t, 2);
    mpz_mul_2exp(t, t, (Lcom + lrel + l1 + l2));
    mpz_mul(Rd, Rd, t);
    divgcd(t, Rn, Rd);

    if (f - f1 > 0)
    {
        mpz_ui_pow_ui(t, m1w1, f - f1);
        mpz_mul(Rn, Rn, t);
    }
    else if (f - f1 < 0)
    {
        mpz_ui_pow_ui(t, m1w1, f1 - f);
        mpz_mul(Rd, Rd, t);
    }
    mpz_ui_pow_ui(t, m2w2, f1 + F);
    mpz_mul(Rn, Rn, t);
    mpz_ui_pow_ui(t, m1w1 + m2w2, f + F);
    mpz_mul(Rd, Rd, t);
    divgcd(t, Rn, Rd);

    mpz_set_ui(sum_n, 0);
    mpz_set_ui(sum_d, 1);
    for (int fa = 0; fa <= minint(f1, F); ++fa)
    {
        const int fb = f1 - fa;
        const int fc = F - fa;
        const int fd = f2 - fc;
        if (fd < 0)
            continue;
        bin(FAn, f1 + 2, fa + 1);
        mpz_mul(FAn, FAn, bin(t, f2 + 2, fc + 1));
        mpz_ui_pow_ui(t, m1w1, fa);
        mpz_mul(FAn, FAn, t);
        mpz_ui_pow_ui(FAd, m2w2, fa);
        divgcd(t, FAn, FAd);
        for (int la = fa & 1; la <= fa; la += 2)
        {
            const int na = (fa - la) / 2;
            const int nla = na + la;
            mpz_mul_ui(An, FAn, 2 * la + 1);
            mpz_mul(An, An, bin(t, fa + 1, na));
            mpz_mul_2exp(An, An, la);
            mpz_mul(Ad, FAd, bin(t, 2 * nla + 1, nla));
            divgcd(t, An, Ad);
            for (int lb = abs(l1 - la); lb <= minint(l1 + la, fb); lb += 2)
            {
                const int nb = (fb - lb) / 2;
                const int nlb = nb + lb;
                mpz_mul_ui(Bn, An, 2 * lb + 1);
                mpz_mul(Bn, Bn, bin(t, fb + 1, nb));
                mpz_mul_2exp(Bn, Bn, lb);
                mpz_mul(Bd, Ad, bin(t, 2 * nlb + 1, nlb));
                mpz_mul(Bn, Bn, omega(tx, t, la, l1, (la + lb + l1) / 2));
                // Δ(lalbl1)
                mpz_mul(Bd, Bd, bin(t, la + lb + l1 + 1, 2 * l1 + 1));
                mpz_mul(Bd, Bd, bin(t, 2 * l1, l1 + la - lb));
                divgcd(t, Bn, Bd);
                for (int lc = abs(Lcom - la); lc <= minint(Lcom + la, fc); lc += 2)
                {
                    const int nc = (fc - lc) / 2;
                    const int nlc = nc + lc;
                    mpz_mul_ui(Cn, Bn, 2 * lc + 1);
                    mpz_mul(Cn, Cn, bin(t, fc + 1, nc));
                    mpz_mul_2exp(Cn, Cn, lc);
                    mpz_mul(Cd, Bd, bin(t, 2 * nlc + 1, nlc));
                    mpz_mul(Cn, Cn, omega(tx, t, la, Lcom, (la + lc + Lcom) / 2));
                    // Δ(lalcl3)
                    mpz_mul(Cd, Cd, bin(t, la + lc + Lcom + 1, 2 * Lcom + 1));
                    mpz_mul(Cd, Cd, bin(t, 2 * Lcom, Lcom + la - lc));
                    divgcd(t, Cn, Cd);
                    const int ldmin = maxint(abs(l2 - lc), abs(lrel - lb));
                    const int ldmax = minint(fd, minint(l2 + lc, lrel + lb));
                    for (int ld = ldmin; ld <= ldmax; ld += 2)
                    {
                        const int nd = (fd - ld) / 2;
                        const int nld = nd + ld;
                        mpz_mul_ui(Dn, Cn, 2 * ld + 1);
                        mpz_mul(Dn, Dn, bin(t, fd + 1, nd));
                        mpz_mul_2exp(Dn, Dn, ld);
                        mpz_mul(Dd, Cd, bin(t, 2 * nld + 1, nld));
                        mpz_mul(Dn, Dn, omega(tx, t, lc, l2, (lc + ld + l2) / 2));
                        mpz_mul(Dn, Dn, omega(tx, t, lb, lrel, (lb + ld + lrel) / 2));
                        // Δ(lbldl4)
                        mpz_mul(Dd, Dd, bin(t, lb + ld + lrel + 1, 2 * lrel + 1));
                        mpz_mul(Dd, Dd, bin(t, 2 * lrel, lrel + lb - ld));
                        // Δ(lcldl2)
                        mpz_mul(Dd, Dd, bin(t, lc + ld + l2 + 1, 2 * l2 + 1));
                        mpz_mul(Dd, Dd, bin(t, 2 * l2, l2 + lc - ld));
                        divgcd(t, Dn, Dd);
                        _m9j(M9j, t, tx, Pt, ABC, la, lb, l1, lc, ld, l2, Lcom, lrel, lam);
                        mpz_mul(Dn, Dn, M9j);
                        divgcd(t, Dn, Dd);
                        if (isodd(ld))
                        {
                            mpz_neg(Dn, Dn);
                        }
                        q_add(t, sum_n, sum_d, Dn, Dd);
                    }
                }
            }
        }
    }
    if (mpz_sgn(sum_n) == 0)
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    int hint = maxint(chi + 2, maxint(Lcom + lrel, l1 + l2) + lam + 1);
    simplify4(t, sum_n, sum_d, Rn, Rd, hint);
    mpz_clear(t);
    mpz_clear(tx);
    mpz_clear(Pt);
    mpz_clear(ABC);
    mpz_clear(M9j);
    mpz_clear(FAn);
    mpz_clear(FAd);
    mpz_clear(An);
    mpz_clear(Ad);
    mpz_clear(Bn);
    mpz_clear(Bd);
    mpz_clear(Cn);
    mpz_clear(Cd);
    mpz_clear(Dn);
    mpz_clear(Dd);
    return hint;
}

int exact_CG(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    if (!check_CG(dj1, dj2, dj3, dm1, dm2, dm3))
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    const int hint = impl_CG(ans, dj1, dj2, dj3, dm1, dm2, dm3);
    if (hint == 0)
    {
        return 0;
    }
    simplify3(ans->sn, ans->sd, ans->rn, ans->rd, hint);
    return hint;
}

int exact_3j(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    if (!check_3j(dj1, dj2, dj3, dm1, dm2, dm3))
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    const int hint = impl_3j(ans, dj1, dj2, dj3, dm1, dm2, dm3);
    if (hint == 0)
    {
        return 0;
    }
    simplify3(ans->sn, ans->sd, ans->rn, ans->rd, hint);
    return hint;
}

int exact_CG0(qsqrt_ptr ans, int j1, int j2, int j3)
{
    if (!check_CG0(j1, j2, j3))
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    const int hint = impl_CG0(ans, j1, j2, j3);
    simplify2(ans->sn, ans->sd, ans->rn, ans->rd, hint);
    return hint;
}

int exact_3j0(qsqrt_ptr ans, int j1, int j2, int j3)
{
    if (!check_CG0(j1, j2, j3))
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    const int hint = impl_3j0(ans, j1, j2, j3);
    simplify2(ans->sn, ans->sd, ans->rn, ans->rd, hint);
    return hint;
}

int exact_6j(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dj4, int dj5, int dj6)
{
    if (!check_6j(dj1, dj2, dj3, dj4, dj5, dj6))
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    int hint = impl_6j(ans, dj1, dj2, dj3, dj4, dj5, dj6);
    if (hint == 0)
    {
        return 0;
    }
    simplify3(ans->sn, ans->sd, ans->rn, ans->rd, hint);
    return hint;
}

int exact_Racah(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dj4, int dj5, int dj6)
{
    const int hint = exact_6j(ans, dj1, dj2, dj5, dj4, dj3, dj6);
    if (hint == 0)
    {
        return 0;
    }
    if (isodd((dj1 + dj2 + dj3 + dj4) / 2))
    {
        mpz_neg(ans->sn, ans->sn);
    }
    return hint;
}

int exact_9j(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9)
{
    if (!check_9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9))
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    const int hint = impl_norm_9j(ans, dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
    if (hint == 0)
    {
        return 0;
    }
    mpz_mul_ui(ans->rd, ans->rd, (unsigned long)(dj3 + 1));
    mpz_mul_ui(ans->rd, ans->rd, (unsigned long)(dj6 + 1));
    mpz_mul_ui(ans->rd, ans->rd, (unsigned long)(dj7 + 1));
    mpz_mul_ui(ans->rd, ans->rd, (unsigned long)(dj8 + 1));
    simplify2(ans->sn, ans->sd, ans->rn, ans->rd, hint);
    return hint;
}

int exact_norm9j(qsqrt_ptr ans, int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9)
{
    if (!check_9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9))
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    const int hint = impl_norm_9j(ans, dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
    if (hint == 0)
    {
        return 0;
    }
    simplify2(ans->sn, ans->sd, ans->rn, ans->rd, hint);
    return hint;
}

int exact_Moshinsky(qsqrt_ptr ans, int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda)
{
    if (!check_Moshinsky(N, L, n, l, n1, l1, n2, l2, lambda))
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    return impl_Moshinsky(ans, N, L, n, l, n1, l1, n2, l2, lambda);
}

int exact_Moshinsky_d(qsqrt_ptr ans, int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda, int m1w1,
                      int m2w2)
{
    if (!check_Moshinsky(N, L, n, l, n1, l1, n2, l2, lambda))
    {
        qsqrt_set_ui(ans, 0);
        return 0;
    }
    return impl_Moshinsky_d(ans, N, L, n, l, n1, l1, n2, l2, lambda, m1w1, m2w2);
}

double ef_CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    if (!check_CG(dj1, dj2, dj3, dm1, dm2, dm3))
    {
        return 0.0;
    }
    qsqrt_t ans;
    qsqrt_init(ans);
    const int hint = impl_CG(ans, dj1, dj2, dj3, dm1, dm2, dm3);
    if (hint == 0)
    {
        qsqrt_clear(ans);
        return 0.0;
    }
    double ret = qsqrt_get_d3(ans);
    qsqrt_clear(ans);
    return ret;
}

double ef_3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    if (!check_3j(dj1, dj2, dj3, dm1, dm2, dm3))
    {
        return 0.0;
    }
    qsqrt_t ans;
    qsqrt_init(ans);
    const int hint = impl_3j(ans, dj1, dj2, dj3, dm1, dm2, dm3);
    if (hint == 0)
    {
        qsqrt_clear(ans);
        return 0.0;
    }
    double ret = qsqrt_get_d3(ans);
    qsqrt_clear(ans);
    return ret;
}

double ef_CG0(int j1, int j2, int j3)
{
    if (!check_CG0(j1, j2, j3))
    {
        return 0.0;
    }
    qsqrt_t ans;
    qsqrt_init(ans);
    const int hint = impl_CG0(ans, j1, j2, j3);
    if (hint == 0)
    {
        qsqrt_clear(ans);
        return 0.0;
    }
    double ret = qsqrt_get_d2(ans);
    qsqrt_clear(ans);
    return ret;
}

double ef_3j0(int j1, int j2, int j3)
{
    if (!check_CG0(j1, j2, j3))
    {
        return 0.0;
    }
    qsqrt_t ans;
    qsqrt_init(ans);
    const int hint = impl_3j0(ans, j1, j2, j3);
    if (hint == 0)
    {
        qsqrt_clear(ans);
        return 0.0;
    }
    double ret = qsqrt_get_d2(ans);
    qsqrt_clear(ans);
    return ret;
}

double ef_6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6)
{
    if (!check_6j(dj1, dj2, dj3, dj4, dj5, dj6))
    {
        return 0.0;
    }
    qsqrt_t ans;
    qsqrt_init(ans);
    const int hint = impl_6j(ans, dj1, dj2, dj3, dj4, dj5, dj6);
    if (hint == 0)
    {
        qsqrt_clear(ans);
        return 0.0;
    }
    double ret = qsqrt_get_d3(ans);
    qsqrt_clear(ans);
    return ret;
}

double ef_Racah(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6)
{
    double ret = ef_6j(dj1, dj2, dj5, dj4, dj3, dj6);
    if (isodd((dj1 + dj2 + dj3 + dj4) / 2))
    {
        ret = -ret;
    }
    return ret;
}

double ef_9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9)
{
    double ret = ef_norm9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
    ret /= sqrt((double)(dj3 + 1) * (dj6 + 1) * (dj7 + 1) * (dj8 + 1));
    return ret;
}

double ef_norm9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9)
{
    if (!check_9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9))
    {
        return 0.0;
    }
    qsqrt_t ans;
    qsqrt_init(ans);
    const int hint = impl_norm_9j(ans, dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
    if (hint == 0)
    {
        qsqrt_clear(ans);
        return 0.0;
    }
    double ret = qsqrt_get_d2(ans);
    qsqrt_clear(ans);
    return ret;
}

double ef_Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda)
{
    if (!check_Moshinsky(N, L, n, l, n1, l1, n2, l2, lambda))
    {
        return 0.0;
    }
    qsqrt_t ans;
    qsqrt_init(ans);
    int hint = impl_Moshinsky(ans, N, L, n, l, n1, l1, n2, l2, lambda);
    if (hint == 0)
    {
        qsqrt_clear(ans);
        return 0.0;
    }
    double ret = qsqrt_get_d(ans);
    qsqrt_clear(ans);
    return ret;
}

double ef_Moshinsky_d(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda, int m1w1, int m2w2)
{
    if (!check_Moshinsky(N, L, n, l, n1, l1, n2, l2, lambda))
    {
        return 0.0;
    }
    qsqrt_t ans;
    qsqrt_init(ans);
    int hint = impl_Moshinsky_d(ans, N, L, n, l, n1, l1, n2, l2, lambda, m1w1, m2w2);
    if (hint == 0)
    {
        qsqrt_clear(ans);
        return 0.0;
    }
    double ret = qsqrt_get_d(ans);
    qsqrt_clear(ans);
    return ret;
}
