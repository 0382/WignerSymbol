/* Reference:
[1] D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, Quantum Theory of Angular Momentum, (World Scientific,
1988).
[2] Buck et al. Nuc. Phys. A 600 (1996) 387-402.
*/

#include "exactWigner.h"
#include <stdio.h>

void test_special_CG(int djmax);
void test_CG0(int jmax);
void test_special_6j(int djmax);
void test_special_Racah(int djmax);
void test_special_9j(int djmax);
void test_Moshinsky();

int main()
{
    test_special_CG(50);
    test_CG0(20);
    test_special_6j(20);
    test_special_Racah(20);
    test_special_9j(10);
    test_Moshinsky();
    return 0;
}

// Ref[1], P248, Sec 8.5, Formula(1)
void test_special_CG(int djmax)
{
    qsqrt_t cg, quick;
    qsqrt_init(cg);
    qsqrt_init(quick);
    int ok = 1;
    int count = 0;
    for (int dj = 0; dj <= djmax; ++dj)
    {
        for (int dm = -dj; dm <= dj; dm += 2)
        {
            int hint = exact_CG(cg, dj, dj, 0, dm, -dm, 0);
            int phase = ((dj - dm) / 2) % 2 == 0 ? 1 : -1;
            qsqrt_set_si4(quick, phase, 1, 1, dj + 1);
            qsqrt_simplify(quick, hint);
            if (!qsqrt_eq(cg, quick))
            {
                printf("test special CG failed for j=%d/2, m=%d/2\n", dj, dm);
                qsqrt_print(cg);
                qsqrt_print(quick);
                ok = 0;
            }
            ++count;
        }
    }
    if (ok)
    {
        printf("%d tests passed for special CG with Jmax=%d/2\n", count, djmax);
    }
    qsqrt_clear(cg);
    qsqrt_clear(quick);
}

void test_CG0(int jmax)
{
    qsqrt_t cg, cg0;
    qsqrt_init(cg);
    qsqrt_init(cg0);
    int ok = 1;
    int count = 0;
    for (int j1 = 0; j1 <= jmax; ++j1)
    {
        for (int j2 = 0; j2 <= jmax; ++j2)
        {
            for (int j3 = 0; j3 <= jmax; ++j3)
            {
                exact_CG(cg, 2 * j1, 2 * j2, 2 * j3, 0, 0, 0);
                exact_CG0(cg0, j1, j2, j3);
                if (!qsqrt_eq(cg, cg0))
                {
                    printf("test CG0 failed for j1=%d, j2=%d, j3=%d\n", j1, j2, j3);
                    qsqrt_print(cg);
                    qsqrt_print(cg0);
                    ok = 0;
                }
                exact_3j(cg, 2 * j1, 2 * j2, 2 * j3, 0, 0, 0);
                exact_3j0(cg0, j1, j2, j3);
                if (!qsqrt_eq(cg, cg0))
                {
                    printf("test 3j0 failed for j1=%d, j2=%d, j3=%d\n", j1, j2, j3);
                    qsqrt_print(cg);
                    qsqrt_print(cg0);
                    ok = 0;
                }
                count += 2;
            }
        }
    }
    if (ok)
    {
        printf("%d tests passed for CG0 with Jmax=%d\n", count, jmax);
    }
    qsqrt_clear(cg);
    qsqrt_clear(cg0);
}

// Ref[1], P299, Sec 9.5, Formula (1)
void test_special_6j(int djmax)
{
    qsqrt_t e6j, quick;
    qsqrt_init(e6j);
    qsqrt_init(quick);
    int ok = 1;
    int count = 0;
    for (int dj1 = 0; dj1 <= djmax; ++dj1)
    {
        for (int dj2 = 0; dj2 <= djmax; ++dj2)
        {
            for (int dj3 = 0; dj3 <= djmax; ++dj3)
            {
                int hint = exact_6j(e6j, dj1, dj2, dj3, dj2, dj1, 0);
                if (hint == 0)
                {
                    continue;
                }
                int phase = ((dj1 + dj2 + dj3) / 2) % 2 == 0 ? 1 : -1;
                qsqrt_set_si4(quick, phase, 1, 1, (dj1 + 1) * (dj2 + 1));
                qsqrt_simplify(quick, hint);
                if (!qsqrt_eq(e6j, quick))
                {
                    printf("test special 6j failed for j1=%d/2, j2=%d/2, j3=%d/2\n", dj1, dj2, dj3);
                    qsqrt_print(e6j);
                    qsqrt_print(quick);
                    ok = 0;
                }
                ++count;
            }
        }
    }
    if (ok)
    {
        printf("%d tests passed for special 6j with Jmax=%d/2\n", count, djmax);
    }
    qsqrt_clear(e6j);
    qsqrt_clear(quick);
}

// Ref[1], P300, Sec 9.5, Formula (2)
void test_special_Racah(int djmax)
{
    qsqrt_t racah, quick;
    qsqrt_init(racah);
    qsqrt_init(quick);
    int ok = 1;
    int count = 0;
    for (int dj1 = 0; dj1 <= djmax; ++dj1)
    {
        for (int dj2 = 0; dj2 <= djmax; ++dj2)
        {
            for (int dj3 = 0; dj3 <= djmax; ++dj3)
            {
                int hint = exact_Racah(racah, 0, dj1, dj2, dj3, dj1, dj2);
                if (hint == 0)
                {
                    continue;
                }
                qsqrt_set_ui4(quick, 1, 1, 1, (dj1 + 1) * (dj2 + 1));
                qsqrt_simplify(quick, hint);
                if (!qsqrt_eq(racah, quick))
                {
                    printf("test special Racah failed for j1=%d/2, j2=%d/2, j3=%d/2\n", dj1, dj2, dj3);
                    qsqrt_print(racah);
                    qsqrt_print(quick);
                    ok = 0;
                }
                ++count;
            }
        }
    }
    if (ok)
    {
        printf("%d tests passed for special Racah with Jmax=%d/2\n", count, djmax);
    }
    qsqrt_clear(racah);
    qsqrt_clear(quick);
}

// Ref[1], P357, Sec 10.9, Formula (1)
void test_special_9j(int djmax)
{
    qsqrt_t e9j, quick;
    qsqrt_init(e9j);
    qsqrt_init(quick);
    int ok = 1;
    int count = 0;
    for (int dj1 = 0; dj1 <= djmax; ++dj1)
    {
        for (int dj2 = 0; dj2 <= djmax; ++dj2)
        {
            for (int dj3 = 0; dj3 <= djmax; ++dj3)
            {
                for (int dj4 = 0; dj4 <= djmax; ++dj4)
                {
                    for (int dj5 = 0; dj5 <= djmax; ++dj5)
                    {
                        for (int dj7 = 0; dj7 <= djmax; ++dj7)
                        {
                            int hint = exact_9j(e9j, dj1, dj2, dj3, dj4, dj5, dj3, dj7, dj7, 0);
                            exact_6j(quick, dj1, dj2, dj3, dj5, dj4, dj7);
                            if (hint == 0)
                            {
                                if (!qsqrt_is_zero(quick))
                                {
                                    printf("test special 9j failed for j1=%d/2, j2=%d/2, j3=%d/2, j4=%d/2, j5=%d/2, "
                                           "j7=%d/2\n",
                                           dj1, dj2, dj3, dj4, dj5, dj7);
                                    qsqrt_print(e9j);
                                    qsqrt_print(quick);
                                    ok = 0;
                                }
                                continue;
                            }
                            int phase = ((dj2 + dj3 + dj4 + dj7) / 2) % 2 == 0 ? 1 : -1;
                            if (phase == -1)
                            {
                                mpz_neg(quick->sn, quick->sn);
                            }
                            mpz_mul_ui(quick->rd, quick->rd, (dj3 + 1) * (dj7 + 1));
                            qsqrt_simplify(quick, hint);
                            if (!qsqrt_eq(e9j, quick))
                            {
                                printf("test special 9j failed for j1=%d/2, j2=%d/2, j3=%d/2, j4=%d/2, j5=%d/2, "
                                       "j7=%d/2\n",
                                       dj1, dj2, dj3, dj4, dj5, dj7);
                                qsqrt_print(e9j);
                                qsqrt_print(quick);
                                ok = 0;
                            }
                            ++count;
                        }
                    }
                }
            }
        }
    }
    if (ok)
    {
        printf("%d tests passed for special 9j with Jmax=%d/2\n", count, djmax);
    }
    qsqrt_clear(e9j);
    qsqrt_clear(quick);
}

// Ref [2]
void test_Moshinsky()
{
    static const int test_set[9][9] = {
        {0, 2, 1, 0, 0, 1, 0, 3, 2}, {0, 1, 0, 5, 0, 1, 0, 5, 6}, {0, 1, 0, 3, 0, 2, 0, 2, 4},
        {1, 3, 0, 1, 0, 2, 0, 4, 3}, {0, 5, 0, 2, 0, 2, 0, 5, 4}, {0, 3, 1, 6, 2, 2, 1, 3, 4},
        {1, 0, 2, 5, 2, 2, 1, 3, 5}, {0, 2, 4, 2, 2, 2, 1, 4, 2}, {3, 2, 0, 4, 2, 2, 1, 4, 4}};
    static const int test_ans[9][4] = {{-1, 2, 7, 10},  {1, 2, 1, 1},      {0, 1, 1, 1},
                                       {-1, 2, 5, 14},  {-1, 6, 1, 2},     {1, 24, 65, 3},
                                       {-5, 96, 13, 7}, {-1, 28, 195, 14}, {4463, 25872, 1, 3}};
    qsqrt_t x, y;
    qsqrt_init(x);
    qsqrt_init(y);
    int ok = 1;
    for (int i = 0; i < 9; ++i)
    {
        int Ncom = test_set[i][0];
        int Lcom = test_set[i][1];
        int nrel = test_set[i][2];
        int lrel = test_set[i][3];
        int n1 = test_set[i][4];
        int l1 = test_set[i][5];
        int n2 = test_set[i][6];
        int l2 = test_set[i][7];
        int lambda = test_set[i][8];
        exact_Moshinsky(x, Ncom, Lcom, nrel, lrel, n1, l1, n2, l2, lambda);
        int sn = test_ans[i][0];
        int sd = test_ans[i][1];
        int rn = test_ans[i][2];
        int rd = test_ans[i][3];
        qsqrt_set_si4(y, sn, sd, rn, rd);
        if (!qsqrt_eq(x, y))
        {
            printf(
                "test Moshinsky failed for Ncom=%d, Lcom=%d, nrel=%d, lrel=%d, n1=%d, l1=%d, n2=%d, l2=%d, lambda=%d\n",
                Ncom, Lcom, nrel, lrel, n1, l1, n2, l2, lambda);
            qsqrt_print(x);
            qsqrt_print(y);
            ok = 0;
        }
    }
    if (ok)
    {
        printf("9 test Moshinsky passed\n");
    }
    qsqrt_clear(x);
    qsqrt_clear(y);
}
