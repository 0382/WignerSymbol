// test the orthogonality of the Moshinsky transformation

#include "WignerSymbol.hpp"
#include <chrono>
#include <iostream>

using namespace std;
using namespace util;

double test_orth(int Emax);
double test_orth2(int Emax);
double bench_mosh(int Emax);

int main()
{
    int Emax = 12;
    // double orth = bench_mosh(Emax);
    double orth = test_orth2(Emax);
    cout << "Orthogonality of Moshinsky transformation: " << orth << endl;
    return 0;
}

double bench_mosh(int Emax)
{
    double sum = 0.0;
    int count = 0;
    wigner_init(Emax, "Moshinsky", 0);
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int E = 0; E <= Emax; ++E)
    {
        const int e = Emax - E;
        for (int e1 = 0; e1 <= Emax; ++e1)
        {
            const int e2 = E + e - e1;
            for (int L = E & 1; L <= E; L += 2)
            {
                const int N = (E - L) / 2;
                for (int l = e & 1; l <= e; l += 2)
                {
                    const int n = (e - l) / 2;
                    for (int l1 = e1 & 1; l1 <= e1; l1 += 2)
                    {
                        const int n1 = (e1 - l1) / 2;
                        for (int l2 = e2 & 1; l2 <= e2; l2 += 2)
                        {
                            const int n2 = (e2 - l2) / 2;
                            const int Lam_max = std::min(L + l, l1 + l2);
                            const int Lam_min = std::max(std::abs(L - l), std::abs(l1 - l2));
                            for (int Lam = Lam_min; Lam <= Lam_max; ++Lam)
                            {
                                sum += Moshinsky(N, L, n, l, n1, l1, n2, l2, Lam);
                                ++count;
                            }
                        }
                    }
                }
            }
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Number of counts: " << count << std::endl;
    auto tms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "Time: " << tms << " ms, " << sum / count * 1e6 << " us per call" << std::endl;
    return sum;
}

double test_orth(int Emax)
{
    int count = 0;
    double delta = 0.0;
    wigner_init(2 * Emax, "Moshinsky", 0);
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int E = 0; E <= Emax; ++E)
    {
        for (int e = 0; e <= Emax; ++e)
        {
            for (int Ep = 0; Ep <= Emax; ++Ep)
            {
                int ep = E + e - Ep;
                if (ep < 0 || ep > Emax)
                    continue;
                for (int L = E & 1; L <= E; L += 2)
                {
                    int N = (E - L) / 2;
                    for (int l = e & 1; l <= e; l += 2)
                    {
                        int n = (e - l) / 2;
                        for (int Lp = Ep & 1; Lp <= Ep; Lp += 2)
                        {
                            int Np = (Ep - Lp) / 2;
                            for (int lp = ep & 1; lp <= ep; lp += 2)
                            {
                                int np = (ep - lp) / 2;
                                double diag = (E == Ep && e == ep && L == Lp && l == lp) ? 1.0 : 0.0;
                                int Lam_max = std::min(L + l, Lp + lp);
                                int Lam_min = std::max(std::abs(L - l), std::abs(Lp - lp));
                                for (int Lam = Lam_min; Lam <= Lam_max; ++Lam)
                                {
                                    double sum = 0.0;
                                    for (int e1 = 0; e1 <= E + e; ++e1)
                                    {
                                        int e2 = E + e - e1;
                                        for (int l1 = e1 & 1; l1 <= e1; ++l1)
                                        {
                                            int n1 = (e1 - l1) / 2;
                                            int l2_min = std::max(e2 & 1, abs(Lam - l1));
                                            int l2_max = std::min(e2, Lam + l1);
                                            for (int l2 = l2_min; l2 <= l2_max; ++l2)
                                            {
                                                int n2 = (e2 - l2) / 2;
                                                double mosh = Moshinsky(N, L, n, l, n1, l1, n2, l2, Lam);
                                                double moshp = Moshinsky(Np, Lp, np, lp, n1, l1, n2, l2, Lam);
                                                sum += mosh * moshp;
                                                ++count;
                                            }
                                        }
                                    }
                                    delta += std::abs(sum - diag);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Number of tests: " << count << std::endl;
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms"
              << std::endl;
    return delta;
}

double test_orth2(int Emax)
{
    int count = 0;
    double delta = 0.0;
    wigner_init(2 * Emax, "Moshinsky", 0);
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int E = 0; E <= Emax; ++E)
    {
        int e = Emax - E;
        for (int Ep = 0; Ep <= Emax; ++Ep)
        {
            int ep = Emax - Ep;
            for (int L = E & 1; L <= E; L += 2)
            {
                int N = (E - L) / 2;
                for (int l = e & 1; l <= e; l += 2)
                {
                    int n = (e - l) / 2;
                    for (int Lp = Ep & 1; Lp <= Ep; Lp += 2)
                    {
                        int Np = (Ep - Lp) / 2;
                        for (int lp = ep & 1; lp <= ep; lp += 2)
                        {
                            int np = (ep - lp) / 2;
                            double diag = (E == Ep && e == ep && L == Lp && l == lp) ? 1.0 : 0.0;
                            int Lam_max = std::min(L + l, Lp + lp);
                            int Lam_min = std::max(std::abs(L - l), std::abs(Lp - lp));
                            for (int Lam = Lam_min; Lam <= Lam_max; ++Lam)
                            {
                                double sum = 0.0;
                                for (int e1 = 0; e1 <= E + e; ++e1)
                                {
                                    int e2 = E + e - e1;
                                    for (int l1 = e1 & 1; l1 <= e1; ++l1)
                                    {
                                        int n1 = (e1 - l1) / 2;
                                        int l2_min = std::max(e2 & 1, abs(Lam - l1));
                                        int l2_max = std::min(e2, Lam + l1);
                                        for (int l2 = l2_min; l2 <= l2_max; ++l2)
                                        {
                                            int n2 = (e2 - l2) / 2;
                                            double mosh = Moshinsky(N, L, n, l, n1, l1, n2, l2, Lam);
                                            double moshp = Moshinsky(Np, Lp, np, lp, n1, l1, n2, l2, Lam);
                                            sum += mosh * moshp;
                                            ++count;
                                        }
                                    }
                                }
                                delta += std::abs(sum - diag);
                            }
                        }
                    }
                }
            }
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Number of tests: " << count << std::endl;
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms"
              << std::endl;
    return delta;
}