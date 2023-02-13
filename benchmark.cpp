#include "WignerSymbol.hpp"
#include <chrono>
#include <gsl/gsl_specfunc.h>

using namespace util;

void time_3j()
{
    using timer_clock = std::chrono::high_resolution_clock;
    auto t1 = timer_clock::now();
    WignerSymbols wigner;
    int N = 40;
    wigner.reserve(N, "2bjmax", 3);
    double x = 0;
    for (int dj1 = 0; dj1 <= N; ++dj1)
    {
        for (int dj2 = 0; dj2 <= N; ++dj2)
        {
            for (int dj3 = std::abs(dj1 - dj2); dj3 <= dj1 + dj2; dj3 += 2)
            {
                for (int dm1 = -dj1; dm1 <= dj1; dm1 += 2)
                {
                    for (int dm2 = -dj2; dm2 <= dj2; dm2 += 2)
                    {
                        for (int dm3 = -dj3; dm3 <= dj3; dm3 += 2)
                        {
                            x += wigner.f3j(dj1, dj2, dj3, dm1, dm2, dm3);
                        }
                    }
                }
            }
        }
    }
    auto t2 = timer_clock::now();
    double y = 0;
    for (int dj1 = 0; dj1 <= N; ++dj1)
    {
        for (int dj2 = 0; dj2 <= N; ++dj2)
        {
            for (int dj3 = std::abs(dj1 - dj2); dj3 <= dj1 + dj2; dj3 += 2)
            {
                for (int dm1 = -dj1; dm1 <= dj1; dm1 += 2)
                {
                    for (int dm2 = -dj2; dm2 <= dj2; dm2 += 2)
                    {
                        for (int dm3 = -dj3; dm3 <= dj3; dm3 += 2)
                        {
                            y += gsl_sf_coupling_3j(dj1, dj2, dj3, dm1, dm2, dm3);
                        }
                    }
                }
            }
        }
    }
    auto t3 = timer_clock::now();
    std::cout << "time 3j, diff = " << x - y << std::endl;
    std::cout << "this code time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms"
              << std::endl;
    std::cout << "gsl library time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
              << " ms" << std::endl;
}

// in time_3j, most Wigner 3j symbol is zero
// here, we test the performance of Wigner 3j symbol which is always non-zero
void time_3j_always_valid()
{
    using timer_clock = std::chrono::high_resolution_clock;
    auto t1 = timer_clock::now();
    WignerSymbols wigner;
    int N = 50;
    wigner.reserve(N, "2bjmax", 3);
    double x = 0;
    for (int dj1 = 0; dj1 <= N; ++dj1)
    {
        for (int dj2 = 0; dj2 <= N; ++dj2)
        {
            for (int dj3 = std::abs(dj1 - dj2); dj3 <= dj1 + dj2; dj3 += 2)
            {
                for (int dm1 = -dj1; dm1 <= dj1; dm1 += 2)
                {
                    for (int dm2 = -dj2; dm2 <= dj2; dm2 += 2)
                    {
                        x += wigner.f3j(dj1, dj2, dj3, dm1, dm2, -dm1 - dm2);
                    }
                }
            }
        }
    }
    auto t2 = timer_clock::now();
    double y = 0;
    for (int dj1 = 0; dj1 <= N; ++dj1)
    {
        for (int dj2 = 0; dj2 <= N; ++dj2)
        {
            for (int dj3 = std::abs(dj1 - dj2); dj3 <= dj1 + dj2; dj3 += 2)
            {
                for (int dm1 = -dj1; dm1 <= dj1; dm1 += 2)
                {
                    for (int dm2 = -dj2; dm2 <= dj2; dm2 += 2)
                    {
                        y += gsl_sf_coupling_3j(dj1, dj2, dj3, dm1, dm2, -dm1 - dm2);
                    }
                }
            }
        }
    }
    auto t3 = timer_clock::now();
    std::cout << "time 3j, diff = " << x - y << std::endl;
    std::cout << "this code time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms"
              << std::endl;
    std::cout << "gsl library time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
              << " ms" << std::endl;
}

void time_6j()
{
    using timer_clock = std::chrono::high_resolution_clock;
    auto t1 = timer_clock::now();
    WignerSymbols wigner;
    int N = 30;
    wigner.reserve(N, "Jmax", 6);
    double x = 0;
    for (int dj1 = N; dj1 <= 2 * N; ++dj1)
    {
        for (int dj2 = N; dj2 <= 2 * N; ++dj2)
        {
            for (int dj3 = N; dj3 <= 2 * N; ++dj3)
            {
                for (int dj4 = N; dj4 <= 2 * N; ++dj4)
                {
                    for (int dj5 = N; dj5 <= 2 * N; ++dj5)
                    {
                        for (int dj6 = N; dj6 <= 2 * N; ++dj6)
                        {
                            x += wigner.f6j(dj1, dj2, dj3, dj4, dj5, dj6);
                        }
                    }
                }
            }
        }
    }
    auto t2 = timer_clock::now();
    double y = 0;
    for (int dj1 = N; dj1 <= 2 * N; ++dj1)
    {
        for (int dj2 = N; dj2 <= 2 * N; ++dj2)
        {
            for (int dj3 = N; dj3 <= 2 * N; ++dj3)
            {
                for (int dj4 = N; dj4 <= 2 * N; ++dj4)
                {
                    for (int dj5 = N; dj5 <= 2 * N; ++dj5)
                    {
                        for (int dj6 = N; dj6 <= 2 * N; ++dj6)
                        {
                            y += gsl_sf_coupling_6j(dj1, dj2, dj3, dj4, dj5, dj6);
                        }
                    }
                }
            }
        }
    }
    auto t3 = timer_clock::now();
    std::cout << "time 6j, diff = " << x - y << std::endl;
    std::cout << "this code time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms"
              << std::endl;
    std::cout << "gsl library time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
              << " ms" << std::endl;
}

void time_6j_always_valid()
{
    using timer_clock = std::chrono::high_resolution_clock;
    auto t1 = timer_clock::now();
    WignerSymbols wigner;
    int N = 50;
    wigner.reserve(N, "2bjmax", 6);
    double x = 0;
    for (int dj1 = 0; dj1 <= N; ++dj1)
    {
        for (int dj2 = 0; dj2 <= N; ++dj2)
        {
            for (int dj4 = 0; dj4 <= N; ++dj4)
            {
                for (int dj5 = wigner.isodd(dj1 + dj2 + dj4); dj5 <= N; dj5 += 2)
                {
                    int dj3_min = std::max(std::abs(dj1 - dj2), std::abs(dj4 - dj5));
                    int dj3_max = std::min(dj1 + dj2, dj4 + dj5);
                    int dj6_min = std::max(std::abs(dj1 - dj5), std::abs(dj2 - dj4));
                    int dj6_max = std::min(dj1 + dj5, dj2 + dj4);
                    for (int dj3 = dj3_min; dj3 <= dj3_max; dj3 += 2)
                    {
                        for (int dj6 = dj6_min; dj6 <= dj6_max; dj6 += 2)
                        {
                            x += wigner.f6j(dj1, dj2, dj3, dj4, dj5, dj6);
                        }
                    }
                }
            }
        }
    }
    auto t2 = timer_clock::now();
    double y = 0;
    for (int dj1 = 0; dj1 <= N; ++dj1)
    {
        for (int dj2 = 0; dj2 <= N; ++dj2)
        {
            for (int dj4 = 0; dj4 <= N; ++dj4)
            {
                for (int dj5 = wigner.isodd(dj1 + dj2 + dj4); dj5 <= N; dj5 += 2)
                {
                    int dj3_min = std::max(std::abs(dj1 - dj2), std::abs(dj4 - dj5));
                    int dj3_max = std::min(dj1 + dj2, dj4 + dj5);
                    int dj6_min = std::max(std::abs(dj1 - dj5), std::abs(dj2 - dj4));
                    int dj6_max = std::min(dj1 + dj5, dj2 + dj4);
                    for (int dj3 = dj3_min; dj3 <= dj3_max; dj3 += 2)
                    {
                        for (int dj6 = dj6_min; dj6 <= dj6_max; dj6 += 2)
                        {
                            y += gsl_sf_coupling_6j(dj1, dj2, dj3, dj4, dj5, dj6);
                        }
                    }
                }
            }
        }
    }
    auto t3 = timer_clock::now();
    std::cout << "time 6j, diff = " << x - y << std::endl;
    std::cout << "this code time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms"
              << std::endl;
    std::cout << "gsl library time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
              << " ms" << std::endl;
}

void time_9j()
{
    using timer_clock = std::chrono::high_resolution_clock;
    auto t1 = timer_clock::now();
    WignerSymbols wigner;
    int N = 8;
    wigner.reserve(N, "Jmax", 9);
    double x = 0;
    for (int dj1 = N; dj1 <= 2 * N; ++dj1)
    {
        for (int dj2 = N; dj2 <= 2 * N; ++dj2)
        {
            for (int dj3 = N; dj3 <= 2 * N; ++dj3)
            {
                for (int dj4 = N; dj4 <= 2 * N; ++dj4)
                {
                    for (int dj5 = N; dj5 <= 2 * N; ++dj5)
                    {
                        for (int dj6 = N; dj6 <= 2 * N; ++dj6)
                        {
                            for (int dj7 = N; dj7 <= 2 * N; ++dj7)
                            {
                                for (int dj8 = N; dj8 <= 2 * N; ++dj8)
                                {
                                    for (int dj9 = N; dj9 <= 2 * N; ++dj9)
                                    {
                                        x += wigner.f9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    auto t2 = timer_clock::now();
    double y = 0;
    for (int dj1 = N; dj1 <= 2 * N; ++dj1)
    {
        for (int dj2 = N; dj2 <= 2 * N; ++dj2)
        {
            for (int dj3 = N; dj3 <= 2 * N; ++dj3)
            {
                for (int dj4 = N; dj4 <= 2 * N; ++dj4)
                {
                    for (int dj5 = N; dj5 <= 2 * N; ++dj5)
                    {
                        for (int dj6 = N; dj6 <= 2 * N; ++dj6)
                        {
                            for (int dj7 = N; dj7 <= 2 * N; ++dj7)
                            {
                                for (int dj8 = N; dj8 <= 2 * N; ++dj8)
                                {
                                    for (int dj9 = N; dj9 <= 2 * N; ++dj9)
                                    {
                                        y += gsl_sf_coupling_9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    auto t3 = timer_clock::now();
    std::cout << "time 9j, diff = " << x - y << std::endl;
    std::cout << "this code time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms"
              << std::endl;
    std::cout << "gsl library time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
              << " ms" << std::endl;
}

void time_9j_always_valid()
{
    using timer_clock = std::chrono::high_resolution_clock;
    auto t1 = timer_clock::now();
    WignerSymbols wigner;
    int N = 10;
    wigner.reserve(N, "Jmax", 9);
    double x = 0;
    for (int dj1 = 0; dj1 <= N; ++dj1)
    {
        for (int dj2 = 0; dj2 <= N; ++dj2)
        {
            for (int dj4 = 0; dj4 <= N; ++dj4)
            {
                for (int dj5 = 0; dj5 <= N; ++dj5)
                {
                    for (int dj3 = std::abs(dj1 - dj2); dj3 <= dj1 + dj2; dj3 += 2)
                    {
                        for (int dj6 = std::abs(dj4 - dj5); dj6 <= dj4 + dj5; dj6 += 2)
                        {
                            for (int dj7 = std::abs(dj1 - dj4); dj7 <= dj1 + dj4; dj7 += 2)
                            {
                                for (int dj8 = std::abs(dj2 - dj5); dj8 <= dj2 + dj5; dj8 += 2)
                                {
                                    int dj9_min = std::max(std::abs(dj3 - dj6), std::abs(dj7 - dj8));
                                    for (int dj9 = dj9_min; dj9 <= 2 * N; dj9 += 2)
                                    {
                                        x += wigner.f9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    auto t2 = timer_clock::now();
    double y = 0;
    for (int dj1 = 0; dj1 <= N; ++dj1)
    {
        for (int dj2 = 0; dj2 <= N; ++dj2)
        {
            for (int dj4 = 0; dj4 <= N; ++dj4)
            {
                for (int dj5 = 0; dj5 <= N; ++dj5)
                {
                    for (int dj3 = std::abs(dj1 - dj2); dj3 <= dj1 + dj2; dj3 += 2)
                    {
                        for (int dj6 = std::abs(dj4 - dj5); dj6 <= dj4 + dj5; dj6 += 2)
                        {
                            for (int dj7 = std::abs(dj1 - dj4); dj7 <= dj1 + dj4; dj7 += 2)
                            {
                                for (int dj8 = std::abs(dj2 - dj5); dj8 <= dj2 + dj5; dj8 += 2)
                                {
                                    int dj9_min = std::max(std::abs(dj3 - dj6), std::abs(dj7 - dj8));
                                    for (int dj9 = dj9_min; dj9 <= 2 * N; dj9 += 2)
                                    {
                                        y += gsl_sf_coupling_9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    auto t3 = timer_clock::now();
    std::cout << "time 9j, diff = " << x - y << std::endl;
    std::cout << "this code time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms"
              << std::endl;
    std::cout << "gsl library time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
              << " ms" << std::endl;
}

int main()
{
    std::cout << "----- test where most results are zeros -----" << std::endl;
    time_3j();
    time_6j();
    time_9j();
    std::cout << "----- test where arguements are always valid -----" << std::endl;
    time_3j_always_valid();
    time_6j_always_valid();
    time_9j_always_valid();
    return 0;
}
