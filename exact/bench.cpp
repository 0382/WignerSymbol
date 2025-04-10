#include "exactWigner.h"
#include <chrono>
#include <iostream>

using pfun6_t = double (*)(int, int, int, int, int, int);
using pfun9_t = double (*)(int, int, int, int, int, int, int, int, int);

double bench_3j(int N, pfun6_t func, const char *name)
{
    auto t1 = std::chrono::high_resolution_clock::now();
    double sum = 0;
    for (int dj1 = 1; dj1 <= N; ++dj1)
    {
        std::cout << dj1 << '\r' << std::flush;
        for (int dj2 = 0; dj2 <= dj1; ++dj2)
        {
            for (int dj3 = std::abs(dj1 - dj2); dj3 <= dj2; dj3 += 2)
            {
                for (int dm1 = -dj1; dm1 <= 0; dm1 += 2)
                {
                    for (int dm2 = -dj2; dm2 <= dj2; dm2 += 2)
                    {
                        const int dm3 = -dm1 - dm2;
                        double w = func(dj1, dj2, dj3, dm1, dm2, dm3);
                        sum += w;
                    }
                }
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << name << " time: " << duration << " ms\n";
    return sum;
}

double bench_6j(int N, pfun6_t func, const char *name)
{
    auto t1 = std::chrono::high_resolution_clock::now();
    double sum = 0;
    for (int dj1 = 1; dj1 <= N; ++dj1)
    {
        std::cout << dj1 << '\r' << std::flush;
        for (int dj2 = 0; dj2 <= dj1; ++dj2)
        {
            for (int dj3 = std::abs(dj1 - dj2); dj3 <= dj2; dj3 += 2)
            {
                for (int dj4 = 0; dj4 <= dj1; ++dj4)
                {
                    for (int dj5 = std::abs(dj3 - dj4); dj5 <= std::min(dj1, dj3 + dj4); dj5 += 2)
                    {
                        const int dj6_min = std::max(std::abs(dj1 - dj5), std::abs(dj2 - dj4));
                        const int dj6_max = std::min(dj1, std::abs(dj2 + dj4));
                        for (int dj6 = dj6_min; dj6 <= dj6_max; dj6 += 2)
                        {
                            double w = func(dj1, dj2, dj3, dj4, dj5, dj6);
                            sum += w;
                        }
                    }
                }
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << name << " time: " << duration << " ms\n";
    return sum;
}

double bench_9j(int N, pfun9_t func, const char *name)
{
    auto t1 = std::chrono::high_resolution_clock::now();
    double sum = 0;
    for (int dj1 = 1; dj1 <= N; ++dj1)
    {
        std::cout << dj1 << '\r' << std::flush;
        for (int dj2 = 0; dj2 <= dj1; ++dj2)
        {
            for (int dj3 = 0; dj3 <= dj2; ++dj3)
            {
                for (int dj4 = 0; dj4 <= dj1; ++dj4)
                {
                    for (int dj12 = std::abs(dj1 - dj2); dj12 <= dj2; dj12 += 2)
                    {
                        for (int dj34 = std::abs(dj3 - dj4); dj34 <= std::min(dj1, dj3 + dj4); dj34 += 2)
                        {
                            for (int dj13 = std::abs(dj1 - dj3); dj13 <= dj3; dj13 += 2)
                            {
                                for (int dj24 = std::abs(dj2 - dj4); dj24 <= std::min(dj1, dj2 + dj4); dj24 += 2)
                                {
                                    int Jmin = std::max(std::abs(dj12 - dj34), std::abs(dj13 - dj24));
                                    int Jmax = std::min(dj1, std::min(dj12 + dj34, dj13 + dj24));
                                    for (int J = Jmin; J <= Jmax; J += 2)
                                    {
                                        double w = func(dj1, dj2, dj12, dj3, dj4, dj34, dj13, dj24, J);
                                        sum += w;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << name << " time: " << duration << " ms\n";
    return sum;
}

int main()
{
    double sum = 0;
    sum = bench_3j(30, ef_3j, "3j");
    std::cout << "3j sum: " << sum << '\n';
    sum = bench_6j(20, ef_6j, "6j");
    std::cout << "6j sum: " << sum << '\n';
    sum = bench_9j(12, ef_9j, "9j");
    std::cout << "9j sum: " << sum << '\n';
    return 0;
}