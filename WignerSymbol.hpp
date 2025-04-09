// MIT License

// Copyright (c) 2022 0382

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once
#ifndef JSHL_WIGNERSYMBOL_HPP
#define JSHL_WIGNERSYMBOL_HPP

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace util
{

#ifdef __SIZEOF_INT128__
using __ubin_t = __uint128_t;
constexpr int _bin_nmax = 131;
#else
using __ubin_t = uint64_t;
constexpr int _bin_nmax = 67;
#endif

class WignerSymbols
{
  public:
    WignerSymbols() : _nmax(_bin_nmax)
    {
        // initialize the data
        _binomial_data.resize(_binomial_data_size(_nmax), 1);
        __ubin_t temp[_bin_nmax / 2 + 1];
        std::fill(temp, temp + _nmax / 2 + 1, __ubin_t(1));
        std::size_t pos = 1;
        for (int n = 1; n <= _nmax; ++n)
        {
            ++pos;
            __ubin_t tk1 = 1;
            for (int k = 1; k <= n / 2; ++k)
            {
                __ubin_t x = (2 * k == n) ? 2 * tk1 : temp[k] + tk1;
                tk1 = temp[k];
                temp[k] = x;
                _binomial_data[pos++] = static_cast<double>(x);
            }
        }
    }

    // judge if a number is a odd number
    static bool isodd(int x) { return x % 2 != 0; }
    // judge if a number is a even number
    static bool iseven(int x) { return x % 2 == 0; }
    // judge if two number are same odd or same even
    static bool is_same_parity(int x, int y) { return iseven(x ^ y); }
    // return (-1)^n
    static int iphase(int x) { return iseven(x) - isodd(x); }
    // check if m-quantum number if one of the components of a the j-quantum number
    static bool check_jm(int dj, int dm) { return is_same_parity(dj, dm) && (std::abs(dm) <= dj); }
    // judge if three angular momentum can couple
    static bool check_couple(int dj1, int dj2, int dj3)
    {
        return dj1 >= 0 && dj2 >= 0 && is_same_parity(dj1 + dj2, dj3) && (dj3 <= (dj1 + dj2)) &&
               (dj3 >= std::abs(dj1 - dj2));
    }
    static bool check_couple_int(int j1, int j2, int j3)
    {
        return j1 >= 0 && j2 >= 0 && (j3 <= (j1 + j2)) && (j3 >= std::abs(j1 - j2));
    }
    // only works for positive n
    static double quick_pow(double x, int n)
    {
        double ans = 1;
        while (n)
        {
            if (n & 1)
                ans = ans * x;
            n = n >> 1;
            x = x * x;
        }
        return ans;
    }

    double binomial(int n, int k) const
    {
        if (unsigned(n) > unsigned(_nmax) || unsigned(k) > unsigned(n))
            return 0;
        else
        {
            k = std::min(k, n - k);
            return _binomial_data[_binomial_index(n, k)];
        }
    }

    double unsafe_binomial(int n, int k) const
    {
        k = std::min(k, n - k);
        return _binomial_data[_binomial_index(n, k)];
    }

    double CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3) const
    {
        if (!(check_jm(dj1, dm1) && check_jm(dj2, dm2) && check_jm(dj3, dm3)))
            return 0;
        if (!check_couple(dj1, dj2, dj3))
            return 0;
        if (dm1 + dm2 != dm3)
            return 0;
        const int J = (dj1 + dj2 + dj3) / 2;
        const int jm1 = J - dj1;
        const int jm2 = J - dj2;
        const int jm3 = J - dj3;
        const int j1mm1 = (dj1 - dm1) / 2;
        const int j2mm2 = (dj2 - dm2) / 2;
        const int j3mm3 = (dj3 - dm3) / 2;
        const int j2pm2 = (dj2 + dm2) / 2;
        const double A = std::sqrt(unsafe_binomial(dj1, jm2) * unsafe_binomial(dj2, jm3) /
                                   (unsafe_binomial(J + 1, jm3) * unsafe_binomial(dj1, j1mm1) *
                                    unsafe_binomial(dj2, j2mm2) * unsafe_binomial(dj3, j3mm3)));
        double B = 0;
        const int low = std::max(0, std::max(j1mm1 - jm2, j2pm2 - jm1));
        const int high = std::min(jm3, std::min(j1mm1, j2pm2));
        for (auto z = low; z <= high; ++z)
        {
            B = -B + unsafe_binomial(jm3, z) * unsafe_binomial(jm2, j1mm1 - z) * unsafe_binomial(jm1, j2pm2 - z);
        }
        return iphase(high) * A * B;
    }

    double CG0(int j1, int j2, int j3) const
    {
        if (!check_couple_int(j1, j2, j3))
            return 0;
        const int J = j1 + j2 + j3;
        if (isodd(J))
            return 0;
        const int g = J / 2;
        return iphase(g - j3) * unsafe_binomial(g, j3) * unsafe_binomial(j3, g - j1) /
               std::sqrt(unsafe_binomial(J + 1, 2 * j3 + 1) * unsafe_binomial(2 * j3, J - 2 * j1));
    }

    double f3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3) const
    {
        if (!(check_jm(dj1, dm1) && check_jm(dj2, dm2) && check_jm(dj3, dm3)))
            return 0;
        if (!check_couple(dj1, dj2, dj3))
            return 0;
        if (dm1 + dm2 + dm3 != 0)
            return 0;
        const int J = (dj1 + dj2 + dj3) / 2;
        const int jm1 = J - dj1;
        const int jm2 = J - dj2;
        const int jm3 = J - dj3;
        const int j1mm1 = (dj1 - dm1) / 2;
        const int j2mm2 = (dj2 - dm2) / 2;
        const int j3mm3 = (dj3 - dm3) / 2;
        const int j1pm1 = (dj1 + dm1) / 2;
        const double A = std::sqrt(unsafe_binomial(dj1, jm2) * unsafe_binomial(dj2, jm1) /
                                   ((J + 1) * unsafe_binomial(J, jm3) * unsafe_binomial(dj1, j1mm1) *
                                    unsafe_binomial(dj2, j2mm2) * unsafe_binomial(dj3, j3mm3)));
        double B = 0;
        const int low = std::max(0, std::max(j1pm1 - jm2, j2mm2 - jm1));
        const int high = std::min(jm3, std::min(j1pm1, j2mm2));
        for (auto z = low; z <= high; ++z)
        {
            B = -B + unsafe_binomial(jm3, z) * unsafe_binomial(jm2, j1pm1 - z) * unsafe_binomial(jm1, j2mm2 - z);
        }
        return iphase(dj1 + (dj3 + dm3) / 2 + high) * A * B;
    }

    double f6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6) const
    {
        if (!(check_couple(dj1, dj2, dj3) && check_couple(dj1, dj5, dj6) && check_couple(dj4, dj2, dj6) &&
              check_couple(dj4, dj5, dj3)))
            return 0;
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
        const double A = std::sqrt(unsafe_binomial(j123 + 1, dj1 + 1) * unsafe_binomial(dj1, jpm123) /
                                   (unsafe_binomial(j156 + 1, dj1 + 1) * unsafe_binomial(dj1, jpm156) *
                                    unsafe_binomial(j453 + 1, dj4 + 1) * unsafe_binomial(dj4, jpm453) *
                                    unsafe_binomial(j426 + 1, dj4 + 1) * unsafe_binomial(dj4, jpm426)));
        double B = 0;
        const int low = std::max(j123, std::max(j156, std::max(j426, j453)));
        const int high = std::min(jpm123 + j453, std::min(jpm132 + j426, jpm231 + j156));
        for (auto x = low; x <= high; ++x)
        {
            B = -B + unsafe_binomial(x + 1, j123 + 1) * unsafe_binomial(jpm123, x - j453) *
                         unsafe_binomial(jpm132, x - j426) * unsafe_binomial(jpm231, x - j156);
        }
        return iphase(high) * A * B / (dj4 + 1);
    }

    double Racah(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6) const
    {
        return iphase((dj1 + dj2 + dj3 + dj4) / 2) * f6j(dj1, dj2, dj5, dj4, dj3, dj6);
    }

    double f9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9) const
    {
        if (!(check_couple(dj1, dj2, dj3) && check_couple(dj4, dj5, dj6) && check_couple(dj7, dj8, dj9) &&
              check_couple(dj1, dj4, dj7) && check_couple(dj2, dj5, dj8) && check_couple(dj3, dj6, dj9)))
            return 0;
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
        const double P0_nu = unsafe_binomial(j123 + 1, dj1 + 1) * unsafe_binomial(dj1, pm123) * //
                             unsafe_binomial(j456 + 1, dj5 + 1) * unsafe_binomial(dj5, pm456) * //
                             unsafe_binomial(j789 + 1, dj9 + 1) * unsafe_binomial(dj9, pm798);
        const double P0_de = unsafe_binomial(j147 + 1, dj1 + 1) * unsafe_binomial(dj1, (dj1 + dj4 - dj7) / 2) *
                             unsafe_binomial(j258 + 1, dj5 + 1) * unsafe_binomial(dj5, (dj2 + dj5 - dj8) / 2) *
                             unsafe_binomial(j369 + 1, dj9 + 1) * unsafe_binomial(dj9, (dj3 + dj9 - dj6) / 2);
        const double P0 = std::sqrt(P0_nu / P0_de);
        const int dtl = std::max(std::abs(dj2 - dj6), std::max(std::abs(dj4 - dj8), std::abs(dj1 - dj9)));
        const int dth = std::min(dj2 + dj6, std::min(dj4 + dj8, dj1 + dj9));
        double PABC = 0;
        for (auto dt = dtl; dt <= dth; dt += 2)
        {
            const int j19t = (dj1 + dj9 + dt) / 2;
            const int j26t = (dj2 + dj6 + dt) / 2;
            const int j48t = (dj4 + dj8 + dt) / 2;
            double Pt_de = unsafe_binomial(j19t + 1, dt + 1) * unsafe_binomial(dt, (dj1 + dt - dj9) / 2) *
                           unsafe_binomial(j26t + 1, dt + 1) * unsafe_binomial(dt, (dj2 + dt - dj6) / 2) *
                           unsafe_binomial(j48t + 1, dt + 1) * unsafe_binomial(dt, (dj4 + dt - dj8) / 2);
            Pt_de *= (dt + 1) * (dt + 1);
            const int xl = std::max(j123, std::max(j369, std::max(j26t, j19t)));
            const int xh = std::min(pm123 + j369, std::min(pm132 + j26t, pm231 + j19t));
            double At = 0;
            for (auto x = xl; x <= xh; ++x)
            {
                At = -At + unsafe_binomial(x + 1, j123 + 1) * unsafe_binomial(pm123, x - j369) *
                               unsafe_binomial(pm132, x - j26t) * unsafe_binomial(pm231, x - j19t);
            }
            const int yl = std::max(j456, std::max(j26t, std::max(j258, j48t)));
            const int yh = std::min(pm456 + j26t, std::min(pm465 + j258, pm564 + j48t));
            double Bt = 0;
            for (auto y = yl; y <= yh; ++y)
            {
                Bt = -Bt + unsafe_binomial(y + 1, j456 + 1) * unsafe_binomial(pm456, y - j26t) *
                               unsafe_binomial(pm465, y - j258) * unsafe_binomial(pm564, y - j48t);
            }
            const int zl = std::max(j789, std::max(j19t, std::max(j48t, j147)));
            const int zh = std::min(pm789 + j19t, std::min(pm798 + j48t, pm897 + j147));
            double Ct = 0;
            for (auto z = zl; z <= zh; ++z)
            {
                Ct = -Ct + unsafe_binomial(z + 1, j789 + 1) * unsafe_binomial(pm789, z - j19t) *
                               unsafe_binomial(pm798, z - j48t) * unsafe_binomial(pm897, z - j147);
            }
            PABC += iphase(xh + yh + zh) * At * Bt * Ct / Pt_de;
        }
        return iphase(dth) * P0 * PABC;
    }

    double norm9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9)
    {
        return f9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9) *
               std::sqrt((dj3 + 1.) * (dj6 + 1.) * (dj7 + 1.) * (dj8 + 1.));
    }

    static double lsjj(int l1, int l2, int dj1, int dj2, int L, int S, int J)
    {
        if (!check_couple(2 * l1, 2 * l2, 2 * L))
            return 0.0;
        if (!check_couple(dj1, dj2, 2 * J))
            return 0.0;

        int LSJcase = 0;
        if (S == 0)
        {
            LSJcase = 1;
        }
        else if (S == 1)
        {
            if (L == J - 1)
                LSJcase = 2;
            else if (J == L)
                LSJcase = 3;
            else if (L == J + 1)
                LSJcase = 4;
        }
        if (LSJcase == 0)
            return 0.0;

        const int m = l1 - l2;
        const int p = l1 + l2;
        const int d1 = dj1 - 2 * l1;
        const int d2 = dj2 - 2 * l2;
        const double den = 2 * (2 * l1 + 1) * (2 * l2 + 1);
        double r = 0.0;

        if (d1 == 1 && d2 == 1)
        {
            if (LSJcase == 1)
            {
                r = (J + p + 2) * (p + 1 - J) / den;
            }
            else if (LSJcase == 2)
            {
                r = (J + m) * (J - m) / static_cast<double>(J * (2 * J + 1));
                r *= (J + p + 1) * (J + p + 2) / den;
            }
            else if (LSJcase == 3)
            {
                if (m == 0)
                    return 0.0;
                r = m * std::abs(m) / static_cast<double>(J * (J + 1));
                r *= (J + p + 2) * (p + 1 - J) / den;
            }
            else if (LSJcase == 4)
            {
                r = -(L + m) * (L - m) / static_cast<double>(L * (2 * J + 1));
                r *= (p - J) * (p + 1 - J) / den;
            }
        }
        else if (d1 == 1 && d2 == -1)
        {
            if (LSJcase == 1)
            {
                r = (J + m + 1) * (J - m) / den;
            }
            else if (LSJcase == 2)
            {
                r = -(p + 1 + J) * (p + 1 - J) / static_cast<double>(J * (2 * J + 1));
                r *= (J + m + 1) * (J + m) / den;
            }
            else if (LSJcase == 3)
            {
                r = (p + 1) * (p + 1) / static_cast<double>(J * (J + 1));
                r *= (J + m + 1) * (J - m) / den;
            }
            else if (LSJcase == 4)
            {
                r = -(p + 1 + L) * (p + 1 - L) / static_cast<double>(L * (2 * J + 1));
                r *= (J - m + 1) * (J - m) / den;
            }
        }
        else if (d1 == -1 && d2 == 1)
        {
            if (LSJcase == 1)
            {
                r = -(J + m) * (J - m + 1) / den;
            }
            else if (LSJcase == 2)
            {
                r = (p + 1 + J) * (p + 1 - J) / static_cast<double>(J * (2 * J + 1));
                r *= (J - m) * (J - m + 1) / den;
            }
            else if (LSJcase == 3)
            {
                r = (p + 1) * (p + 1) / static_cast<double>(J * (J + 1));
                r *= (J + m) * (J - m + 1) / den;
            }
            else if (LSJcase == 4)
            {
                r = (p + 1 + L) * (p + 1 - L) / static_cast<double>(L * (2 * J + 1));
                r *= (J + m + 1) * (J + m) / den;
            }
        }
        else if (d1 == -1 && d2 == -1)
        {
            if (LSJcase == 1)
            {
                r = (J + p + 1) * (p - J) / den;
            }
            else if (LSJcase == 2)
            {
                r = -(J + m) * (J - m) / static_cast<double>(J * (2 * J + 1));
                r *= (p - J) * (p + 1 - J) / den;
            }
            else if (LSJcase == 3)
            {
                if (m == 0)
                    return 0.0;
                r = -m * std::abs(m) / static_cast<double>(J * (J + 1));
                r *= (J + p + 1) * (p - J) / den;
            }
            else if (LSJcase == 4)
            {
                r = (L + m) * (L - m) / static_cast<double>(L * (2 * J + 1));
                r *= (J + p + 1) * (J + p + 2) / den;
            }
        }

        return std::copysign(std::sqrt(std::abs(r)), r);
    }

    double _m9j(int j1, int j2, int j3, int j4, int j5, int j6, int j7, int j8, int j9) const
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
        double sum = 0.0;
        const int tl = std::max(std::abs(j2 - j6), std::max(std::abs(j4 - j8), std::abs(j1 - j9)));
        const int th = std::min(j2 + j6, std::min(j4 + j8, j1 + j9));
        for (int t = tl; t <= th; ++t)
        {
            const int j19t = j1 + j9 + t;
            const int j26t = j2 + j6 + t;
            const int j48t = j4 + j8 + t;
            const int dt = 2 * t;
            double Pt_de = unsafe_binomial(j19t + 1, dt + 1) * unsafe_binomial(dt, j1 + t - j9) *
                           unsafe_binomial(j26t + 1, dt + 1) * unsafe_binomial(dt, j2 + t - j6) *
                           unsafe_binomial(j48t + 1, dt + 1) * unsafe_binomial(dt, j4 + t - j8);
            Pt_de *= (dt + 1) * (dt + 1);
            const int xl = std::max(j123, std::max(j369, std::max(j26t, j19t)));
            const int xh = std::min(pm123 + j369, std::min(pm132 + j26t, pm231 + j19t));
            double At = 0.0;
            for (int x = xl; x <= xh; ++x)
            {
                At = -At + unsafe_binomial(x + 1, j123 + 1) * unsafe_binomial(pm123, x - j369) *
                               unsafe_binomial(pm132, x - j26t) * unsafe_binomial(pm231, x - j19t);
            }
            const int yl = std::max(j456, std::max(j26t, std::max(j258, j48t)));
            const int yh = std::min(pm456 + j26t, std::min(pm465 + j258, pm564 + j48t));
            double Bt = 0.0;
            for (int y = yl; y <= yh; ++y)
            {
                Bt = -Bt + unsafe_binomial(y + 1, j456 + 1) * unsafe_binomial(pm456, y - j26t) *
                               unsafe_binomial(pm465, y - j258) * unsafe_binomial(pm564, y - j48t);
            }
            const int zl = std::max(j789, std::max(j19t, std::max(j48t, j147)));
            const int zh = std::min(pm789 + j19t, std::min(pm798 + j48t, pm897 + j147));
            double Ct = 0.0;
            for (int z = zl; z <= zh; ++z)
            {
                Ct = -Ct + unsafe_binomial(z + 1, j789 + 1) * unsafe_binomial(pm789, z - j19t) *
                               unsafe_binomial(pm798, z - j48t) * unsafe_binomial(pm897, z - j147);
            }
            sum += iphase(xh + yh + zh) * At * Bt * Ct / Pt_de;
        }
        return sum;
    }

    // Buck et al. Nuc. Phys. A 600 (1996) 387-402
    double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda,
                     double tan_beta = 1.0) const
    {
        if (!check_couple_int(L, l, lambda) || !check_couple_int(l1, l2, lambda))
            return 0.0;
        // check energy conservation
        const int f1 = 2 * n1 + l1;
        const int f2 = 2 * n2 + l2;
        const int F = 2 * N + L;
        const int f = 2 * n + l;
        if (f1 + f2 != f + F)
            return 0;

        const int nl1 = n1 + l1;
        const int nl2 = n2 + l2;
        const int NL = N + L;
        const int nl = n + l;
        const int chi = f1 + f2;

        const double cos_beta = 1.0 / std::sqrt(1.0 + tan_beta * tan_beta);
        const double sin_beta = tan_beta * cos_beta;
        double pre = unsafe_binomial(chi + 2, f1 + 1) / unsafe_binomial(chi + 2, F + 1);
        pre *= unsafe_binomial(L + l + lambda + 1, 2 * lambda + 1) * unsafe_binomial(2 * lambda, lambda + L - l);
        pre /= unsafe_binomial(l1 + l2 + lambda + 1, 2 * lambda + 1) * unsafe_binomial(2 * lambda, lambda + l1 - l2);

        pre *= (2 * l1 + 1) * unsafe_binomial(2 * nl1 + 1, nl1) / (unsafe_binomial(f1 + 1, n1) * quick_pow(2.0, l1));
        pre *= (2 * l2 + 1) * unsafe_binomial(2 * nl2 + 1, nl2) / (unsafe_binomial(f2 + 1, n2) * quick_pow(2.0, l2));
        pre *= (2 * L + 1) * unsafe_binomial(2 * NL + 1, NL) / (unsafe_binomial(F + 1, N) * quick_pow(2.0, L));
        pre *= (2 * l + 1) * unsafe_binomial(2 * nl + 1, nl) / (unsafe_binomial(f + 1, n) * quick_pow(2.0, l));
        pre = std::sqrt(pre) / ((f1 + 2) * (f2 + 2));

        double sum = 0.0;
        for (int fa = 0; fa <= std::min(f1, F); ++fa)
        {
            const int fb = f1 - fa;
            const int fc = F - fa;
            const int fd = f2 - fc;
            if (fd < 0)
                continue;
            const double tfa = quick_pow(sin_beta, fa + fd) * quick_pow(cos_beta, fb + fc) *
                               unsafe_binomial(f1 + 2, fa + 1) * unsafe_binomial(f2 + 2, fc + 1);
            for (int la = fa & 0x01; la <= fa; la += 2)
            {
                const int na = (fa - la) / 2;
                const int nla = na + la;
                const double _ta =
                    quick_pow(2, la) * (2 * la + 1) * unsafe_binomial(fa + 1, na) / unsafe_binomial(2 * nla + 1, nla);
                const double ta = tfa * _ta;
                for (int lb = std::abs(l1 - la); lb <= std::min(l1 + la, fb); lb += 2)
                {
                    const int nb = (fb - lb) / 2;
                    const int nlb = nb + lb;
                    const double _tb = quick_pow(2, lb) * (2 * lb + 1) * unsafe_binomial(fb + 1, nb) /
                                       unsafe_binomial(2 * nlb + 1, nlb);
                    const int g1 = (la + lb + l1) / 2;
                    const double t1 = unsafe_binomial(g1, l1) * unsafe_binomial(l1, g1 - la);
                    const double tb = ta * _tb * t1;
                    for (int lc = std::abs(L - la); lc <= std::min(L + la, fc); lc += 2)
                    {
                        const int nc = (fc - lc) / 2;
                        const int nlc = nc + lc;
                        const double _tc = quick_pow(2, lc) * (2 * lc + 1) * unsafe_binomial(fc + 1, nc) /
                                           unsafe_binomial(2 * nlc + 1, nlc);
                        const int g3 = (la + lc + L) / 2;
                        const double t3 = unsafe_binomial(g3, L) * unsafe_binomial(L, g3 - la);
                        const double d3 = (2 * L + 1) * unsafe_binomial(la + lc + L + 1, 2 * L + 1) *
                                          unsafe_binomial(2 * L, L + la - lc);
                        const double tc = tb * _tc * t3 / d3;
                        const int ldmin = std::max(std::abs(l2 - lc), std::abs(l - lb));
                        const int ldmax = std::min(fd, std::min(l2 + lc, l + lb));
                        for (int ld = ldmin; ld <= ldmax; ld += 2)
                        {
                            const int nd = (fd - ld) / 2;
                            const int nld = nd + ld;
                            const double _td = quick_pow(2, ld) * (2 * ld + 1) * unsafe_binomial(fd + 1, nd) /
                                               unsafe_binomial(2 * nld + 1, nld);
                            const int g2 = (lc + ld + l2) / 2;
                            const double t2 = unsafe_binomial(g2, l2) * unsafe_binomial(l2, g2 - lc);
                            const int g4 = (lb + ld + l) / 2;
                            const double t4 = unsafe_binomial(g4, l) * unsafe_binomial(l, g4 - lb);
                            const double d4 = (2 * l + 1) * unsafe_binomial(lb + ld + l + 1, 2 * l + 1) *
                                              unsafe_binomial(2 * l, l + lb - ld);
                            const double td = tc * _td * t2 * t4 / d4;
                            const double m9j = _m9j(la, lb, l1, lc, ld, l2, L, l, lambda);
                            sum += iphase(ld) * td * m9j;
                        }
                    }
                }
            }
        }
        return pre * sum;
    }

    double dfunc(int dj, int dm1, int dm2, double beta) const
    {
        if (!(check_jm(dj, dm1) && check_jm(dj, dm2)))
            return 0.;
        const int jm1 = (dj - dm1) / 2;
        const int jp1 = (dj + dm1) / 2;
        const int jm2 = (dj - dm2) / 2;
        const int mm = (dm1 + dm2) / 2;
        const double c = std::cos(beta / 2);
        const double s = std::sin(beta / 2);
        const int kmin = std::max(0, -mm);
        const int kmax = std::min(jm1, jm2);
        double sum = 0.;
        for (int k = kmin; k <= kmax; ++k)
        {
            sum = -sum + unsafe_binomial(jm1, k) * unsafe_binomial(jp1, mm + k) * quick_pow(c, mm + 2 * k) *
                             quick_pow(s, jm1 + jm2 - 2 * k);
        }
        sum = iphase(jm2 + kmax) * sum;
        sum = sum * std::sqrt(unsafe_binomial(dj, jm1) / unsafe_binomial(dj, jm2));
        return sum;
    }

    void reserve(int num, std::string type, int rank)
    {
        if (type == "Jmax")
        {
            switch (rank)
            {
                case 3: fill_binomial_data(3 * num + 1); break;
                case 6: fill_binomial_data(4 * num + 1); break;
                case 9: fill_binomial_data(5 * num + 1); break;
                default:
                {
                    std::cerr << "Error: rank must be 3, 6, or 9" << std::endl;
                    std::exit(-1);
                }
            }
        }
        else if (type == "2bjmax")
        {
            switch (rank)
            {
                case 3: fill_binomial_data(2 * num + 1); break;
                case 6: fill_binomial_data(3 * num + 1); break;
                case 9: fill_binomial_data(4 * num + 1); break;
                default:
                {
                    std::cerr << "Error: rank must be 3, 6, or 9" << std::endl;
                    std::exit(-1);
                }
            }
        }
        else if (type == "Moshinsky")
        {
            fill_binomial_data(6 * num + 1);
        }
        else if (type == "nmax")
        {
            fill_binomial_data(num);
        }
        else
        {
            std::cerr << "Error: type must be Jmax, 2bjmax, Moshinsky or nmax" << std::endl;
            std::exit(-1);
        }
    }

  private:
    std::vector<double> _binomial_data;
    int _nmax;
    static std::size_t _binomial_data_size(int n)
    {
        std::size_t x = n / 2 + 1;
        return x * (x + isodd(n));
    }
    static std::size_t _binomial_index(int n, int k)
    {
        std::size_t x = n / 2 + 1;
        return x * (x - iseven(n)) + k;
    }
    void fill_binomial_data(int nmax)
    {
        if (nmax <= _nmax)
            return;
        std::vector<double> old_data = _binomial_data;
        std::size_t reserve_size = _binomial_data_size(nmax);
        if (reserve_size > std::numeric_limits<int>::max())
        {
            std::cerr << "Error: nmax too large" << std::endl;
            std::exit(-1);
        }
        _binomial_data.resize(reserve_size);
        std::copy(std::begin(old_data), std::end(old_data), _binomial_data.begin());
        for (int n = _nmax + 1; n <= nmax; ++n)
        {
            for (int k = 0; k <= n / 2; ++k)
            {
                _binomial_data[_binomial_index(n, k)] = binomial(n - 1, k) + binomial(n - 1, k - 1);
            }
            ++_nmax;
        }
        _nmax = nmax;
    }
};

inline WignerSymbols wigner;

inline void wigner_init(int num, std::string type, int rank) { wigner.reserve(num, type, rank); }

inline double fast_binomial(int n, int k) { return wigner.binomial(n, k); }

// CG coefficient for two spin-1/2
inline double CGspin(int dm1, int dm2, int S)
{
    static constexpr double inv_sqrt_2 = 0.7071067811865476;
    static constexpr double values[2][2][2] = {0.0, 1.0, -inv_sqrt_2, inv_sqrt_2, inv_sqrt_2, inv_sqrt_2, 0.0, 1.0};
    if (unsigned(S) > 1)
        return 0;
    if (std::abs(dm1) != 1 || std::abs(dm2) != 1)
        return 0;
    return values[dm1 > 0][dm2 > 0][S];
}

// <S12,M12|1/2,m1;1/2,m2><S,M|S12,M12;1/2,m3>
inline double CG3spin(int dm1, int dm2, int dm3, int S12, int dS)
{
    static constexpr double values[2][2][2][3] = {
        0.0,                 // (-1/2, -1/2, -1/2,   0,  1/2) -> 0
        0.0,                 // (-1/2, -1/2, -1/2,   1,  1/2) -> 0
        1.0,                 // (-1/2, -1/2, -1/2,   1,  3/2) -> 1
        0.0,                 // (-1/2, -1/2,  1/2,   0,  1/2) -> 0
        -0.816496580927726,  // (-1/2, -1/2,  1/2,   1,  1/2) -> -sqrt(2/3)
        0.5773502691896257,  // (-1/2, -1/2,  1/2,   1,  3/2) -> sqrt(1/3)
        -0.7071067811865476, // (-1/2,  1/2, -1/2,   0,  1/2) -> -sqrt(1/2)
        0.408248290463863,   // (-1/2,  1/2, -1/2,   1,  1/2) -> sqrt(1/6)
        0.5773502691896257,  // (-1/2,  1/2, -1/2,   1,  3/2) -> sqrt(1/3)
        -0.7071067811865476, // (-1/2,  1/2,  1/2,   0,  1/2) -> -sqrt(1/2)
        -0.408248290463863,  // (-1/2,  1/2,  1/2,   1,  1/2) -> -sqrt(1/6)
        0.5773502691896257,  // (-1/2,  1/2,  1/2,   1,  3/2) -> sqrt(1/3)
        0.7071067811865476,  // ( 1/2, -1/2, -1/2,   0,  1/2) -> sqrt(1/2)
        0.408248290463863,   // ( 1/2, -1/2, -1/2,   1,  1/2) -> sqrt(1/6)
        0.5773502691896257,  // ( 1/2, -1/2, -1/2,   1,  3/2) -> sqrt(1/3)
        0.7071067811865476,  // ( 1/2, -1/2,  1/2,   0,  1/2) -> sqrt(1/2)
        -0.408248290463863,  // ( 1/2, -1/2,  1/2,   1,  1/2) -> -sqrt(1/6)
        0.5773502691896257,  // ( 1/2, -1/2,  1/2,   1,  3/2) -> sqrt(1/3)
        0.0,                 // ( 1/2,  1/2, -1/2,   0,  1/2) -> 0
        0.816496580927726,   // ( 1/2,  1/2, -1/2,   1,  1/2) -> sqrt(2/3)
        0.5773502691896257,  // ( 1/2,  1/2, -1/2,   1,  3/2) -> sqrt(1/3)
        0.0,                 // ( 1/2,  1/2,  1/2,   0,  1/2) -> 0
        0.0,                 // ( 1/2,  1/2,  1/2,   1,  1/2) -> 0
        1.0                  // ( 1/2,  1/2,  1/2,   1,  3/2) -> 1
    };
    if (unsigned(S12) > 1)
        return 0;
    if (S12 == 0 && dS != 1)
        return 0;
    if (S12 == 1 && (dS != 1 && dS != 3))
        return 0;
    if (std::abs(dm1) != 1 || std::abs(dm2) != 1 || std::abs(dm3) != 1)
        return 0;
    return values[dm1 > 0][dm2 > 0][dm3 > 0][(S12 + dS) / 2];
}

inline double CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    return wigner.CG(dj1, dj2, dj3, dm1, dm2, dm3);
}

inline double CG0(int j1, int j2, int j3) { return wigner.CG0(j1, j2, j3); }

inline double wigner_3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3)
{
    return wigner.f3j(dj1, dj2, dj3, dm1, dm2, dm3);
}

inline double wigner_6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6)
{
    return wigner.f6j(dj1, dj2, dj3, dj4, dj5, dj6);
}

inline double Racah(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6)
{
    return wigner.Racah(dj1, dj2, dj3, dj4, dj5, dj6);
}

inline double wigner_9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9)
{
    return wigner.f9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
}

inline double wigner_norm9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9)
{
    return wigner.norm9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
}

inline double lsjj(int l1, int l2, int dj1, int dj2, int L, int S, int J)
{
    return WignerSymbols::lsjj(l1, l2, dj1, dj2, L, S, J);
}

inline double dfunc(int dj, int dm1, int dm2, double beta) { return wigner.dfunc(dj, dm1, dm2, beta); }

inline double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda, double tan_beta = 1.0)
{
    return wigner.Moshinsky(N, L, n, l, n1, l1, n2, l2, lambda, tan_beta);
}

} // end namespace util

#endif // JSHL_WIGNERSYMBOL_HPP
