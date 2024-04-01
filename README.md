# WignerSymbol

[![test](https://github.com/0382/WignerSymbol/actions/workflows/test.yml/badge.svg)](https://github.com/0382/WignerSymbol/actions/workflows/test.yml)

[中文](README_zh.md)

Calculate CG coefficient, Racah coefficient, and Wigner 3j, 6j, 9j coefficient. Calculation formula please see [CGcoefficient.jl](https://github.com/0382/CGcoefficient.jl).

## Usage

This is a head-only library, and only one file `WignerSymbol.hpp` is needed. 

```cpp
#include "WignerSymbol.hpp"
using namespace util;

int djmax = 21;
wigner_init(djmax, "2bjmax", 6);
double x = wigner_6j(dj1, dj2, dj3, dj4, dj5, dj6);
```

### Limitation 

For quite large quantum number, the package will give wrong answer, since it use float number arithmetic.
However, it is trustworthy for most real real world numerical calculation system. For `f9j`, it works at least about `Jmax = 25`.

## API

```cpp
// reserve binomial table
void wigner_init(int num, std::string type, int rank);
// fast access binomial table, it may return 0 for very large `n`
double fast_binomial(int n, int k);
// CG coefficient
double CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
// CG coefficient for two spin-1/2, equivalent to `CG(1, 1, 2*S, ds1, ds2, ds1+ds2)`, and faster
double CGspin(int ds1, int ds2, int S);
// CG coefficient with m1 == m2 == m3 == 0
double CG0(int j1, int j2, int j3);
// Wigner 3j symbol
double wigner_3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
// Wigner 6j symbol
double wigner_6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
// Racah coefficient
double Racah(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
// Wigner 9j symbol
double wigner_9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9);
// normalized Wigner 9j symbol
double wigner_norm9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9);
// LS-coupling to jj-coupling transformation coefficient
double lsjj(int l1, int l2, int dj1, int dj2, int L, int S, int J);
// Wigner d-function <j,m1|exp(i*beta*jy)|j,m2>
double dfunc(int dj, int dm1, int dm2, double beta);
// Moshinsky bracket，Ref: Buck et al. Nuc. Phys. A 600 (1996) 387-402
double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda, double tan_beta = 1.0);
```

Becase the angular momentum qunatum number can be half integers, people often use double of the exact quantum number as arguments. In this library, we also use the same convention. However, this library contains some other functions like `Moshinsky`, which only needs orbital quantum number, using doubled arguments is not needed.

So we take such convention: the arguments starts with `d`-character means double of the quantum number, while other arguments represent just the equivalent quantum number.

So if you want to calculate `<10|1/2,1/2;1/2,-1/2>`, you should call like this,

```cpp
double x = CG(1, 1, 2, 1, -1, 0);
double y = CGspin(1, -1, 1); // we use `S`, not `dS`, because `S` can only be `0, 1`
```

## The `wigner_init` function

We calculate the Wigner Symbols with `binomial`s, and we will store some binomials first, then when we need one binomial, we just load it. In this package, the `binomial` function is only valid in the stored range. If you call a `binomial` function out of the range, it just gives you `0`.

The `binomial` table is stored in the `WignerSymbols` class。 We define a globle variable `wigner` of the class to serve for all the functions. (The `inline` variable needs the `c++17` feature.)

```cpp
inline WignerSymbols wigner;
inline void wigner_init(int num, std::string type, int rank) { wigner.reserve(num, type, rank); }
```

When constructing the `WignerSymbols` object, it will store binomials from `binomial(0, 0)` to `binomial(67, 33)` (All can be exactly represented with a `uint64_t`). And one can use `wigner_init` function to extend the `binomial` table.

The `wigner_init` actually extends the maximum `n` for `binomial(n, k)`. To help users to make sure which `nmax` is safe for all of the following calculations, it defines several modes.

### `"2bjmax"`

```cpp
wigner_init(21, "2bjmax", 6);
```

This means the maximum single particle angular momentum is `21/2`, and thus the maximum two-body coupled angular momentum is `21`, and the `rank = 6` means you only need to calculate CG and/or 6j symbols, you don't need to calculate 9j symbols.

The `rank` can only be `3, 6, 9`, which respectively means `wigner_3j & CG`, `wigner_6j & Racah` and `wigner_9j` level calculation.

The `"2bjmax"` mode means your calculation only consider two-body coupling, and no three-body coupling. This mode assumes that in all these coefficients, at least one of the angular momentun is a single particle angular momentum, thus in this example no larger than `21/2`. With this assumption, `"2bjmax"` mode will use less memory than `"Jmax"` mode.

### `"Jmax"`

`"Jmax"` means the global maximum angular momentum, for every parameters.

```cpp
wigner_init(21, "Jmax", 9);
```

This means in all the calculations, `Jmax = 21`, and we calculate upto 9j symbols.

In the `"Jmax"` mode, it is always safe with out any assumption. Even having three-body coupling, you just need to use the maximum three body coupled angular momentum as `Jmax`, although it will cost more memory.

Actually, the memory used for store the binomials is not very large. For example, in the `"Jmax"` mode, and `Jmax = 200` for 9j calculations, the memory cost is just 2MB.

### `"nmax"`

The `"nmax"` mode directly set `nmax` of the `binomial` table, and the `rank` parameter is ignored. This maybe useful when you only want to calculate `binomial`s using this library.

### Exact `nmax`

The following table shows the exact `nmax` setted in different condition. See [Estimate-the-capacity](https://0382.github.io/CGcoefficient.jl/stable/formula/#Estimate-the-capacity) for details.

|                                       | Calculate range |  CG & 3j   | 6j & Racah |     9j     |
| :-----------------------------------: | :-------------: | :--------: | :--------: | :--------: |
|           meaning of `type`           | `type`\\`rank`  |     3      |     6      |     9      |
|         max angular momentum          |    `"Jmax"`     | `3*Jmax+1` | `4*Jmax+1` | `5*Jmax+1` |
| max two-body coupled angular momentum |   `"2bjmax"`    | `2*jmax+1` | `3*jmax+1` | `4*jmax+1` |
|             max binomial              |    `"nmax"`     |   `nmax`   |   `namx`   |   `nmax`   |

### Thread safety

The `wigner_init` function is **not** thread safe. So you shuld not call `winger_init` function dymanically in a multi-threading program. The correct way to use this package is find the maximum angular momentum quantum number in you system, and call `wigner_init` at the beginning of the code, and then don't call it any more.
