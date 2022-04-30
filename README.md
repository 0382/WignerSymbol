# WignerSymbol

[中文](README_zh.md)

Calculate CG coefficient, Racah coefficient, and Wigner 3j, 6j, 9j coefficient. Calculation formula please see [CGcoefficient.jl](https://github.com/0382/CGcoefficient.jl).

## Usage

Copy the `WignerSymbol.hpp` file to you code, and include it.

```cpp
using namespace jshl;
WignerSymbols wigner;
int djmax = 21;
wigner.reserve(djmax, "2bjmax", 6);
double x = wigner.f6j(dj1, dj2, dj3, dj4, dj5, dj6);
```

## API
```cpp
// binomial
double WignerSymbols::binomial(int n, int k);
// CG coefficient
double WignerSymbols::CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
// Wigner 3j symbol
double WignerSymbols::f3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
// Wigner 6j symbol
double WignerSymbols::f6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
// Racah coefficient
double WignerSymbols::Racah(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
// Wigner 9j symbol
double WignerSymbols::f9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9);
```
Apart from `binomial`, all the functions use double of the real angular momentum quantum number to avoid half integers. So if you want to calculate `<10|1/2,1/2;1/2,-1/2>`, you should call like this,
```cpp
WignerSymbols wigner;
double x = wigner.CG(1,1,2,1,-1,0);
```

## The `reserve` function

We calculate the Wigner Symbols with `binomial`s, and we will store some binomials first, then when we need one binomial, we just load it. In this package, the `binomial` function is only valid in the stored range. If you call a `binomial` function out of the range, it just gives you `0`.

When constructing the `WignerSymbols` object, it will store binomials from `binomial(0, 0)` to `binomial(67, 33)`. You can use `reserve` function to extent the range. The `reserve` function is
```cpp
void WignerSymbols::reserve(int num, std::string type, int rank)
```
and its parameters means

|                          |    Calculate range    |   CG & 3j   | 6j & Racah  |     9j      |
| :----------------------: | :------------: | :---------: | :---------: | :---------: |
|       meaning of `type`           | `type`\\`rank` |      3      |      6      |      9      |
|      max angular momentum mode    |    `"Jmax"`    | `3*Jmax+1`  | `4*Jmax+1`  | `5*Jmax+1`  |
| max two-body coupled angular momentum mode |   `"2bjmax"`    | `2*jmax+1` | `3*jmax+1` | `4*jmax+1` |
|    max binomial mode             |    `"nmax"`    |   `nmax`    |   `namx`    |   `nmax`    |

The value in the table means the minimum binomial range to guarantee the Wigner Symbol calculation. You do not need to rememmber those values. You just need to find the maximum angular momentum in you canculation, wen call the `reserve` function.

### `"2bjmax"`

For example

```cpp
WignerSymbols wigner;
int djmax = 21;
wigner.reserve(djmax, "2bjmax", 6);
```

This means the maximum single particle angular momentum is `21/2`, and thus the maximum two-body coupled angular momentum is `21`, and the `rank = 6` means you only need to calculate CG and/or 6j symbols, you don't need to calculate 9j symbol. In this example, the `n` of the stored maximum binomial coefficient is

```cpp
nmax = 3*jmax + 1 = 3 * 21 + 1 = 64
```

Because `64` is less than `67`, the `reserve` function will do nothing.

The `"2bjmax"` mode means your calculation only consider two-body coupling, and no three-body coupling. This mode assumes that in all these coefficients, at least one of the angular momentun is just a single particle angular momentum, thus in this example no larger than `21/2`. With this assumption, `"2bjmax"` mode will use less memory than `"Jmax"` mode.

### `"Jmax"`

`"Jmax"` means the global maximum angular momentum, for every parameters. For another example,

```cpp
WignerSymbols wigner;
int Jmax = 21;
wigner.reserve(Jmax, "Jmax", 6);
```

This means in this system, `Jmax = 21`, and we calculate to 6j symbols. Here
```cpp
nmax = 4*Jmax + 1 = 85
```
So `reserve` function will extent the stored binomials range. For quantum many body calculation with only two-body coupling, if single particle `jmax = 21/2`, then `Jmax = 21`. The `"2bjmax"` mode will cost less storage.

In the `"Jmax"` mode, it is always safe with out any assumption. Even having three-body coupling, you just need to use the maximum three body coupled angular momentum as `Jmax`, although it will cost more memory.

Actually, the memory used for store the binomials is not very large. A simple estimate is
```cpp
2 * nmax * nmax // Byte
```
Even for `nmax = 1000` (`Jmax = 200` in 9j calculation, which is absolutly enough for most calculation), the memory cost is just 2MB, which is not very large.

### `"nmax"`

The `"nmax"` mode directly set `nmax`, and the `rank` parameter is ignored. This maybe useful when you only want to calculate `binomial`s using this package.

### Thread safety

The `reserve` is **not** thread safe. So you shuld not call `reserve` function dymanically in a multi-threading program. The correct way to use this package is find the maximum angular momentum quantum number in you system, and call `reserve` at the beginning of the code, and then don't call it any more.
