# WignerSymbol

[中文](README.md)

Calculate CG coefficient, Racah coefficient, and Wigner 3j, 6j, 9j coefficient. Calculation formula please see [CGcoefficient.jl](https://github.com/0382/CGcoefficient.jl).

## Usage

Copy the `WignerSymbol.hpp` file to you code, and include it.

```cpp
using namespace jshl;
WignerSymbols wigner;
int djmax = 21;
wigner.reserve(djmax, "djmax", 6);
double x = wigner.f6j(dj1, dj2, dj3, dj4, dj5, dj6);
```

## API
```cpp
// 二项式系数
double WignerSymbols::binominal(int n, int k);
// CG系数
double WignerSymbols::CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
// 3j系数
double WignerSymbols::f3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
// 6j系数
double WignerSymbols::f6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
// Racah系数
double WignerSymbols::Racah(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
// 9j系数
double WignerSymbols::f9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9);
```
Apart from `binominal`, all the functions use double of the real angular momentum quantum number to avoid half integers. So if you want to calculate `<10|1/2,1/2;1/2,-1/2>`, you should call like this,
```cpp
WignerSymbols wigner;
double x = wigner.CG(1,1,2,1,-1,0);
```

## The `reserve` function

We calculate the Wigner Symbols with `binominal`s, and we will store some binominals first, then when we need one binominal, we just load it. In this package, the `binominal` is only valid in the stored range. If you call a `binominal` out of the range, it just give you `0`.

When constructing the `WignerSymbols` object, it will store binominals from `binominal(0, 0)` to `binominal(67, 33)`. You can use `reserve` function to extent the range. The `reserve` function is
```cpp
void WignerSymbols::reserve(int num, std::string type, int rank)
```
and its parameters are

|                          |    Calculate range    |   CG & 3j   | 6j & Racah  |     9j      |
| :----------------------: | :------------: | :---------: | :---------: | :---------: |
|       meaning of `type`           | `type`\\`rank` |      3      |      6      |      9      |
|      max angular momentum mode    |    `"Jmax"`    | `3*Jmax+1`  | `4*Jmax+1`  | `5*Jmax+1`  |
| max single particle momentum mode |   `"djmax"`    | `2*djmax+1` | `3*djmax+1` | `4*djmax+1` |
|    max binominal mode             |    `"nmax"`    |   `nmax`    |   `namx`    |   `nmax`    |

The value in the table means the minimum binominal range to guarantee the Wigner Symbol calculation. You do not need to rememmber those values. You just need to find the maximum angular momentum in you canculation, wen call the `reserve` function.

### `"djmax"`

For example

```cpp
WignerSymbols wigner;
int djmax = 21;
wigner.reserve(djmax, "djmax", 6);
```

This means the maximum single particle angular momentum is `21/2`, and you only need to calculate CG and/or 6j symbols, you don't need to calculate 9j symbol. In this example, the `n` of the stored maximum binominal coefficient is

```cpp
nmax = 3*djmax + 1 = 3*21 + 1 = 64
```

Because `64` is less than `67`, the `reserve` function will do nothing.

Warning: This range is deduced assuming only two body coupled angular momentum. If in you system there are three body coupling, this range may give wrong result.

### `"Jmax"`

`"Jmax"` means the maximum angular momentum, for every parameters. For another example,

```cpp
WignerSymbols wigner;
int Jmax = 21;
wigner.reserve(Jmax, "Jmax", 6);
```

This means in this system, `Jmax = 21`, and we calculate to 6j symbols. Here
```cpp
nmax = 4*Jmax + 1 = 85
```
So `reserve` function will extents the stored binominals range. For quantum many body calculation with only two body coupling, if single particle `jmax = 21/2`, then `Jmax = 21`. The `"djmax"` mode will use less storage.

In the `"Jmax"` mode, it is always safe even having three body coupling. You just need to use the maximum three body coupled angular momentum as `Jmax`, although it will cost more memory.

Actually, the memory used for store the binominals is not very large. A simple estimate is
```cpp
2 * nmax * nmax // Byte
```
Even for `nmax = 1000` (`Jmax = 200` in 9j calculation, which is absolutly enough for most calculation), the memory cost is just 2MB, which is not very large.

### `"nmax"`

The `"nmax"` mode directly set `nmax`, and the `rand` parameter is ignored. This maybe useful when you only want to calculate `binominal`s using this package.

### Thread safety

The `reserve` is **not** thread safe. So you shuld not call `reserve` function dymanically in a multi-threading program. The correct way to use this package is find the maximum angular momentum quantum number in you system, and call `reserve` at the beginning of the code, and then don't call it any more.
