# WignerSymbol

[![test](https://github.com/0382/WignerSymbol/actions/workflows/test.yml/badge.svg)](https://github.com/0382/WignerSymbol/actions/workflows/test.yml)

[English](README.md)

计算CG系数，Racah系数，Wigner 3j, 6j, 9j 系数，Moshinsky括号等。其中一些公式请看[CGcoefficient.jl](https://0382.github.io/CGcoefficient.jl/stable/wigner/)。

## 使用方法

这是一个单头文件库，直接把文件复制到你的项目里面。

```cpp
#include "WignerSymbol.hpp"
using namespace util;

int djmax = 21;
wigner_init(djmax, "2bjmax", 6);
double x = wigner_6j(dj1, dj2, dj3, dj4, dj5, dj6);
```

### 不足

由于使用了浮点数计算，这个包在较大的角动量时会给出错误的结果。
详细的误差和性能测试，请看：[wigner-benchmark](https://github.com/0382/wigner-benchmark)。

## 提供的函数
```cpp
// 预计算二项式系数表
void wigner_init(int num, std::string type, int rank);
// 快速访问二项式系数表，在`n`很大时它可能失效（返回零）
double fast_binomial(int n, int k);
// CG系数
double CG(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
// 两个 1/2 自旋的CG系数
double CGspin(int ds1, int ds2, int S);
// 三个 1/2 自旋两次耦合 <S12,M12|1/2,m1;1/2,m2><S,M|S12,M12;1/2,m3>
double CG3spin(int dm1, int dm2, int dm3, int S12, int dS);
// CG 系数特殊情况 m1 == m2 == m3 == 0
double CG0(int j1, int j2, int j3);
// 3j系数
double wigner_3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
// 6j系数
double wigner_6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
// Racah系数
double Racah(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);
// 9j系数
double wigner_9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9);
// normalized 9j系数
double wigner_norm9j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6, int dj7, int dj8, int dj9);
// LS耦合到jj耦合的转换系数
double lsjj(int l1, int l2, int dj1, int dj2, int L, int S, int J);
// Wigner d函数 <j,m1|exp(i*beta*jy)|j,m2>
double dfunc(int dj, int dm1, int dm2, double beta);
// Moshinsky 括号，参考: Buck et al. Nuc. Phys. A 600 (1996) 387-402
double Moshinsky(int N, int L, int n, int l, int n1, int l1, int n2, int l2, int lambda, double tan_beta = 1.0);
```

角动量量子数可以是半整数，通常的库会使用实际量子数的两倍作为函数参数。本函数库也采取类似的做法。不过，考虑到有一些函数，比如`Moshinsky`系数只和轨道量子数有关，不存在半整数问题，没必要多此一举用两倍参数。

因此，我们约定，参数名以`d`开头，表示实际量子数的两倍，否者仅表示实际量子数。

所以要计算`<10|1/2,1/2;1/2,-1/2>`这个CG系数，你可以这样计算
```cpp
WignerSymbols wigner;
double x = CG(1, 1, 2, 1, -1, 0);
double y = CGspin(1, -1, 1); // 注意最后一个参数是`S`而不是`dS`，因为在此情况下，`S`只能是`0,1`没必要用两倍参数
```

## `wigner_init`函数

所有的系数都是使用二项式系数`binomial`来计算的，这个库会先存储下部分二项式系数，之后用到的时候直接取。所以这个库调用的`binomial`函数只有参数的范围在存储范围之内才获得正确的结果，否则直接返回零。

二项式系数表存储在一个`WignerSymbols`类中。我们定义了一个该类型的全局变量`wigner`来服务所有的函数。（`inline`全局变量需要`c++17`的特性。）

```cpp
inline WignerSymbols wigner;
inline void wigner_init(int num, std::string type, int rank) { wigner.reserve(num, type, rank); }
```

默认构造函数会存下`binomial(0, 0)`到`binomial(67, 33)`所有的二项式系数（这些正好都可以被`uint64_t`精确表示）。你可以使用`wigner_init`函数扩充这个范围。

`wigner_init`实际上做的事情是扩充二项式表，不过为了更明确如何设置对于后面的计算是安全的，我们定义了几种设置模式。

### `"2bjmax"`

```cpp
wigner_init(21, "2bjmax", 6);
```

表示的意义是，体系中最大的单粒子轨道的角动量为`21/2`，于是最大可能的两体耦合角动量为`21`，同时代码中仅计算CG系数和6j系数，而不会计算9j系数。

`"2bjmax"`意味着你的程序仅需要计算两体耦合而不需要计算三体耦合或更高的耦合。这个模式假定在所有的Wigner系数计算中，至少有一个角动量是单粒子角动量，也即在这个例子中不超过`21/2`。在这个假定下，`"2bjmax"`模式能够少存一些二项式系数，占用内存更少。

### `"Jmax"`

`"Jmax"`最大角动量模式则不区分角动量的来源，表示全局的最大角动量。比如如果上述例子改成

```cpp
wigner_init(21, "Jmax", 9);
```

表示我们体系中最大可能的角动量为`Jmax = 21`，计算到9j系数。

`"Jmax"`模式没有任何假定，它总是安全的。即使考虑到三体耦合，只要使用三体的总`Jmax`来`wigner_init`，就能够保证不会溢出，尽管这可能会浪费一些空间。

不过实际上存下二项式系数占用的内存并不多，算到9j系数，`Jmax = 200`，需要的内存也只有2MB，相对于科学计算的其他内存占用来说实在不算什么。

### `"nmax"`

`"nmax"`模式就是直接设置`nmax`，此时`rank`参数不起作用。可能只有你不想算Wigner系数而只想用这个库来计算二项式系数时会有用。

### 实际的`nmax`

下表给出各种情况下实际设置的二项式表的`nmax`。参考[Estimate-the-capacity](https://0382.github.io/CGcoefficient.jl/stable/formula/#Estimate-the-capacity)查看详细推导过程。

|                    |    计算范围    |  CG & 3j   | 6j & Racah |     9j     |
| :----------------: | :------------: | :--------: | :--------: | :--------: |
|    `type`的意义    | `type`\\`rank` |     3      |     6      |     9      |
|   最大角动量模式   |    `"Jmax"`    | `3*Jmax+1` | `4*Jmax+1` | `5*Jmax+1` |
| 最大两体角动量模式 |   `"2bjmax"`   | `2*jmax+1` | `3*jmax+1` | `4*jmax+1` |
| 最大二项式系数模式 |    `"nmax"`    |   `nmax`   |   `namx`   |   `nmax`   |


### 线程安全

注意：`wigner_init`函数**不是**线程安全的，如果你的程序是并行的，不要动态地调用`wigner_init`函数。正确是使用本库的方法是，先计算出体系最大角动量，然后在程序开始时调用一次`wigner_init`函数，之后就不应该继续调用这个函数了。其他函数只是读取二项式系数表，并行地调用完全没有问题。

## 参考资料

1. T. Engeland and M. Hjorth-Jensen, the Oslo-FCI code. https://github.com/ManyBodyPhysics/CENS.
2. A. N. Moskalev D. A. Varshalovich and V. K. Khersonskii, *Quantum theory of angular momentum*.