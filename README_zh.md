# WignerSymbol

[English](README.md)

计算CG系数，Racah系数，Wigner 3j, 6j, 9j 系数。公式参考自[CGcoefficient.jl](https://github.com/0382/CGcoefficient.jl)。

## 使用方法

这是一个单头文件库，直接把文件复制到你的项目里面。

```cpp
using namespace jshl;
WignerSymbols wigner;
int djmax = 21;
wigner.reserve(djmax, "2bjmax", 6);
double x = wigner.f6j(dj1, dj2, dj3, dj4, dj5, dj6);
```

## 提供的函数
```cpp
// 二项式系数
double WignerSymbols::binomial(int n, int k);
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
// Wigner d函数 <j,m1|exp(i*beta*jy)|j,m2>
double WignerSymbols::dfunc(int dj, int dm1, int dm2, double beta);
```
其中，除了`binomial`函数之外，其余函数均使用真实角动量量子数的两倍作为参数，这是为了处理半整数角动量的情况。所以要计算`<10|1/2,1/2;1/2,-1/2>`这个CG系数，你需要调用的是
```cpp
WignerSymbols wigner;
double x = wigner.CG(1,1,2,1,-1,0);
```

## `reserve`函数

所有的系数都是使用二项式系数`binomial`来计算的，这个库会先存储下部分二项式系数，之后用到的时候直接取。所以这个库调用的`binomial`函数只有参数的范围在存储范围之内才获得正确的结果，否则直接返回零。

在`WignerSymbols`对象初始化时，会先存下`binomial(0, 0)`到`binomial(67, 33)`所有的二项式系数。你可以使用`reserve`函数扩充这个范围。`reserve`函数签名为
```cpp
void WignerSymbols::reserve(int num, std::string type, int rank)
```
其含义如下

|                       |    计算范围    |   CG & 3j   | 6j & Racah  |     9j      |
| :-------------------: | :------------: | :---------: | :---------: | :---------: |
|       `type`的意义    | `type`\\`rank` |      3      |      6      |      9      |
|      最大角动量模式     |    `"Jmax"`    | `3*Jmax+1`  | `4*Jmax+1`  | `5*Jmax+1`  |
|    最大两体角动量模式    |   `"2bjmax"`    | `2*jmax+1` | `3*jmax+1` | `4*jmax+1` |
|    最大二项式系数模式    |    `"nmax"`    |   `nmax`    |   `namx`    |   `nmax`    |

表格中的数据是最极端条件下保证不溢出最少要存多少二项式系数。不过使用的时候你不需要记住这些数值，只需要根据你的程序中出现的最大角动量和计算范围调用一下`reserve`函数就好了。

### `"2bjmax"`

例如

```cpp
WignerSymbols wigner;
int djmax = 21;
wigner.reserve(djmax, "2bjmax", 6);
```

表示的意义是，我们体系中最大的单粒子轨道的角动量为`21/2`，于是最大可能的两体耦合角动量为`21`，同时代码中仅计算CG系数和6j系数，而不会计算9j系数。在这个例子中，保存的`binomial(n, k)`系数中，最大的`n`为

```cpp
nmax = 3*djmax + 1 = 3 * 21 + 1 = 64
```

由于`64`小于初始化时预存的`67`，所以`reserve`函数不会做任何事情。

`"2bjmax"`意味着你的程序仅需要计算两体耦合而不需要计算三体耦合或更高的耦合。这个模式假定在所有的Wigner系数计算中，至少有一个角动量是单粒子角动量，也即在这个例子中不超过`21/2`。在这个假定下，`"2bjmax"`模式能够少存一些二项式系数，占用内存更少。

### `"Jmax"`

`"Jmax"`最大角动量模式则不区分角动量的来源，表示全局的最大角动量。比如如果上述例子改成

```cpp
WignerSymbols wigner;
int Jmax = 21;
wigner.reserve(Jmax, "Jmax", 6);
```

表示我们体系中最大可能的角动量为`Jmax = 21`，只计算到6j系数。那么对应的

```cpp
nmax = 4*Jmax + 1 = 85
```

此时将会增加存储的二项式系数范围。对于量子多体计算来说，如果单粒子轨道最大角动量为`21/2`，那么两体耦合的最大角动量为`21`。采用`"2bjmax"`模式能够更加节省内存。

`"Jmax"`模式没有任何假定，它总是安全的。即使考虑到三体耦合，只要使用三体的总`Jmax`来`reserve`，就能够保证不会溢出，尽管这可能会浪费一些空间。

不过实际上存下二项式系数占用的内存并不多，简单的估算其占用内存约为
```cpp
2 * nmax * nmax // Byte
```
即使`nmax = 1000`（对应计算到9j系数，`Jmax = 200`，对大多数计算绰绰有余），需要的内存也只有2MB，相对于科学计算的其他内存占用来说实在不算什么。

### `"nmax"`

`"nmax"`模式就是直接设置`nmax`，此时`rank`参数不起作用。可能只有你不想算Wigner系数而只想用这个库来计算二项式系数时会有用。

### 线程安全

注意：`reserve`函数**不是**线程安全的，如果你的程序是并行的，不要动态地调用`reserve`函数。正确是使用本库的方法是，先计算出体系最大角动量，然后在程序开始时调用一次`reserve`函数，之后就不应该继续调用这个函数了。

不过，你可以为每个线程创建一个`WignerSymbols`对象。由于存储`binomial`系数所占的内存实在很小，我没有用静态成员存储它们，这也意味着各线程之间的对象是完全互不干扰的。同样由于占用内存其实很小，所以同时创建好几个`WignerSymbols`并没有什么代价，是完全可以的。