# Arageli

This is a fork of [Arageli project](http://www.arageli.org) imported from SourceForge SVN.

Arageli is C++ library for computations in **ar**ithmetic, **a**lgebra, **ge**ometry, **l**inear and **i**nteger linear programming. Arageli provides routines supporting precise, i.e. symbolic or algebraic, computations. It contains definitions of basic algebraic structures such as integer numbers with arbitrary precision, rational numbers, vectors, matrices, polynomials etc.

See project description in original [README](/README)


## Goals of this fork

1. Compile library by recent versions of GCC and Clang  
1. Improve perfomance of arbitrary precision computations: classes *big_int* and *big_float*   

Supporting MSVC *is not* a goal. However, it will be great if there is a volunteer who can with MSVC.

## Building with CMake (Linux)

###  Setting up the build

Create a directory for a build

```
mkdir mybuild
cd mybuild
```

Optionaly override default compiler 

```
export CC=/usr/bin/clang-15
export CXX=/usr/bin/clang++
```

Prepare makefiles 

```
cmake ..
```

### Building the library
``` 
make
```

### Running tests
```
make runtests
```

### Installing (may need root access) 
```    
make install
```

### Building Hello, World application

hello_arageli.cpp:

```
#include <iostream>
#include <arageli/arageli.hpp>

int main()
{
    using Arageli::big_int;

    big_int a = "222222222222222222222222222222";
    big_int b = "111111111111111111111111111111";

    std::cout << a  << " + " << b << " = " << a + b << "\n";
    std::cout << a  << " - " << b << " = " << a - b << "\n";
    std::cout << a  << " * " << b << " = " << a * b << "\n";
    std::cout << a  << " / " << b << " = " << a / b << "\n";
}
```

```
c++ -O2 hello_arageli.cpp -o hello_arageli -larageli
./hello_arageli
```
