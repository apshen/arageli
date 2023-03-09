# Arageli

This is a fork of [Arageli project](http://www.arageli.org) imported from SourceForge SVN.

See project description in original [README](/README)

## Building with cmake (Linux)

###  Setting up the build
```
mkdir mybuild
cd mybuild
cmake ..
```

### Building the library
``` 
make
```
    
### Installing (may need root access) 
```    
make install
```

### Running tests
```
make runtests
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
g++ -O2 hello_arageli.cpp -o hello_arageli -larageli
./hello_arageli
```
    


    
