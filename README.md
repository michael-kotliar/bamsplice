# geep
## Installation
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make .
```

### Note

if after `make install`
you got `error while loading shared libraries: libbamtools.so.2.4.1: cannot open shared object file: No such file or directory`,
make sure that `LD_LIBRARY_PATH` includes path to `usr/local/lib/bamtools`
