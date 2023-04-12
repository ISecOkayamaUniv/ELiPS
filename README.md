# ELiPS
This document describes how to make existing ELiPS library working in Linux environment.
This is expected that it will work any 64bit Unix distribution Ububtu.
Autotools installation may vary for Linux.
Please keep in mind that it is still in developing phase. 
If found any bug related to installation, please infrom in `daichi_hattori@s.okayama-u.ac.jp`

## Preparation
ELiPS needs `gmp` and `autoconf` and `libtoopl`.
In Ubuntu, it can install like following command.
```
sudo apt install libgmp-dev autoconf libtool
```

## Build
You can build using following step:
```
git clone https://github.com/ISecOkayamaUniv/ELiPS.git
cd ELiPS
autoreconf -i
./configure
make
sudo make install
```
Still there is no single header. Therefore please `/usr/local/include/ELiPS` directory to get the header declaration. 

If you face `cannot open shared object file: No such file or directory` while running then follow this steps:
```
sudo ldconfig
LD_LIBRARY_PATH=/usr/local/lib
```
Check `echo $LD_LIBRARY_PATH`. If path is set then run again.


## Example
`ELiPS/test/example` is example of ELiPS. This can runs following steps.
It confirms pairing on bls12.
```
gcc example.c -lgmp -lelips
./a.out
```

## Uninstall
```
sudo make uninstall
