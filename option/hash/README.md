# how to use hash in ELiPS

1. install libssl-dev.
```
sudo apt install libssl-dev
```
2. compile with elips_sha512.h.
```
gcc ???.c -lgmp -lelips -lcrypto
```

Please refer to "hash_test.c".