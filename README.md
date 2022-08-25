# MAYO
Implementation of the MAYO signature scheme with the parameter set n=62,m=60,o=6,k=10,q=31. 

The E_{ij} matrices are chosen as multiplication by 1,x,x^2,...,x^54 in the extension field F_31[x]/(x^60 + x^2 + 3x+ 27).

For the description of the MAYO scheme, read the paper "MAYO: Practical Post-Quantum Signatures from Oil-and-Vinegar Maps" https://ia.cr/2021/1144

to run the benchmarks: 

```
  git submodule update --init
  cd MAYO
  make test
  ./test
```

NOTE: the key and signature sizes are a bit larger than advertized in the paper, because the keys are elements mod 31, and the implementation just stores one element per Byte. Compressing the keys and signatures to store 8 field elements per 5 Bytes should be straightforward, but hasn't been done yet. 
