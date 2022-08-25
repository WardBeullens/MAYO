# MAYO
Implementation of the MAYO signature scheme with the parameter set n=62,m=60,o=6,k=10,q=31. 

For the description of the MAYO scheme, read the paper "MAYO: Practical Post-Quantum Signatures from Oil-and-Vinegar Maps" https://ia.cr/2021/1144

to run the benchmarks: 

```
  git submodule update --init
  cd MAYO
  make test
  ./test
```
