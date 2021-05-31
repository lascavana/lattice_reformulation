### Requirements
[NTL library](https://cs.uwaterloo.ca/~echrzano/tour-unix.html)

### Instructions
run
```shell
g++ -g -O2 -std=c++11 -pthread reduce.cpp -o reduce -lntl -lgmp -lm
```
then
```shell
./reduce <inputfile> <outputfile>
```

this will generate three lp files:
```shell
<outputfile>_orig.lp
<outputfile>_ahl.lp
<outputfile>_cea.lp
```