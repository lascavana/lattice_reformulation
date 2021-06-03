### Requirements
[NTL library](https://cs.uwaterloo.ca/~echrzano/tour-unix.html)

### Input
The input file must be a ```.txt``` file in the following form:
```shell
n m
A -b
c
u
```
where n is the number of variables, m is the number of constraints, A is the constraint matrix (m x n), b is the right-hand side (m), c is the objective function (n), and u is the vector containing the variable upper bounds (n).

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

### Optional arguments
```shell
./reduce <inputfile> <outputfile> --translate
```
Will print text files ```translate_W.txt``` and ```translate_x0.txt``` containing the necessary information to translate between the original basis and the AHL basis.
 

#### Final note
If NTL is not installed in the standard directory (```/usr/local```), but rather in some directory ```INSTALL_DIR``` you need to specify so by instead running
```shell
g++ -g -O2 -std=c++11  -I$INSTALL_DIR/include -L$INSTALL_DIR/lib -pthread reduce.cpp -o reduce -lntl -lgmp -lm
```