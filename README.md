# Lattice reformulations

This repository provides code to reproduce the experiments of [Aardal et al. 2023](https://www.sciencedirect.com/science/article/pii/S0167637723000652?ref=pdf_download&fr=RR-2&rr=7e5876097f890e78).

## Requirements
For the reformulation step:
* [NTL](https://libntl.org/)

For running the experiments:
* [SCIP](https://scipopt.org/#scipoptsuite)
* A Python environment with [PySCIPOpt](https://github.com/scipopt/PySCIPOpt)

## Instructions
### Step 1: generate the (original) instances
Code to generate the multi-row instances can be found under `benchmarks/<name>`. For the single-row instances, we provide the already generated instances under the same directory. The procedure to generate them from scratch is as follows: first, generate the constraint coefficients by running `benchmarks/single_row/generate_coefficients.py` (with the appropriate parameters). Then use the code provided in [this other repo](https://github.com/lascavana/FrobeniusNum) to generate the right-hand-side.

### Step 2: reformulate the instances
Run `sh 01_reduce.sh` to compile the reformulation code and generate the reformulations of all instances. Make sure that the environment variable `SCIPOPTDIR` is correctly set to your SCIP installation path. Other than the instances, this process generates log files that are stored in the `log/` directory.

### Step 3: solve
To reproduce the experiments, run
```
python 02_solve.py <name>
```
where `<name>` is in `{struct_s, struct_b, nostruct_s, nostruct_b, ms, gap, ca}`. This will generate a csv file in the `results` directory.

### Step 4: analysis
To analyze the performance of the different formulations, run
```
python 03_analyze.py <name> <output>
```
with `<output>` being either `per_instance`, `aggregated` or `boxplot`. The average and maximum reformulation time for the multi-row instances can also be obtained by running
```
python 03_analyze.py <name> reftime
```
