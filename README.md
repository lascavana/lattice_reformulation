# Lattice reformulations

This repository provides code to reproduce the experiments of [Aardal et al. 2023](https://www.sciencedirect.com/science/article/pii/S0167637723000652?ref=pdf_download&fr=RR-2&rr=7e5876097f890e78).

## Step 1: generate the (original) instances
Code to generate the multi-row instances can be found under `benchmarks/<name>`. For the single-row instances, we provide the already generated instances under the same directory. The procedure to generate them from scratch is as follows: first, generate the constraint coefficients by running `benchmarks/single_row/generate_coefficients.py` (with the appropriate parameters). Then use the code provided in [this other repo](https://github.com/lascavana/FrobeniusNum) to generate the right-hand-side.

## Usage

Modify file `01_reduce.sh` to compile and run the reformulation code. This will generate a new lp file with name `ahl_<inputfilename>.lp`.

## Requirements

* [SCIP](https://scipopt.org/#scipoptsuite)
* [NTL](https://libntl.org/)
