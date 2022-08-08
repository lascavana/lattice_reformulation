# Kernel-lattice reformulation

Generates a reformulation of a MILP based on the [AHL method](https://pubsonline.informs.org/doi/pdf/10.1287/moor.25.3.427.12219).

## Usage

Modify file `01_reduce.sh` to compile and run the reformulation code. This will generate a new lp file with name `ahl_<inputfilename>.lp`.

## Requirements

* [SCIP](https://scipopt.org/#scipoptsuite)
* [NTL](https://libntl.org/)
