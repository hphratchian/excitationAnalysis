# excitationAnalysis

## Overview and Execution
This program carries out excitation analysis for Delta-SCF calculations. As input to the program, one provides two Gaussian-style matrix files that correspond to the converged initial and final electronic states. If the initial and final state matrix files are called `initial.mat` and `final.mat`, respectively, the program can be executed using the command

```
excitationAnalysis.exe initial.mat 001final.mat > log
```

where the file `log` will have the output of the program after its execution.

---

## Compilation
The repo includes a makefile that has been built for use with the fortran compilers included in the `nvhpc` package. The package must be available and properly loaded into the user's environment to compile the program.

To compile the program, simply execute `make` at the command line from within the repo's main directory.

---

## Exmample Calculations
The `GTest` sub-directory includes Gaussian input files that can be run to generate matrix files that can be used with excitationAnalysis.
