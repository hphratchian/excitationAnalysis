# excitationAnalysis
This program carries out excitation analysis for Delta-SCF calculations. As input to the program, one provides two Gaussian-style matrix files that correspond to the converged initial and final electronic states. If the initial and final state matrix files are called `initial.mat` and `final.mat`, respectively, the program can be executed using the command

```
excitationAnalysis.exe initial.mat 001final.mat > log
```

where the file `log` will have the output of the program after its execution.
