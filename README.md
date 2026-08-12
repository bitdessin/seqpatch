[![build](https://github.com/bitdessin/seqpatch/actions/workflows/build_pkg.yaml/badge.svg?style=flat-square)](https://github.com/bitdessin/seqpatch/actions/workflows/build_pkg.yaml)

# seqpatch

`r pkg("seqpatch")` is a compact R toolkit for analyzing bulk RNA-seq count data.
A core feature is normalization based on the DEGES framework,
the differentially expressed gene (DEG) elimination strategy [@TbT; @TCC].
DEGES iteratively identifies candidate DEGs and reduces their influence
when estimating sample specific scaling factors.
This approach is intended for datasets affected by compositional imbalance
or strongly asymmetric DEG patterns, for which standard normalization assumptions may not hold.
The package also provides visualization functions and utilities for simulating count data
to support methodological evaluation.


## Documentation

- https://bitdessin.github.io/seqpatch/


## Citation

When using the DEGES normalization workflow, please cite:

Sun J, Nishiyama T, Shimizu K, Kadota K.
TCC: an R package for comparing tag count data with robust normalization strategies.
BMC Bioinformatics. 14:219, 2013.
[10.1186/1471-2105-14-219](https://doi.org/10.1186/1471-2105-14-219)
