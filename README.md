If you are using banana in a publication or presentation, please cite the following paper describing automatic bubble detection with banana: https://ui.adsabs.harvard.edu/abs/2021A%26A...653A..16C/abstract

This repository has been archived at zenodo: 

[![DOI](https://zenodo.org/badge/256713274.svg)](https://zenodo.org/badge/latestdoi/256713274)

# banana
Tools for creating Minkowski maps, astronomical bubble detection, and related functions for FITS-files

Irreducible Minkowski Tensors are described in https://morphometry.org/theory/anisotropy-analysis-by-imt/

This work is based on an older version of papaya2 (https://morphometry.org/software/papaya2/)

For further details, see doc.pdf

# clone and build

```
git clone https://github.com/ccollischon/banana.git
mkdir banana-build
cd banana-build
cmake ../banana
make
```

If a library (e.g. CCFits) is installed in a non-standard path you can add this path when calling cmake:

```
cmake -DCMAKE_CXX_FLAGS="-L /path/to/lib" ../banana
```
