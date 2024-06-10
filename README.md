
[![ci](https://github.com/kmlefran/compchem-si/actions/workflows/ci.yml/badge.svg)](https://github.com/kmlefran/compchem-si/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/compchem-si/badge/?version=latest)](https://compchem-si.readthedocs.io/en/latest/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/github/kmlefran/compchem-si/badge.svg?branch=main)](https://coveralls.io/github/kmlefran/compchem-si?branch=main)
[![PyPI version](https://badge.fury.io/py/compchem-si.svg)](https://badge.fury.io/py/compchem-si)
# Molecule Reorientation and .sum file handling

This package provides utilities to take Gaussian output files and generate a page containing information on the molecule. The page currently contains a visual representation of the molecule, its associated energy in Hartree, and its molecular geometry.

# Authors
Kevin Lefrancois-Gagnon
Robert C. Mawhinney

# Installation
```
pip install compchem-si
```

# Example Usage:
If you are in a working directory with all of your .log files that you want in the SI, you can run the following to generate the SI as SupplementaryInformation.pdf
```python
constuct_si()
```

Alternatively, if you only want some of the files to be included in the SI, run:
```python
constuct_si(log_file_list = ['log1.log', 'log2.log',...],out_name="my_si_name.pdf")
```

# Features to Develop
1. Handle non-optimization jobs
2. Include more features in SI file
3. Optimize formatting of SI

# Note
This is an initial release, the SI page format is pretty barebones for now, but will be improved over time
