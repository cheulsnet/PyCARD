# PyCARD

PyCARD is a Python implementation of CARD for spatially informed cell-type deconvolution of spatial transcriptomics data.

PyCARD provides a Python-based workflow for:

* Reference-based CARD deconvolution
* Reference-free CARDfree deconvolution
* Single-cell resolution mapping (`scMapping`)
* Spatial resolution enhancement / imputation
* Visualization of CARD outputs
* C++-accelerated computation through `pybind11`
* Configurable CARD optimization parameters, including user-defined `phi` settings

This release has been validated in a clean Conda environment using Python 3.10.19.

---

## Requirements

### Tested environment

The current release has been tested with:

* Linux / WSL Ubuntu
* Python 3.10.19
* GNU C/C++ compiler
* CMake
* Eigen3
* Armadillo
* R 4.1.2
* `pybind11` 2.11.1

R support is used by the current CARDfree implementation through `rpy2`, including the R-based NMF backend options.

---

## Installation

### 1. Create a Conda environment

```bash
conda create -n pycard python=3.10.19 pip -y
conda activate pycard
```

For a clean installation, make sure that packages from a user-level Python installation are not being injected into the environment:

```bash
unset PYTHONPATH
export PYTHONNOUSERSITE=1
```

You can verify the active Python installation with:

```bash
python --version
which python
```

The tested Python version is:

```text
Python 3.10.19
```

---

### 2. Install system dependencies

On Ubuntu / WSL Ubuntu:

```bash
sudo apt update
sudo apt install -y build-essential cmake libeigen3-dev libarmadillo-dev
```

Check that R is available:

```bash
R --version
```

The release-validation environment used R 4.1.2.

---

### 3. Install Python dependencies

From the PyCARD project root:

```bash
python -m pip install -r requirements.txt
```

The validated direct dependencies are:

```text
numpy==1.26.4
pandas==2.2.3
scipy==1.11.4
scikit-learn==1.3.2

anndata==0.10.9

matplotlib==3.8.4
seaborn==0.13.1

nimfa==1.4.0
mlpack==4.3.0.post2

geopandas==0.14.2
alphashape==1.3.1
shapely==2.0.2

pyreadr==0.5.0
rpy2==3.5.16

pybind11==2.11.1
```

---

## Build the C++ modules

PyCARD uses two C++ extension modules:

* `CARDref_module`
* `CARDfree_module`

Both modules are built using the single top-level `CMakeLists.txt`.

From the project root:

```bash
cmake -S . -B build \
  -Dpybind11_DIR="$(python -m pybind11 --cmakedir)"
```

Compile both modules:

```bash
cmake --build build -j
```

A successful build should generate files similar to:

```text
build/CARDref_module.cpython-310-x86_64-linux-gnu.so
build/CARDfree_module.cpython-310-x86_64-linux-gnu.so
```

Verify the reference-based module:

```bash
PYTHONPATH=./build python -c \
"import CARDref_module; print('CARDref OK:', CARDref_module.__file__)"
```

Verify the CARDfree module:

```bash
PYTHONPATH=./build python -c \
"import CARDfree_module; print('CARDfree OK:', CARDfree_module.__file__)"
```

---

## Input data

The example workflow uses the following input files:

```text
data/
├── spatial_count.csv
├── spatial_location.csv
├── sc_count.csv
├── sc_meta.csv
└── markerList.csv
```

The large example / benchmark datasets are intentionally excluded from this source repository.

Benchmark data used with CARD are available from the original CARD repository:

https://github.com/YMa-lab/CARD/tree/master/data

Place the required input files in a local `data/` directory.

By default, PyCARD reads data from:

```text
./data
```

A different input directory can be selected with:

```bash
PYTHONPATH=./build python main.py \
  --data-dir /path/to/data
```

---

# Quick Start

## Reference-based CARD deconvolution

Standard reference-based CARD deconvolution runs by default:

```bash
PYTHONPATH=./build python main.py
```

For additional result and profiling information:

```bash
PYTHONPATH=./build python main.py --verbose
```

The default deconvolution settings include:

```text
phi-grid = 0.9
max-iter = 200
epsilon = 1e-4
```

The `phi` values can be configured using:

```text
--phi-grid
```

For example:

```bash
PYTHONPATH=./build python main.py \
  --phi-grid 0.5,0.7,0.9 \
  --verbose
```

Warm-start behavior across `phi` values can optionally be enabled with:

```bash
PYTHONPATH=./build python main.py \
  --phi-grid 0.5,0.7,0.9 \
  --warm-start-phi \
  --verbose
```

---

## Single-cell resolution mapping

Single-cell resolution mapping is performed after standard CARD deconvolution.

```bash
PYTHONPATH=./build python main.py \
  --run-scmapping
```

Relevant parameters include:

```text
--scmapping-num-cell
--scmapping-ncore
```

Example:

```bash
PYTHONPATH=./build python main.py \
  --run-scmapping \
  --scmapping-num-cell 20 \
  --scmapping-ncore 10
```

---

## Spatial resolution enhancement / imputation

Spatial resolution enhancement / imputation is also performed after standard CARD deconvolution.

```bash
PYTHONPATH=./build python main.py \
  --run-imputation
```

Relevant parameters include:

```text
--imputation-num-grids
--imputation-in-neighbor
--imputation-show-grid
```

Example:

```bash
PYTHONPATH=./build python main.py \
  --run-imputation \
  --imputation-num-grids 2000 \
  --imputation-in-neighbor 10
```

---

## Reference-free CARDfree deconvolution

CARDfree can be enabled with:

```bash
PYTHONPATH=./build python main.py \
  --run-cardfree
```

For verbose output:

```bash
PYTHONPATH=./build python main.py \
  --run-cardfree \
  --verbose
```

In the current command-line workflow, standard reference-based CARD deconvolution runs by default; `--run-cardfree` additionally executes the CARDfree workflow.

The NMF backend for CARDfree is selected using:

```text
--cardfree-nmf-select
```

Available selections are:

```text
1  scikit-learn NMF
2  nimfa NMF
3  mlpack NMF
4  R NMF package
5  RcppML-based R NMF
```

The current default is:

```text
--cardfree-nmf-select 4
```

For R-based execution, an R library location can optionally be specified:

```bash
PYTHONPATH=./build python main.py \
  --run-cardfree \
  --set-r-libs-user ~/R/library
```

---

# Python API

The command-line workflow in `main.py` uses the same PyCARD functions that can also be called from Python.

## Load data

```python
from CARD.dataload import _load_data

spatial_count = _load_data("./data/spatial_count.csv", "csv")
spatial_location = _load_data("./data/spatial_location.csv", "csv")
sc_count = _load_data("./data/sc_count.csv", "csv")
sc_meta = _load_data("./data/sc_meta.csv", "csv")
```

The data loader supports multiple formats, including CSV, Feather, Excel, and RData where supported by the corresponding dependencies.

---

## Reference-based deconvolution

```python
from CARD.CARDrun import _run_CARD_deconvolution

CARD_obj = _run_CARD_deconvolution(
    sc_count=sc_count,
    sc_meta=sc_meta,
    spatial_count=spatial_count,
    spatial_location=spatial_location
)
```

---

## Single-cell mapping

```python
from CARD.CARDrun import _run_CARD_scMapping

scMapping = _run_CARD_scMapping(
    CARD_obj,
    shapeSpot="Square",
    numCell=20,
    ncore=10
)
```

---

## Imputation

```python
from CARD.CARDrun import _run_CARD_imputation

CARDimput_obj = _run_CARD_imputation(
    CARD_obj,
    num_grids=2000,
    in_neighbor=10,
    exclude=None,
    showGrid=True
)
```

---

## CARDfree

```python
from CARD.CARDrun import _run_CARDfree_deconvolution

CARDfree_obj = _run_CARDfree_deconvolution(
    spatial_count=spatial_count,
    spatial_location=spatial_location,
    nmfSelect=4
)
```

---

# Visualization

PyCARD provides visualization functionality through `_run_CARD_visualization`.

Available visualization modes include:

```text
visual_pie
visual_prop
visual_prop_2CT
visual_cor
visual_imput_loc
visual_gene
visual_scmapping
```

Example:

```python
from CARD.CARDrun import _run_CARD_visualization

visual_functions = [
    (
        "visual_pie",
        {
            "proportion": prop_data,
            "spatial_location": spatial_data
        }
    ),
    (
        "visual_prop",
        {
            "proportion": prop_data,
            "spatial_location": spatial_data,
            "ct_visualize": ct_data
        }
    )
]

_run_CARD_visualization(visual_functions)
```

Additional examples are available in `main.py` and the PyCARD source modules.

---

# Command-Line Options

Display the current command-line options with:

```bash
PYTHONPATH=./build python main.py --help
```

Principal options include:

```text
--data-dir
--set-r-libs-user

--run-deconvolution
--run-scmapping
--run-imputation
--run-cardfree

--phi-grid
--warm-start-phi
--max-iter
--epsilon

--scmapping-num-cell
--scmapping-ncore

--imputation-num-grids
--imputation-in-neighbor
--imputation-show-grid

--cardfree-nmf-select
--verbose
```

---

# Release Validation

The release was validated using a newly created Conda environment rather than relying on the original development environment.

The validation procedure included:

1. Creating a clean Conda environment using Python 3.10.19.
2. Removing access to user-level Python packages during validation.
3. Installing the Python dependencies from `requirements.txt`.
4. Configuring the project with the unified `CMakeLists.txt`.
5. Building both `CARDref_module` and `CARDfree_module`.
6. Successfully importing both C++ extension modules.
7. Running reference-based CARD deconvolution on the example dataset.
8. Running CARDfree deconvolution on the example dataset.
9. Confirming generation of cell-type proportion outputs from both workflows.

The clean validation environment successfully completed both reference-based CARD and CARDfree workflows.

---

# Project Structure

```text
PyCARD/
├── CARD/
│   ├── __init__.py
│   ├── CARDSCMapping.py
│   ├── CARDimputation.py
│   ├── CARDprop.py
│   ├── CARDrefFree.py
│   ├── CARDrun.py
│   ├── dataload.py
│   ├── pyNMFpackage.py
│   ├── utilities.py
│   └── visualization.py
│
├── src/
│   ├── CARDref.cpp
│   ├── CARDref.hpp
│   ├── CARDref_wrapper.cpp
│   ├── CARDfree.cpp
│   ├── CARDfree.hpp
│   └── CARDfree_wrapper.cpp
│
├── CMakeLists.txt
├── main.py
├── requirements.txt
├── .gitignore
├── README.md
└── LICENSE
```

Build outputs, compiled extension modules, debug artifacts, editor metadata, and large datasets are excluded from version control.

---

# Original CARD

The original CARD implementation and documentation are available at:

https://yma-lab.github.io/CARD/

Original CARD repository:

https://github.com/YMa-lab/CARD

Original CARD benchmark data:

https://github.com/YMa-lab/CARD/tree/master/data

---

# PyCARD Repository

https://github.com/faishim/PyCARD

---

# Citation

If you use PyCARD, please cite the original CARD publication:

Ying Ma and Xiang Zhou. Spatially informed cell-type deconvolution for spatial transcriptomics. *Nature Biotechnology* 40, 1349–1359 (2022).

Please also cite the PyCARD protocol publication once its final citation information is available.

---

# License

PyCARD is released under the MIT License.

See the `LICENSE` file for details.
