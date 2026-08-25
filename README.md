<!--
<p align="center">
  <img src="https://github.com/pnnl-predictive-phenomics/crop/raw/main/docs/source/logo.png" height="150">
</p>
-->

<h1 align="center">
  CROP
</h1>

<p align="center">
    <a href="https://github.com/pnnl-predictive-phenomics/crop/actions/workflows/tests.yml">
        <img alt="Tests" src="https://github.com/pnnl-predictive-phenomics/crop/workflows/Tests/badge.svg" />
    </a>
    <a href="https://pypi.org/project/crop">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/crop" />
    </a>
    <a href="https://pypi.org/project/crop">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/crop" />
    </a>
    <a href="https://github.com/pnnl-predictive-phenomics/crop/blob/main/LICENSE">
        <img alt="PyPI - License" src="https://img.shields.io/pypi/l/crop" />
    </a>
    <a href='https://crop.readthedocs.io/en/latest/?badge=latest'>
        <img src='https://readthedocs.org/projects/crop/badge/?version=latest' alt='Documentation Status' />
    </a>
    <a href="https://codecov.io/gh/pnnl-predictive-phenomics/crop/branch/main">
        <img src="https://codecov.io/gh/pnnl-predictive-phenomics/crop/branch/main/graph/badge.svg" alt="Codecov status" />
    </a>  
    <a href="https://github.com/cthoyt/cookiecutter-python-package">
        <img alt="Cookiecutter template from @cthoyt" src="https://img.shields.io/badge/Cookiecutter-snekpack-blue" /> 
    </a>
    <a href='https://github.com/psf/black'>
        <img src='https://img.shields.io/badge/code%20style-black-000000.svg' alt='Code style: black' />
    </a>
    <a href="https://github.com/pnnl-predictive-phenomics/crop/blob/main/.github/CODE_OF_CONDUCT.md">
        <img src="https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg" alt="Contributor Covenant"/>
    </a>
</p>

Consistent Reproduction of Phenotype (CROP) is an mixed integer linear programming (MILP) algorithm for finding reactions to remove from a model to improve the prediction of no-growth phenotypes while ensuring that correctly predicted growth phenotypes are still preserved. 

## 💪 Getting Started

See the [CROP Notebook](https://github.com/pnnl-predictive-phenomics/crop/blob/main/docs/notebook/CROP.ipynb).

## Verify an Applied Correction

After applying CROP's suggested reaction removals to a model, verify every observed
phenotype with fresh flux balance analysis:

```bash
crop verify \
    --model corrected-model.json \
    --phenotypes phenotypes.json \
    --media media.json
```

The corrected model may be COBRA JSON (`.json`) or SBML (`.xml` or `.sbml`). The
phenotype file is keyed by condition:

```json
{
    "glucose": {"observed": "growth", "predicted": "growth"},
    "lactose": {"observed": "no_growth", "predicted": "growth"}
}
```

The media file uses COBRA's convention in which negative lower bounds allow uptake:

```json
{
    "glucose": {"EX_glc": -10.0, "EX_lac": 0.0},
    "lactose": {"EX_glc": 0.0, "EX_lac": -10.0}
}
```

By default, observed growth passes when biomass is at least `2.0`, and observed
no-growth passes when biomass is at most `1.0` or FBA is infeasible. Configure these
limits with `--minimum-growth` and `--maximum-nogrowth`. Use `--biomass-rxn` to
override the model's objective reaction.

For automation, request JSON and optionally write it to a file:

```bash
crop verify \
    --model corrected-model.xml \
    --phenotypes phenotypes.json \
    --media media.json \
    --format json \
    --output verification.json
```

The command exits with status `0` when all observations are reproduced, `1` when
verification completes with phenotype failures, and `2` for invalid inputs.

## Compare Single- and Multi-Phenotype Reconciliation

Generate a deterministic toy comparison of independent one-phenotype-at-a-time
reconciliation and CROP's joint multi-phenotype constraints:

```bash
uv run python scripts/compare_growmatch_crop.py --output-dir comparison_output
```

The independent arm is GrowMatch-style: each false-growth phenotype is solved with
only its local growth control, and the two deletion sets are then combined. It uses
CROP's optimizer to isolate the effect of condition scope; it does not invoke the
separate GrowMatch implementation or claim numerical parity with its genome-scale
model.

The script applies every suggested deletion set and reruns FBA under all four media.
It writes a biomass-flux plot, a phenotype match matrix, the underlying CSV, and the
reaction-removal sets as JSON. In this cross-coupled model, independent fixes suppress
both false-growth phenotypes but clobber both omitted growth controls. Joint CROP
preserves both growth phenotypes while suppressing both false-growth phenotypes.

## 🚀 Installation

<!-- Uncomment this section after your first ``tox -e finish``
The most recent release can be installed from
[PyPI](https://pypi.org/project/crop/) with:

```shell
$ pip install crop
```
-->

The most recent code and data can be installed directly from GitHub with:

```bash
$ pip install git+https://github.com/pnnl-predictive-phenomics/crop.git
```

## 👐 Contributing

Contributions, whether filing an issue, making a pull request, or forking, are appreciated. See
[CONTRIBUTING.md](https://github.com/pnnl-predictive-phenomics/crop/blob/master/.github/CONTRIBUTING.md) for more information on getting involved.

## 👋 Attribution

### ⚖️ License

The code in this package is licensed under the MIT License.

<!--
### 📖 Citation

Citation goes here!
-->

<!--
### 🎁 Support

This project has been supported by the following organizations (in alphabetical order):

- [Harvard Program in Therapeutic Science - Laboratory of Systems Pharmacology](https://hits.harvard.edu/the-program/laboratory-of-systems-pharmacology/)

-->

<!--
### 💰 Funding

This project has been supported by the following grants:

| Funding Body                                             | Program                                                                                                                       | Grant           |
|----------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------|-----------------|
| DARPA                                                    | [Automating Scientific Knowledge Extraction (ASKE)](https://www.darpa.mil/program/automating-scientific-knowledge-extraction) | HR00111990009   |
-->

### 🍪 Cookiecutter

This package was created with [@audreyfeldroy](https://github.com/audreyfeldroy)'s
[cookiecutter](https://github.com/cookiecutter/cookiecutter) package using [@cthoyt](https://github.com/cthoyt)'s
[cookiecutter-snekpack](https://github.com/cthoyt/cookiecutter-snekpack) template.

## 🛠️ For Developers

<details>
  <summary>See developer instructions</summary>

The final section of the README is for if you want to get involved by making a code contribution.

### Development Installation

To install in development mode, use the following:

```bash
$ git clone git+https://github.com/pnnl-predictive-phenomics/crop.git
$ cd crop
$ pip install -e .
```

### 🥼 Testing

After cloning the repository and installing `tox` with `pip install tox`, the unit tests in the `tests/` folder can be
run reproducibly with:

```shell
$ tox
```

Additionally, these tests are automatically re-run with each commit in a [GitHub Action](https://github.com/pnnl-predictive-phenomics/crop/actions?query=workflow%3ATests).

### 📖 Building the Documentation

The documentation can be built locally using the following:

```shell
$ git clone git+https://github.com/pnnl-predictive-phenomics/crop.git
$ cd crop
$ tox -e docs
$ open docs/build/html/index.html
``` 

The documentation automatically installs the package as well as the `docs`
extra specified in the [`setup.cfg`](setup.cfg). `sphinx` plugins
like `texext` can be added there. Additionally, they need to be added to the
`extensions` list in [`docs/source/conf.py`](docs/source/conf.py).

### 📦 Making a Release

After installing the package in development mode and installing
`tox` with `pip install tox`, the commands for making a new release are contained within the `finish` environment
in `tox.ini`. Run the following from the shell:

```shell
$ tox -e finish
```

This script does the following:

1. Uses [Bump2Version](https://github.com/c4urself/bump2version) to switch the version number in the `setup.cfg`,
   `src/crop/version.py`, and [`docs/source/conf.py`](docs/source/conf.py) to not have the `-dev` suffix
2. Packages the code in both a tar archive and a wheel using [`build`](https://github.com/pypa/build)
3. Uploads to PyPI using [`twine`](https://github.com/pypa/twine). Be sure to have a `.pypirc` file configured to avoid the need for manual input at this
   step
4. Push to GitHub. You'll need to make a release going with the commit where the version was bumped.
5. Bump the version to the next patch. If you made big changes and want to bump the version by minor, you can
   use `tox -e bumpversion -- minor` after.
</details>
