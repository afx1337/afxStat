# afxStat

**afxStat** is a MATLAB toolbox for permutation-based statistical inference in neuroimaging.
It provides flexible random permutation testing for voxel-wise and cluster-wise inference within the general linear model (GLM) framework, with support for nuisance covariates via the Freedman–Lane procedure.

The toolbox is designed for researchers who want a lightweight, scriptable alternative for permutation-based neuroimaging statistics and lesion-symptom mapping workflows.

---

## Features

* Random permutation testing for GLM-based neuroimaging analyses
* Voxel-wise and cluster-wise family-wise error (FWE) corrected inference
* Freedman–Lane permutation procedure for designs with nuisance covariates
* Support for binary and continuous outcome variables
* Script-based workflows for reproducible analyses
* Optional GUI for users preferring interactive setup

---

## Methodological Background

Permutation inference is implemented using nonparametric randomization / permutation testing in the GLM framework.

For methodological details, see:

* Nichols TE, Holmes AP (2002). *Nonparametric permutation tests for functional neuroimaging: a primer with examples.* Human Brain Mapping.
* Winkler AM et al. (2014). *Permutation inference for the general linear model.* NeuroImage.

---

## Requirements

* MATLAB **R2017b or newer** (tested)
* **SPM12**
* Tested on:

  * Linux
  * Windows

> Octave compatibility has not been tested.

---

## Installation

Clone the repository and add the `scripts` folder to your MATLAB path:

```matlab
addpath('scripts')
```

Make sure SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) is available on your MATLAB path before running analyses.

---

## Repository Structure

```text
afxStat/
├── doc/        Documentation / manual
├── masks/      Example analysis masks
├── results/    Default output directory
└── scripts/    MATLAB source code
```

### Included Masks

* `brainmask.nii` – Whole-brain analysis mask
* `gmmask_20.nii` – Gray matter mask (SPM12 GM TPM > 20%)
* `gmmask_20_ext.nii` – Extended gray matter mask including subcortical / infratentorial regions

---

## Quick Start

The recommended high-level entry point for most analyses is:

```matlab
destFolder = afxStatFiles( ...
    imgFiles, ...
    [], ...                    % flipLR
    [], ...                    % FWHM smoothing
    [], ...                    % threshold / binarization
    X, ...
    {[1 -1]}, ...             % contrast(s)
    'masks/brainmask.nii', ...
    5000, ...
    'voxel', ...
    true, ...
    0.001, ...
    0.05, ...
    'results/example');
```

---

## Core Workflow

### Input Arguments (`afxStatFiles`)

| Argument      | Description                                       |
| ------------- | ------------------------------------------------- |
| `imgFiles`    | Cell array of input image filenames               |
| `flipLR`      | Optional left/right flip vector                   |
| `FWHM`        | Optional smoothing kernel                         |
| `thr`         | Optional threshold / binarization value           |
| `X`           | Design matrix                                     |
| `contrasts`   | Cell array of contrast vectors                    |
| `maskFile`    | Analysis mask                                     |
| `nPerms`      | Number of permutations                            |
| `inference`   | `'voxel'` or `'cluster'`                          |
| `FWE`         | Apply family-wise error correction (`true/false`) |
| `threshVox`   | Voxel-level threshold                             |
| `threshClust` | Cluster-level threshold                           |
| `destFolder`  | Output directory                                  |

---

### Output Files

Each analysis generates output in the specified destination folder (or in `results/` by default).

| File                    | Description                                                |
| ----------------------- | ---------------------------------------------------------- |
| `TMap_XXX.nii`          | Raw statistical map (t-values)                             |
| `TMap_XXX_filtered.nii` | Thresholded significant result after permutation inference |
| `mask.nii`              | Analysis mask used                                         |
| `sumMap.nii`            | Input overlap / summary map (VLSM only)                    |
| `info_XXX.mat`          | MATLAB metadata / analysis information                     |
| `info_XXX.txt`          | Human-readable analysis summary                            |
| `info_XXX.json`         | Machine-readable analysis summary                          |

---

## Other Main Interfaces

### `afxStatXlsx`

Excel-based interface for analyses defined via design spreadsheets.

Convenient for workflows where design matrices and file mappings are maintained in tabular form.

---

### `afxVlsmFiles` / `afxVlsmXlsx`

Specialized wrappers for voxel-based lesion-symptom mapping (VLSM).

Includes lesion-analysis specific preprocessing such as overlap thresholding and lesion volume regression.

---

### `afxSDSMFiles`

Wrapper for structural disconnection-symptom mapping workflows.

---

### `afxStat`

Low-level interface operating directly on matrix data (`Y`, `X`) rather than image files.

Useful when preprocessing / masking has already been performed externally.

---

## Reference

```matlab
afxStatFiles(imgFiles, flipLR, FWHM, thr, X, contrasts, maskFile, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment)
afxStatXlsx(designFile, contrasts, maskFile, nPerms, inference, FWE, threshVox, threshClust)
afxVlsmFiles(imgFiles, flipLR, FWHM, X, minOverlap, regressLesion, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment)
afxVlsmXlsx(designFile, minOverlap, regressLesion, nPerms, inference, FWE, threshVox, threshClust)
afxSDSMFiles(imgNetwork, flipLR, FWHM, thr, imgLesion, X, minOverlap, contrasts, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment)
afxStat(Y, X, contrasts, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment)
```

---

## Methodological Recommendations

Recommended defaults / best practices:

* **≥ 5000 permutations**
* **Enable Freedman–Lane** whenever nuisance covariates are present

  * In afxStat this is done by passing a **negative** permutation count
    (e.g. `nPerms = -5000`)
* **Use FWE correction**
* **Voxel-wise inference** is standard in many lesion-symptom mapping applications
* **Cluster-wise inference** may provide increased sensitivity but should be interpreted carefully
* For lesion analyses, a **minimum overlap of 5 subjects** is commonly recommended

---

## Freedman–Lane Permutation Procedure

When nuisance regressors / covariates are present in the design matrix, standard row-wise permutation is generally invalid.

afxStat supports the **Freedman–Lane procedure** for valid permutation inference in such designs.

Activate it by specifying a **negative** number of permutations:

```matlab
nPerms = -5000;
```

See Winkler et al. (2014) for details.

---

## GUI Usage

Although afxStat is primarily intended for script-based workflows, a basic graphical user interface is available:

```matlab
afxStatGUI
```

The GUI provides interactive access to common workflows but exposes only a subset of the toolbox’s flexibility.

Additional GUI documentation can be found in the `/doc` folder (in German).

---

## Citation

If you use this toolbox in scientific work, please cite the repository:

```text
afxStat. Available at: https://github.com/afx1337/afxStat
```

Please also cite the methodological references listed above where appropriate.

---

## Support / Contributions

Bug reports, feature requests, and contributions are welcome via GitHub Issues / Pull Requests.

---

## Documentation Note

Parts of this documentation were drafted with the assistance of a large language model and reviewed by the authors.
