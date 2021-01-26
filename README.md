[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/MolSSI/mmic_mda/workflows/CI/badge.svg)](https://github.com/MolSSI/mmic_mda/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/MolSSI/mmic_mda/branch/master/graph/badge.svg)](https://codecov.io/gh/MolSSI/mmic_mda/branch/master)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/MolSSI/mmic_mda.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/MolSSI/mmic_mda/context:python)

MDAnalysis translator for MMSchema
==============================
This is part of the [MolSSI](http://molssi.org) Molecular Mechanics Interoperable Components ([MMIC](https://github.com/MolSSI/mmic)) project. This package provides translators between MMSchema and [MDAnalysis](https://github.com/MDAnalysis/mdanalysis).

![image](mmic_mda/data/imgs/component.png)

# Basic usage
**mmic_mda** provides 3 classes of translators for: molecules, trajectories, and forcefields.

## Molecules
```python
from mmic_mda.models import MdaMol

# Convert MMSchema to MDAnalysis molecule
mda_mol = MdaMol.from_schema(mm_mol) -> MDAnalysis.Universe

# Convert MDAnalysis to MMSchema molecule
mm_mol = MdaMol.to_schema(mda_mol) -> mmelemental.models.molecule.Mol

```
# Under the hood
## Molecules
The `from_schema` and `to_schema` methods in the `MdaMol` model use translation components provided by **mmic_mda** and **MMElemental** to convert between MMSchema and MDAnalysis.

```python
from mmic_mda.components import MdaToMolComponent, MolToMdaComponent
from mmic_mda.models.import MdaMol
from mmelemental.models.molecule import Mol
```

### MMSchema to MDAnalysis molecule
```python
# Creating MMSchema molecule
mm_mol = Mol.from_file(path_to_file)

# Preparing translation input for a molecule
trans_input = TransInput(mol=mda_mol)

# Running compute
mda_mol = MolToMdaComponent.compute(trans_input)
```

### MDAnalysis to MMSchema molecule
```python
# Creating MDAnalysis input
mda_uni = mda.Universe(path_to_file)
mda_mol = mmic_mda.models.MdaMol(mol=mda_uni)

# Preparing translation input for a molecule
trans_input = TransInput(mol=mda_mol)

# Running compute
mm_mol = Translator.compute(trans_input)
```


### Copyright
Copyright (c) 2021, MolSSI


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
