[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/MolSSI/mmic_parmed/workflows/CI/badge.svg)](https://github.com/MolSSI/mmic_parmed/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/MolSSI/mmic_parmed/branch/master/graph/badge.svg)](https://codecov.io/gh/MolSSI/mmic_parmed/branch/master)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/MolSSI/mmic_parmed.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/MolSSI/mmic_parmed/context:python)

MDAnalysis translator for MMSchema
==============================
This is part of the [MolSSI](http://molssi.org) Molecular Mechanics Interoperable Components ([MMIC](https://github.com/MolSSI/mmic)) project. This package provides translators between MMSchema and [MDAnalysis](https://github.com/MDAnalysis/mdanalysis).

![image](mmic_parmed/data/imgs/component.png)

# Basic API
**mmic_parmed** provides 3 classes of translators for: molecules, trajectories, and forcefields.

## Molecules
```python
from mmic_parmed.models import MdaMol

# Convert MMSchema to MDAnalysis molecule
mda_mol = MdaMol.from_schema(mm_mol) -> MDAnalysis.Universe

# Convert MDAnalysis to MMSchema molecule
mm_mol = MdaMol.to_schema(mda_mol) -> mmelemental.models.molecule.Mol

```
# Under the hood
## Molecules
The `from_schema` and `to_schema` methods in the `MdaMol` model use translation components provided by **mmic_parmed** and **MMElemental** to convert between MMSchema and MDAnalysis.

```python
from mmic_parmed.components import MdaToMolComponent, MolToMdaComponent
from mmic_parmed.models.import MdaMol
from mmelemental.models.molecule import Mol
```

### MMSchema to MDAnalysis molecule
```python
# Creating MMSchema molecule
mm_mol = Mol.from_file(path_to_file)

# Running translator compute
mda_mol = MolToMdaComponent.compute(mm_mol)
```

### MDAnalysis to MMSchema molecule
```python
# Creating MDAnalysis input
mda_uni = mda.Universe(path_to_file)
mda_mol = mmic_parmed.models.MdaMol(mol=mda_uni)

# Running translator compute
mm_mol = Translator.compute(mda_mol)
```


### Copyright
Copyright (c) 2021, MolSSI


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
