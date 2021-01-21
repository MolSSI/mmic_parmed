[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/MolSSI/mmic_mda/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/mmic_mda/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/mmic_mda/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/mmic_mda/branch/master)

MDAnalysis translator for MMSchema
==============================
This is part of the [MolSSI](http://molssi.org) Molecular Mechanics Interoperable Components ([MMIC](https://github.com/MolSSI/mmic)) project. This package provides translators between MMSchema and [MDAnalysis](https://github.com/MDAnalysis/mdanalysis).

![image](mmic_mda/data/imgs/component.png)

# MDAnalysis translator in action

```python
from mmic_mda.components import Translator
from mmic_mda.models.import TransInput, MdaMol
from mmelemental.models import Mol
```

### MMSchema to MDAnalysis molecule
```python
# Creating MMSchema molecule
mm_mol = Mol.from_file(path_to_file)

# Preparing translation input for a molecule
trans_input = TransInput(mol=mm_mol)

# Running compute
mda_mol = Translator.compute(trans_input)
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
