from mmelemental.models import ForceField
import parmed

ff = ForceField.from_file("alanine.top")
pff = ff.to_data("parmed").data

gro_pff = parmed.gromacs.GromacsTopologyFile.from_structure(pff)

ff.to_file("tmp.top")
