#!/usr/bin/python3

from rdkit import Chem
from rdkit.Chem import Draw
import sys

inf=sys.argv[1];
mol=Chem.MolFromPDBFile(inf,sanitize=True,proximityBonding=False)
smi=Chem.MolToSmiles(mol)
print(smi)
