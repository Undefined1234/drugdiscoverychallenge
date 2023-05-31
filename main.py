#this will be the main script for our drugdiscoverychallenge code. Please only merge working code here. 
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
IPythonConsole.ipython_useSVG=True


m = Chem.MolFromSmiles("c1ccccc1-C-c1ccccc1")
m

import pandas as pd
df = pd.read_csv('tested_molecules-1.csv',index_col=0)
df.info()
#Draw.MolsToGridImage([Chem.MolFromSmiles(s) for s in pu], molsPerRow=4, subImgSize=(200,200), legends=pu)