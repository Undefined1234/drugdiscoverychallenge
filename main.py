#this will be the main script for our drugdiscoverychallenge code. Please only merge working code here. 
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
IPythonConsole.ipython_useSVG=True

m = Chem.MolFromSmiles("c1ccccc1-C-c1ccccc1")
m