#this will be the main script for our drugdiscoverychallenge code. Please only merge working code here. 
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
IPythonConsole.ipython_useSVG=True

m = Chem.MolFromSmiles('Cc1ccccc1')
print(Descriptors.RingCount(m))
