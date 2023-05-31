import os
import pandas as pd
from rdkit.Chem import PandasTools, Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

from IPython.display import display



TESTED_MOLECULE_PATH = "data\\tested_molecules-1.csv"

if __name__ =="__main__":
    tested_molecules = pd.read_csv(os.path.join(os.getcwd(), TESTED_MOLECULE_PATH))
    PandasTools.AddMoleculeColumnToFrame(tested_molecules, smilesCol='SMILES')



    # r = PandasTools.FrameToGridImage(tested_molecules, molsPerRow=32, maxMols=1000)

    # display(r)

    desc_list = [n[0] for n in Descriptors._descList]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_list)

    rdkit_desc = [calc.CalcDescriptors(m) for m in tested_molecules["SMOL"]]

    # calc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_list)