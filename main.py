#this will be the main script for our drugdiscoverychallenge code. Please only merge working code here. 

import pandas as pd
df = pd.read_csv('tested_molecules-1.csv')
df.info()
#Draw.MolsToGridImage([Chem.MolFromSmiles(s) for s in pu], molsPerRow=4, subImgSize=(200,200), legends=pu)