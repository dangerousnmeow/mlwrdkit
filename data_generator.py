from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Draw import SimilarityMaps
import pprint
from rdkit.Chem import PandasTools
import collections
import pandas as pd
import numpy as np


def data_generator(positive_hits,data):
    #load dataset
    data = PandasTools.LoadSDF(data,smilesName='SMILES',molColName='mol', includeFingerprints=True,embedProps=True)
    positive_hits = PandasTools.LoadSDF(positive_hits,smilesName='SMILES',molColName='mol', includeFingerprints=True,embedProps=True)
        
    #generate hit names
    hitname = list(positive_hits.NAME)
        
    #test if the compound in dataset is unique
    if len(np.unique(hitname)) != len(hitname):
        print([item for item, count in collections.Counter(hitname).items() if count > 1])
    else:
    #get the label 0-negative,1-positive
        y = pd.to_numeric(data["NAME"].isin(hitname).astype('uint8'))
        data['mol'] = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in data['mol']]
        data_name = [f'Bit_{i}' for i in range(2048)]
        data_bits = [list(x) for x in data['mol']]
        X = pd.DataFrame(data_bits, columns=data_name)
    return X,y
