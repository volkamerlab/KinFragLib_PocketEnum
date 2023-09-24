import pandas as pd
from rdkit import Chem
from rdkit.Chem import (Descriptors,
    Lipinski, QED, Draw, AllChem)
import json
from kinfraglib import utils
import matplotlib.pyplot as plt
import seaborn as sns

def read_mols(path_to_mols):
    """
    Read ligands from result file.

    Parameters
    ----------
    path_to_lib : str
        Path to results .sdf file.

    Returns
    -------
    pandas.DataFrame
        ligands details details, i.e. SMILES, and RDKit molecules.
    """
    to_float_probs = ['BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_LOWER_BOUNDARY [nM]', 
                        'BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_UPPER_BOUNDARY [nM]', 
                        'BIOSOLVEIT.LOGP', 'BIOSOLVEIT.MOLECULAR_WEIGHT', 
                        'BIOSOLVEIT.TPSA', 'BIOSOLVEIT.DOCKING_SCORE']
    to_dict_probs = ['fragment_ids', 'smiles_fragments_dummy', 'smiles_fragments']

    data = [[
                mol # 3D conformation of ligand
            ] + [    
                (float(mol.GetProp(prop_name)) if prop_name in to_float_probs else mol.GetProp(prop_name).replace("\'", "\""))
                    for prop_name in mol.GetPropNames() 
            ] for mol in Chem.SDMolSupplier(str(path_to_mols), removeHs=False)]

    data_df = pd.DataFrame(
        data,
        columns = ['ROMol'] + list(data[0][0].GetPropNames())
    ).apply(lambda x: x.apply(float) if x.name in to_float_probs else x).apply(lambda x: x.apply(json.loads) if x.name in to_dict_probs else x)
    
    data_df['num_fragments'] = [len(x) for x in data_df['fragment_ids']]
    return data_df