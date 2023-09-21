import pandas as pd
from rdkit import (Chem, DataStructs, Geometry)
from rdkit.Chem import (Descriptors,
    Lipinski, QED, Draw, AllChem, PandasTools, rdFingerprintGenerator)
import json
from kinfraglib import utils
import matplotlib.pyplot as plt
import seaborn as sns
import utils_eval
import math
from copy import deepcopy
from rdkit.ML.Cluster import Butina
import numpy as np
from chembl_webresource_client.new_client import new_client

def most_similar_chembl_ligand(ligand_inchi, chembl):
    """
    Get the most similar ChEMBL ligand (ChEMBL compound ID and Tanimoto similarity) to the query ligand.

    Parameters
    ----------
    ligand_inchi : str
        Recombined ligand (InChI)
    kinodata : pandas.DataFrame
        kinodata ligands, column fingerprint necessary.

    Returns
    -------
    tuple of (str, str, str, float)
        ChEMBL assay ID, ChEMBL target ID, ChEMBL compound ID and Tanimoto similarity of kinodata ligand most similar to the query ligand.
    """
    try:

        # get ROMol from recombined ligand InChI
        ligand = Chem.MolFromInchi(ligand_inchi)

        # generate query ligand fingerprint
        rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
        query_fingerprint = rdkit_gen.GetFingerprint(ligand)

        # get ChEMBL fingerprints as list
        chembl_fingerprints = chembl.fingerprint.to_list()

        # get pairwise similarities
        chembl['similarity'] = DataStructs.BulkTanimotoSimilarity(query_fingerprint, chembl_fingerprints)

        # get ligand with maximal similarity
        chembl_most_similar_ix = chembl.similarity.idxmax()

        return [
            chembl.loc[chembl_most_similar_ix].chembl_id,
            round(chembl.loc[chembl_most_similar_ix].similarity, 2)
        ]

    except Exception as e:
        
        print(f'Most similar ChEMBL ligand search problem for {ligand_inchi}: {e}')
        return [None, None]

USE_MORGAN = False
CHEMBL_PATH = 'evaluation/chembl_33.sdf'
COMPOUNDS_PATH = 'final_results/3amb/results.sdf'
OUTPATH = 'evaluation/chembl_most_similar_rdkit.csv'

# read chembl data
chembl_data = PandasTools.LoadSDF(
   CHEMBL_PATH,
   embedProps=True
)

print(chembl_data.shape)

print('added inchi')

# generate
if USE_MORGAN:
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
    chembl_data['fingerprint'] = chembl_data['ROMol'].map(lambda x: morgan_gen.GetFingerprint(x))
else:
    rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
    chembl_data['fingerprint'] = chembl_data['ROMol'].map(lambda x: rdkit_gen.GetFingerprint(x))

print('calculated fingerprints')

# read results data
data = utils_eval.read_mols(COMPOUNDS_PATH)
data['binding_affinity'] = data.apply(lambda x: (x['BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_UPPER_BOUNDARY [nM]'] + x['BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_LOWER_BOUNDARY [nM]'])/2, axis=1)
# post filtering
data_post_filtered = data[data['binding_affinity'] <= 1000].copy().reset_index(drop=True) 

data_post_filtered['inchi'] = data_post_filtered.apply(lambda x: Chem.MolToInchi(utils.standardize_mol(x.ROMol)), axis=1)

# calculated most similar kinodata ligand
most_similar_chembl_ligands = [most_similar_chembl_ligand(ligand_inchi, chembl_data) for ligand_inchi in data_post_filtered.inchi]
data_post_filtered['most_similar_chembl_ligand.compound_id'] = [res[0] for res in most_similar_chembl_ligands]
data_post_filtered['most_similar_chembl_ligand.similarity'] = [res[1] for res in most_similar_chembl_ligands]

# save to csv
data_post_filtered.to_csv(OUTPATH)