import collections
import logging
import pathlib
import time
import warnings
import datetime

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

import pandas as pd
from rdkit.Chem import PandasTools, rdFingerprintGenerator
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
import requests
from rdkit.Chem import PandasTools
from tqdm.auto import tqdm
from kinfraglib.utils import standardize_mol
from src.evaluation.utils import (read_mols)
import redo
import requests_cache
import nglview
import pypdb
import biotite.database.rcsb as rcsb
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools

from opencadd.structure.superposition.api import align, METHODS
from opencadd.structure.core import Structure

# Disable some unneeded warnings
logger = logging.getLogger("opencadd")
logger.setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

# Cache requests -- this will speed up repeated queries to PDB
requests_cache.install_cache("rcsb_pdb", backend="memory")


def _fetch_pdb_nonpolymer_info(pdb_id):
    """
    Fetch nonpolymer data from rcsb.org.
    Thanks @BJWiley233 and Rachel Green for this GraphQL solution.
    """
    query = (
        """{
          entry(entry_id: "%s") {
            nonpolymer_entities {
              pdbx_entity_nonpoly {
                comp_id
                name
                rcsb_prd_id
              }
            }
          }
        }"""
        % pdb_id
    )

    query_url = f"https://data.rcsb.org/graphql?query={query}"
    response = requests.get(query_url)
    response.raise_for_status()
    info = response.json()
    return info

def _fetch_ligand_expo_info(ligand_id):
  url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"
  r = requests.get(url)
  data = r.json()
  chem = data.get("chem_comp", {})
  desc = data.get("rcsb_chem_comp_descriptor", {})
  return {
      "@structureId": ligand_id,
      "@chemicalID": chem.get("id"),
      "@type": chem.get("type"),
      "@molecularWeight": chem.get("formula_weight"),
      "chemicalName": chem.get("name"),
      "formula": chem.get("formula"),
      "InChI": desc.get("InChI"),
      "InChIKey": desc.get("InChIKey"),
      "smiles": desc.get("smiles"),
  }
  
  
# TODO rewrite for chembl

def get_ligands(pdb_id):
    """
    RCSB has not provided a new endpoint for ligand information yet. As a
    workaround we are obtaining extra information from ligand-expo.rcsb.org,
    using HTML parsing. Check Talktorial T011 for more info on this technique!
    """
    pdb_info = _fetch_pdb_nonpolymer_info(pdb_id)
    ligand_expo_ids = [
        nonpolymer_entities["pdbx_entity_nonpoly"]["comp_id"]
        for nonpolymer_entities in pdb_info["data"]["entry"]["nonpolymer_entities"]
    ]

    ligands = {}
    for ligand_expo_id in ligand_expo_ids:
        ligand_expo_info = _fetch_ligand_expo_info(ligand_expo_id)
        ligands[ligand_expo_id] = ligand_expo_info

    return ligands

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
        chembl["similarity"] = DataStructs.BulkTanimotoSimilarity(
            query_fingerprint, chembl_fingerprints
        )

        # get ligand with maximal similarity
        chembl_most_similar_ix = chembl.similarity.idxmax()

        return [
            chembl.loc[chembl_most_similar_ix]['@chemicalID'],
            round(chembl.loc[chembl_most_similar_ix].similarity, 2),
        ]

    except Exception as e:

        print(f"Most similar ChEMBL ligand search problem for {ligand_inchi}: {e}")
        return [None, None]
    
    
def get_data_for_ids(pdb_ids):
    rows = []
    for pdb_id in pdb_ids:
        ligands = get_ligands(pdb_id)

        ligands = {k: l for k, l in ligands.items() if l['@type'] == 'non-polymer'}
        
        ligand_id, properties = max(
            ((k, l) for k, l in ligands.items() if l['@type'] == 'non-polymer' and l['@molecularWeight'] != None),
            key=lambda kv: kv[1].get("@molecularWeight", 0), default=[None, None]
        )

        if ligand_id and properties:
            properties['@structureId'] = ligand_id
            rows.append(properties)
        else:
            print(f"{pdb_id}: no non-polymer entry")
            
    return pd.DataFrame(rows)

def generate_and_save(ligands, data, path):
    PandasTools.AddMoleculeColumnToFrame(ligands,'smiles')
    ligands = ligands.dropna(ignore_index=True, subset="ROMol")
    ligands["ROMol_standardized"] = ligands["ROMol"].apply(standardize_mol)
    rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
    ligands["fingerprint"] = ligands["ROMol_standardized"].map(
            lambda x: rdkit_gen.GetFingerprint(x)
        )
    x = [most_similar_chembl_ligand(ligand_inchi, ligands)
                for ligand_inchi in data.inchi]
    data = data.copy()
    data["most_similar_chembl_ligand.ligand_id"] = [
        res[0] for res in x
    ]
    data["most_similar_chembl_ligand.similarity"] = [
        res[1] for res in x
    ]
    data.to_csv(path)

# ===== constants =====

uniprot_id = "P25321" #hamster PKA
experimental_method = "X-RAY DIFFRACTION"
min_ligand_molecular_weight = 100.0

query_by_uniprot_id = rcsb.FieldQuery(
    "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
    exact_match=uniprot_id,
)
query_by_experimental_method = rcsb.FieldQuery("exptl.method", exact_match=experimental_method)
query_by_ligand_mw = rcsb.FieldQuery(
    "chem_comp.formula_weight", molecular_definition=True, greater=min_ligand_molecular_weight
)
query_by_non_polymer = rcsb.FieldQuery(
    "rcsb_entry_info.nonpolymer_entity_count", greater=0
)
query_by_title = rcsb.FieldQuery(
    "struct.title", 
    contains_phrase="kinase"
)

# === load data ===

data = read_mols("../results/5n1f/results.sdf")

# post filtering
data = data[data["binding_affinity"] <= 1000].reset_index(drop=True)

# ===== hamster ====== 
print("==== Hamster =====")

hamster_query = rcsb.CompositeQuery(
    [
        query_by_uniprot_id,
        query_by_experimental_method,
        query_by_ligand_mw,
        query_by_non_polymer
    ],
    "and",
)
hamster_pdb_ids = rcsb.search(hamster_query)

print(f"Number of matches: {len(hamster_pdb_ids)}")

hamster_ligands = get_data_for_ids(hamster_pdb_ids)

print(f"Number of structures: {len(hamster_ligands)}")

generate_and_save(hamster_ligands, data, "hamster_pka_compare.csv")

print("==== Kinase =====")

kinase_query = rcsb.CompositeQuery(
    [
        query_by_title,
        query_by_experimental_method,
        query_by_ligand_mw,
        query_by_non_polymer
    ],
    "and",
)
kinase_pdb_ids = rcsb.search(kinase_query)

print(f"Number of matches: {len(kinase_pdb_ids)}")

kinase_ligands = get_data_for_ids(kinase_pdb_ids)

print(f"Number of structures: {len(kinase_ligands)}")

generate_and_save(kinase_ligands, data, "kinase_compare.csv")
