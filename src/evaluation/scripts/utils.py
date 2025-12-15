import requests

import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, PandasTools
from kinfraglib.utils import standardize_mol


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
            chembl.loc[chembl_most_similar_ix].chembl_id,
            round(chembl.loc[chembl_most_similar_ix].similarity, 2),
        ]

    except Exception as e:

        print(f"Most similar ChEMBL ligand search problem for {ligand_inchi}: {e}")
        return [None, None]


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


def most_similar_ligand(query_inchi, targets):
    """
    Get the most similar target ligand (compound ID and Tanimoto similarity) to the query ligand.

    Parameters
    ----------
    ligand_inchi : str
        Recombined ligand (InChI)
    targets : Dataframe with fingerprint column

    Returns
    -------
    tuple of (str, str, str, float)
        assay ID, ChEMBL target ID, ChEMBL compound ID and Tanimoto similarity of kinodata ligand most similar to the query ligand.
    """
    try:

        # get ROMol from recombined ligand InChI
        ligand = Chem.MolFromInchi(query_inchi)

        # generate query ligand fingerprint
        rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
        query_fingerprint = rdkit_gen.GetFingerprint(ligand)

        # get target fingerprints as list
        target_fingerprints = targets.fingerprint.to_list()

        # get pairwise similarities
        targets["similarity"] = DataStructs.BulkTanimotoSimilarity(
            query_fingerprint, target_fingerprints
        )

        # get ligand with maximal similarity
        chembl_most_similar_ix = targets.similarity.idxmax()
        chembl_least_similar_ix = targets.similarity.idxmin()
        chembl_sim_mean = targets.similarity.mean()

        return (
            [
                targets.loc[chembl_most_similar_ix]["@chemicalID"],
                round(targets.loc[chembl_most_similar_ix].similarity, 2),
            ],
            [
                targets.loc[chembl_least_similar_ix]["@chemicalID"],
                round(targets.loc[chembl_least_similar_ix].similarity, 2),
            ],
            chembl_sim_mean,
        )

    except Exception as e:

        print(f"Most similar ChEMBL ligand search problem for {query_inchi}: {e}")
        return [None, None]


def get_data_for_ids(pdb_ids):
    """
    Redrives ligand information for pdb_ids

    Parameters
    ----------
    pdb_ids : list(str)
        PDB ids of structures from which the ligands should be determined

    Returns
    -------
    Dataframe
        Redrived ligands
    """
    rows = []
    for pdb_id in pdb_ids:
        ligands = get_ligands(pdb_id)

        ligands = {k: l for k, l in ligands.items() if l["@type"] == "non-polymer"}

        ligand_id, properties = max(
            (
                (k, l)
                for k, l in ligands.items()
                if l["@type"] == "non-polymer" and l["@molecularWeight"] != None
            ),
            key=lambda kv: kv[1].get("@molecularWeight", 0),
            default=[None, None],
        )

        if ligand_id and properties:
            properties["@structureId"] = ligand_id
            rows.append(properties)
        else:
            print(f"{pdb_id}: no non-polymer entry")

    return pd.DataFrame(rows)


def generate_and_save(target_ligands, query_ligands, path):
    """
    Determines the most similar compounds from the given target data for each of the query ligands
    and saves to information to a CSV file

    Parameters
    ----------
    query_ligands : Dataframe
        Query ligands
    target_lignads : Dataframe
        Target ligands
    path : str
        CSV file wher results are saved

    """

    PandasTools.AddMoleculeColumnToFrame(target_ligands, "smiles")
    target_ligands = target_ligands.dropna(ignore_index=True, subset="ROMol")
    target_ligands["ROMol_standardized"] = target_ligands["ROMol"].apply(
        standardize_mol
    )
    rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
    target_ligands["fingerprint"] = target_ligands["ROMol_standardized"].map(
        lambda x: rdkit_gen.GetFingerprint(x)
    )
    x = [
        most_similar_ligand(ligand_inchi, target_ligands)
        for ligand_inchi in query_ligands.inchi
    ]
    query_ligands = query_ligands.copy()
    query_ligands["most_similar.ligand_id"] = [res[0][0] for res in x]
    query_ligands["most_similar.similarity"] = [res[0][1] for res in x]
    query_ligands["least_similar.ligand_id"] = [res[1][0] for res in x]
    query_ligands["least_similar.similarity"] = [res[1][1] for res in x]
    query_ligands["mean_similarity"] = [res[2] for res in x]
    query_ligands.to_csv(path)
