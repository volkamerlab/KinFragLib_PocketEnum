from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


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
