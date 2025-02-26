from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


def most_similar_database_ligand(ligand_inchi, database_ligand, id_columns):
    """
    Get the most similar database ligand (ID and Tanimoto similarity) to the query ligand.

    Parameters
    ----------
    ligand_inchi : str
        Recombined ligand (InChI)
    database_ligand : pandas.DataFrame
        database_ligands ligands, column fingerprint necessary.
    id_columns: list(str)
        Column names of identifier columns (colums that should be copied)

    Returns
    -------
    tuple of (list[str], float)
        database ligand ID, Tanimoto similarity of database ligand most similar to the query ligand.
    """
    try:

        # get ROMol from recombined ligand InChI
        ligand = Chem.MolFromInchi(ligand_inchi)

        # generate query ligand fingerprint
        rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
        query_fingerprint = rdkit_gen.GetFingerprint(ligand)

        # get database fingerprints as list
        database_fingerprints = database_ligand.fingerprint.to_list()

        # get pairwise similarities
        database_ligand["similarity"] = DataStructs.BulkTanimotoSimilarity(
            query_fingerprint, database_fingerprints
        )

        # get ligand with maximal similarity
        database_most_similar_ix = database_ligand.similarity.idxmax()

        return [[database_ligand.loc[database_most_similar_ix][idx] for idx in id_columns],
            round(database_ligand.loc[database_most_similar_ix].similarity, 2),
        ]

    except Exception as e:

        print(f"Most similar database ligand search problem for {ligand_inchi}: {e}")
        return [None, None]
