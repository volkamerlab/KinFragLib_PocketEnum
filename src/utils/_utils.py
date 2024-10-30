from rdkit.Chem import rdMolAlign

import random
import os


def calc_distance_matrix(molecules: list) -> list:
    """
    Calculates the distance matrix (based on RMSD) that can be used for clustering

    Parameters
    ----------
    molecules: List(Mol)
        List of molecules for which the matrix is calculated
    """
    atom_mapping = [[j, j] for j in range(molecules[0].GetNumAtoms())]

    # for each combination calculate the RMSD (without considering symmetry)
    return [
        rdMolAlign.CalcRMS(molecules[i], molecules[j], map=[atom_mapping])
        for i in range(len(molecules))
        for j in range(i)
    ]

def seed_everything(seed: int = 41):
    """ Code adapted from TeachOpenCADD
    Set seeds
    """
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)