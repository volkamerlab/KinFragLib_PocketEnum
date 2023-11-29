from kinfraglib import utils

class Recombination:
    def __init__(self, fragment_ids, bonds, smiles, smiles_dummy):
        self.ligand = None
        self.bonds = bonds
        self.fragments = fragment_ids
        self.smiles = smiles
        self.smiles_dummy = smiles_dummy
    def add_fragment(self, fragment_id, bonds, smiles_dummy, smiles):
        """
        Adds a fragment given by it's id and bonds

        Parameters
        ----------
        fragment_id: str
            Format: <subpocket>_<id>
        bonds: List(str)
            Format of a bond [<subpocket>_<atom_id>, <subpocket>_<atom_id>]
        smiles_dummy: str
            SMILES with dummy atom of fragment that should be added 
        smiles: str
            SMILES (without dummy atom) atom of fragment that should be added
        """
        if fragment_id not in self.fragments:
            self.fragments.append(fragment_id)
        self.smiles[fragment_id[:2]] = smiles
        self.smiles_dummy[fragment_id[:2]] = smiles_dummy
        self.bonds += bonds
    def construct(self, fragment_library):
        """
        Constructs a Mol-object from the recombination and stores it

        Parameters
        ----------
        fragment_library: Dict
            Library containing all fragments where the index should match to the fragment ids
        """
        self.ligand = utils.construct_ligand(self.fragments, self.bonds, fragment_library)
    def copy(self):
        """
        Copies itself

        Returns
        ----------
        New instance of the given recombination

        Parameters
        ----------
        fragment_library: Dict
            Library containing all fragments where the index should match to the fragment ids
        """
        return Recombination(self.fragments.copy(), self.bonds.copy(), self.smiles_dummy.copy(), self.smiles.copy())