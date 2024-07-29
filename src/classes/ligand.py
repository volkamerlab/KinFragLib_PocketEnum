from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.ML.Cluster import Butina

import logging

from utils.brics_rules import is_brics_bond
from classes.recombination import Recombination
from utils._utils import calc_distance_matrix

SUBPOCKETS = ["AP", "FP", "SE", "GA", "B1", "B2"]


class Pose:
    def __init__(self, ROMol, docking_score):
        self.docking_score = docking_score
        self.ROMol = ROMol
        self.binding_affinity_upper = None
        self.binding_affinity_lower = None

    def _is_valid_operand(self, other):
        return hasattr(other, "docking_score")


class Ligand:
    def __init__(
        self,
        ROMol,
        fragment_ids: dict,
        recombination: Recombination,
        smiles_dummy: dict,
        smiles: dict,
    ):
        self.ROMol = ROMol
        self.poses = []
        self.fragment_ids = fragment_ids
        self.dummy_atoms = {}
        self.min_docking_score = None
        self.min_binding_affinity = None
        self.recombinations = []
        self.recombination: Recombination = recombination
        self.smiles_dummy = smiles_dummy
        self.smiles = smiles
        self.num_hyde_violations = 0
        self.mean_hyde_displacement_undropped = 0  # mean of displacement of all poses that were not dropped due to displacement violation
        self.mean_hyde_displacement = 0  # mean of displacement of all poses (including violations)
        self.poses_pre_filtered = None

    def to_sdf(self, sdf_path):
        """
        Writes the current object to a sdf file.

        Parameters
        ----------
        sdf_path: str or pathlib.path
            Path to output sdf file.
        pH: float
            Protonation at given pH.
        """

        # protonate
        molecule = Chem.AddHs(self.ROMol)

        # 3D generation & optimization of the ligand itself
        if AllChem.EmbedMolecule(molecule, randomSeed=0xF00D) < 0:
            # molecule is too big
            if Chem.AllChem.EmbedMolecule(molecule, randomSeed=0xF00D, useRandomCoords=True) != 0:
                # embedding wasn't succesful
                logging.error("Could not embed molecule")
                return False
        try:
            Chem.AllChem.UFFOptimizeMolecule(molecule)
            with Chem.SDWriter(str(sdf_path)) as w:
                w.write(molecule)
            return True
        except ValueError:
            logging.error("Could not optimize molecule")
            return False

    def compute_unique_id(self) -> str:
        """
        Computes a unique id for the ligand

        Returns
        --------
        str: a unique id
        """
        # scheme: [AP_index][FP_index][SE_index][GA_index][B1_index][B2_index] each 4 digits
        return "L" + "".join([f"{(self.fragment_ids.get(sp) or 0):04d}" for sp in SUBPOCKETS])

    def add_pose(self, pose: Pose):
        """
        Adds a pose to the given ligand.

        Parameters
        ----------
        pose: Pose
            Pose object that should be added
        """
        if self.min_docking_score == None or pose.docking_score < self.min_docking_score:
            self.min_docking_score = pose.docking_score
        if pose.binding_affinity_lower and (
            self.min_binding_affinity == None
            or pose.binding_affinity_lower < self.min_binding_affinity
        ):
            self.min_binding_affinity = pose.binding_affinity_lower
        self.poses.append(pose)
        return

    def get_best_pose(self) -> Pose:
        """
        Get Pose with the best docking score

        Resturns
        --------
        Pose
            pose with lowest docking score
        """
        return (
            min(self.poses, key=lambda p: p.binding_affinity_lower)
            if self.min_binding_affinity
            else min(self.poses, key=lambda p: p.docking_score)
        )
    
    def get_min_score(self) -> float:
        """
        Get the minimum score (docking or hyde score)

        Resturns
        --------
        float
            minimum score
        """
        return self.min_binding_affinity or self.min_docking_score

    def choose_template_poses(self, num_templates=None):
        """
        Chooses *num_templates*-best poses according to docking result and diversity and removes the non-top-scored or invalid poses

        Parameters
        ----------
        num_templates: int
            max. templates to choose, if None, num_templates == len(path)
        """
        self.poses_pre_filtered = self.poses.copy()

        if not len(self.poses):
            # nothing to choose from
            return
        if not num_templates:
            # default of num_poses = amount of poses
            return

        # always choose the pose with the best docking score
        choosen_poses = [min(self.poses, key=lambda p: p.docking_score)]
        self.poses.remove(choosen_poses[0])

        # define mapping of atoms => avoid a maximal common substructure search (not needed because the structures should be the same)
        atom_mapping = [[j, j] for j in range(choosen_poses[0].ROMol.GetNumAtoms())]

        # append num_templates - 1 more templates according to highest mean of RMSD (to already choosen poses) and the lowest docking-score
        for _ in range(num_templates - 1):
            if not len(self.poses):
                break
            # choose minimal pose according to -2 * mean(RMSD to choosen poses) + dockin_score
            pose = min(
                self.poses,
                key=lambda pose: -2
                * sum(
                    rdMolAlign.CalcRMS(pose.ROMol, x.ROMol, map=[atom_mapping])
                    for x in choosen_poses
                )
                / len(choosen_poses)
                + pose.docking_score,
            )
            choosen_poses.append(pose)
            self.poses.remove(pose)
        # overwrite poses with choosen poses
        self.poses = choosen_poses

    def choose_template_poses_cluster_based(self, num_templates=None, distance_threshold=1.5):
        """
        Chooses up to *num_templates*-best poses according to docking result and diversity and
        removes the non-top-scored or invalid poses.
        Poses are clustered according to RMSD first. Then at most one pose per cluster is choosen.

        Parameters
        ----------
        num_templates: int
            max. templates to choose, if None, num_templates == number of poses that have rmsd >= 1.5 A
        """

        self.poses_pre_filtered = self.poses.copy()

        if not len(self.poses):
            # nothing to choose from
            return

        # calculate the distance matrix according to RMSD
        dists_RMS = calc_distance_matrix([pose.ROMol for pose in self.poses])

        # cluster poses according to the distance matrix
        clusters = Butina.ClusterData(
            dists_RMS, len(self.poses), distance_threshold, isDistData=True, reordering=True
        )

        # num_templates = min(amount clusters, num_templates) => ensures that at most one pose per cluster is choosen
        num_templates = (
            num_templates if num_templates and num_templates <= len(clusters) else len(clusters)
        )

        # only use the best pose (according to docking score) per cluster
        if self.min_binding_affinity:
            clustered_pose_gen = (
                min([self.poses[idx] for idx in cluster], key=lambda p: p.binding_affinity_lower)
                for cluster in clusters
            )
            self.poses = sorted(clustered_pose_gen, key=lambda p: p.binding_affinity_lower)[
                :num_templates
            ]
        else:
            clustered_pose_gen = (
                min([self.poses[idx] for idx in cluster], key=lambda p: p.docking_score)
                for cluster in clusters
            )
            self.poses = sorted(clustered_pose_gen, key=lambda p: p.docking_score)[:num_templates]

    def calculate_missing_dummy_atoms(self, fragment_library):
        """
        Determines dummy atoms that aren't calculated

        Parameters
        ----------
        fragment_library: dict
            Library containing all fragments where the index should match to the fragment ids
        """
        for subpocket, id in self.fragment_ids.items():
            if subpocket not in self.dummy_atoms.keys():
                # calculate the dummy atoms if there are no dummy atoms for the subpocket
                fragment = Chem.RemoveHs(fragment_library[subpocket]["ROMol_original"][id])
                self.dummy_atoms[subpocket] = [
                    (
                        f"{subpocket}_{i}",
                        a.GetNeighbors()[0].GetProp("environment"),
                        a.GetProp("subpocket"),
                    )
                    for i, a in enumerate(fragment.GetAtoms())
                    if a.GetSymbol() == "*"
                ]

    def recombine(self, fragment_id, subpocket, fragment_library):
        """
        Recombines the ligand with the fragment (subpocket_fragment_id)

        Returns
        ----------
        True if recombination could be found else False

        Parameters
        ----------
        fragment_library: dict
            Library containing all fragments where the index should match to the fragment ids
        """
        # make sure that all dummy atoms are added to self.dummy_atoms
        self.calculate_missing_dummy_atoms(fragment_library)
        # this is need to detemine the correct atom-id
        fragment = Chem.RemoveHs(fragment_library[subpocket]["ROMol_original"][fragment_id])
        # dummy atoms of the fragment that should to be recombined
        dummy_atoms = [
            (
                f"{subpocket}_{i}",
                a.GetNeighbors()[0].GetProp("environment"),
                a.GetProp("subpocket"),
            )
            for i, a in enumerate(fragment.GetAtoms())
            if a.GetSymbol() == "*"
        ]
        counter_recombinations = 0
        counter_unambiguous_bonds = 0
        for sp in self.dummy_atoms.keys():
            # for every subpocket (used by ligand): add a connection if valid
            matching_dummies = [
                (id, env) for id, env, con in dummy_atoms if con == sp
            ]  # dummy atoms of fragment that have a connection to the current subpocket
            matching_dummies_2 = [
                (id, env) for id, env, con in self.dummy_atoms[sp] if con == subpocket
            ]  # dummy atoms of ligand (in subpocket sb) that have a connection to the subpocket of the fragment
            if len(matching_dummies) != 1 or len(matching_dummies_2) != 1:
                # if there are more than 1 connection to one subpocket (for now we only allow one connection between subpockets) OR no connection
                # TODO maybe add all possibities
                counter_unambiguous_bonds += 1
                continue
            id, env = matching_dummies[0]
            id_2, env_2 = matching_dummies_2[0]
            if not is_brics_bond(env, env_2):
                continue
            new_rec = self.recombination.copy()
            new_rec.add_fragment(
                subpocket + "_" + str(fragment_id),
                [[id, id_2]],
                fragment_library[subpocket]["smiles_dummy"][fragment_id],
                fragment_library[subpocket]["smiles"][fragment_id],
            )
            new_rec.construct(fragment_library)
            if new_rec.ligand != None:
                counter_recombinations += 1
                self.recombinations.append(new_rec)
        return counter_recombinations, counter_unambiguous_bonds


def from_recombination(recombination) -> Ligand:
    """
    Converts a given Recombination-object to a Ligand

    Returns
    ----------
    Ligand

    Parameters
    ----------
    recombination: Recombination
        Recombination that should be converted
    """
    fragment_ids = {x[:2]: int(x[3:]) for x in recombination.fragments}
    smiles = recombination.smiles.copy()
    smiles_dummy = recombination.smiles_dummy.copy()
    return Ligand(recombination.ligand, fragment_ids, recombination, smiles_dummy, smiles)
