from pathlib import Path
import subprocess
import os
from typing import Any

from rdkit import Chem
import sys
from rdkit.ML.Cluster import Butina

from rdkit.Chem import AllChem, Draw

from rdkit.Chem import rdMolAlign

from rdkit.Chem.PropertyMol import PropertyMol
from functools import reduce
import pandas as pd
from kinfraglib import utils, filters
from brics_rules import is_brics_bond
import logging

def core_docking(path_fragment, path_config, path_output, path_flexx, print_output=False):
    """
    runs FlexX docking

    Returns
    ----------
    List of all poses (docking result) as Mol (Properties: docking-score, pose)

    Parameters
    ----------
    path_fragment: pathlib.path
        Path to core-fragment sdf-file
    path_config: pathlib.path
        Path to FlexX-config file.
    path_flexx: pathlib.path
        Path to FlexX
    path_output: pathlib.path
        Path to output file
    """
    # core docking
    output_text = subprocess.run(
        [
            str(Path('.') / path_flexx),
            "-i",
            str(path_fragment),
            "--docking-definition",
            str(path_config),
            "-o",
            str(path_output),
        ],
        capture_output=True
    )
    if print_output:
        print(output_text.stdout)

    # read core fragments poses from sdf
    docked_core_fragments = []
    if os.stat(str(path_output)).st_size:   # only acces reult file if at least one pose was generated
        for i, molecule in enumerate(Chem.SDMolSupplier(str(path_output))):
            molecule.SetProp('pose', str(i))
            docked_core_fragments.append(molecule)
    
    return docked_core_fragments

def remove_files(*path_files: Path):
    """
    Deletes the given files if they exist (this function should be used after docking)

    Parameters
    ----------
    path_files: **pathlib.path
        Paths to files that should be deleted
    """
    for path_file in path_files:
        if os.path.exists(path_file):
            os.remove(str(path_file))

def hyde_scoring(path_docking_results, path_config, path_output, print_output=False):
    """
    runs hydescoring

    Returns
    ----------
    List of all poses (docking result) as Mol (Properties: BIOSOLVEIT.DOCKING_SCORE, pose)

    Parameters
    ----------
    path_fragment: pathlib.path
        Path to fragment sdf-file
    path_config: pathlib.path
        Path to hyde-config file.
    path_output: pathlib.path
        Path to output file
    """
    output_text = subprocess.run(
        [
            str('./Hydescorer.app/Contents/MacOS/hydescorer'),
            "-i",
            str(path_docking_results),
            "--binding-site-definition",
            str(path_config),
            "-o",
            str(path_output),
        ],
        capture_output=True
    )
    if print_output:
        print(output_text.stdout)

    # read results from sdf
    opt_fragments = []
    if os.stat(str(path_output)).st_size and os.stat(str(path_docking_results)).st_size:  # only acces reult file if at least one pose was generated
        for i, molecules in enumerate(zip(Chem.SDMolSupplier(str(path_docking_results)), Chem.SDMolSupplier(str(path_output)))):
            molecule_docking, molecule_opt = molecules
            molecule_opt.SetProp('pose', str(i))
            molecule_opt.SetProp('BIOSOLVEIT.DOCKING_SCORE', molecule_docking.GetProp('BIOSOLVEIT.DOCKING_SCORE'))
            opt_fragments.append(molecule_opt)
    
    return opt_fragments

def calc_distance_matrix(molecules):
    """
    Calculates the distance matrix (based on RMSD) that can be used for clustering

    Parameters
    ----------
    molecules: List(Mol)
        List of molecules for which the matrix is calculated
    """
    atom_mapping = [[j, j] for j in range(molecules[0].GetNumAtoms())]

    # for each combination calculate the RMSD (without considering symmetry)
    return [rdMolAlign.CalcRMS(molecules[i], molecules[j], map = [atom_mapping]) for i in range(len(molecules)) for j in range(i)]



def template_docking(path_fragment, path_template, path_config, path_output, path_flexx, print_output=False):
    """
    runs FlexX template-docking

    Returns
    ----------
    List of all poses (docking result) as Mol (Properties: BIOSOLVEIT.DOCKING_SCORE, pose)

    Parameters
    ----------
    path_fragment: pathlib.path
        Path to fragment sdf-file
    path_template: pathlib.path
        Path to template-fragment sdf-file
    path_config: pathlib.path
        Path to FlexX-config file.
    path_flexx: pathlib.path
        Path to FlexX
    path_output: pathlib.path
        Path to output file
    """
    # template docking
    output_text = subprocess.run(
            [
                str(path_flexx),
                "-i",
                str(path_fragment),
                "--docking-definition",
                str(path_config),
                "-o",
                str(path_output),
                "-t",
                str(path_template),
                "--docking-type",
                "1",
                #"--stereo-mode"
                #,"3" 
            ],
            capture_output=True,  # needed to capture output text
        )
    if print_output:
        print(output_text)

    # read fragments poses from sdf
    docked_fragments = []
    if not os.path.exists(str(path_output)): # If docking failed TODO docking shouldn't fail
        return docked_fragments
    if os.stat(str(path_output)).st_size:
        for i, molecule in enumerate(Chem.SDMolSupplier(str(path_output))):
            molecule.SetProp('pose', str(i))
            docked_fragments.append(molecule)
    
    return docked_fragments

#   ===================  CLASS DEFINITIONS ================

class Pose:
    def __init__(self, ROMol, docking_score):
        self.docking_score = docking_score
        self.ROMol = ROMol
    def _is_valid_operand(self, other):
        return hasattr(other, 'docking_score')

class Ligand:
    def __init__(self, ROMol, fragment_ids, recombination):
        self.ROMol = ROMol
        self.poses = []
        self.fragment_ids = fragment_ids
        self.dummy_atoms = {}
        self.min_docking_score = None
        self.recombinations = []
        self.recombination : Recombination = recombination
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
        if AllChem.EmbedMolecule(molecule, randomSeed=0xf00d) < 0:
            # molecule is too big
            if Chem.AllChem.EmbedMolecule(molecule , randomSeed=0xf00d, useRandomCoords=True) != 0:
                # embedding wasn't succesful
                logging.error('Could not embed molecule')
                return False
        try :
            Chem.AllChem.UFFOptimizeMolecule(molecule)
            with Chem.SDWriter(str(sdf_path)) as w:
                w.write(molecule)
            return True
        except ValueError :
            logging.error('Could not optimize molecule')
            return False
        
    def add_pose(self, pose : Pose):
        """
        Adds a pose to the given ligand.

        Parameters
        ----------
        pose: Pose
            Pose object that should be added
        """
        if self.min_docking_score == None or pose.docking_score < self.min_docking_score:
            self.min_docking_score = pose.docking_score
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
        return min(self.poses, key=lambda p: p.docking_score)
        
    def choose_template_poses(self, num_templates=None):
        """
        Chooses *num_templates*-best poses according to docking result and diversity and removes the non-top-scored or invalid poses

        Parameters
        ----------
        num_templates: int
            max. templates to choose, if None, num_templates == len(path)
        """
        if not len(self.poses):
            # nothing to choose from
            return
        if not num_templates:
            # default of num_poses = amount of poses
            num_templates = len(self.poses)
        
        # always choose the pose with the best docking score
        choosen_poses = [min(self.poses, key=lambda p: p.docking_score)]
        self.poses.remove(choosen_poses[0])

        # -----> TODO REMOVE IF FIXED
        for pose in self.poses:
                if choosen_poses[0].ROMol.GetNumAtoms() != pose.ROMol.GetNumAtoms():
                    # this should't happen TODO should be fixed by cleaning files
                    print("Couldn't match pose")
                    pose.ROMol = choosen_poses[0].ROMol
        # -------|

        # define mapping of atoms => avoid a maximal common substructure search (not needed because the structures should be the same)
        atom_mapping = [[j, j] for j in range(choosen_poses[0].ROMol.GetNumAtoms())]

        # append num_templates - 1 more templates according to highest mean of RMSD (to already choosen poses) and the lowest docking-score
        for _ in range(num_templates - 1):
            if not len(self.poses):
                break
            # choose minimal pose according to -2 * mean(RMSD to choosen poses) + dockin_score
            pose = min(self.poses, key=lambda pose : -2 * sum(rdMolAlign.CalcRMS(pose.ROMol, x.ROMol, map=[atom_mapping]) for x in choosen_poses) / len(choosen_poses) + pose.docking_score)
            choosen_poses.append(pose)
            self.poses.remove(pose)
        # overwrite poses with choosen poses
        self.poses = choosen_poses

    def choose_template_poses_cluster_based(self, num_templates=None, distance_threshold = 1.5):
        """
        Chooses up to *num_templates*-best poses according to docking result and diversity and removes the non-top-scored or invalid poses.
        Poses are clustered according to RMSD first. Then at most one pose per cluster is choosen. 

        Parameters
        ----------
        num_templates: int
            max. templates to choose, if None, num_templates == number of poses that have rmsd >= 1.5 A
        """
        if not len(self.poses):
            # nothing to choose from
            return

        # calculate the distance matrix according to RMSD
        dists_RMS = calc_distance_matrix([pose.ROMol for pose in self.poses])

        # cluster poses according to the distance matrix 
        clusters = Butina.ClusterData(dists_RMS, len(self.poses), distance_threshold, isDistData=True, reordering=True)

        if (not num_templates) or num_templates >= len(clusters):
            # default of num_poses = amount of poses
            num_templates = len(clusters)

        # only use the best pose (according to docking score) per cluster
        clustered_pose_gen = (min([self.poses[idx] for idx in cluster], key = lambda p : p.docking_score) for cluster in clusters)

        # num_templates = min(amount clusters, num_templates) => ensures that at most one pose per cluster is choosen
        num_templates = num_templates if num_templates and num_templates <= len(clusters) else len(clusters)

        self.poses = sorted(clustered_pose_gen, key = lambda p : p.docking_score)[:num_templates]

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
                fragment = Chem.RemoveHs(fragment_library[subpocket]['ROMol_original'][id])
                self.dummy_atoms[subpocket] = [(f"{subpocket}_{i}", a.GetNeighbors()[0].GetProp('environment'), a.GetProp('subpocket')) for i, a in enumerate(fragment.GetAtoms()) if a.GetSymbol() == '*']
    
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
        fragment = Chem.RemoveHs(fragment_library[subpocket]['ROMol_original'][fragment_id])
        # dummy atoms of the fragment that should to be recombined
        dummy_atoms = [(f"{subpocket}_{i}", a.GetNeighbors()[0].GetProp('environment'), a.GetProp('subpocket')) for i, a in enumerate(fragment.GetAtoms()) if a.GetSymbol() == '*']
        new_rec = None
        for sp in self.dummy_atoms.keys():
            # for every subpocket (used by ligand): add a connection if valid
            matching_dummies = [(id, env) for id, env, con in dummy_atoms if con == sp] # dummy atoms of fragment that have a connection to the current subpocket
            matching_dummies_2 = [(id, env) for id, env, con in self.dummy_atoms[sp] if con == subpocket] # dummy atoms of ligand (in subpocket sb) that have a connection to the subpocket of the fragment
            if len(matching_dummies) != 1 or len(matching_dummies_2) != 1:
                # if there are more than 1 connection to one subpocket (for now we only allow one connection between subpockets) OR no connection
                # TODO maybe add all possibities 
                continue
            id, env = matching_dummies[0]
            id_2, env_2 = matching_dummies_2[0]
            if not is_brics_bond(env, env_2):
                continue
            if new_rec == None:
                new_rec = self.recombination.copy()
            new_rec.add_fragment(subpocket + "_" + str(fragment_id), [[id, id_2]])
        if new_rec != None:
            new_rec.construct(fragment_library)
            if new_rec.ligand != None:
                self.recombinations.append(new_rec)
                return True
        return False
    
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
    fragments = {x[:2]: int(x[3:]) for x in recombination.fragments}
    return Ligand(recombination.ligand, fragments, recombination)

class Recombination:
    def __init__(self,fragment_ids, bonds):
        self.ligand = None
        self.bonds = bonds
        self.fragments = fragment_ids
    def add_fragment(self, fragment_id, bonds):
        """
        Adds a fragment given by it's id and bonds

        Parameters
        ----------
        fragment_id: str
            Format: <subpocket>_<id>
        bonds: List(str)
            Format of a bond [<subpocket>_<atom_id>, <subpocket>_<atom_id>]
        """
        if fragment_id not in self.fragments:
            self.fragments.append(fragment_id)
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
        return Recombination(self.fragments.copy(), self.bonds.copy())
    
class Filter:
    def __init__(self, name, params: dict):
        self.name = name.lower()
        self.params = params
    def get_param(self, param_name):
        """
        Get a filter-parameter by it's name

        Returns
        ----------
        Paramter "param_name" if exists else None

        Parameters
        ----------
        param_name: str
            Name of the parameter
        """
        if param_name in self.params.keys():
            return self.params[param_name]
        return None
    def apply_filter(self, fragment_library):
        """
        Applies the filter to the given fragment_library

        Returns
        ----------
        Filtered library

        Parameters
        ----------
        fragment_library: Dict
            Library containing all fragments where the index should match to the fragment ids
        """
        if self.name == 'pains':
            fragment_library, _ = filters.unwanted_substructures.get_pains(fragment_library)
        elif self.name == 'brenk':
            fragment_library, _ = filters.unwanted_substructures.get_brenk(fragment_library, Path(self.get_param('path_data')))
        elif self.name == 'ro3':
            fragment_library = filters.drug_likeness.get_ro3_frags(fragment_library)
        elif self.name == 'qed':
            fragment_library = filters.drug_likeness.get_qed(fragment_library, cutoff_val=self.get_param('cutoff_val'))
        elif self.name == 'bb':
            fragment_library = filters.synthesizability.check_building_blocks(
                fragment_library,
                str(str(self.get_param('path_data')) + "/Enamine_Building_Blocks.sdf"),
            )
        elif self.name == 'syba':
            fragment_library = filters.synthesizability.calc_syba(fragment_library, cutoff=self.get_param('cutoff_val'))
            
        # remove all filtered fragments from the fragment library
        for sp in fragment_library.keys():
            t = fragment_library[sp]['bool_' + self.name]
            fragment_library[sp] = fragment_library[sp][fragment_library[sp]['bool_' + self.name] == 1]
        
        return fragment_library
    
def append_ligands_to_file(ligands : 'list[Ligand]', path_output, condition = lambda l: True):
    """
    Appends the best poses of the given ligand to the output file if they fullfill condition and sets the fragment_ids and bonds as properties

    Parameters
    ----------
    ligands: list[Ligand]
        List containing all ligands
    path_output: pathlib.path
        Path to output file
    condition: Ligand -> bool
        Function that should return true iff the given ligand should be added to the given file, defaults to True
    """

    with open(str(path_output), 'at')  as outf:  # needed to be able to APPEND molecules to the result file
        with Chem.SDWriter(outf) as w:
            ligand: Ligand
            for ligand in ligands:
                if not condition(ligand):
                    continue
                mol = ligand.get_best_pose().ROMol
                mol.SetProp('Fragments', str(ligand.recombination.fragments))
                mol.SetProp('Bonds', str(ligand.recombination.bonds))
                w.write(mol)
