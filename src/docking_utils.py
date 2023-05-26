from pathlib import Path
import subprocess
import os
from typing import Any

from rdkit import Chem
import sys

from rdkit.Chem import AllChem

from rdkit.Chem import rdMolAlign

from rdkit.Chem.PropertyMol import PropertyMol
from functools import reduce
import pandas as pd
from kinfraglib import utils
from brics_rules import is_brics_bond

PATH_FLEXX = Path('./FlexX.app/Contents/MacOS/FlexX')

def core_docking(path_fragment, path_config, path_output, print_output=False):
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
    path_output: pathlib.path
        Path to output file
    """
    # core docking
    output_text = subprocess.run(
        [
            str(PATH_FLEXX),
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
    if os.stat(str(path_output)).st_size:
        for i, molecule in enumerate(Chem.SDMolSupplier(str(path_output))):
            molecule.SetProp('pose', str(i))
            docked_core_fragments.append(molecule)
    
    return docked_core_fragments

def template_docking(path_fragment, path_template, path_config, path_output, print_output=False):
    """
    runs FlexX template-docking

    Returns
    ----------
    List of all poses (docking result) as Mol (Properties: docking-score, pose)

    Parameters
    ----------
    path_fragment: pathlib.path
        Path to core-fragment sdf-file
    path_template: pathlib.path
        Path to template-fragment sdf-file
    path_config: pathlib.path
        Path to FlexX-config file.
    path_output: pathlib.path
        Path to output file
    template_pose: int
        ID of template pose  
    """
    # template docking
    output_text = subprocess.run(
            [
                str(PATH_FLEXX),
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
    if not os.path.exists(str(path_output)):
        return docked_fragments
    if os.stat(str(path_output)).st_size:
        for i, molecule in enumerate(Chem.SDMolSupplier(str(path_output))):
            molecule.SetProp('pose', str(i))
            docked_fragments.append(molecule)
    
    return docked_fragments

#   ===================  Class-based ================

class Pose:
    def __init__(self, ROMol, docking_score, hyde_score = None):
        self.docking_score = docking_score
        self.ROMol = ROMol
        self.hyde_score = hyde_score
    def _is_valid_operand(self, other):
        return hasattr(other, 'docking_score')
    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.hyde_score != None and other.hyde_score != None:
            return self.hyde_score < other.hyde_score
        else:
            return self.docking_score < other.docking_score
    def __le__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.hyde_score != None and other.hyde_score != None:
            return self.hyde_score <= other.hyde_score
        else:
            return self.docking_score <= other.docking_score
    def __gt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.hyde_score != None and other.hyde_score != None:
            return self.hyde_score > other.hyde_score
        else:
            return self.docking_score > other.docking_score
    def __ge__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.hyde_score != None and other.hyde_score != None:
            return self.hyde_score > other.hyde_score
        else:
            return self.docking_score > other.docking_score
    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.hyde_score != None and other.hyde_score != None:
            return self.hyde_score == other.hyde_score
        else:
            return self.docking_score == other.docking_score
    def __ne__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.hyde_score != None and other.hyde_score != None:
            return self.hyde_score != other.hyde_score
        else:
            return self.docking_score != other.docking_score

class Ligand:
    def __init__(self, ROMol, fragment_ids, recombination):
        self.ROMol = ROMol
        self.poses = []
        self.fragment_ids = fragment_ids
        self.dummy_atoms = {}
        self.min_docking_score = None
        self.recombinations = []
        self.recombination = recombination
    def to_sdf(self, sdf_path, pH=7.4):
        """
        Convert a SMILES string to a sdf file needed by docking programs of the AutoDock family.

        Parameters
        ----------
        sdf_path: str or pathlib.path
            Path to output sdf file.
        pH: float
            Protonation at given pH.
        """
        # TODO whats about PH-protonation
        # TODO is ROMol_original with dummy?
        molecule = Chem.AddHs(self.ROMol)
        status = AllChem.EmbedMolecule(molecule)
        status = AllChem.UFFOptimizeMolecule(molecule)
        w = Chem.SDWriter(str(sdf_path))
        w.write(molecule)
        return
    def add_pose(self, pose : Pose):
        if self.min_docking_score == None or pose.docking_score < self.min_docking_score:
            self.min_docking_score = pose.docking_score
        self.poses.append(pose)
        return
    def choose_template_poses(self, num_templates=None):
        """
        Chooses *num_templates*-best poses according to docking result and removes the non-top-scored or invalid poses

        Parameters
        ----------
        num_templates: int
            max. templates to choose, if None, num_templates == len(path)
            Path to output file
        """
        if not len(self.poses):
            return
        if not num_templates:
            num_templates = len(self.poses)

        self.poses.sort()
        
        # Always choose pose with the best docking score
        choosen_poses = [self.poses[0]]
        self.poses = self.poses[1:]
        for pose in self.poses:
                if choosen_poses[0].ROMol.GetNumAtoms() != pose.ROMol.GetNumAtoms():
                    print("Couldn't match pose")
                    pose.ROMol = choosen_poses[0].ROMol
        atom_mapping = [[j, j] for j in range(choosen_poses[0].ROMol.GetNumAtoms())]
        for pose in self.poses:
            if choosen_poses[0].ROMol.GetNumAtoms() != pose.ROMol.GetNumAtoms():
                print("Couldn't match pose")
                pose.ROMol = choosen_poses[0].ROMol
        # Append num_templates - 1 more templates according to highest RMSD to already choosen poses and lowest docking-score
        for _ in range(num_templates - 1):
            # TODO make this more efficient
            if not len(self.poses):
                break
            self.poses.sort(key= lambda pose : -2 * sum(rdMolAlign.CalcRMS(pose.ROMol, x.ROMol, map=[atom_mapping]) for x in choosen_poses) / len(choosen_poses) + pose.docking_score)
            choosen_poses.append(self.poses[0])
            self.poses = self.poses[1:]
        self.poses = choosen_poses
    def calculate_missing_dummy_atoms(self, fragment_library):
        for subpocket, id in self.fragment_ids.items():
            if subpocket not in self.dummy_atoms.keys():
                fragment = Chem.RemoveHs(fragment_library[subpocket]['ROMol_original'][id])
                self.dummy_atoms[subpocket] = [(f"{subpocket}_{i}", a.GetNeighbors()[0].GetProp('environment'), a.GetProp('subpocket')) for i, a in enumerate(fragment.GetAtoms()) if a.GetSymbol() == '*']
    def recombine(self, fragment_id, subpocket, fragment_library):
        self.calculate_missing_dummy_atoms(fragment_library)
        fragment = Chem.RemoveHs(fragment_library[subpocket]['ROMol_original'][fragment_id])
        dummy_atoms = [(f"{subpocket}_{i}", a.GetNeighbors()[0].GetProp('environment'), a.GetProp('subpocket')) for i, a in enumerate(fragment.GetAtoms()) if a.GetSymbol() == '*']
        new_rec = None
        for sp in self.dummy_atoms.keys():
            matching_dummies = [(id, env, con) for id, env, con in dummy_atoms if con == sp]
            matching_dummies_2 = [(id, env, con) for id, env, con in self.dummy_atoms[sp] if con == subpocket]
            if len(matching_dummies) != 1 or len(matching_dummies_2) != 1:
                continue
            id, env, con = matching_dummies[0]
            id_2, env_2, con_2 = matching_dummies_2[0]
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
def from_recombination(recombination):
    fragments = {x[:2]: int(x[3:]) for x in recombination.fragments}
    return Ligand(recombination.ligand, fragments, recombination)

class Recombination:
    def __init__(self,fragment_ids, bonds):
        self.ligand = None
        self.bonds = bonds
        self.fragments = fragment_ids
    def add_fragment(self, fragment_id, bonds):
        if fragment_id not in self.fragments:
            self.fragments.append(fragment_id)
        self.bonds += bonds
    def construct(self, fragment_library):
        print("===> recombining:", self.fragments, self.bonds)
        self.ligand = utils.construct_ligand(self.fragments, self.bonds, fragment_library)
    def copy(self):
        return Recombination(self.fragments.copy(), self.bonds.copy())
    