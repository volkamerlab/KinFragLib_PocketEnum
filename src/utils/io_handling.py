from pathlib import Path
import os

from rdkit import Chem

from classes.ligand import Ligand


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


def write_violations_to_file(violations: list, path_output: Path):
    """
    Writes all violations (tuple of mols) to the output file

    Parameters
    ----------
    molecules: list[tuple(mol, mol)]
        List containing all violations
    path_output: pathlib.path
        Path to output file
    """
    with Chem.SDWriter(str(path_output)) as w:
        for i, violation in enumerate(violations):
            violation[0].SetProp("_Name", f"V{i:08}_pre")
            violation[1].SetProp("_Name", f"V{i:08}_opt")
            w.write(violation[0])
            w.write(violation[1])


def append_ligands_to_file(
    ligands: "list[Ligand]",
    ligands_filtered: "list[Ligand]",
    path_output: Path,
    condition=lambda l: True,
):
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

    with open(
        str(path_output), "at"
    ) as outf:  # needed to be able to APPEND molecules to the result file
        with Chem.SDWriter(outf) as w:
            ligand: Ligand
            for ligand in ligands:
                if not condition(ligand):
                    continue
                mol = ligand.get_best_pose().ROMol
                mol.SetProp("_Name", ligand.compute_unique_id())
                if ligand in ligands_filtered:
                    mol.SetProp("filtered", "0")
                else:
                    mol.SetProp("filtered", "1")
                w.write(mol)


def write_all_poses_to_file(
    molecules: "list[Ligand]", ligands_filtered: "list[Ligand]", path_output
):
    """
    Writes all poses of the given molecules to the output file

    Parameters
    ----------
    molecules: list[Ligand]
        List containing all molecules
    path_output: pathlib.path
        Path to output file
    """
    with Chem.SDWriter(str(path_output)) as w:
        ligand: Ligand
        for ligand in molecules:
            # set name
            ligand_id = ligand.compute_unique_id()
            for i, pose in enumerate(ligand.poses_pre_filtered or ligand.poses):
                pose.ROMol.SetProp("_Name", ligand_id + f"P{i:04d}")
                if ligand in ligands_filtered:
                    # molecule was chosen
                    if pose in ligand.poses:
                        pose.ROMol.SetProp("filtered", "0")  # pose was chosen
                    else:
                        pose.ROMol.SetProp("filtered", "2")  # pose was not chosen
                else:
                    # molecule was not chosen
                    pose.ROMol.SetProp("filtered", "1")  # pose was chosen
                w.write(pose.ROMol)
