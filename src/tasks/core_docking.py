from pathlib import Path
import threading
import logging
from rdkit import Chem
import os
import subprocess

from classes.config import Config
from classes.ligand import Pose, Ligand
from utils.io_handling import remove_files
from tasks._utils import hyde_scoring
from utils._utils import calc_distance_matrix


def core_docking_task(config: Config, core_fragment):
    """
    runs a FlexX core docking task. It should be used for multithreading.

    Returns
    ----------
    List of all poses (docking result) as Ligand (Properties: docking-score, pose)

    Parameters
    ----------
    config: Config
        Object, containing all neccessary information
    core_fragment: Ligand
        Fragment that should be docked in this task

    """
    thread_id = threading.get_ident()
    logging.debug(
        "Docking of "
        + config.core_subpocket
        + "-Fragment: "
        + str(core_fragment.fragment_ids[config.core_subpocket])
    )
    # safe fragment as sdf file
    if not core_fragment.to_sdf(
        config.path_temp / ("cand_thread_" + str(thread_id) + "_core_fragment.sdf")
    ):
        # could not generate 3D-conformation
        logging.error(
            "Could not write Fragemnt: "
            + str(core_fragment.fragment_ids)
            + " to files due to 3d-generation-error"
        )
        return [core_fragment]

    res_docking = _core_docking(
        config.path_temp / ("cand_thread_" + str(thread_id) + "_core_fragment.sdf"),
        config.path_structure_config / (config.core_subpocket + ".flexx"),
        config.path_temp / ("dock_thread_" + str(thread_id) + "_core_fragment.sdf"),
        config.path_flexx,
        core_fragment,
    )

    violations = []  # track violations
    sum_hyde_displacement_violations = 0  # used to calculate the mean displacement of hyde
    sum_hyde_displacement = 0  # used to calculate the mean displacement of hyde

    if config.use_hyde and len(res_docking):
        # if path to hyde is given: perform hyde_scoring and opt.
        res_hyde = hyde_scoring(
            config.path_temp / ("dock_thread_" + str(thread_id) + "_core_fragment.sdf"),
            config.path_structure_config / (config.core_subpocket + ".hydescorer"),
            config.path_temp / ("hyde_thread_" + str(thread_id) + "_core_fragment.sdf"),
            config.path_hyde,
            core_fragment,
        )

        remove_files(
            config.path_temp / ("dock_thread_" + str(thread_id) + "_core_fragment.sdf"),
            config.path_temp / ("cand_thread_" + str(thread_id) + "_core_fragment.sdf"),
            config.path_temp / ("hyde_thread_" + str(thread_id) + "_core_fragment.sdf"),
        )

        for i, (conformer_docking, conformer_hyde) in enumerate(
            zip(res_docking, res_hyde)
        ):  # safe every resulting pose within the fragment
            displacement = calc_distance_matrix([conformer_hyde, conformer_docking])[0]
            if displacement > config.hyde_displacement_cutoff:
                logging.warning(
                    f"""VIOLATION: RMSD between HYDE and docking pose {core_fragment.fragment_ids} {i}: 
                                {displacement}"""
                )
                violations.append((conformer_docking, conformer_hyde))
                sum_hyde_displacement_violations += displacement
                # drop pose if rmds > hyde_cutoff
                continue
            sum_hyde_displacement += displacement
            pose = Pose(conformer_hyde, float(conformer_hyde.GetProp("BIOSOLVEIT.DOCKING_SCORE")))
            pose.binding_affinity_upper = float(
                conformer_hyde.GetProp("BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_UPPER_BOUNDARY [nM]")
            )
            pose.binding_affinity_lower = float(
                conformer_hyde.GetProp("BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_LOWER_BOUNDARY [nM]")
            )
            core_fragment.add_pose(pose)
    else:
        remove_files(
            config.path_temp / ("dock_thread_" + str(thread_id) + "_core_fragment.sdf"),
            config.path_temp / ("cand_thread_" + str(thread_id) + "_core_fragment.sdf"),
        )

        for conformer_docking in res_docking:  # safe every resulting pose within the fragment
            pose = Pose(
                conformer_docking, float(conformer_docking.GetProp("BIOSOLVEIT.DOCKING_SCORE"))
            )
            core_fragment.add_pose(pose)

    logging.debug(str(len(core_fragment.poses)) + " poses have been generated")

    if len(core_fragment.poses):
        # if fragment could be docked, save fragment (including it's poses)
        logging.debug(
            "Best score: "
            + str(core_fragment.min_binding_affinity or core_fragment.min_docking_score)
        )
        core_fragment.mean_hyde_displacement_undropped = sum_hyde_displacement / len(
            core_fragment.poses
        )
        core_fragment.mean_hyde_displacement = (
            sum_hyde_displacement + sum_hyde_displacement_violations
        ) / (len(core_fragment.poses) + len(violations))
        core_fragment.num_hyde_violations = len(violations)
        return [core_fragment, violations]
    # return empty list if no pose was found
    return []


def _core_docking(
    path_fragment: Path,
    path_config: Path,
    path_output: Path,
    path_flexx: Path,
    ligand: Ligand,
    print_output=False,
) -> list:
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
    ligand: Ligand
        Ligand object of molecule to dock
    """
    # core docking
    output_text = subprocess.run(
        [
            str(Path(".") / path_flexx),
            "-i",
            str(path_fragment),
            "--docking-definition",
            str(path_config),
            "-o",
            str(path_output),
        ],
        capture_output=True,
    )
    if print_output:
        print(output_text.stderr)

    # read core fragments poses from sdf
    docked_core_fragments = []
    if os.stat(
        str(path_output)
    ).st_size:  # only acces reult file if at least one pose was generated
        for molecule in Chem.SDMolSupplier(str(path_output)):
            molecule.SetProp("fragment_ids", str(ligand.fragment_ids))
            molecule.SetProp("smiles_ligand", Chem.MolToSmiles(molecule))
            molecule.SetProp("smiles_fragments_dummy", str(ligand.smiles_dummy))
            molecule.SetProp("smiles_fragments", str(ligand.smiles))
            docked_core_fragments.append(molecule)

    return docked_core_fragments
