import threading
import os
import subprocess
import logging
from pathlib import Path
from rdkit import Chem

from classes.config import Config
from classes.recombination import Recombination
from classes.ligand import Pose, Ligand
from classes import ligand
from utils._utils import calc_distance_matrix
from utils.io_handling import remove_files
from tasks._utils import hyde_scoring


def template_docking_task(
    config: Config, subpocket: str, recombination: Recombination, poses: list
):
    """
    runs a FlexX template docking task. It should be used for multithreading.

    Returns
    ----------
    List of all poses (docking result) as Ligand (Properties: docking-score, pose)

    Parameters
    ----------
    config: Config
        Object, containing all neccessary information
    subpocket: str
        Subpocket where to place the fragment
    recombination: Recombination
        Recombination that should be docked in this task
    poses: List(Pose)
        Poses that are used as templates

    """
    thread_id = threading.get_ident()
    logging.debug("Template docking of Recombination: " + str(recombination.fragments))
    fragment = ligand.from_recombination(recombination)

    # write recombination that should be docked to file
    if not fragment.to_sdf(
        config.path_temp / ("cand_thread_" + str(thread_id) + "_" + subpocket + "_fragment.sdf")
    ):
        logging.error(
            "Could not write fragemnt: "
            + str(fragment.fragment_ids)
            + " to files due to 3d-generation-error"
        )
        return [fragment]

    sum_hyde_displacement_violations = 0  # used to calculate the mean displacement of hyde
    sum_hyde_displacement = 0  # used to calculate the mean displacement of hyde
    violations = []
    succesfully_docked = 0

    # for every choosen pose: perform template docking
    for i, pose in enumerate(poses):
        logging.debug("Template " + str(i + 1) + " / " + str(len(poses)))
        # write template to sdf file
        with Chem.SDWriter(
            str(
                config.path_temp
                / ("temp_thread_" + str(thread_id) + "_" + subpocket + "_fragment.sdf")
            )
        ) as w:
            w.write(pose.ROMol)
        # template docking (FlexX)
        res_docking = _template_docking(
            config.path_temp
            / ("cand_thread_" + str(thread_id) + "_" + subpocket + "_fragment.sdf"),
            config.path_temp
            / ("temp_thread_" + str(thread_id) + "_" + subpocket + "_fragment.sdf"),
            config.path_structure_config / (subpocket + ".flexx"),
            config.path_temp / ("dock_thread_" + str(thread_id) + "_" + "fragments.sdf"),
            config.path_flexx,
            fragment,
        )

        succesfully_docked += 1 if len(res_docking) else 0

        if len(res_docking) and config.use_hyde:
            res_hyde = hyde_scoring(
                config.path_temp / ("dock_thread_" + str(thread_id) + "_" + "fragments.sdf"),
                config.path_structure_config / (subpocket + ".hydescorer"),
                config.path_temp / ("hyde_thread_" + str(thread_id) + "_fragment.sdf"),
                config.path_hyde,
                fragment,
            )
            # remove files containg docking results and template
            remove_files(
                config.path_temp / ("dock_thread_" + str(thread_id) + "_" + "fragments.sdf"),
                config.path_temp
                / ("temp_thread_" + str(thread_id) + "_" + subpocket + "_fragment.sdf"),
                config.path_temp / ("hyde_thread_" + str(thread_id) + "_fragment.sdf"),
            )

            # safe resulting poses within fragment
            for i, (conformer_hyde, conformer_docking) in enumerate(
                zip(res_hyde, res_docking)
            ):  # safe every resulting pose within the fragment
                pose = Pose(
                    conformer_hyde, float(conformer_hyde.GetProp("BIOSOLVEIT.DOCKING_SCORE"))
                )
                pose.binding_affinity_upper = float(
                    conformer_hyde.GetProp(
                        "BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_UPPER_BOUNDARY [nM]"
                    )
                )
                pose.binding_affinity_lower = float(
                    conformer_hyde.GetProp(
                        "BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_LOWER_BOUNDARY [nM]"
                    )
                )
                rmsd = calc_distance_matrix([conformer_hyde, conformer_docking])[0]

                if rmsd > config.hyde_displacement_cutoff:
                    violations.append((conformer_docking, conformer_hyde))
                    sum_hyde_displacement_violations += rmsd
                    logging.debug(
                        f"VIOLATION: RMSD between HYDE and docking pose {fragment.fragment_ids} {i}: {rmsd}"
                    )
                else:
                    sum_hyde_displacement += rmsd
                    fragment.add_pose(pose)
        else:
            # remove files containg docking results and template
            remove_files(
                config.path_temp / ("dock_thread_" + str(thread_id) + "_" + "fragments.sdf"),
                config.path_temp
                / ("temp_thread_" + str(thread_id) + "_" + subpocket + "_fragment.sdf"),
            )
            for conformer in res_docking:  # safe every resulting pose within the fragment
                pose = Pose(conformer, float(conformer.GetProp("BIOSOLVEIT.DOCKING_SCORE")))
                fragment.add_pose(pose)

        logging.debug(str(len(fragment.poses)) + " poses have been generated")

        if len(fragment.poses):
            logging.debug(
                "Best  score: " + str(fragment.min_binding_affinity or fragment.min_docking_score)
            )

    # remove molecule file
    remove_files(
        config.path_temp / ("cand_thread_" + str(thread_id) + "_" + subpocket + "_fragment.sdf")
    )

    if len(fragment.poses):
        fragment.mean_hyde_displacement_undropped = sum_hyde_displacement / len(fragment.poses)
        fragment.mean_hyde_displacement = (
            sum_hyde_displacement + sum_hyde_displacement_violations
        ) / (len(fragment.poses) + len(violations))
        fragment.num_hyde_violations = len(violations)
        # safe recombination as result only if at least one pose was generated
        return [fragment, violations, succesfully_docked]
    return []


def _template_docking(
    path_fragment: Path,
    path_template: Path,
    path_config: Path,
    path_output: Path,
    path_flexx: Path,
    ligand: Ligand,
    print_output=False,
) -> list:
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
        Path to Flexx
    path_output: pathlib.path
        Path to output file
    ligand: Ligand
        Ligand object of molecule to dock
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
            # "--stereo-mode"
            # ,"3"
        ],
        capture_output=True,  # needed to capture output text
    )
    if print_output:
        print(output_text.stderr)

    # read fragments poses from sdf
    docked_fragments = []
    if os.stat(str(path_output)).st_size:
        for molecule in Chem.SDMolSupplier(str(path_output)):
            molecule.SetProp("fragment_ids", str(ligand.fragment_ids))
            molecule.SetProp("smiles_ligand", Chem.MolToSmiles(molecule))
            molecule.SetProp("smiles_fragments_dummy", str(ligand.smiles_dummy))
            molecule.SetProp("smiles_fragments", str(ligand.smiles))
            docked_fragments.append(molecule)

    return docked_fragments
