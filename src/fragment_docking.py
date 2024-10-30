from kinfraglib import utils, filters
from rdkit import Chem
from rdkit.Chem import rdFMCS
import logging
import time
import sys
import wandb
import argparse

import json
from concurrent.futures import ThreadPoolExecutor, as_completed

from classes.config import Config
from tasks.core_docking import core_docking_task
from tasks.template_docking import template_docking_task
from tasks._utils import prepare_core_fragments
from tasks.compound_filtering import cluster_based_compound_filtering
from utils._utils import seed_everything

from utils.io_handling import (
    write_all_poses_to_file,
    write_violations_to_file,
    append_ligands_to_file,
)

if __name__ == "__main__":
    # parse command line arguments
    parser = argparse.ArgumentParser(
        prog=sys.argv[0], description="Generates compounds for a given kinase"
    )
    parser.add_argument(
        "-d",
        "--definitions",
        default="definitions.json",
        help="JSON file with program configuration",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="output.json",
        help="Name of output JSON file with program statistics",
    )
    parser.add_argument(
        "-r", "--results", default="results", help="Folder, where results are placed"
    )
    parser.add_argument(
        "-log",
        "--loglevel",
        default="INFO",
        help="Logging level (error, warning, info, or debug). Example --loglevel debug, default=info",
    )

    args = parser.parse_args()

    # parse configs
    config = Config()
    config.parse(args.definitions)
    config.initialize_folders(args.results)

    # init logging
    numeric_level = getattr(logging, args.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.loglevel)

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=numeric_level,
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    wandb.init(
        # Set the project where this run will be logged
        name=config.pdb_code,
        project="subpocket_based_docking_kinases",
    )

    start_time_all = time.time()

    # prepare directory storing all inforamtion for output json log file
    output_logs = {
        "pdbCode": config.pdb_code,
        "CoreSubpocket": config.core_subpocket,
        "Subpockets": config.subpockets,
        "UsedHyde": config.use_hyde,
        "HydeDisplacementCutoff": (
            config.hyde_displacement_cutoff if config.use_hyde else None
        ),
        "NumberPosesPerFragment": config.poses_per_fragment,
        "NumberFragmentsPerIterations": config.fragments_per_iteration,
        "UseClusterBasedPoseFiltering": config.cluster_based,
        "DistanceThresholdClustering": (
            config.cluster_threshold if config.cluster_based else None
        ),
        "Filters": [filter.name for filter in config.filters],
        "GeneratedPoses": {},
        "ChosenPoses": {},
        "ChosenMolecules": {},
        "UndockableMolecules": {},
        "Recombinations": {},
        "RunTimeTotal": 0,
        "RunTimeCoreDocking": 0,
        "RunTimeTemplateDocking": {},
        "DockedMolecules": {},
        "MeanHydeDisplacement": {},
        "MeanHydeDisplacementIncludingViolations": {},
        "NumDisplacementViolations": {},
        "DockingRuns": {},
        "SuccesfullDockingRuns": {},
        "NumMultipleBondsBetweenFragmentAndLigand": {},
        "NumMultipleBondsBetweenTwoFragments": {},
    }

    # set seeds
    seed_everything(config.seed)

    # ==== PREPROCESSING =======

    logging.info("Preprocessing started")

    # Access fragment library
    fragment_library_original = utils.read_fragment_library(
        config.path_kinfraglib / "fragment_library"
    )

    logging.debug(
        "KinFragLib has been loaded: "
        + str(
            [
                sp + ": " + str(len(fragment_library_original[sp]))
                for sp in fragment_library_original.keys()
            ]
        )
    )

    logging.debug("Start prefiltering")
    # removing fragments in pool X
    # removing duplicates  # TODO according to smiles with dummy
    # removing fragments without dummy atoms (unfragmented ligands)
    # removing fragments only connecting to pool X
    fragment_library = filters.prefilters.pre_filters(fragment_library_original)
    logging.debug(
        "Finished prefiltering: "
        + str(
            [
                sp + ": " + str(len(fragment_library[sp]))
                for sp in fragment_library.keys()
            ]
        )
    )

    # apply filters
    for filter in config.filters:
        logging.debug("Applying filter: " + filter.name)
        l_before = [
            sp + ": " + str(len(fragment_library[sp])) for sp in fragment_library.keys()
        ]
        fragment_library = filter.apply_filter(fragment_library)
        logging.debug(
            filter.name
            + " applied: "
            + str(
                [
                    sp + ": " + str(len(fragment_library[sp]))
                    for sp in fragment_library.keys()
                ]
            )
        )

    # determine size of fragment library and log it
    library_size = {sp: len(fragment_library[sp]) for sp in fragment_library.keys()}
    output_logs["FragmentLibrarySize"] = library_size
    logging.info("Preprocessing finished")
    logging.info("Size of fragment library" + str(library_size))

    # ===== CORE DOCKING =======

    logging.info(
        "Core docking of "
        + str(len(fragment_library[config.core_subpocket]))
        + " "
        + config.core_subpocket
        + "-Fragments"
    )

    if config.core_subpocket not in fragment_library:
        logging.error("Fragment library must at least contain 1 core fragment")
        raise Exception(
            "Fragment library must at least contain 1 core fragment but was 0"
        )

    core_fragments = prepare_core_fragments(fragment_library, config)

    # core docking
    docking_results = []

    violations = []
    unsuccesfull_3d_generations = []
    start_time = time.time()

    wandb.log({"total_frags_SP0": len(core_fragments)})
    docking_run_counter = 0

    with ThreadPoolExecutor(config.num_threads) as executor:
        # submit core docking tasks
        features = [
            executor.submit(core_docking_task, config, core_fragment)
            for core_fragment in core_fragments
        ]
        # iterate over all completetd tasks
        for feature in as_completed(features):
            try:
                # get result
                result = feature.result()
            except Exception as exc:
                logging.error("Generated an exception during core_docking: %s" % (exc))
            else:
                if len(result) == 2:
                    docking_results.append(result[0])
                    violations += result[1]
                elif len(result) == 1:
                    unsuccesfull_3d_generations += result
                docking_run_counter += 1
                wandb.log({"docked_frags_SP0": docking_run_counter})

    # core docking logs
    output_logs["RunTimeCoreDocking"] = time.strftime(
        "%H:%M:%S", time.gmtime(round(time.time() - start_time, 2))
    )
    output_logs["GeneratedPoses"]["SP0"] = sum(
        len(mol.poses) for mol in docking_results
    )
    output_logs["UndockableMolecules"]["SP0"] = len(core_fragments) - len(
        docking_results
    )
    output_logs["DockedMolecules"]["SP0"] = len(docking_results)
    output_logs["MeanHydeDisplacement"]["SP0"] = sum(
        mol.mean_hyde_displacement_undropped for mol in docking_results
    ) / (len(docking_results) or 1)
    output_logs["MeanHydeDisplacementIncludingViolations"]["SP0"] = sum(
        mol.mean_hyde_displacement for mol in docking_results
    ) / (len(docking_results) or 1)
    output_logs["NumDisplacementViolations"]["SP0"] = sum(
        mol.num_hyde_violations for mol in docking_results
    )
    output_logs["DockingRuns"]["SP0"] = len(core_fragments)
    output_logs["SuccesfullDockingRuns"]["SP0"] = len(docking_results)

    # write violations so sdf file
    write_violations_to_file(violations, config.path_results / "violations_SP0.sdf")

    logging.info("Runtime: %s" % (time.time() - start_time))

    # ===== TEMPLATE DOCKING ======

    for subpocket in config.subpockets:
        if subpocket not in fragment_library:
            logging.warning(
                f"Fragment library contains 0 {subpocket}-fragments, continue with next subpocket"
            )
            continue

        logging.info(
            "Template docking of "
            + str(len(fragment_library[subpocket]))
            + " "
            + subpocket
            + "-Fragments"
        )

        # filtering
        if config.use_hyde:
            docking_results.sort(key=lambda l: l.min_binding_affinity)
        else:
            docking_results.sort(key=lambda l: l.min_docking_score)

        # save state before filtering
        docking_results_pre_filtering = docking_results.copy()

        # choose n best fragments
        docking_results = cluster_based_compound_filtering(
            docking_results, config.fragments_per_iteration, config, subpocket
        )

        # choose k best conformers (per fragment)
        for ligand in docking_results:
            if config.cluster_based:
                ligand.choose_template_poses_cluster_based(
                    config.poses_per_fragment, config.cluster_threshold
                )
            else:
                ligand.choose_template_poses(config.poses_per_fragment)

        # log chosen poses and molecules
        output_logs["ChosenMolecules"][
            "SP" + str(config.subpockets.index(subpocket))
        ] = len(docking_results)
        output_logs["ChosenPoses"]["SP" + str(config.subpockets.index(subpocket))] = (
            sum(len(mol.poses) for mol in docking_results)
        )

        # store the best pose of a docked molecule if it consists of more than 1 fragment in the overall result file
        append_ligands_to_file(
            docking_results_pre_filtering,
            docking_results,
            config.path_results / "results.sdf",
            lambda l: len(l.fragment_ids) > 1,
        )

        # store all poses of all docked molecules in one file for the current subpocket iteration
        write_all_poses_to_file(
            docking_results_pre_filtering,
            docking_results,
            config.path_results
            / ("SP" + str(config.subpockets.index(subpocket)) + ".sdf"),
        )

        # recombining

        candidates = []

        num_recombinations = 0

        counter_num_multpl_bonds = 0
        counter_num_unambigious_bonds = 0
        # try recombine every ligand (comb. of fragmnets) with every fragment of the current subpocket
        for ligand in docking_results:
            for fragment_idx in fragment_library[subpocket].index:
                num_multpl_bonds, num_unambigious_bonds = ligand.recombine(
                    fragment_idx, subpocket, fragment_library
                )  # possible recombinations are stored within ligand
                counter_num_multpl_bonds += num_multpl_bonds
                counter_num_unambigious_bonds += num_unambigious_bonds
            if len(ligand.recombinations):
                num_recombinations += len(ligand.recombinations)
                # use a ligand for template docking iff at least one recombination was produced for a ligand
                candidates.append(ligand)

        logging.debug("Generated " + str(num_recombinations) + " recombinations")
        output_logs["Recombinations"][
            "SP" + str(config.subpockets.index(subpocket) + 1)
        ] = num_recombinations
        output_logs["NumMultipleBondsBetweenFragmentAndLigand"][
            "SP" + str(config.subpockets.index(subpocket) + 1)
        ] = counter_num_multpl_bonds
        output_logs["NumMultipleBondsBetweenTwoFragments"][
            "SP" + str(config.subpockets.index(subpocket) + 1)
        ] = counter_num_unambigious_bonds
        # template docking

        docking_results = []
        violations = []

        num_docking_runs = sum(len(c.poses) * len(c.recombinations) for c in candidates)
        num_succ_docking_runs = 0

        start_time = time.time()

        wandb.log(
            {
                "total_frags_SP"
                + str(config.subpockets.index(subpocket) + 1): num_recombinations
            }
        )
        docking_run_counter = 0

        with ThreadPoolExecutor(config.num_threads) as executor:
            # submit template docking tasks
            features = [
                executor.submit(
                    template_docking_task,
                    config,
                    subpocket,
                    recombination,
                    ligand.poses,
                )
                for ligand in candidates
                for recombination in ligand.recombinations
            ]
            # iterate over all completetd tasks
            for count, feature in enumerate(as_completed(features)):
                try:
                    # get docking result
                    result = feature.result()
                except Exception as exc:
                    logging.error(
                        "Generated an exception during template_docking: %s" % (exc)
                    )
                else:
                    docking_run_counter += 1
                    wandb.log(
                        {
                            "docked_frags_SP"
                            + str(
                                config.subpockets.index(subpocket) + 1
                            ): docking_run_counter
                        }
                    )
                    if len(result) == 3:
                        # docking was succesfull
                        docking_results.append(result[0])
                        violations += result[1]
                        num_succ_docking_runs += result[2]
                    elif len(result) == 1:
                        # could not genrate the 3d conformation
                        unsuccesfull_3d_generations += result
                        num_docking_runs -= len(result[0].poses) * len(
                            result[0].recombinations
                        )

        logging.info(
            f"Runtime template-docking ({subpocket}): {(time.time() - start_time)}"
        )

        # template docking logs
        output_logs["RunTimeTemplateDocking"][
            "SP" + str(config.subpockets.index(subpocket) + 1)
        ] = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
        output_logs["GeneratedPoses"][
            "SP" + str(config.subpockets.index(subpocket) + 1)
        ] = sum(len(mol.poses) for mol in docking_results)
        output_logs["UndockableMolecules"][
            "SP" + str(config.subpockets.index(subpocket) + 1)
        ] = num_recombinations - len(docking_results)
        output_logs["DockedMolecules"][
            "SP" + str(config.subpockets.index(subpocket) + 1)
        ] = len(docking_results)
        output_logs["MeanHydeDisplacement"][
            "SP" + str(config.subpockets.index(subpocket) + 1)
        ] = sum(mol.mean_hyde_displacement_undropped for mol in docking_results) / (
            len(docking_results) or 1
        )
        output_logs["MeanHydeDisplacementIncludingViolations"][
            "SP" + str(config.subpockets.index(subpocket) + 1)
        ] = sum(mol.mean_hyde_displacement for mol in docking_results) / (
            len(docking_results) or 1
        )
        output_logs["NumDisplacementViolations"][
            "SP" + str(config.subpockets.index(subpocket) + 1)
        ] = sum(mol.num_hyde_violations for mol in docking_results)
        output_logs["DockingRuns"][
            "SP" + str(config.subpockets.index(subpocket) + 1)
        ] = num_docking_runs
        output_logs["SuccesfullDockingRuns"][
            "SP" + str(config.subpockets.index(subpocket) + 1)
        ] = num_succ_docking_runs
        # write violations so sdf file
        write_violations_to_file(
            violations,
            config.path_results
            / ("violations_SP" + str(config.subpockets.index(subpocket) + 1) + ".sdf"),
        )

    # store all poses of all docked molecules in one file for the current subpocket iteration
    write_all_poses_to_file(
        docking_results,
        docking_results,
        config.path_results / ("SP" + str(len(config.subpockets)) + ".sdf"),
    )

    output_logs["Unsuccesful3DGenerations"] = [
        Chem.MolToSmiles(mol.ROMol) for mol in unsuccesfull_3d_generations
    ]

    # log total run time
    output_logs["RunTimeTotal"] = time.strftime(
        "%H:%M:%S", time.gmtime(round(time.time() - start_time_all, 2))
    )
    logging.info("Runtime: %s" % (time.time() - start_time_all))

    # ===== EVALUATION ======

    # write logs to json output file
    with open(args.output, "w") as json_file:
        json.dump(output_logs, json_file, indent=4)

    if len(docking_results):
        if config.use_hyde:
            docking_results.sort(key=lambda l: l.min_binding_affinity)
        else:
            docking_results.sort(key=lambda l: l.min_docking_score)

        logging.debug(
            "Best recombination: "
            + str(list(docking_results[0].fragment_ids.items()))
            + " Score: "
            + str(
                docking_results[0].min_binding_affinity
                or docking_results[0].min_docking_score
            )
        )

        # store generated ligands
        append_ligands_to_file(
            docking_results,
            docking_results,
            config.path_results / "results.sdf",
            lambda l: len(l.fragment_ids) > 1,
        )
