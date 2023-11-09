from pathlib import Path
from kinfraglib import utils, filters
import docking_utils
from rdkit import Chem
from pathlib import Path
import logging
import time
import os,sys
import wandb 
import argparse

import json
import threading_docking
from functools import partial
from concurrent.futures import ThreadPoolExecutor, as_completed

if __name__ == '__main__':


    # parse command line arguments
    parser = argparse.ArgumentParser(
                    prog=sys.argv[0],
                    description='Generates compounds for a given kinase')
    parser.add_argument('-d', '--definitions', default='definitions.json', help='JSON file with program configuration')
    parser.add_argument('-o', '--output', default='output.json', help='Name of output JSON file with program statistics')
    parser.add_argument('-r', '--results', default='results', help='Folder, where results are placed')

    args = parser.parse_args()

    # load program definitions
    with open(args.definitions, 'rt') as json_file:
        definitions = json.load(json_file)

    # Definitions
    HERE = Path().resolve()
    PATH_TEMP = HERE / 'temp' / definitions['pdbCode']
    PATH_DATA = Path(definitions['KinFragLib'])
    PATH_FLEXX = Path(definitions['FlexX'])
    PATH_TO_CONFIGS = Path(definitions['Config']) / definitions['pdbCode']
    PATH_TO_RESULTS = HERE / args.results / definitions['pdbCode']
    PATH_HYDE = Path(definitions['Hyde']) if definitions['UseHyde'] else None
    use_hyde = definitions['UseHyde']
    cutoff_hyde_displacement = definitions.get('HydeDisplacementCutoff') or 2.5     # rmds cutoff when a optimized pose can be droped 
    num_conformers = definitions['NumberPosesPerFragment']                      # amount of conformers to choose per docked fragment  (according to docking score and diversity)
    num_fragments_per_iterations = definitions['NumberFragmentsPerIterations']  # amount of fragments to choose per docking iteration (according to docking score)

    num_threads = definitions.get('NumberThreads')                              # number of threads to use, if not specified min(32, os.cpu_count() + 4) are used. If multithreading isn't wanted, NumberThreads should be set to 1
    cluster_based_pose_filtering = definitions['UseClusterBasedPoseFiltering']  
    clusted_pose_filter_dist_threshold = definitions.get('DistanceThresholdClustering') or 1.5 # threshold that should be used for pose clustering

    # init logging
    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    wandb.init(
        # Set the project where this run will be logged
        name=definitions['pdbCode'], 
        project="subpocket_based_docking_kinases")

    # create temp folder if it does not exists
    if not os.path.exists(HERE / 'temp'):
        os.mkdir(HERE / 'temp')
    if not os.path.exists(PATH_TEMP):
        os.mkdir(PATH_TEMP)

    if not os.path.exists(HERE / args.results):
        raise OSError(f"Result folder {str(HERE / args.results)} does not exists")
    
    if not os.path.exists(PATH_TO_RESULTS):
        os.mkdir(PATH_TO_RESULTS)
    
    with open(PATH_TO_RESULTS/ 'results.sdf', 'wt') as file:
        pass

    # define filters
    filters_ = [docking_utils.Filter(name, values) for name, values in definitions['Filters'].items()]

    start_time_all = time.time()

    # define pockets
    core_subpocket = definitions['CoreSubpocket']

    subpockets : list = definitions['Subpockets']

    # prepare directory storing all inforamtion for output json log file
    output_logs = {'pdbCode': definitions['pdbCode'],
                   'CoreSubpocket': core_subpocket,
                   'Subpockets': subpockets,
                   'UsedHyde': use_hyde, 'HydeDisplacementCutoff': cutoff_hyde_displacement if use_hyde else None, 
                   'NumberPosesPerFragment': num_conformers, 'NumberFragmentsPerIterations': num_fragments_per_iterations,
                   'UseClusterBasedPoseFiltering': cluster_based_pose_filtering, 'DistanceThresholdClustering': clusted_pose_filter_dist_threshold if cluster_based_pose_filtering else None, 
                   'Filters': definitions['Filters'],
                   'GeneratedPoses': {},
                   'ChosenPoses': {},
                   'ChosenMolecules': {},
                   'UndockableMolecules': {},
                   'Recombinations': {},
                   'RunTimeTotal': 0,
                   'RunTimeCoreDocking': 0,
                   'RunTimeTemplateDocking': {},
                   'DockedMolecules': {},
                   'MeanHydeDisplacement': {},
                   'MeanHydeDisplacementIncludingViolations': {},
                   'NumDisplacementViolations': {},
                   'DockingRuns': {},
                   'SuccesfullDockingRuns':{},
                   'NumMultipleBondsBetweenFragmentAndLigand': {},
                   'NumMultipleBondsBetweenTwoFragments': {}}

    # ==== PREPROCESSING =======

    logging.info('Preprocessing started')

    # Access fragment library
    fragment_library_original = utils.read_fragment_library(PATH_DATA / "fragment_library")

    logging.debug('KinFragLib has been loaded: ' + str([sp + ': ' + str(len(fragment_library_original[sp])) for sp in fragment_library_original.keys()]))

    logging.debug('Start prefiltering')
    # removing fragments in pool X
    # removing duplicates  # TODO according to smiles with dummy
    # removing fragments without dummy atoms (unfragmented ligands)
    # removing fragments only connecting to pool X
    fragment_library = filters.prefilters.pre_filters(fragment_library_original)
    logging.debug('Finished prefiltering: ' + str([sp + ': ' + str(len(fragment_library[sp])) for sp in fragment_library.keys()]))

    # apply filters
    for filter in filters_:
        logging.debug('Applying filter: ' + filter.name)
        l_before = [sp + ': ' + str(len(fragment_library[sp])) for sp in fragment_library.keys()]
        fragment_library = filter.apply_filter(fragment_library)
        logging.debug(filter.name + ' applied: ' + str([sp + ': ' + str(len(fragment_library[sp])) for sp in fragment_library.keys()]))

    # determine size of fragment library and log it
    library_size = {sp: len(fragment_library[sp]) for sp in fragment_library.keys()}
    output_logs['FragmentLibrarySize'] = library_size
    logging.info('Preprocessing finished\n Size of fragment library' + str(library_size))

    # ===== CORE DOCKING =======

    logging.info('Core docking of ' + str(len(fragment_library[core_subpocket])) + ' ' + core_subpocket + '-Fragments')

    core_fragments = []

    # prepare all core fragments
    for i in fragment_library[core_subpocket].index:
        smiles = fragment_library[core_subpocket]['smiles'][i]
        smiles_dummy = fragment_library[core_subpocket]['smiles_dummy'][i]
        core_fragments.append(docking_utils.Ligand(fragment_library[core_subpocket]['ROMol'][i], {core_subpocket: i}, 
                                                   docking_utils.Recombination([core_subpocket + "_" + str(i)], [], {core_subpocket: smiles}, {core_subpocket: smiles_dummy}), 
                                                   {core_subpocket: smiles_dummy}, {core_subpocket: smiles}))

    # core docking 
    docking_results = []

    violations = []
    unsuccesfull_3d_generations = []
    start_time = time.time()

    # create partial docking task to avoid function calls with many of 
    core_docking_task = partial(threading_docking.core_docking_task, PATH_TEMP,  PATH_TO_CONFIGS, PATH_TO_CONFIGS, PATH_FLEXX, PATH_HYDE, cutoff_hyde_displacement)

    wandb.log({'total_frags_SP0': len(core_fragments)})
    docking_run_counter = 0

    with ThreadPoolExecutor(num_threads) as executor:
        # submit core docking tasks
        features = [executor.submit(core_docking_task, core_subpocket, core_fragment) for core_fragment in core_fragments]
        # iterate over all completetd tasks
        for feature in as_completed(features):
            try:
               # get result
               result = feature.result()
            except Exception as exc:
                logging.error('Generated an exception during core_docking: %s' % (exc)) 
            else:
                if len(result) == 2:
                    docking_results.append(result[0])
                    violations += result[1]
                elif len(result) == 1:
                    unsuccesfull_3d_generations += result
                docking_run_counter += 1
                wandb.log({"docked_frags_SP0": docking_run_counter})


    # core docking logs
    output_logs['RunTimeCoreDocking'] =   time.strftime('%H:%M:%S', time.gmtime(round(time.time() - start_time, 2)))
    output_logs['GeneratedPoses']['SP0'] = sum(len(mol.poses) for mol in docking_results)
    output_logs['UndockableMolecules']['SP0'] = len(core_fragments) - len(docking_results)
    output_logs['DockedMolecules']['SP0'] = len(docking_results)
    output_logs['MeanHydeDisplacement']['SP0'] = sum(mol.mean_hyde_displacement_undropped for mol in docking_results) / (len(docking_results) or 1)
    output_logs['MeanHydeDisplacementIncludingViolations']['SP0'] = sum(mol.mean_hyde_displacement for mol in docking_results) / (len(docking_results) or 1)
    output_logs['NumDisplacementViolations']['SP0'] = sum(mol.num_hyde_violations for mol in docking_results)
    output_logs['DockingRuns']['SP0'] = len(core_fragments)
    output_logs['SuccesfullDockingRuns']['SP0'] = len(docking_results)

    # write violations so sdf file
    docking_utils.write_violations_to_file(violations, PATH_TO_RESULTS / 'violations_SP0.sdf')

    logging.info("Runtime: %s" % (time.time() - start_time))

    # ===== TEMPLATE DOCKING ======

    for subpocket in subpockets:
        logging.info('Template docking of ' + str(len(fragment_library[subpocket])) + ' ' + subpocket + '-Fragments')

        # filtering
        if use_hyde:
            docking_results.sort(key=lambda l: l.min_binding_affinity)
        else:
            docking_results.sort(key=lambda l: l.min_docking_score)

        # save state before filtering
        docking_results_pre_filtering = docking_results.copy()

        # choose n best fragments
        docking_results = docking_results[:min(len(docking_results), num_fragments_per_iterations)]

        # choose k best conformers (per fragment)
        for ligand in docking_results:
            if cluster_based_pose_filtering:
              ligand.choose_template_poses_cluster_based(num_conformers, clusted_pose_filter_dist_threshold)
            else:
              ligand.choose_template_poses(num_conformers)

        # log chosen poses and molecules
        output_logs['ChosenMolecules']['SP' + str(subpockets.index(subpocket))] = len(docking_results)
        output_logs['ChosenPoses']['SP' + str(subpockets.index(subpocket))] = sum(len(mol.poses) for mol in docking_results)

        # store the best pose of a docked molecule if it consists of more than 1 fragment in the overall result file
        docking_utils.append_ligands_to_file(docking_results_pre_filtering, docking_results, PATH_TO_RESULTS/ 'results.sdf', lambda l: len(l.fragment_ids) > 1)

        # store all poses of all docked molecules in one file for the current subpocket iteration
        docking_utils.write_all_poses_to_file(docking_results_pre_filtering, docking_results, PATH_TO_RESULTS/ ('SP' + str(subpockets.index(subpocket)) + '.sdf'))


        # recombining

        candidates = []

        num_recombinations = 0

        counter_num_multpl_bonds = 0
        counter_num_unambigious_bonds = 0
        # try recombine every ligand (comb. of fragmnets) with every fragment of the current subpocket
        for ligand in docking_results:
            for fragment_idx in fragment_library[subpocket].index:
                num_multpl_bonds, num_unambigious_bonds = ligand.recombine(fragment_idx, subpocket, fragment_library) # possible recombinations are stored within ligand
                counter_num_multpl_bonds += num_multpl_bonds
                counter_num_unambigious_bonds += num_unambigious_bonds
            if len(ligand.recombinations):
                num_recombinations += len(ligand.recombinations)
                # use a ligand for template docking iff at least one recombination was produced for a ligand
                candidates.append(ligand)

        logging.debug('Generated ' + str(num_recombinations) + ' recombinations')
        output_logs['Recombinations']['SP' + str(subpockets.index(subpocket) + 1)] = num_recombinations
        output_logs['NumMultipleBondsBetweenFragmentAndLigand']['SP' + str(subpockets.index(subpocket) + 1)] = counter_num_multpl_bonds
        output_logs['NumMultipleBondsBetweenTwoFragments']['SP' + str(subpockets.index(subpocket) + 1)] = counter_num_unambigious_bonds
        # template docking

        docking_results = []
        violations = []

        num_docking_runs = sum(len(c.poses) * len(c.recombinations) for c in candidates)
        num_succ_docking_runs = 0
        
        start_time = time.time()
        
        wandb.log({'total_frags_SP' + str(subpockets.index(subpocket) + 1): num_recombinations})
        docking_run_counter = 0

        with ThreadPoolExecutor(num_threads) as executor:
            # create partial template docking function due to huge amount of arguments
            task = partial(threading_docking.template_docking_task, PATH_TEMP,  PATH_TO_CONFIGS, PATH_TO_CONFIGS, PATH_FLEXX,PATH_HYDE, cutoff_hyde_displacement)
            # submit template docking tasks
            features = [executor.submit(task, subpocket, recombination, ligand.poses) for ligand in candidates for recombination in ligand.recombinations]
            # iterate over all completetd tasks
            for count, feature in enumerate(as_completed(features)):
                try:
                    # get docking result
                    result = feature.result()
                except Exception as exc:
                    logging.error('Generated an exception during template_docking: %s' % (exc)) 
                else:
                    docking_run_counter += 1
                    wandb.log({"docked_frags_SP" + str(subpockets.index(subpocket) + 1): docking_run_counter})
                    if len(result) == 3:
                        # docking was succesfull
                        docking_results.append(result[0])
                        violations += result[1]
                        num_succ_docking_runs += result[2]
                    elif len(result) == 1:
                        # could not genrate the 3d conformation
                        unsuccesfull_3d_generations += result
                        num_docking_runs -= len(result[0].poses) * len(result[0].recombinations)
  
        logging.info(f"Runtime template-docking ({subpocket}): {(time.time() - start_time)}")

        # template docking logs
        output_logs['RunTimeTemplateDocking']['SP' + str(subpockets.index(subpocket) + 1)] = time.strftime('%H:%M:%S', time.gmtime(time.time() - start_time))
        output_logs['GeneratedPoses']['SP' + str(subpockets.index(subpocket) + 1)] = sum(len(mol.poses) for mol in docking_results)
        output_logs['UndockableMolecules']['SP' + str(subpockets.index(subpocket) + 1)] = num_recombinations - len(docking_results)
        output_logs['DockedMolecules']['SP' + str(subpockets.index(subpocket) + 1)] = len(docking_results)
        output_logs['MeanHydeDisplacement']['SP' + str(subpockets.index(subpocket) + 1)] = sum(mol.mean_hyde_displacement_undropped for mol in docking_results) / (len(docking_results) or 1)
        output_logs['MeanHydeDisplacementIncludingViolations']['SP' + str(subpockets.index(subpocket) + 1)] = sum(mol.mean_hyde_displacement for mol in docking_results) / (len(docking_results) or 1) 
        output_logs['NumDisplacementViolations']['SP' + str(subpockets.index(subpocket) + 1)] = sum(mol.num_hyde_violations for mol in docking_results)
        output_logs['DockingRuns']['SP' + str(subpockets.index(subpocket) + 1)] = num_docking_runs
        output_logs['SuccesfullDockingRuns']['SP' + str(subpockets.index(subpocket) + 1)] = num_succ_docking_runs
        # write violations so sdf file
        docking_utils.write_violations_to_file(violations, PATH_TO_RESULTS / ('violations_SP' + str(subpockets.index(subpocket) + 1) + '.sdf'))

    # store all poses of all docked molecules in one file for the current subpocket iteration
    docking_utils.write_all_poses_to_file(docking_results, docking_results, PATH_TO_RESULTS/ ('SP' + str(len(subpockets)) + '.sdf'))

    output_logs['Unsuccesful3DGenerations'] = [Chem.MolToSmiles(mol.ROMol) for mol in unsuccesfull_3d_generations]

    # log total run time
    output_logs['RunTimeTotal'] = time.strftime('%H:%M:%S', time.gmtime(round(time.time() - start_time_all, 2)))
    logging.info("Runtime: %s" % (time.time() - start_time_all))

    # ===== EVALUATION ======

    # write logs to json output file
    with open(args.output, 'w') as json_file:
        json.dump(output_logs, json_file, indent=4)

    if len(docking_results):
        if use_hyde:
            docking_results.sort(key=lambda l: l.min_binding_affinity)
        else:
            docking_results.sort(key=lambda l: l.min_docking_score)

        logging.debug("Best recombination: " + str(list(docking_results[0].fragment_ids.items())) + " Score: " + str(docking_results[0].min_binding_affinity or docking_results[0].min_docking_score))

        # store generated ligands
        docking_utils.append_ligands_to_file(docking_results, docking_results, PATH_TO_RESULTS/ 'results.sdf', lambda l: len(l.fragment_ids) > 1)

        

