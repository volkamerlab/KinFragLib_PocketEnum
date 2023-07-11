from pathlib import Path
from kinfraglib import utils, filters
import docking_utils
from rdkit import Chem
from pathlib import Path
import logging
import time
import os


import json
import threading_docking
from functools import partial
from concurrent.futures import ThreadPoolExecutor, as_completed

if __name__ == '__main__':

    with open('definitions.json', 'rt') as json_file:
        definitions = json.load(json_file)

    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')

    # Definitions
    HERE = Path().resolve()
    PATH_DATA = Path(definitions['KinFragLib'])
    PATH_FLEXX = Path(definitions['FlexX'])
    PATH_TO_DOCKING_CONFIGS = Path(definitions['DockingConfig']) / definitions['pdbCode']
    PATH_TO_HYDE_CONFIGS = HERE / 'hyde_config' / definitions['pdbCode']
    PATH_TO_SDF_FRAGMENTS = HERE / 'data/fragments' / definitions['pdbCode']
    PATH_TO_DOCKING_RESULTS = HERE / 'data/docking' / definitions['pdbCode']
    PATH_TO_HYDE_RESULTS = HERE / 'data/scoring' / definitions['pdbCode']
    PATH_TO_TEMPLATES =  HERE / 'data/templates' / definitions['pdbCode']
    PATH_TO_RESULTS = HERE / 'data/results' / definitions['pdbCode']
    PATH_HYDE = Path(definitions['Hyde']) if definitions['UseHyde'] else None
    use_hyde = definitions['UseHyde']
    cutoff_hyde_displacement = definitions.get('HydeDisplacementCutoff') or 1.5     # rmds cutoff when a optimized pose can be droped 

    # clear results file if it exists
    if os.path.exists(PATH_TO_RESULTS/ 'results.sdf'):
        with open(PATH_TO_RESULTS/ 'results.sdf', 'wt') as file:
            pass

    num_fragments = 1                                                          # number of fragments to use 
    num_conformers = definitions['NumberPosesPerFragment']                      # amount of conformers to choose per docked fragment  (according to docking score and diversity)
    num_fragments_per_iterations = definitions['NumberFragmentsPerIterations']  # amount of fragments to choose per docking iteration (according to docking score)

    num_threads = definitions.get('NumberThreads')                              # number of threads to use, if not specified min(32, os.cpu_count() + 4) are used. If multithreading isn't wanted, NumberThreads should be set to 1
    cluster_based_pose_filtering = definitions['UseClusterBasedPoseFiltering']  
    clusted_pose_filter_dist_threshold = definitions.get('DistanceThresholdClustering') # threshold that shuld be used for pose clustering

    # define filters
    filters_ = [docking_utils.Filter(name, values) for name, values in definitions['Filters'].items()]

    start_time_all = time.time()

    # define pockets
    core_subpocket = definitions['CoreSubpocket']

    subpockets = definitions['Subpockets']

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
    for filter in []:#filters_:
        logging.debug('Applying filter: ' + filter.name)
        l_before = [sp + ': ' + str(len(fragment_library[sp])) for sp in fragment_library.keys()]
        fragment_library = filter.apply_filter(fragment_library)
        logging.debug(filter.name + ' applied: ' + str([sp + ': ' + str(len(fragment_library[sp])) for sp in fragment_library.keys()]))

    logging.info('Preprocessing finished\n Size of fragment library' + str([sp + ': ' + str(len(fragment_library[sp])) for sp in fragment_library.keys()]))

    # ===== CORE DOCKING =======

    logging.info('Core docking of ' + str(len(fragment_library[core_subpocket])) + ' ' + core_subpocket + '-Fragments')

    core_fragments = []

    # prepare all core fragments
    for i in fragment_library[core_subpocket].index:
        core_fragments.append(docking_utils.Ligand(fragment_library[core_subpocket]['ROMol'][i], {core_subpocket: i}, docking_utils.Recombination([core_subpocket + "_" + str(i)], [])))

    # core docking 
    docking_results = []
    
    start_time = time.time()

    # create partial docking task to avoid function calls with many of 
    core_docking_task = partial(threading_docking.core_docking_task, PATH_TO_SDF_FRAGMENTS,  PATH_TO_DOCKING_CONFIGS, 
                                PATH_TO_DOCKING_RESULTS, PATH_TO_HYDE_RESULTS, PATH_TO_HYDE_CONFIGS, PATH_FLEXX, PATH_HYDE, cutoff_hyde_displacement)

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
                docking_results += result

    logging.info("Runtime: %s" % (time.time() - start_time))

    # ===== TEMPLATE DOCKING ======

    for subpocket in subpockets:
        logging.info('Template docking of ' + str(len(fragment_library[subpocket])) + ' ' + subpocket + '-Fragments')

        # filtering
        if use_hyde:
            docking_results.sort(key=lambda l: l.min_binding_affinity)
        else:
            docking_results.sort(key=lambda l: l.min_docking_score)

        # store intermediate results if docking result is negativ and consists of more than 1 fragment
        docking_utils.append_ligands_to_file(docking_results, PATH_TO_RESULTS/ 'results.sdf', lambda l: len(l.fragment_ids) > 1)

        # choose n best fragments
        docking_results = docking_results[:min(len(docking_results), num_fragments_per_iterations)]

        # choose k best conformers (per fragment)
        for ligand in docking_results:
            if cluster_based_pose_filtering:
              ligand.choose_template_poses_cluster_based(num_conformers, clusted_pose_filter_dist_threshold)
            else:
              ligand.choose_template_poses(num_conformers)

        # only for evaluation (statistics) purpose 
        with open("statistics.txt", "a") as f:
            for ligand in docking_results:
                f.write(str(list(ligand.fragment_ids.items())) + ": " + str(ligand.min_binding_affinity or ligand.min_docking_score) + "\n")

        # recombining

        candidates = []

        num_recombinations = 0

        # try recombine every ligand (comb. of fragmnets) with every fragment of the current subpocket
        for ligand in docking_results:
            for fragment_idx in fragment_library[subpocket].index:
                ligand.recombine(fragment_idx, subpocket, fragment_library) # possible recombinations are stored within ligand
            if len(ligand.recombinations):
                num_recombinations += len(ligand.recombinations)
                # use a ligand for template docking iff at least one recombination was produced for a ligand
                candidates.append(ligand)

        logging.debug('Generated ' + str(num_recombinations) + ' recombinations')
        # template docking

        docking_results = []
        
        start_time = time.time()

        with ThreadPoolExecutor(num_threads) as executor:
            # create partial template docking function due to huge amount of arguments
            task = partial(threading_docking.template_docking_task, PATH_TO_SDF_FRAGMENTS,  PATH_TO_DOCKING_CONFIGS, PATH_TO_DOCKING_RESULTS, 
                           PATH_TO_HYDE_RESULTS, PATH_TO_HYDE_CONFIGS, PATH_FLEXX,PATH_HYDE, PATH_TO_TEMPLATES, cutoff_hyde_displacement)
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
                    docking_results += result
  
        logging.info(f"Runtime template-docking ({subpocket}): {(time.time() - start_time)}")

    logging.info("Runtime: %s" % (time.time() - start_time_all))

    # ===== EVALUATION ======
    if len(docking_results):
        if use_hyde:
            docking_results.sort(key=lambda l: l.min_binding_affinity)
        else:
            docking_results.sort(key=lambda l: l.min_docking_score)

        with open("statistics.txt", "a") as f:
            for ligand in docking_results:
                f.write(str(list(ligand.fragment_ids.items())) + ": " + str(ligand.min_docking_score) + "\n")

        logging.debug("Best recombination: " + str(list(docking_results[0].fragment_ids.items())) + " Score: " + str(docking_results[0].min_binding_affinity or docking_results[0].min_docking_score))

        # store intermediate results if docking score <= 50 and consists of more than 1 fragment
        docking_utils.append_ligands_to_file(docking_results, PATH_TO_RESULTS/ 'results.sdf', lambda l: len(l.fragment_ids) > 1)

        min_pose = min(docking_results[0].poses, key=lambda p: p.binding_affinity_lower or p.docking_score)
        with Chem.SDWriter(str(PATH_TO_DOCKING_RESULTS / ('final_fragment.sdf'))) as w:
            w.write(min_pose.ROMol)
