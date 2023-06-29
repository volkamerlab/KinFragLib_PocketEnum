from pathlib import Path
from kinfraglib import utils, filters
import docking_utils
from rdkit import Chem
from pathlib import Path
import logging
from queue import Queue
import time

import json
import threading_docking
import itertools

if __name__ == '__main__':

    with open('definitions.json', 'rt') as json_file:
        definitions = json.load(json_file)

    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')

    # Definitions
    HERE = Path().resolve()
    PATH_DATA = Path(definitions['KinFragLib'])
    PATH_FLEXX = Path(definitions['FlexX'])
    PATH_TO_DOCKING_CONFIGS = Path(definitions['DockingConfig']) / definitions['pdbCode']
    PATH_TO_SDF_FRAGMENTS = HERE / 'data/fragments' / definitions['pdbCode']
    PATH_TO_DOCKING_RESULTS = HERE / 'data/docking' / definitions['pdbCode']
    PATH_TO_HYDE_RESULTS = HERE / 'data/scoring' / definitions['pdbCode']
    PATH_TO_TEMPLATES =  HERE / 'data/templates' / definitions['pdbCode']

    num_fragments = 4                                                           # number of fragments to use 
    num_conformers = definitions['NumberPosesPerFragment']                      # amount of conformers to choose per docked fragment  (according to docking score and diversity)
    num_fragments_per_iterations = definitions['NumberFragmentsPerIterations']  # amount of fragments to choose per docking iteration (according to docking score)

    # define filters
    filters_ = [docking_utils.Filter(name, values) for name, values in definitions['Filters'].items()]

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
    for filter in filters_:
        logging.debug('Applying filter: ' + filter.name)
        l_before = [sp + ': ' + str(len(fragment_library[sp])) for sp in fragment_library.keys()]
        fragment_library = filter.apply_filter(fragment_library)
        logging.debug(filter.name + ' applied: ' + str([sp + ': ' + str(len(fragment_library[sp])) for sp in fragment_library.keys()]))

    logging.info('Preprocessing finished\n Size of fragment library' + str([sp + ': ' + str(len(fragment_library[sp])) for sp in fragment_library.keys()]))

    # ===== CORE DOCKING =======

    logging.info('Core docking of ' + str(len(fragment_library[core_subpocket])) + ' ' + core_subpocket + '-Fragments')

    core_fragments = []

    # prepare all core fragments
    for i in fragment_library[core_subpocket][:num_fragments].index:
        core_fragments.append(docking_utils.Ligand(fragment_library[core_subpocket]['ROMol'][i], {core_subpocket: i}, docking_utils.Recombination([core_subpocket + "_" + str(i)], [])))

    # core docking 
    # define result table per job
    docking_results = []

    queue = Queue()

    for x in range(8):
        worker = threading_docking.DockingWorker(queue, x, PATH_TO_SDF_FRAGMENTS, PATH_TO_DOCKING_CONFIGS, PATH_TO_DOCKING_RESULTS, PATH_TO_TEMPLATES ,PATH_FLEXX)
        worker.daemon = True
        worker.start()
    
    start_time = time.time()

    for core_fragment in core_fragments:
        _results = []
        docking_results.append(_results)
        queue.put(('core_docking', _results, core_fragment, core_subpocket))
    
    queue.join()

    logging.info("Runtime: %s" % (time.time() - start_time))
    docking_results = list(itertools.chain.from_iterable(docking_results))
    
    
    

    # ===== TEMPLATE DOCKING ======

    for subpocket in subpockets:
        logging.info('Template docking of ' + str(len(fragment_library[subpocket])) + ' ' + subpocket + '-Fragments')

        # filtering

        docking_results.sort(key=lambda l: l.min_docking_score)

        # choose n best fragments
        docking_results = docking_results[:min(len(docking_results), num_fragments_per_iterations)]

        # choose k best conformers (per fragment)
        for ligand in docking_results:
            ligand.choose_template_poses(num_conformers)

        # only for evaluation (statistics) purpose 
        with open("statistics.txt", "a") as f:
            for ligand in docking_results:
                f.write(str(list(ligand.fragment_ids.items())) + ": " + str(ligand.min_docking_score) + "\n")

        # recombining

        candidates = []

        num_recombinations = 0

        # try recombine every ligand (comb. of fragmnets) with every fragment of the current subpocket
        for ligand in docking_results:
            for fragment_idx in fragment_library[subpocket][:num_fragments].index:
                ligand.recombine(fragment_idx, subpocket, fragment_library) # possible recombinations are stored within ligand
            if len(ligand.recombinations):
                num_recombinations += len(ligand.recombinations)
                # use a ligand for template docking only if at least one recombination was produced for a ligand
                candidates.append(ligand)

        logging.debug('Generated ' + str(num_recombinations) + ' recombinations')
        # template docking

        docking_results = []
        
        start_time = time.time()

        ligand : docking_utils.Ligand
        for ligand in candidates:
            for recombination in ligand.recombinations:
                _results = []
                docking_results.append(_results)
                queue.put(('template_docking', _results, (ligand.poses, recombination), subpocket))
        
        queue.join()   

        logging.info("Runtime: %s" % (time.time() - start_time))

        docking_results = list(itertools.chain.from_iterable(docking_results))
                

    # ===== EVALUATION ======
    if len(docking_results):
        docking_results.sort(key=lambda l: l.min_docking_score)
        with open("statistics.txt", "a") as f:
            for ligand in docking_results:
                f.write(str(list(ligand.fragment_ids.items())) + ": " + str(ligand.min_docking_score) + "\n")
        logging.debug("Best recombination: " + str(list(docking_results[0].fragment_ids.items())) + " Score: " + str(docking_results[0].min_docking_score))

        min_pose = min(docking_results[0].poses, key=lambda p: p.docking_score)
        with Chem.SDWriter(str(PATH_TO_DOCKING_RESULTS / ('final_fragment.sdf'))) as w:
            w.write(min_pose.ROMol)