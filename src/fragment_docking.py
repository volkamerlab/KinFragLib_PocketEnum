from pathlib import Path
from kinfraglib import utils, filters
import docking_utils
from rdkit import Chem
from pathlib import Path
import logging

import json

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
PATH_TO_RESULTS = HERE / 'data/results' / definitions['pdbCode']

num_fragments = 2                                                           # number of fragments to use 
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
docking_results = []

for core_fragment in core_fragments:
    logging.debug('Docking of ' + core_subpocket + "-Fragment: " + str(core_fragment.fragment_ids[core_subpocket]))
    # safe fragment as sdf file
    if not core_fragment.to_sdf(PATH_TO_SDF_FRAGMENTS / 'core_fragment.sdf'):
        # could not generate 3D-conformation
        logging.error('Could not write Fragemnt: ' + str(core_fragment.fragment_ids) + ' to files due to 3d-generation-error')
        continue

    res = docking_utils.core_docking(PATH_TO_SDF_FRAGMENTS / 'core_fragment.sdf', PATH_TO_DOCKING_CONFIGS / (core_subpocket + '.flexx'), PATH_TO_DOCKING_RESULTS / 'core_fragments.sdf', PATH_FLEXX)
    docking_utils.remove_files(PATH_TO_DOCKING_RESULTS / 'core_fragments.sdf', PATH_TO_SDF_FRAGMENTS / 'core_fragment.sdf')

    for conformer in res:   # safe every resulting pose within the fragment
        pose = docking_utils.Pose(conformer, float(conformer.GetProp('BIOSOLVEIT.DOCKING_SCORE')))
        core_fragment.add_pose(pose)
    logging.debug(str(len(res)) + ' poses have been generated')
    if len(res):
        # if fragment could be docked, save fragment (including it's poses)
        logging.debug("Best docking score: " + str(core_fragment.min_docking_score))
        docking_results.append(core_fragment)

# ===== TEMPLATE DOCKING ======

for subpocket in subpockets:
    logging.info('Template docking of ' + str(len(fragment_library[subpocket])) + ' ' + subpocket + '-Fragments')

    # filtering

    docking_results.sort(key=lambda l: l.min_docking_score)

    # store intermediate results if docking result is negativ and consists of more than 1 fragment
    docking_utils.append_ligands_to_file(docking_results, PATH_TO_RESULTS/ 'results.sdf', lambda l: l.min_docking_score <= 0 and len(l.fragment_ids) > 1)

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

    ligand : docking_utils.Ligand
    for ligand in candidates:
        for recombination in ligand.recombinations:
            logging.debug('Template docking of Recombination: ' + str(recombination.fragments))
            fragment = docking_utils.from_recombination(recombination)

            # write recombination that should be docked to file
            if not fragment.to_sdf(PATH_TO_SDF_FRAGMENTS /(subpocket + '_fragment.sdf')):
                logging.error('Could not write fragemnt: ' + str(fragment.fragment_ids) + ' to files due to 3d-generation-error')
                continue

            # for every choosen pose: perform template docking
            for i, pose in enumerate(ligand.poses):
                logging.debug('Template ' + str(i + 1) + ' / ' + str(len(ligand.poses)))
                # write template to sdf file
                with Chem.SDWriter(str(PATH_TO_TEMPLATES / (subpocket + '_fragment.sdf'))) as w:
                    w.write(pose.ROMol)
                # template docking (FlexX)
                res = docking_utils.template_docking(PATH_TO_SDF_FRAGMENTS / (subpocket + '_fragment.sdf'), PATH_TO_TEMPLATES / (subpocket + '_fragment.sdf'), PATH_TO_DOCKING_CONFIGS / (subpocket + '.flexx'), PATH_TO_DOCKING_RESULTS / ('fragments.sdf'), PATH_FLEXX)

                # remove files containg docking results and template
                docking_utils.remove_files(PATH_TO_DOCKING_RESULTS / ('fragments.sdf'))

                # safe resulting poses within fragment
                logging.debug(str(len(res)) + ' poses have been generated')
                for conformer in res:
                    pose = docking_utils.Pose(conformer, float(conformer.GetProp('BIOSOLVEIT.DOCKING_SCORE')))
                    fragment.add_pose(pose)
                if len(res):
                    logging.debug("Best docking score: " + str(fragment.min_docking_score))
            if len(fragment.poses):
                # safe recombination as result only if at least one pose was generated
                docking_results.append(fragment)

# ===== EVALUATION ======
if len(docking_results):
    docking_results.sort(key=lambda l: l.min_docking_score)

    docking_utils.append_ligands_to_file(docking_results, PATH_TO_RESULTS/ 'results.sdf')

    with open("statistics.txt", "a") as f:
        for ligand in docking_results:
            f.write(str(list(ligand.fragment_ids.items())) + ": " + str(ligand.min_docking_score) + "\n")
    logging.debug("Best recombination: " + str(list(docking_results[0].fragment_ids.items())) + " Score: " + str(docking_results[0].min_docking_score))

    min_pose = min(docking_results[0].poses, key=lambda p: p.docking_score)
    with Chem.SDWriter(str(PATH_TO_DOCKING_RESULTS / ('final_fragment.sdf'))) as w:
        w.write(min_pose.ROMol)