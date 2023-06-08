from pathlib import Path
from kinfraglib import utils, filters
import docking_utils
from rdkit import Chem
from pathlib import Path

import pandas as pd



# Definitions
HERE = Path().resolve()
PATH_DATA = HERE / "KinFragLib/data"
PATH_DATA_BRENK = HERE / "KinFragLib/data/filters/Brenk"
PATH_DATA_ENAMINE = HERE / "KinFragLib/data/filters/DataWarrior"
PATH_COMBINATORIAL_LIBRARY = HERE /'KinFragLib/data/combinatorial_library/combinatorial_library_deduplicated.json'
PATH_TO_DOCKING_CONFIGS = HERE / 'docking_config/5l4q'
PATH_TO_HYDE_CONFIGS = HERE / 'hyde_config/5l4q'
PATH_TO_SDF_FRAGMENTS = HERE / 'data/fragments/5l4q'
PATH_TO_DOCKING_RESULTS = HERE / 'data/docking/5l4q'
PATH_TO_HYDE_RESULTS = HERE / 'data/scoring/5l4q'
PATH_TO_TEMPLATES =  HERE / 'data/templates/5l4q'

num_fragments = 5                   # number of fragments to use 
num_conformers = 5                  # amount of conformers to choose per docked fragment  (according to docking score and diversity)
num_fragments_per_iterations = 5  # amount of fragments to choose per docking iteration (according to docking score)

# define filters
filters_ = []

"""filters_.append(docking_utils.Filter('pains', {}))
filters_.append(docking_utils.Filter('brenk', {'path_data': PATH_DATA_BRENK}))
filters_.append(docking_utils.Filter('ro3', {})) # TODO add params
filters_.append(docking_utils.Filter('qed', {'cutoff_val': 0.492}))
filters_.append(docking_utils.Filter('bb', {'path_data': PATH_DATA_ENAMINE}))
filters_.append(docking_utils.Filter('syba', {'cutoff_val': 0}))"""

# define pockets
core_subpocket = 'AP'

subpockets = ['FP', 'SE']

# ==== PREPROCESSING =======

# Access fragment library
fragment_library_original = utils.read_fragment_library(PATH_DATA / "fragment_library")

print(str([sp + ': ' + str(len(fragment_library_original[sp])) for sp in fragment_library_original.keys()]))

print('Prefiltering ===>')
# removing fragments in pool X
# removing duplicates
# removing fragments without dummy atoms (unfragmented ligands)
# removing fragments only connecting to pool X
fragment_library = filters.prefilters.pre_filters(fragment_library_original)
print(str([sp + ': ' + str(len(fragment_library[sp])) for sp in fragment_library.keys()]))

# apply filters
for filter in filters_:
    print('Filtering: ' + filter.name + ' ===>')
    l_before = [sp + ': ' + str(len(fragment_library[sp])) for sp in fragment_library.keys()]
    fragment_library = filter.apply_filter(fragment_library)
    print(str([sp + ': ' + str(len(fragment_library[sp])) for sp in fragment_library.keys()]))

# ===== CORE DOCKING =======

core_fragments = []

# prepare all core fragments
for i in fragment_library[core_subpocket][:num_fragments].index:
    core_fragments.append(docking_utils.Ligand(fragment_library[core_subpocket]['ROMol'][i], {core_subpocket: i}, docking_utils.Recombination([core_subpocket + "_" + str(i)], [])))

# core docking 
docking_results = []

for core_fragment in core_fragments:
    print("=== Docking of " + core_subpocket + "-Fragment: " + str(core_fragment.fragment_ids[core_subpocket]) + " ====")
    # safe fragment as sdf file
    core_fragment.to_sdf(PATH_TO_SDF_FRAGMENTS / 'core_fragment.sdf')

    res = docking_utils.core_docking(PATH_TO_SDF_FRAGMENTS / 'core_fragment.sdf', PATH_TO_DOCKING_CONFIGS / (core_subpocket + '.flexx'), PATH_TO_DOCKING_RESULTS / 'core_fragments.sdf')
    docking_utils.remove_files(PATH_TO_DOCKING_RESULTS / 'core_fragments.sdf', PATH_TO_SDF_FRAGMENTS / 'core_fragment.sdf')

    for conformer in res:   # safe every resulting pose within the fragment
        pose = docking_utils.Pose(conformer, float(conformer.GetProp('BIOSOLVEIT.DOCKING_SCORE')))
        core_fragment.add_pose(pose)
    if len(res):
        # if fragment could be docked, save fragment (including it's poses)
        print("  best docking score: " + str(core_fragment.min_docking_score))
        docking_results.append(core_fragment)

# ===== TEMPLATE DOCKING ======

for subpocket in subpockets:

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

    # try recombine every ligand (comb. of fragmnets) with every fragment of the current subpocket
    for ligand in docking_results:
        for fragment_idx in fragment_library[subpocket][:num_fragments].index:
            ligand.recombine(fragment_idx, subpocket, fragment_library) # possible recombinations are stored within ligand
        if len(ligand.recombinations):
            # use a ligand for template docking only if at least one recombination was produced for a ligand
            candidates.append(ligand)

    # template docking

    docking_results = []

    ligand : docking_utils.Ligand
    for ligand in candidates:
        for recombination in ligand.recombinations:
            fragment = docking_utils.from_recombination(recombination)

            # write recombination that should be docked to file
            fragment.to_sdf(PATH_TO_SDF_FRAGMENTS /(subpocket + '_fragment.sdf'))

            # for every choosen pose: perform template docking
            for i, pose in enumerate(ligand.poses):
                # write template to sdf file
                with Chem.SDWriter(str(PATH_TO_TEMPLATES / (subpocket + '_fragment.sdf'))) as w:
                    w.write(pose.ROMol)
                # template docking (FlexX)
                print("=== Docking of " + str(list(fragment.fragment_ids.items())) + " Pose: " + str(i) +  " ====")
                res = docking_utils.template_docking(PATH_TO_SDF_FRAGMENTS / (subpocket + '_fragment.sdf'), PATH_TO_TEMPLATES / (subpocket + '_fragment.sdf'), PATH_TO_DOCKING_CONFIGS / (subpocket + '.flexx'), PATH_TO_DOCKING_RESULTS / ('fragments.sdf'))
                
                # remove files containg docking results and template
                docking_utils.remove_files(PATH_TO_DOCKING_RESULTS / ('fragments.sdf'))

                # safe resulting poses within fragment
                for conformer in res:
                    pose = docking_utils.Pose(conformer, float(conformer.GetProp('BIOSOLVEIT.DOCKING_SCORE')))
                    fragment.add_pose(pose)
                if len(res):
                    print("  best docking score: " + str(fragment.min_docking_score))
            if len(fragment.poses):
                # safe recombination as result only if at least one pose was generated
                docking_results.append(fragment)

# ===== EVALUATION ======
if len(docking_results):
    docking_results.sort(key=lambda l: l.min_docking_score)
    with open("statistics.txt", "a") as f:
        for ligand in docking_results:
            f.write(str(list(ligand.fragment_ids.items())) + ": " + str(ligand.min_docking_score) + "\n")
    print("Best recombination: " + str(list(docking_results[0].fragment_ids.items())) + " Score: " + str(docking_results[0].min_docking_score))

    min_pose = min(docking_results[0].poses, key=lambda p: p.docking_score)
    with Chem.SDWriter(str(PATH_TO_DOCKING_RESULTS / ('final_fragment.sdf'))) as w: 
        w.write(min_pose.ROMol)