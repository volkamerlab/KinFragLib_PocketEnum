from pathlib import Path
from kinfraglib import utils
import docking_utils
from rdkit import Chem
from pandas import DataFrame
import pandas as pd
from brics_rules import is_brics_bond
from itertools import permutations

# Definitions
PATH_TO_LIB = Path('KinFragLib/data/fragment_library')
PATH_COMBINATORIAL_LIBRARY = 'KinFragLib/data/combinatorial_library/combinatorial_library_deduplicated.json'
PATH_TO_DOCKING_CONFIGS = Path('docking_config/5l4q')
PATH_TO_SDF_FRAGMENTS = Path('data/fragments/5l4q/')
PATH_TO_DOCKING_RESULTS = Path('data/docking/5l4q/')
PATH_TO_TEMPLATES =  Path('data/templates/5l4q/')

num_ap_fragments = 10
num_conformers = 5

core_subpocket = 'AP'

subpockets = ['FP', 'SE']

# Access fragment library
fragment_library = utils.read_fragment_library(PATH_TO_LIB)

core_fragments = []

# pre selection  has to be applied here
for i in fragment_library[core_subpocket][:20].index:
    core_fragments.append(docking_utils.Ligand(fragment_library[core_subpocket]['ROMol'][i], {core_subpocket: i}, docking_utils.Recombination([core_subpocket + "_" + str(i)], [])))

docking_results = []

# core docking
for core_fragment in core_fragments:
    # safe fragment as sdf file
    print("=== Docking of " + core_subpocket + "-Fragment: " + str(core_fragment.fragment_ids[core_subpocket]) + " ====")
    core_fragment.to_sdf(PATH_TO_SDF_FRAGMENTS / 'core_fragment.sdf')
    res = docking_utils.core_docking(PATH_TO_SDF_FRAGMENTS / 'core_fragment.sdf', PATH_TO_DOCKING_CONFIGS / (core_subpocket + '.flexx'), PATH_TO_DOCKING_RESULTS / 'core_fragments.sdf')
    for conformer in res:
        pose = docking_utils.Pose(conformer, float(conformer.GetProp('BIOSOLVEIT.DOCKING_SCORE')))
        core_fragment.add_pose(pose)
    if len(res):
        print("  best docking score: " + str(core_fragment.min_docking_score))
        docking_results.append(core_fragment)

# Evaluation
# HYDE TODO


used_subpockets = [core_subpocket]
subpocket = subpockets[0]

for subpocket in subpockets:

    for ligand in docking_results:
        ligand.choose_template_poses(num_conformers)

    docking_results.sort(key=lambda l: l.min_docking_score)

    with open("statistics.txt", "a") as f:
        for ligand in docking_results:
            f.write(str(list(ligand.fragment_ids.items())) + ": " + str(ligand.min_docking_score) + "\n")
       

    candidates = []
    for ligand in docking_results:
        for fragment_idx in fragment_library[subpocket][:20].index:
            ligand.recombine(fragment_idx, subpocket, fragment_library)
        if len(ligand.recombinations):
            candidates.append(ligand)

    docking_results = []

    for ligand in candidates:
        for recombination in ligand.recombinations:
            fragment = docking_utils.from_recombination(recombination)
            fragment.to_sdf(PATH_TO_SDF_FRAGMENTS /(subpocket + '_fragment.sdf'))
            for i, pose in enumerate(ligand.poses):
                w = Chem.SDWriter(str(PATH_TO_TEMPLATES / (subpocket + '_fragment.sdf')))  # TODO: sometimes it doesn't write anything to the file
                w.write(pose.ROMol)

                print("=== Docking of " + str(list(fragment.fragment_ids.items())) + " Pose: " + str(i) +  " ====")
                res = docking_utils.template_docking(PATH_TO_SDF_FRAGMENTS / (subpocket + '_fragment.sdf'), PATH_TO_TEMPLATES / (subpocket + '_fragment.sdf'), PATH_TO_DOCKING_CONFIGS / (subpocket + '.flexx'), PATH_TO_DOCKING_RESULTS / ('fragments.sdf'))
                for conformer in res:
                    pose = docking_utils.Pose(conformer, float(conformer.GetProp('BIOSOLVEIT.DOCKING_SCORE')))
                    fragment.add_pose(pose)
                if len(res):
                    print("  best docking score: " + str(fragment.min_docking_score))
            if len(fragment.poses):
                docking_results.append(fragment)

if len(docking_results):
    docking_results.sort(key=lambda l: l.min_docking_score)
    with open("statistics.txt", "a") as f:
        for ligand in docking_results:
            f.write(str(list(ligand.fragment_ids.items())) + ": " + str(ligand.min_docking_score) + "\n")
    print("Best recombination: " + str(list(docking_results[0].fragment_ids.items())) + " Score: " + str(docking_results[0].min_docking_score))

    min_pose = min(docking_results[0].poses)
    w = Chem.SDWriter(str(PATH_TO_DOCKING_RESULTS / ('final_fragment.sdf')))  # TODO: sometimes it doesn't write anything to the file
    w.write(min_pose.ROMol)