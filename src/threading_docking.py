import threading
import docking_utils 
import logging
from rdkit import Chem

def core_docking_task(PATH_TO_SDF_FRAGMENTS,  PATH_TO_DOCKING_CONFIGS, PATH_TO_DOCKING_RESULTS, PATH_FLEXX, core_subpocket: str, core_fragment):
    """
    runs a FlexX core docking task. It should be used for multithreading.

    Returns
    ----------
    List of all poses (docking result) as Ligand (Properties: docking-score, pose)

    Parameters
    ----------
    PATH_TO_SDF_FRAGMENTS: pathlib.path
        Path to directory where sdf files of the core fragments are stored intermediately
    path_config: pathlib.path
        Path to FlexX-config file.
    path_output: pathlib.path
        Path to diectory of output files
    path_flexx: pathlib.path
        Path to FlexX
    core_subpocket: str
        Subpocket where to place the core fragment
    core_fragment: Ligand
        Fragment that should be docked in this task
    
    """
    thread_id = threading.get_ident()
    logging.debug('Docking of ' + core_subpocket + "-Fragment: " + str(core_fragment.fragment_ids[core_subpocket]))
    # safe fragment as sdf file
    if not core_fragment.to_sdf(PATH_TO_SDF_FRAGMENTS / ('thread_' + str(thread_id) + '_core_fragment.sdf')):
        # could not generate 3D-conformation
        logging.error('Could not write Fragemnt: ' + str(core_fragment.fragment_ids) + ' to files due to 3d-generation-error')
        return []

    res = docking_utils.core_docking(PATH_TO_SDF_FRAGMENTS / ('thread_' + str(thread_id) + '_core_fragment.sdf'), PATH_TO_DOCKING_CONFIGS / (core_subpocket + '.flexx'), PATH_TO_DOCKING_RESULTS / ('thread_' + str(thread_id) + '_core_fragment.sdf'), PATH_FLEXX)
    docking_utils.remove_files(PATH_TO_DOCKING_RESULTS / ('thread_' + str(thread_id) + '_core_fragment.sdf'), PATH_TO_SDF_FRAGMENTS / ('thread_' + str(thread_id) + '_core_fragment.sdf'))

    for conformer in res:   # safe every resulting pose within the fragment
        pose = docking_utils.Pose(conformer, float(conformer.GetProp('BIOSOLVEIT.DOCKING_SCORE')))
        core_fragment.add_pose(pose)
    logging.debug(str(len(res)) + ' poses have been generated')
    if len(res):
        # if fragment could be docked, save fragment (including it's poses)
        logging.debug("Best docking score: " + str(core_fragment.min_docking_score))
        return [core_fragment]
    # if no pose was found: return empty list
    return []

def template_docking_task(PATH_TO_SDF_FRAGMENTS,  PATH_TO_DOCKING_CONFIGS, PATH_TO_DOCKING_RESULTS, PATH_FLEXX, PATH_TO_TEMPLATES, subpocket, recombination: docking_utils.Recombination, poses):
    """
    runs a FlexX template docking task. It should be used for multithreading.

    Returns
    ----------
    List of all poses (docking result) as Ligand (Properties: docking-score, pose)

    Parameters
    ----------
    PATH_TO_SDF_FRAGMENTS: pathlib.path
        Path to directory where sdf files of thefragments are stored intermediately
    PATH_TO_DOCKING_CONFIGS: pathlib.path
        Path to FlexX-config file.
    PATH_TO_DOCKING_RESULTS: pathlib.path
        Path to diectory of output files
    PATH_FLEXX: pathlib.path
        Path to FlexX
    PATH_TO_TEMPLATES: pathlib.path
        Path to directory where sdf files of the templates are stored intermediately
    subpocket: str
        Subpocket where to place the fragment
    recombination: Recombination
        Recombination that should be docked in this task
    poses: List(Pose)
        Poses that are used as templates
    
    """
    thread_id = threading.get_ident()
    logging.debug('Template docking of Recombination: ' + str(recombination.fragments))
    fragment = docking_utils.from_recombination(recombination)

    # write recombination that should be docked to file
    if not fragment.to_sdf(PATH_TO_SDF_FRAGMENTS /('thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf')):
        logging.error('Could not write fragemnt: ' + str(fragment.fragment_ids) + ' to files due to 3d-generation-error')
        return []

    # for every choosen pose: perform template docking
    for i, pose in enumerate(poses):
        logging.debug('Template ' + str(i + 1) + ' / ' + str(len(poses)))
        # write template to sdf file
        with Chem.SDWriter(str(PATH_TO_TEMPLATES / ('thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'))) as w:
            w.write(pose.ROMol)
        # template docking (FlexX)
        res = docking_utils.template_docking(PATH_TO_SDF_FRAGMENTS / ('thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'), PATH_TO_TEMPLATES / ('thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'), PATH_TO_DOCKING_CONFIGS / (subpocket + '.flexx'), PATH_TO_DOCKING_RESULTS / ('thread_' + str(thread_id) + '_' + 'fragments.sdf'), PATH_FLEXX)
        docking_utils.remove_files(PATH_TO_DOCKING_RESULTS / ('thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'), PATH_TO_SDF_FRAGMENTS / ('thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'))

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
        return [fragment]
    return []

