import threading
import docking_utils 
import logging
from rdkit import Chem

from classes.config import Config
from classes.recombination import Recombination
from classes.ligand import Pose
from classes import ligand

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
    logging.debug('Docking of ' + config.core_subpocket + "-Fragment: " + str(core_fragment.fragment_ids[config.core_subpocket]))
    # safe fragment as sdf file
    if not core_fragment.to_sdf(config.path_temp/ ('cand_thread_' + str(thread_id) + '_core_fragment.sdf')):
        # could not generate 3D-conformation
        logging.error('Could not write Fragemnt: ' + str(core_fragment.fragment_ids) + ' to files due to 3d-generation-error')
        return [core_fragment]

    res_docking = docking_utils.core_docking(config.path_temp/ ('cand_thread_' + str(thread_id) + '_core_fragment.sdf'), config.path_temp / (config.core_subpocket + '.flexx'), 
                                             config.path_temp / ('dock_thread_' + str(thread_id) + '_core_fragment.sdf'), config.path_flexx, core_fragment.fragment_ids, core_fragment.smiles_dummy, core_fragment.smiles)
    
    violations = []  # track violations
    sum_hyde_displacement_violations = 0   # used to calculate the mean displacement of hyde
    sum_hyde_displacement = 0   # used to calculate the mean displacement of hyde

    if config.use_hyde and len(res_docking):
        # if path to hyde is given: perform hyde_scoring and opt.
        res_hyde = docking_utils.hyde_scoring(config.path_temp / ('dock_thread_' + str(thread_id) + '_core_fragment.sdf'), config.path_structure_config / (config.core_subpocket + '.hydescorer'), 
                                              config.path_temp / ('hyde_thread_' + str(thread_id) + '_core_fragment.sdf'), config.path_hyde,
                                              core_fragment.fragment_ids, core_fragment.smiles_dummy, core_fragment.smiles)

        docking_utils.remove_files(config.path_temp / ('dock_thread_' + str(thread_id) + '_core_fragment.sdf'), config.path_temp/ ('cand_thread_' + str(thread_id) + '_core_fragment.sdf'), 
                                   config.path_temp / ('hyde_thread_' + str(thread_id) + '_core_fragment.sdf'))

        for i, (conformer_docking, conformer_hyde) in enumerate(zip(res_docking, res_hyde)):   # safe every resulting pose within the fragment
            displacement = docking_utils.calc_distance_matrix([conformer_hyde, conformer_docking])[0]
            if displacement > config.hyde_displacement_cutoff:
                logging.warning(f"""VIOLATION: RMSD between HYDE and docking pose {core_fragment.fragment_ids} {i}: 
                                {displacement}""")
                violations.append((conformer_docking, conformer_hyde))
                sum_hyde_displacement_violations += displacement
                # drop pose if rmds > hyde_cutoff
                continue
            sum_hyde_displacement += displacement
            pose = Pose(conformer_hyde, float(conformer_hyde.GetProp('BIOSOLVEIT.DOCKING_SCORE')))
            pose.binding_affinity_upper = float(conformer_hyde.GetProp('BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_UPPER_BOUNDARY [nM]'))
            pose.binding_affinity_lower = float(conformer_hyde.GetProp('BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_LOWER_BOUNDARY [nM]')) 
            core_fragment.add_pose(pose)
    else:
        docking_utils.remove_files(config.path_temp / ('dock_thread_' + str(thread_id) + '_core_fragment.sdf'), config.path_temp/ ('cand_thread_' + str(thread_id) + '_core_fragment.sdf'))

        for conformer_docking in res_docking:   # safe every resulting pose within the fragment
            pose = Pose(conformer_docking, float(conformer_docking.GetProp('BIOSOLVEIT.DOCKING_SCORE'))) 
            core_fragment.add_pose(pose)

    logging.debug(str(len(core_fragment.poses)) + ' poses have been generated')

    if len(core_fragment.poses):
        # if fragment could be docked, save fragment (including it's poses)
        logging.debug("Best score: " + str(core_fragment.min_binding_affinity or core_fragment.min_docking_score))
        core_fragment.mean_hyde_displacement_undropped = sum_hyde_displacement / len(core_fragment.poses)
        core_fragment.mean_hyde_displacement = (sum_hyde_displacement + sum_hyde_displacement_violations) / (len(core_fragment.poses) + len(violations))
        core_fragment.num_hyde_violations = len(violations)
        return [core_fragment, violations]
    # return empty list if no pose was found
    return []

def template_docking_task(config: Config, subpocket: str, recombination: Recombination, poses: list):
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
    logging.debug('Template docking of Recombination: ' + str(recombination.fragments))
    fragment = ligand.from_recombination(recombination)

    # write recombination that should be docked to file
    if not fragment.to_sdf(config.path_temp /('cand_thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf')):
        logging.error('Could not write fragemnt: ' + str(fragment.fragment_ids) + ' to files due to 3d-generation-error')
        return [fragment]
    
    sum_hyde_displacement_violations = 0   # used to calculate the mean displacement of hyde
    sum_hyde_displacement = 0   # used to calculate the mean displacement of hyde
    violations = []
    succesfully_docked = 0
    

    # for every choosen pose: perform template docking
    for i, pose in enumerate(poses):
        logging.debug('Template ' + str(i + 1) + ' / ' + str(len(poses)))
        # write template to sdf file
        with Chem.SDWriter(str(config.path_temp / ('temp_thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'))) as w:
            w.write(pose.ROMol)
        # template docking (FlexX)
        res_docking = docking_utils.template_docking(config.path_temp /('cand_thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'), config.path_temp / ('temp_thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'), 
                                                     config.path_structure_config / (subpocket + '.flexx'), config.path_temp / ('dock_thread_' + str(thread_id) + '_' + 'fragments.sdf'), config.path_flexx,
                                                     fragment.fragment_ids, fragment.smiles_dummy, fragment.smiles)
        
        succesfully_docked += 1 if len(res_docking) else 0 

        if len(res_docking) and config.use_hyde:
            res_hyde = docking_utils.hyde_scoring(config.path_temp / ('dock_thread_' + str(thread_id) + '_' + 'fragments.sdf'), config.path_structure_config / (subpocket + '.hydescorer'), 
                                                  config.path_temp / ('hyde_thread_' + str(thread_id) + '_fragment.sdf'), config.path_hyde,
                                                  fragment.fragment_ids, fragment.smiles_dummy, fragment.smiles)
            # remove files containg docking results and template
            docking_utils.remove_files(config.path_temp / ('dock_thread_' + str(thread_id) + '_' + 'fragments.sdf'), config.path_temp / ('temp_thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'),
                                       config.path_temp / ('hyde_thread_' + str(thread_id) + '_fragment.sdf'))

            # safe resulting poses within fragment
            for i, (conformer_hyde, conformer_docking) in enumerate(zip(res_hyde, res_docking)):   # safe every resulting pose within the fragment
                pose = Pose(conformer_hyde, float(conformer_hyde.GetProp('BIOSOLVEIT.DOCKING_SCORE')))
                pose.binding_affinity_upper = float(conformer_hyde.GetProp('BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_UPPER_BOUNDARY [nM]'))
                pose.binding_affinity_lower = float(conformer_hyde.GetProp('BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_LOWER_BOUNDARY [nM]')) 
                rmsd = docking_utils.calc_distance_matrix([conformer_hyde, conformer_docking])[0]

                if rmsd > config.hyde_displacement_cutoff:
                    violations.append((conformer_docking, conformer_hyde))
                    sum_hyde_displacement_violations += rmsd
                    logging.warning(f"VIOLATION: RMSD between HYDE and docking pose {fragment.fragment_ids} {i}: {rmsd}")
                else:
                    sum_hyde_displacement += rmsd
                    fragment.add_pose(pose)
        else:
            # remove files containg docking results and template
            docking_utils.remove_files(config.path_temp / ('dock_thread_' + str(thread_id) + '_' + 'fragments.sdf'), config.path_temp / ('temp_thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'))
            for conformer in res_docking:   # safe every resulting pose within the fragment
                pose = Pose(conformer, float(conformer.GetProp('BIOSOLVEIT.DOCKING_SCORE'))) 
                fragment.add_pose(pose)

        logging.debug(str(len(fragment.poses)) + ' poses have been generated')

        if len(fragment.poses):
            logging.debug("Best  score: " + str(fragment.min_binding_affinity or fragment.min_docking_score))

    # remove molecule file
    docking_utils.remove_files(config.path_temp /('cand_thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'))

    if len(fragment.poses):
        fragment.mean_hyde_displacement_undropped = sum_hyde_displacement / len(fragment.poses)
        fragment.mean_hyde_displacement = (sum_hyde_displacement + sum_hyde_displacement_violations) / (len(fragment.poses) + len(violations))
        fragment.num_hyde_violations = len(violations)
        # safe recombination as result only if at least one pose was generated
        return [fragment, violations, succesfully_docked]
    return []

