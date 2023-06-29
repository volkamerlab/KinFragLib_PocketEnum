import os
from queue import Queue
from threading import Thread
import docking_utils 
import logging
from rdkit import Chem

def core_docking_task(thread_id, PATH_TO_SDF_FRAGMENTS,  PATH_TO_DOCKING_CONFIGS, PATH_TO_DOCKING_RESULTS, PATH_FLEXX, docking_results, core_fragment, core_subpocket):
    logging.debug('Docking of ' + core_subpocket + "-Fragment: " + str(core_fragment.fragment_ids[core_subpocket]))
    # safe fragment as sdf file
    if not core_fragment.to_sdf(PATH_TO_SDF_FRAGMENTS / ('thread_' + str(thread_id) + '_core_fragment.sdf')):
        # could not generate 3D-conformation
        logging.error('Could not write Fragemnt: ' + str(core_fragment.fragment_ids) + ' to files due to 3d-generation-error')
        return

    res = docking_utils.core_docking(PATH_TO_SDF_FRAGMENTS / ('thread_' + str(thread_id) + '_core_fragment.sdf'), PATH_TO_DOCKING_CONFIGS / (core_subpocket + '.flexx'), PATH_TO_DOCKING_RESULTS / ('thread_' + str(thread_id) + '_core_fragment.sdf'), PATH_FLEXX)
    docking_utils.remove_files(PATH_TO_DOCKING_RESULTS / ('thread_' + str(thread_id) + '_core_fragment.sdf'), PATH_TO_SDF_FRAGMENTS / ('thread_' + str(thread_id) + '_core_fragment.sdf'))

    for conformer in res:   # safe every resulting pose within the fragment
        pose = docking_utils.Pose(conformer, float(conformer.GetProp('BIOSOLVEIT.DOCKING_SCORE')))
        core_fragment.add_pose(pose)
    logging.debug(str(len(res)) + ' poses have been generated')
    if len(res):
        # if fragment could be docked, save fragment (including it's poses)
        logging.debug("Best docking score: " + str(core_fragment.min_docking_score))
        docking_results.append(core_fragment)

def template_docking_task(thread_id, PATH_TO_SDF_FRAGMENTS,  PATH_TO_DOCKING_CONFIGS, PATH_TO_DOCKING_RESULTS, PATH_FLEXX, PATH_TO_TEMPLATES, docking_results, recombination, poses, subpocket):
    
    logging.debug('Template docking of Recombination: ' + str(recombination.fragments))
    fragment = docking_utils.from_recombination(recombination)

    # write recombination that should be docked to file
    if not fragment.to_sdf(PATH_TO_SDF_FRAGMENTS /('thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf')):
        logging.error('Could not write fragemnt: ' + str(fragment.fragment_ids) + ' to files due to 3d-generation-error')
        return

    # for every choosen pose: perform template docking
    for i, pose in enumerate(poses):
        logging.debug('Template ' + str(i + 1) + ' / ' + str(len(poses)))
        # write template to sdf file
        with Chem.SDWriter(str(PATH_TO_TEMPLATES / ('thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'))) as w:
            w.write(pose.ROMol)
        # template docking (FlexX)
        res = docking_utils.template_docking(PATH_TO_SDF_FRAGMENTS / ('thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'), PATH_TO_TEMPLATES / ('thread_' + str(thread_id) + '_' + subpocket + '_fragment.sdf'), PATH_TO_DOCKING_CONFIGS / (subpocket + '.flexx'), PATH_TO_DOCKING_RESULTS / ('thread_' + str(thread_id) + '_' + 'fragments.sdf'), PATH_FLEXX)

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

class DockingWorker(Thread):

    def __init__(self, queue, thead_id, PATH_TO_SDF_FRAGMENTS,  PATH_TO_DOCKING_CONFIGS, PATH_TO_DOCKING_RESULTS, PATH_TO_TEMPLATES, PATH_FLEXX):
        Thread.__init__(self)
        self.queue = queue
        self.thread_id = thead_id
        self.PATH_TO_SDF_FRAGMENTS = PATH_TO_SDF_FRAGMENTS
        self.PATH_TO_DOCKING_CONFIGS = PATH_TO_DOCKING_CONFIGS
        self.PATH_TO_DOCKING_RESULTS = PATH_TO_DOCKING_RESULTS
        self.PATH_FLEXX = PATH_FLEXX
        self.PATH_TO_TEMPLATES = PATH_TO_TEMPLATES

    def run(self):
        while True:
            # Get the work from the queue and expand the tuple
            task, docking_results, ligand, subpocket = self.queue.get()
            try:
                if task == 'core_docking':
                    core_docking_task(self.thread_id, self.PATH_TO_SDF_FRAGMENTS,  self.PATH_TO_DOCKING_CONFIGS, self.PATH_TO_DOCKING_RESULTS, self.PATH_FLEXX, docking_results, ligand, subpocket)
                elif task == 'template_docking':
                    poses, recombination = ligand
                    template_docking_task(self.thread_id, self.PATH_TO_SDF_FRAGMENTS,  self.PATH_TO_DOCKING_CONFIGS, self.PATH_TO_DOCKING_RESULTS, self.PATH_FLEXX, self.PATH_TO_TEMPLATES, docking_results, recombination, poses, subpocket)
            finally:
                self.queue.task_done()


