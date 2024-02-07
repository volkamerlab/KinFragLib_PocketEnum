import logging
import math
import random
import statistics

from rdkit import DataStructs
from rdkit.Chem import Draw, rdFingerprintGenerator
from rdkit.ML.Cluster import Butina

from classes.config import Config
from classes.ligand import Ligand


def cluster_based_compound_filtering(docking_results: list, num_candidates: int, config: Config, SP: str) -> list:
    # TODO set seed
    """
        Chooses *num_candidates* ligands according to docking result and diversity.
        Ligands are clustered according to TODO (Tanimoto similarity?) first and cluster score is calculated 
        based on the mean of the second quartile of the ligand docking scores. 
        Iterativly the best scored ligands are choosen from randomly drawn clusters based on a cluster-score based distribution. 

        Returns
        ----------
        List of selected candidates as Ligand

        Parameters
        ----------
        docking_results: list(Ligand)
            ligands that should be filtered
        num_candidates: int
            ligands/fragments to choose

        
        """
    # TODO later add
    # if len(docking_results) <= num_candidates:
    #   return docking_results
    # cluster compounds 
    clusters = _cluster_compounds(docking_results)

    logging.info(f"Created {len(clusters)} clusters")
    logging.info(f"Number of clusters with 1 compound: {sum(len(c) == 1 for c in clusters)}")
    logging.info(f"Avg. cluster size: {statistics.mean(len(c) for c in clusters)}")
    logging.info(f"Max. cluster size: {max(len(c) for c in clusters)}")

    # TODO plot cluster score against avg. 
    # for evaluation purpose:
    # std_dev_per_cluster = [statistics.stdev(ligand.min_binding_affinity if config.use_hyde else ligand.min_docking_score for ligand in cluster) 
    #                       for cluster in clusters]
    mean_score_per_cluster = [statistics.mean(ligand.min_binding_affinity if config.use_hyde else ligand.min_docking_score for ligand in cluster) 
                           for cluster in clusters]
    max_score_per_cluster = [max(ligand.min_binding_affinity if config.use_hyde else ligand.min_docking_score for ligand in cluster) 
                           for cluster in clusters]
    
    # logging.info(f"Mean standard derivation of ligand scores within clusters")
    #TODO maybe line plot
    
    for i, cluster in enumerate(clusters):
        img = Draw.MolsToGridImage(
            [ligand.ROMol for ligand in cluster[:10]],
            molsPerRow=5,
            returnPNG=False
        )

        img.save(f"cluster_{i}_{SP}.png")


    # sort ligands in clusters (ascending)
    for cluster in clusters:
        if config.use_hyde:
            cluster.sort(key=lambda l: l.min_binding_affinity)
        else:
            cluster.sort(key=lambda l: l.min_docking_score)

    # per cluster: calculate score
    scores = _calc_cluster_scores(clusters, config.use_hyde)

    # select k clusters (?) or all clusters with nm <1000 ?
    # TODO skip this for now -> can be changed

    # sort compounds per cluster 

    # softmax
    probabilities = _draw_distribution(scores)

    candidates = []
    for _ in range(num_candidates):
        # choose cluster
        choosen_cluster_idx = random.choices(range(len(clusters)), probabilities)[0]

        # since ligands are sorted within one clusters, it is sufficient to choose and remove the first ligand
        # in order to choose the best scored one
        selected_compound = clusters[choosen_cluster_idx].pop(0) # select and remove compound
        candidates.append(selected_compound)

        if not len(clusters[choosen_cluster_idx]):
            # if selected cluster is empty now -> remove it and update propabilities
            del clusters[choosen_cluster_idx]
            del scores[choosen_cluster_idx]

            if not len(clusters):
                # if nor cluster (hence compound) left -> early break
                break
            probabilities = _draw_distribution(scores)

    return candidates

def _cluster_compounds(docking_results: list, cutoff: float = 0.2) -> list:
    """
    clusters compounds using Tanimoto similarity and Butina cluster

    Returns
    ----------
    Clusters (list(list(ligand)))

    Parameters
    ----------
    docking_results: list(Ligand)
        ligands that should be clustered
    """

    # code adapted from TeachOpenCADD (T005)
    # https://github.com/volkamerlab/teachopencadd/blob/master/teachopencadd/talktorials/T005_compound_clustering/talktorial.ipynb
    rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
    fingerprints = [rdkit_gen.GetFingerprint(ligand.ROMol) for ligand in docking_results]

    distance_matrix = _tanimoto_distance_matrix(fingerprints)

    clusters_idxs = Butina.ClusterData(distance_matrix, len(fingerprints), cutoff, isDistData=True)
    
    return [[docking_results[idx] for idx in cluster] for cluster in clusters_idxs]

def _calc_cluster_scores(clusters: list, use_hyde_score: bool) -> list:
    """
    Calculates a score for each clusters: mean of all scores within the interquartile range

    Returns
    ----------
    Mapped scores for each cluster

    Parameters
    ----------
    clusters: list(list(Ligand))
        clusters that should be scored
    """

    scores = []
    for cluster in clusters:
        # assuming that ligands clusters are sorted
        first_quratile_idx = len(cluster)//4 
        # do not consider first and last 25% (according to ligand score)
        # to somehow ignore outliers
        scores_within_iqr = (ligand.min_binding_affinity if use_hyde_score else ligand.min_docking_score 
                             for ligand in cluster[first_quratile_idx:len(cluster)-first_quratile_idx])
        scores.append(statistics.mean(scores_within_iqr))

    return scores

def _tanimoto_distance_matrix(fp_list):
    # copied from TeachOpenCADD: 
    # https://github.com/volkamerlab/teachopencadd/blob/master/teachopencadd/talktorials/T005_compound_clustering/talktorial.ipynb
    """Calculate distance matrix for fingerprint list"""
    dissimilarity_matrix = []
    # Notice how we are deliberately skipping the first and last items in the list
    # because we don't need to compare them against themselves
    for i in range(1, len(fp_list)):
        # Compare the current fingerprint against all the previous ones in the list
        similarities = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list[:i])
        # Since we need a distance matrix, calculate 1-x for every element in similarity matrix
        dissimilarity_matrix.extend([1 - x for x in similarities])
    return dissimilarity_matrix

def _draw_distribution(cluster_scores: list, p: int = 100) -> list:
    """
    Draws a probability distribution based on the cluster scores using softmin

    Returns
    ----------
    Mapped propabilities for each cluster

    Parameters
    ----------
    cluster_scores: list(float)
        Scores of clusters
    """
    print(cluster_scores)
    denominator = sum(math.exp(-c/p) for c in cluster_scores)
    print([math.exp(-c/p)/denominator for c in cluster_scores])
    return [math.exp(-c/p)/denominator for c in cluster_scores]
    