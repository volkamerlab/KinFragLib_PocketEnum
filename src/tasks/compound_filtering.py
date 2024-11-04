import json
import logging
import math
import random
import statistics
from pathlib import Path
from typing import Tuple

from rdkit import DataStructs
from rdkit.Chem import Draw, rdFingerprintGenerator
from rdkit.ML.Cluster import Butina

from classes.config import Config


def cluster_based_compound_filtering(
    docking_results: list, num_candidates: int, config: Config, SP: str
) -> list:
    """
    Chooses *num_candidates* ligands according to docking result and diversity.
    Ligands are clustered according to Tanimoto similarity first and a cluster score is calculated
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

    if len(docking_results) <= num_candidates:
        return docking_results

    # cluster compounds
    clusters = _cluster_compounds(docking_results)

    # sort ligands in clusters (ascending)
    for cluster in clusters:
        if config.use_hyde:
            cluster.sort(key=lambda l: l.min_binding_affinity)
        else:
            cluster.sort(key=lambda l: l.min_docking_score)

    # per cluster: calculate score
    scores = _calc_cluster_scores(clusters, config.use_hyde)

    # sort clusters and scores
    clusters, scores = _sort_clusters(clusters, scores)

    logging.info(f"Created {len(clusters)} clusters")
    logging.info(
        f"Number of clusters with 1 compound: {sum(len(c) == 1 for c in clusters)}"
    )
    logging.info(f"Avg. cluster size: {statistics.mean(len(c) for c in clusters)}")
    logging.info(f"Max. cluster size: {max(len(c) for c in clusters)}")

    stats = {
        "NumberClustersPre": len(clusters),
        "SingletonClusters": {sum(len(c) == 1 for c in clusters)},
        "MeanClusterSize": {statistics.mean(len(c) for c in clusters)},
        "MaxClusterSize": {max(len(c) for c in clusters)},
    }
    stats["ScoresPre"] = scores

    # consider only best 90% of cluster
    scores = scores[: math.ceil(len(scores) * 0.9)]
    clusters = clusters[: len(scores)]

    stats["MolScoresWithinCluster"] = [
        [ligand.get_min_score() for ligand in cluster] for cluster in clusters
    ]
    stats["SingletonClusters"] = sum(len(c) == 1 for c in clusters)
    stats["MeanClusterSize"] = statistics.mean(len(c) for c in clusters)
    stats["MaxClusterSize"] = max(len(c) for c in clusters)

    # softmax
    probabilities = _draw_distribution(scores, config.P)

    stats["Distribution"] = probabilities

    # store stats to file
    with open(config.path_results / f"{SP}_filtering_stats.json", "w") as json_file:
        json.dump(stats, json_file, indent=4)

    candidates = []

    for _ in range(num_candidates):
        # choose cluster
        choosen_cluster_idx = random.choices(range(len(clusters)), probabilities)[0]

        # since ligands are sorted within one clusters, it is sufficient to choose and remove the first ligand
        # in order to choose the best scored one
        selected_compound = clusters[choosen_cluster_idx].pop(
            0
        )  # select and remove compound
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
    Clusters compounds using Tanimoto similarity and Butina cluster

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
    fingerprints = [
        rdkit_gen.GetFingerprint(ligand.ROMol) for ligand in docking_results
    ]

    distance_matrix = _tanimoto_distance_matrix(fingerprints)

    clusters_idxs = Butina.ClusterData(
        distance_matrix, len(fingerprints), cutoff, isDistData=True
    )

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
        second_quratile_idx = len(cluster) - len(cluster) // 4
        # do not consider the last 25% (according to ligand score)
        # to somehow ignore negative outliers
        # do not ignore positvie outliers since they will be selected first => hence they should have a rather big influence
        # on the overall score

        scores_within_iqr = (
            ligand.min_binding_affinity if use_hyde_score else ligand.min_docking_score
            for ligand in cluster[:second_quratile_idx]
        )
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


def _draw_distribution(cluster_scores: list, p: int = 1) -> list:
    """
    Draws a probability distribution based on the inverse cluster scores using softmax

    Returns
    ----------
    Mapped propabilities for each cluster

    Parameters
    ----------
    cluster_scores: list(float)
        Scores of clusters
    """
    # min_max_scaling [0,100]Union
    x_min = min(cluster_scores)
    x_max = max(cluster_scores)

    if x_max == x_min:
        # all scores are equal or only one cluster exits
        # -> uniform distribution
        return [1 / len(cluster_scores) for _ in cluster_scores]

    min_max_scaling = lambda c: ((c - x_min) * 100) / (x_max - x_min)

    denominator = sum(math.exp(-min_max_scaling(c) / p) for c in cluster_scores)

    return [math.exp(-min_max_scaling(c) / p) / (denominator) for c in cluster_scores]


def _save_clusters_to_image(folder: Path, subpocket: str, clusters: list):
    """
    Saves mols of clusters with >= 2 molecules as image

    Parameters
    ----------
    folder: Path
        Path to folder, where images should be places
    subpocket: str
        Current subpocket
    Clusters: list()
        Clusters
    """
    for i, cluster in enumerate(clusters):
        if len(cluster) == 1:
            continue

        img = Draw.MolsToGridImage(
            [ligand.ROMol for ligand in cluster[:10]], molsPerRow=5, returnPNG=False
        )

        img.save(f"{folder}/cluster_{i}_{subpocket}.png")


def _sort_clusters(clusters: list, scores: list):
    """
    Sorts clusters and scores synchronous based on scores, such that they still correspont to each other

    Returns
    --------
    (list, list)
        Pair of sorted (clusters, scores)

    Parameters
    ----------
    clusters: list(list(Ligand))
        Clusters
    scores: list(float)
        Scores of clusters

    """
    # indexes that correspond to sorted clusters/scores
    sorted_indexes = sorted(range(len(scores)), key=lambda i: scores[i])

    sorted_clusters = [clusters[idx] for idx in sorted_indexes]
    sorted_scores = [scores[idx] for idx in sorted_indexes]

    return sorted_clusters, sorted_scores
