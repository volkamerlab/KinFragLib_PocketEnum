import pandas as pd
from rdkit import Chem, Geometry
from rdkit.Chem import Descriptors, Lipinski, QED, Draw, AllChem, DataStructs
import json
from kinfraglib import utils
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit.Chem.Draw import rdMolDraw2D, rdDepictor
from statistics import mean

from collections import defaultdict

from PIL import Image as pilImage
from io import BytesIO
import copy

SUBPOCKETS = ["AP", "SE", "FP", "GA", "B1", "B2"]


def get_binding_affinity(mol):
    """
    Extracts and calculates HYDES binding affinity from ROMol (Mean of upper and and lower HYDE affinity estimate).

    Parameters
    ----------
    mol : ROMol
        molecule object with BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_UPPER_BOUNDARY [nM]
            and BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_LOWER_BOUNDARY [nM] as property

    Returns
    -------
    float
        mean binding affinity
    """
    return (
        float(mol.GetProp("BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_UPPER_BOUNDARY [nM]"))
        + float(mol.GetProp("BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_LOWER_BOUNDARY [nM]"))
    ) / 2


def get_number_of_fragments(mol):
    """
    Extracts number of fragemnts in given ligand

    Parameters
    ----------
    mol : ROMol
        molecule object with 'fragment_ids' property

    Returns
    -------
    int
        number of fragments
    """
    return len(json.loads(mol.GetProp("fragment_ids").replace("'", '"')))


def get_fragmnent_ids(mol, subpockets):
    """
    Extracts fragment_ids from given ligand

    Parameters
    ----------
    mol : ROMol
        molecule object with 'fragment_ids' property
    subpockets : list(str)
        list of subpocket names

    Returns
    -------
    list(int)
        fragment ids in order of subpockets list
    """

    fragment_ids = json.loads(mol.GetProp("fragment_ids").replace("'", '"'))
    return [fragment_ids.get(sp) for sp in subpockets]


def get_fragmnent_smiles(mol, subpockets):
    """
    Extracts smiles of fragments from given ligand

    Parameters
    ----------
    mol : ROMol
        molecule object with 'smiles_fragments' property
    subpockets : list(str)
        list of subpocket names

    Returns
    -------
    list(int)
        fragment ids in order of subpockets list
    """

    fragment_smiles = json.loads(mol.GetProp("smiles_fragments").replace("'", '"'))

    return [fragment_smiles.get(sp) for sp in subpockets]


def read_mols(path_to_mols):
    """
    Read ligands from result file.

    Parameters
    ----------
    path_to_lib : str
        Path to results .sdf file.

    Returns
    -------
    pandas.DataFrame
        ligands details details, i.e. SMILES, and RDKit molecules.
    """

    to_dict_probs = ["fragment_ids", "smiles_fragments_dummy", "smiles_fragments"]

    data = [
        [
            mol,
            get_binding_affinity(mol),
            float(mol.GetProp("BIOSOLVEIT.DOCKING_SCORE")),
            get_number_of_fragments(mol),
            Chem.MolToInchi(utils.standardize_mol(mol)),
        ]  # 3D conformation of ligand
        + get_fragmnent_ids(mol, SUBPOCKETS)
        + get_fragmnent_smiles(mol, SUBPOCKETS)
        for mol in Chem.SDMolSupplier(str(path_to_mols), removeHs=False)
    ]

    data_df = pd.DataFrame(
        data,
        columns=["ROMol", "binding_affinity", "docking_score", "num_fragments", "inchi"]
        + SUBPOCKETS
        + [sp + "_smiles" for sp in SUBPOCKETS],
    ).apply(lambda x: x.apply(json.loads) if x.name in to_dict_probs else x)

    return data_df


# from adapted from https://greglandrum.github.io/rdkit-blog/posts/2021-08-07-rgd-and-highlighting.html
def highlight_scaffold(
    mol, patt, color, width=350, height=200, fillRings=True, legend=""
):
    # copy the molecule and core
    mol = Chem.Mol(mol)

    # ------------------
    #  set up our colormap
    #   the three choices here are all "colorblind" colormaps

    # ----------------------
    # Identify and store which atoms, bonds, and rings we'll be highlighting
    highlightatoms = defaultdict(list)
    highlightbonds = defaultdict(list)
    atomrads = {}
    widthmults = {}

    hit_ats = list(mol.GetSubstructMatch(patt))

    rings = []
    rinfo = mol.GetRingInfo()
    for at_idx in list(mol.GetSubstructMatch(patt)):
        highlightatoms[at_idx].append(color)
        atomrads[at_idx] = 0.4
    if fillRings:
        for aring in rinfo.AtomRings():
            tring = []
            allFound = True
            for aid in aring:
                if aid in hit_ats:
                    tring.append(aid)
            if allFound:
                rings.append((tring, color))
    for qbnd in patt.GetBonds():
        batom = hit_ats[qbnd.GetBeginAtomIdx()]
        eatom = hit_ats[qbnd.GetEndAtomIdx()]
        bndIdx = mol.GetBondBetweenAtoms(batom, eatom).GetIdx()
        highlightbonds[bndIdx].append(color)
        widthmults[bndIdx] = 2

    d2d = rdMolDraw2D.MolDraw2DCairo(width, height)
    dos = d2d.drawOptions()
    dos.useBWAtomPalette()

    dos.legendFontSize = 65

    # ----------------------
    # if we are filling rings, go ahead and do that first so that we draw
    # the molecule on top of the filled rings
    if fillRings and rings:
        # a hack to set the molecule scale
        d2d.DrawMoleculeWithHighlights(
            mol,
            legend,
            dict(highlightatoms),
            dict(highlightbonds),
            atomrads,
            widthmults,
        )
        d2d.ClearDrawing()
        conf = mol.GetConformer()
        for aring, color in rings:
            ps = []
            for aidx in aring:
                pos = Geometry.Point2D(conf.GetAtomPosition(aidx))
                ps.append(pos)
            d2d.SetFillPolys(True)
            d2d.SetColour(color)
            d2d.DrawPolygon(ps)
        dos.clearBackground = False

    # ----------------------
    # now draw the molecule, with highlights:
    d2d.DrawMoleculeWithHighlights(
        mol, legend, dict(highlightatoms), dict(highlightbonds), atomrads, widthmults
    )
    d2d.FinishDrawing()
    png = d2d.GetDrawingText()
    return png


def draw_multiple(ms, pattern, legends=None, nPerRow=5, subImageSize=(250, 200)):
    # "Tol" colormap from https://davidmathlogic.com/colorblind
    # colors = [(51,34,136),(17,119,51),(68,170,153),(136,204,238),(221,204,119),(204,102,119),(170,68,153),(136,34,85)]
    # "IBM" colormap from https://davidmathlogic.com/colorblind
    # colors = [(100,143,255),(120,94,240),(220,38,127),(254,97,0),(255,176,0)]
    # Okabe_Ito colormap from https://jfly.uni-koeln.de/color/
    colors = [
        (0, 158, 115),
        (86, 180, 233),
        (204, 121, 167),
        (0, 114, 178),
        (230, 159, 0),
        (213, 94, 0),
        (240, 228, 66),
    ]
    for i, x in enumerate(colors):
        colors[i] = tuple(y / 255 for y in x)

    nRows = len(ms) // nPerRow
    if len(ms) % nPerRow:
        nRows += 1
    nCols = nPerRow
    imgSize = (subImageSize[0] * nCols, subImageSize[1] * nRows)
    res = pilImage.new("RGB", imgSize)

    for i, m in enumerate(ms):
        col = i % nPerRow
        row = i // nPerRow
        if legends:
            legend = legends[i]
        else:
            legend = ""
        png = highlight_scaffold(
            m,
            pattern[i],
            legend=legend,
            width=subImageSize[0],
            height=subImageSize[1],
            color=colors[i // nPerRow],
        )
        bio = BytesIO(png)
        img = pilImage.open(bio)
        res.paste(img, box=(col * subImageSize[0], row * subImageSize[1]))
    bio = BytesIO()
    res.save(bio, format="PNG")
    return bio.getvalue()


def draw_colored_scaffold_ligands(data):
    """
    Draws ligands with colored scaffolds

    Parameters
    ----------
    data : pandas.DataFrame
        Compounds

    """

    # prepare colored ligands
    mols = []
    scaffolds_pattern = []
    legend = []
    for i in range(5):
        scaffold_compounds = (
            data[data["scaffold_id"] == i].copy().drop_duplicates(subset=["inchi"])
        )
        _mols = [copy.deepcopy(mol) for mol in scaffold_compounds["ROMol"][:5]]
        for mol in _mols:
            rdDepictor.Compute2DCoords(mol)
        mols += _mols
        _patt = [
            AllChem.MolFromSmarts(smiles)
            for smiles in scaffold_compounds["Murcko_SMILES"][:5]
        ]
        for mol in _patt:
            rdDepictor.Compute2DCoords(mol)
        for mol in _mols:
            rdDepictor.GenerateDepictionMatching2DStructure(mol, _patt[0])
        scaffolds_pattern += _patt
        legend += [f"Scaffold {i} compound" for _ in range(5)]

    return draw_multiple(mols, scaffolds_pattern, legends=legend)


def tanimoto_distance_matrix(fp_list):
    """Calculate distance matrix for fingerprint list"""
    dissimilarity_matrix = []
    for i in range(1, len(fp_list)):
        similarities = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list[:i])
        dissimilarity_matrix.append([1 - x for x in similarities])
    return dissimilarity_matrix


def mean_tanimoto_similarity(all_fingerprints, target_fingerprint) -> float:
    """
    Calculates mean tanimoto similarity between fingerprint and all_fingerprints (excluding self similarity)

    Parameters
    ----------
    fingerprint : BitVec
        Fingerprint of target ligand
    all_fingerprints : Iterable
        Fingerprints of of all ligans

    Returns
    -------
    float
        Mean tanimoto similarity between fingerprint and all_fingerprints (excluding self similarity)
    """
    similarities = DataStructs.BulkTanimotoSimilarity(
        target_fingerprint, all_fingerprints
    )
    similarities.remove(1)  # self similarity
    return mean(similarities)


def max_tanimoto_similarity(all_fingerprints, target_fingerprint) -> float:
    """
    Calculates the maximum tanimoto similarity between fingerprint and all_fingerprints (excluding self similarity)

    Parameters
    ----------
    fingerprint : BitVec
        Fingerprint of target ligand
    all_fingerprints : Iterable
        Fingerprints of of all ligans

    Returns
    -------
    float
        Maximum tanimoto similarity between fingerprint and all_fingerprints (excluding self similarity)
    """
    similarities = DataStructs.BulkTanimotoSimilarity(
        target_fingerprint, all_fingerprints
    )
    similarities.remove(1)  # self similarity
    return max(similarities)


def tanimoto_similarity_matrix(fingerprints) -> list:
    """
    Calculates the Tanimoto similarity matrix for the given fingerprints

    Parameters
    ----------
    all_fingerprints : Iterable
        Fingerprints of of all ligans

    Returns
    -------
    list
        Similarity matrix
    """
    return [DataStructs.BulkTanimotoSimilarity(fp, fingerprints) for fp in fingerprints]
