import copy
import json
from collections import defaultdict
from copy import deepcopy
from io import BytesIO
import os
from math import pi
from statistics import mean

import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client
from kinfraglib.utils import standardize_mol
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections import register_projection
from matplotlib.projections.polar import PolarAxes
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D
from PIL import Image as pilImage
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, DataStructs, Draw, rdFMCS
from rdkit.Chem.Draw import rdDepictor, rdMolDraw2D
from tqdm.auto import tqdm

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

def pplistdir(path, pre = ""):
    print(pre + " - " + str(path).split("/")[-1])
    if os.path.isfile(path):
        return
    
    for file in sorted(os.listdir(path)):
        pplistdir(path / file, pre + " |\t")

    return

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


def get_fragment_ids(mol, subpockets):
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


def get_fragment_smiles(mol, subpockets, dummy_atoms=False):
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

    fragment_smiles = json.loads(
        mol.GetProp(
            "smiles_fragments_dummy" if dummy_atoms else "smiles_fragments"
        ).replace("'", '"')
    )

    return [fragment_smiles.get(sp) for sp in subpockets]


def read_mols(path_to_mols, remove_dupliactes=True):
    """
    Read ligands from result file.

    Parameters
    ----------
    path_to_lib : str
        Path to results .sdf file.
    remove_duplicates : bool
        If true, molecules are deduplicated such the only the one with the best binding affinity is kept.
        Note: duplicates might appear since similar fragments might have different dummy atoms.

    Returns
    -------
    pandas.DataFrame
        ligands details details, i.e. SMILES, and RDKit molecules.
    """

    data = [
        [
            mol,
            get_binding_affinity(mol),
            float(mol.GetProp("BIOSOLVEIT.DOCKING_SCORE")),
            get_number_of_fragments(mol),
            Chem.MolToInchi(standardize_mol(mol)),
        ]
        + get_fragment_smiles(mol, SUBPOCKETS)
        + get_fragment_smiles(mol, SUBPOCKETS, dummy_atoms=True)
        + get_fragment_ids(mol, SUBPOCKETS)
        for mol in Chem.SDMolSupplier(str(path_to_mols), removeHs=False)
    ]

    data_df = pd.DataFrame(
        data,
        columns=["ROMol", "binding_affinity", "docking_score", "num_fragments", "inchi"]
        + [sp + "_smiles" for sp in SUBPOCKETS]
        + [sp + "_smiles_dummy" for sp in SUBPOCKETS]
        + SUBPOCKETS,
    ).sort_values(by=["binding_affinity"])

    # drop NA columns (columns where all entries are None/NA)
    data_df = data_df.dropna(axis=1, how="all")

    if remove_dupliactes:
        data_df = data_df.drop_duplicates(subset="inchi").reset_index(drop=True)

    return data_df


# code from adapted from https://greglandrum.github.io/rdkit-blog/posts/2021-08-07-rgd-and-highlighting.html
def highlight_scaffold(
    mol, patt, color, width=350, height=200, fillRings=True, legend=""
):
    """
    Highlights pattern in given molecule.

    Parameters
    ----------
    mol : ROMol
        molecule
    patt : ROMol
        pattern to highlight

    Returns
    -------
    PNG
        PNG of highlighted molecule
    """

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


# code from adapted from https://greglandrum.github.io/rdkit-blog/posts/2021-08-07-rgd-and-highlighting.html
def highlight_scaffold_multiple(
    ms, pattern, legends=None, n_per_row=5, sub_image_size=(250, 200)
):
    """
    Highlights pattern in respective molecules.

    Parameters
    ----------
    ms : list(ROMol)
        List of molecules
    patt : ROMol
        List of pattern to highlight in molcecule of ms with the same in index
    legends : list(string)
        Legends for ech molecule
    n_per_row : int
        Number of molecules per row
    sub_image_size : tuple
        Size of subimages

    Returns
    -------
    PNG
        Image with highlighted molecules
    """
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

    nRows = len(ms) // n_per_row
    if len(ms) % n_per_row:
        nRows += 1
    nCols = n_per_row
    imgSize = (sub_image_size[0] * nCols, sub_image_size[1] * nRows)
    res = pilImage.new("RGB", imgSize)

    for i, m in enumerate(ms):
        col = i % n_per_row
        row = i // n_per_row
        if legends:
            legend = legends[i]
        else:
            legend = ""
        png = highlight_scaffold(
            m,
            pattern[i],
            legend=legend,
            width=sub_image_size[0],
            height=sub_image_size[1],
            color=colors[i // n_per_row],
        )
        bio = BytesIO(png)
        img = pilImage.open(bio)
        res.paste(img, box=(col * sub_image_size[0], row * sub_image_size[1]))
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

    return highlight_scaffold_multiple(mols, scaffolds_pattern, legends=legend)


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
        Fingerprints of of all ligands

    Returns
    -------
    list
        Similarity matrix
    """
    return [DataStructs.BulkTanimotoSimilarity(fp, fingerprints) for fp in fingerprints]


def get_chembl_compounds_from_id(chembl_ids: list) -> pd.DataFrame:
    """
    Retrives ChEMBL compounds from given IDs

    Parameters
    ----------
    chembl_ids : list()
        List of ChEMBL IDs

    Returns
    -------
    DatFrame
        ChEMBL compounds (ID, smiles, ROMol)
    """

    # get compounds from client
    compounds_api = new_client.molecule
    compounds_provider = compounds_api.filter(molecule_chembl_id__in=chembl_ids).only(
        "molecule_chembl_id", "molecule_structures"
    )

    compounds = [
        [
            compound["molecule_chembl_id"],
            Chem.rdmolfiles.MolFromMolBlock(compound["molecule_structures"]["molfile"]),
            compound["molecule_structures"]["standard_inchi"],
        ]
        for compound in tqdm(compounds_provider)
    ]

    return pd.DataFrame(
        compounds,
        columns=["chembl_id", "ROMol", "inchi"],
    )


def highlight_molecules(
    molecules, mcs, number, label=True, same_orientation=True, **kwargs
):
    """Highlight the MCS in our query molecules
    Function taken and adapted from https://github.com/volkamerlab/teachopencadd/blob/master/teachopencadd/talktorials/T006_compound_maximum_common_substructures/talktorial.ipynb
    """
    molecules = deepcopy(molecules)
    # convert MCS to molecule
    pattern = Chem.MolFromSmarts(mcs.smartsString)
    rdDepictor.Compute2DCoords(pattern)
    # find the matching atoms in each molecule
    matching = [molecule.GetSubstructMatch(pattern) for molecule in molecules[:number]]

    legends = None
    if label is True:
        legends = [x.GetProp("_Name") for x in molecules]

    # Align by matched substructure so they are depicted in the same orientation
    # Adapted from: https://gist.github.com/greglandrum/82d9a86acb3b00d3bb1df502779a5810
    if same_orientation:
        for mol in molecules:
            rdDepictor.GenerateDepictionMatching2DStructure(mol, pattern)

    return Draw.MolsToGridImage(
        molecules[:number],
        legends=legends,
        molsPerRow=2,
        highlightAtomLists=matching[:number],
        subImgSize=(300, 200),
        returnPNG=False,
        **kwargs,
    )


def save_chemb_mcs_to_file(
    most_similar_chembl_compounds, most_similar_pka_compounds, directory
):
    """Determines the MCS between each ChEMBL compound and the respective given compounds and saves them to png"""

    for id in most_similar_chembl_compounds.index:
        mol_chembl = [most_similar_chembl_compounds["ROMol"][id]]
        chembl_id = most_similar_chembl_compounds["chembl_id"][id]
        mol_chembl[0].SetProp("_Name", "")
        mols = mol_chembl + [
            most_similar_pka_compounds[
                most_similar_pka_compounds["most_similar_chembl_ligand.chembl_id"]
                == chembl_id
            ]["ROMol"].iloc[0]
        ]
        mcs1 = rdFMCS.FindMCS(mols, ringMatchesRingOnly=True, completeRingsOnly=True)
        m1 = Chem.MolFromSmarts(mcs1.smartsString)
        for i, mol in enumerate(mols):
            rdDepictor.Compute2DCoords(mol)
            if i:
                mol.SetProp("_Name", "")
        img = highlight_molecules(mols, mcs1, len(mols), same_orientation=True)
        img.save(directory / f"{chembl_id}.png")

# code is copied and adapted from https://matplotlib.org/stable/gallery/specialty_plots/radar_chart.html
def radar_factory(num_vars, frame='polygon'):
    """
    Returns radar chart with num_vars axes.

    Parameters
    ----------
    num_vars : int
        Number of axis (variables).
    frame : {'circle', 'polygon'}
        Shape of frame surrounding Axes.
    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2 * pi, num_vars, endpoint=False)

    class RadarTransform(PolarAxes.PolarTransform):
            # Paths with non-unit interpolation steps correspond to gridlines,
            # in which case we force interpolation (to defeat PolarTransform's
            # autoconversion to circular arcs).
        def transform_path_non_affine(self, path):
            if path._interpolation_steps > 1:
                path = path.interpolated(num_vars)
            return Path(self.transform(path.vertices), path.codes)

    class RadarAxes(PolarAxes):
        name = "radar"
        PolarTransform = RadarTransform

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # rotate plot so that the first axis is at the top
            self.set_theta_zero_location("N")

        def fill(self, *args, closed=True, **kwargs):
            """Override fill so that line is closed by default"""
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            if x[0] != x[-1]:
                x = np.append(x, x[0])
                y = np.append(y, y[0])
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
            # in axes coordinates.
            if frame == "circle":
                return Circle((0.5, 0.5), 0.5)
            elif frame == "polygon":
                return RegularPolygon((0.5, 0.5), num_vars, radius=0.5, edgecolor="k")
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

        def _gen_axes_spines(self):
            if frame == "circle":
                return super()._gen_axes_spines()
            elif frame == "polygon":
                # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                spine = Spine(
                    axes=self,
                    spine_type="circle",
                    path=Path.unit_regular_polygon(num_vars),
                )
                # unit_regular_polygon gives a polygon of radius 1 centered at
                # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                # 0.5) in axes coordinates.
                spine.set_transform(
                    Affine2D().scale(0.5).translate(0.5, 0.5) + self.transAxes
                )
                return {"polar": spine}
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

    register_projection(RadarAxes)
    return theta
