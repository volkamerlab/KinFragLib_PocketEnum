import pandas as pd
from rdkit import Chem, Geometry
from rdkit.Chem import Descriptors, Lipinski, QED, Draw, AllChem
import json
from kinfraglib import utils
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit.Chem.Draw import rdMolDraw2D

from collections import defaultdict

from PIL import Image as pilImage
from io import BytesIO


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
    to_float_probs = [
        "BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_LOWER_BOUNDARY [nM]",
        "BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_UPPER_BOUNDARY [nM]",
        "BIOSOLVEIT.LOGP",
        "BIOSOLVEIT.MOLECULAR_WEIGHT",
        "BIOSOLVEIT.TPSA",
        "BIOSOLVEIT.DOCKING_SCORE",
    ]
    to_dict_probs = ["fragment_ids", "smiles_fragments_dummy", "smiles_fragments"]

    data = [
        [mol]  # 3D conformation of ligand
        + [
            (
                float(mol.GetProp(prop_name))
                if prop_name in to_float_probs
                else mol.GetProp(prop_name).replace("'", '"')
            )
            for prop_name in mol.GetPropNames()
        ]
        for mol in Chem.SDMolSupplier(str(path_to_mols), removeHs=False)
    ]

    data_df = (
        pd.DataFrame(data, columns=["ROMol"] + list(data[0][0].GetPropNames()))
        .apply(lambda x: x.apply(float) if x.name in to_float_probs else x)
        .apply(lambda x: x.apply(json.loads) if x.name in to_dict_probs else x)
    )

    data_df["num_fragments"] = [len(x) for x in data_df["fragment_ids"]]
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
