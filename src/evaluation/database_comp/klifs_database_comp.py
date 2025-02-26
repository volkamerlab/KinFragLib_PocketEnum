import argparse
import logging
import os
import sys
import pandas as pd

from pathlib import Path

from kinfraglib.utils import standardize_mol
from rdkit import Chem
from rdkit.Chem import PandasTools, rdFingerprintGenerator

from src.evaluation.database_comp.utils import most_similar_database_ligand
from src.evaluation.utils import read_mols

if __name__ == "__main__":
    # parse command line arguments
    parser = argparse.ArgumentParser(
        prog=sys.argv[0], description="Determines the most similar klifs compounds for each target compounds."
    )
    parser.add_argument(
        "-m",
        "--use_morgan",
        action="store_true",
        help="If flag is set, Morgan fingerprints will be used. Otherwise (by default), RDKit fingerprints are used.",
    )
    parser.add_argument(
        "-p",
        "--post_filtering",
        default=1000,
        help="Estimated binding affinity threshold <p> for postfiltering. If > 0, compounds with estimated binding affinity > <P> nM will be removed from given dataset as a post filtering step. Thus, if set to 0, this post filtering step will NOT be applied.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Filename of output CSV file.",
    )
    parser.add_argument(
        "-i",
        "--path_klifs_download",
        help="Path to KLIFS_download",
    )
    parser.add_argument(
        "-g",
        "--organism",
        default=["HUMAN"],
        nargs="*",
        help="Either HUMAN, MOUSE, or both",
    )
    parser.add_argument(
        "-k",
        "--kinase",
        nargs="*",
        default="*",
        help="OPTIONAL: specify kinase (s) (needs to correspond to filname in organism directory!)",
    )
    parser.add_argument(
        "-l",
        "--path_ligands",
        help="Path to SDF file of all compounds (ligands) that should be used as target for comparison.",
    )
    parser.add_argument(
        "-log",
        "--loglevel",
        default="INFO",
        help="Logging level (error, warning, info, or debug). Example --loglevel debug, default=info",
    )

    args = parser.parse_args()

    # init logging
    numeric_level = getattr(logging, args.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.loglevel)

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=numeric_level,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    kinases = args.kinase
    organisms = args.organism
    use_morgan = args.use_morgan
    path_klifs = Path(args.path_klifs_download)
    path_ligands = args.path_ligands
    path_output = args.output
    threshold_post_filtering = args.post_filtering

    klifs_ligands = []
    ignore_s = []

    for organism in organisms:
        if kinases == '*':
            kinase_dirs = [Path(f.path) for f in os.scandir(path_klifs / organism) if f.is_dir()]
        else:
            kinase_dirs = [(path_klifs / organism / kinase) for kinase in kinases]

        for kinase_dir in kinase_dirs:
            for pdb in (Path(f.path) for f in os.scandir(kinase_dir) if f.is_dir() and os.path.isfile(Path(f.path) / "ligand.mol2")):
                # structure with valid ligand
                ligand = Chem.MolFromMol2File(str(pdb / "ligand.mol2"), removeHs=False)
                if ligand != None:
                    klifs_ligands.append({'id': pdb.stem, 'ROMol': ligand, 'kinase': kinase_dir.stem})
                else: 
                    ignore_s.append(pdb.stem)

    with open('e_pars.txt', 'w') as f:
        f.write(str(ignore_s))

    logging.info(f"Loaded {len(klifs_ligands)} KLIFS ligands")

    klifs_data = pd.DataFrame(klifs_ligands, columns=['id','ROMol', 'kinase'])

    # standardize
    klifs_data["ROMol_standardized"] = klifs_data["ROMol"].apply(standardize_mol)

    # remove ligands with none
    r = klifs_data[klifs_data["ROMol_standardized"].isna()]['id']
    with open('e_stand.txt', 'w') as f:
        f.write(str(r))

    klifs_data = klifs_data[klifs_data["ROMol_standardized"].notna()].reset_index()

    logging.info(f"Standardized {klifs_data.shape[0]} klifs ligands")

    # generate fingerprints
    if use_morgan:
        morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
        klifs_data["fingerprint"] = klifs_data["ROMol_standardized"].map(
            lambda x: morgan_gen.GetFingerprint(x)
        )
    else:
        rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
        klifs_data["fingerprint"] = klifs_data["ROMol_standardized"].map(
            lambda x: rdkit_gen.GetFingerprint(x)
        )

    logging.info(f"Calculated {'Morgan' if use_morgan else 'RDKit'} fingerprints of {klifs_data.shape[0]} klifs ligands")

    # read results data
    data = read_mols(path_ligands)

    logging.info(f"Loaded {data.shape[0]} target compounds")

    if threshold_post_filtering > 0:
        num_ligands = data.shape[0]
        # post filtering
        data = data[data["binding_affinity"] <= threshold_post_filtering].reset_index(drop=True)

        logging.info(
            f"Removed {num_ligands - data.shape[0]} with estimated binding affinity > {threshold_post_filtering} nM (post filtering)."
        )

    logging.info(
        f"Start calculating the most similar klifs ligands for {data.shape[0]} compounds ---->"
    )

    # calculated most similar kinodata ligand
    most_similar_klifs_ligands =    [
            most_similar_database_ligand(ligand_inchi, klifs_data, ['id', 'kinase'])
                for ligand_inchi in data.inchi
        ]

    logging.info(
        f"Calculated the most similar chembl ligands for {data.shape[0]} compounds"
    )

    # add id and similarity of most similar chembl ligands
    data["most_similar_klifs_ligand.id"] = [
        res[0][0] for res in most_similar_klifs_ligands
    ]
    data["most_similar_klifs_ligand.kinase"] = [
        res[0][1] for res in most_similar_klifs_ligands
    ]
    data["most_similar_klifs_ligand.similarity"] = [
        res[1] for res in most_similar_klifs_ligands
    ]

    # save to csv
    data.to_csv(path_output)
    logging.info(f"Saved information of the {data.shape[0]} compounds to {path_output}")
