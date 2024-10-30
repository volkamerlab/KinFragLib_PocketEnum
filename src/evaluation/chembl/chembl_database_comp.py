import argparse
import logging
import sys
from concurrent.futures import ProcessPoolExecutor

from kinfraglib.utils import standardize_mol
from rdkit.Chem import PandasTools, rdFingerprintGenerator

from src.evaluation.chembl.utils import most_similar_chembl_ligand
from src.evaluation.utils import read_mols

if __name__ == "__main__":
    # parse command line arguments
    parser = argparse.ArgumentParser(
        prog=sys.argv[0], description="Determines the most similar chembl compounds for each target compounds."
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
        action="store_true",
        help="If flag is set, compounds with estimated bidning affnity > 1,000 nM will be removed from given dataset as a post filtering step.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Filename of output CSV file.",
    )
    parser.add_argument(
        "-c",
        "--path_chembl",
        help="Path to chembl SDF file.",
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

    use_morgan = args.use_morgan
    path_chembl = args.path_chembl
    path_ligands = args.path_ligands
    path_output = args.output
    post_filtering = args.post_filtering

    # read chembl data
    chembl_data = PandasTools.LoadSDF(path_chembl, embedProps=True)

    logging.info(f"Loaded {chembl_data.shape[0]} chembl ligands")

    # # standardize
    # with ProcessPoolExecutor() as executor:
    #     chemb_stand_mols = executor.map(standardize_mol, chembl_data["ROMol"])

    chembl_data["ROMol_standardized"] = chembl_data["ROMol"].apply(standardize_mol)

    logging.info(f"Standardized {chembl_data.shape[0]} chembl ligands")

    # generate fingerprints
    if use_morgan:
        morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
        chembl_data["fingerprint"] = chembl_data["ROMol_standardized"].map(
            lambda x: morgan_gen.GetFingerprint(x)
        )
    else:
        rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
        chembl_data["fingerprint"] = chembl_data["ROMol_standardized"].map(
            lambda x: rdkit_gen.GetFingerprint(x)
        )

    logging.info(f"Calculated {'Morgan' if use_morgan else 'RDKit'} fingerprints of {chembl_data.shape[0]} chembl ligands")

    # read results data
    data = read_mols(path_ligands)

    logging.info(f"Loaded {data.shape[0]} target compounds")

    if post_filtering:
        num_ligands = data.shape[0]
        # post filtering
        data = data[data["binding_affinity"] <= 1000].reset_index(drop=True)

        logging.info(
            f"Removed {num_ligands - data.shape[0]} with estimated binding affinity > 1,000 nM (post filtering)."
        )

    logging.info(
        f"Start calculating the most similar chembl ligands for {data.shape[0]} compounds using {1000000} threads ---->"
    )

    # calculated most similar kinodata ligand
    # with ProcessPoolExecutor() as executor:
    #     most_similar_chembl_ligands = executor.map(lambda ligand_inchi: most_similar_chembl_ligand(ligand_inchi, chembl_data), data['inchi'])
    most_similar_chembl_ligands =    [
            most_similar_chembl_ligand(ligand_inchi, chembl_data)
                for ligand_inchi in data.inchi
        ]

    logging.info(
        f"Calculated the most similar chembl ligands for {data.shape[0]} compounds"
    )

    # add id and similarity of most similar chembl ligands
    data["most_similar_chembl_ligand.chembl_id"] = [
        res[0] for res in most_similar_chembl_ligands
    ]
    data["most_similar_chembl_ligand.similarity"] = [
        res[1] for res in most_similar_chembl_ligands
    ]

    # save to csv
    data.to_csv(path_output)
    logging.info(f"Saved information of the {data.shape[0]} compounds to {path_output}")
