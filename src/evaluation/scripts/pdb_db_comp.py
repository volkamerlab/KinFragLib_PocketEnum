import logging
import warnings
import logging

import pandas as pd
from src.evaluation.utils import read_mols
from src.evaluation.scripts.utils import get_data_for_ids, generate_and_save
import requests_cache
import biotite.database.rcsb as rcsb


# Select which subsets of cocrystallized ligands should be analysed 
HAMSTER_PKA = True # hamster PKA 
PKA = True # PKA ligands (from arbitrary organism)
KINASES = True # all kinase ligands

# If set to true, the analysis is performed only on the ligands, 
# selected for synthesis 
# (on both, the generated and the modified versions of the ligands)
SYNETHSIZED_ONLY = False


if __name__ == "__main__":
    # Disable some unneeded warnings
    logger = logging.getLogger("opencadd")
    logger.setLevel(logging.ERROR)
    warnings.filterwarnings("ignore")

    # Cache requests -- this will speed up repeated queries to PDB
    requests_cache.install_cache("rcsb_pdb", backend="memory")

    uniprot_id = "P25321"  # hamster PKA
    experimental_method = "X-RAY DIFFRACTION"
    min_ligand_molecular_weight = 100.0

    query_by_experimental_method = rcsb.FieldQuery(
        "exptl.method", exact_match=experimental_method
    )
    query_by_ligand_mw = rcsb.FieldQuery(
        "chem_comp.formula_weight",
        molecular_definition=True,
        greater=min_ligand_molecular_weight,
    )
    query_by_non_polymer = rcsb.FieldQuery(
        "rcsb_entry_info.nonpolymer_entity_count", greater=0
    )
    query_by_title = rcsb.FieldQuery("struct.title", contains_phrase="kinase")

    query_by_uniprot_id = rcsb.FieldQuery(
        "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        exact_match=uniprot_id,
    )

    query_by_pka_title = rcsb.FieldQuery("struct.title", contains_phrase="PKA")

    query_by_prot_a_title = rcsb.FieldQuery(
        "struct.title", contains_phrase="Protein kinase A"
    )

    query_by_camp_title = rcsb.FieldQuery(
        "struct.title", contains_phrase="cAMP-dependent protein kinase"
    )

    # === load data ===

    if SYNETHSIZED_ONLY:
        data_proposed = read_mols("../results_5n1f_25_02/5n1f/proposed_ligands.sdf", docking_score=False)
        data_modified = read_mols("../results_5n1f_25_02/5n1f/adapted_mols.sdf", docking_score=False)
        data_proposed["source"] = "proposed"
        data_modified["source"] = "modified"
        data = pd.concat([data_proposed, data_modified], ignore_index=True)
    else:
        data = read_mols("../results_5n1f_25_02/5n1f/results.sdf")

        # post filtering
        data = data[data["binding_affinity"] <= 1000].reset_index(drop=True)

    # ===== hamster ======
    if HAMSTER_PKA:
        print("==== Hamster =====")

        hamster_query = rcsb.CompositeQuery(
            [
                query_by_uniprot_id,
                query_by_experimental_method,
                query_by_ligand_mw,
                query_by_non_polymer,
            ],
            "and",
        )
        hamster_pdb_ids = rcsb.search(hamster_query)

        print(f"Number of matches: {len(hamster_pdb_ids)}")

        hamster_ligands = get_data_for_ids(hamster_pdb_ids)

        print(f"Number of structures: {len(hamster_ligands)}")

        generate_and_save(hamster_ligands, data, "hamster_pka_compare_synth.csv")

    if PKA:
        print("====== PKA =======")

        pkd_title_query = rcsb.CompositeQuery(
            [query_by_camp_title, query_by_pka_title, query_by_prot_a_title],
            "or",
        )

        pka_query = rcsb.CompositeQuery(
            [
                pkd_title_query,
                query_by_experimental_method,
                query_by_ligand_mw,
                query_by_non_polymer,
            ],
            "and",
        )
        pka_pdb_ids = rcsb.search(pka_query)

        print(f"Number of matches: {len(pka_pdb_ids)}")

        pka_ligands = get_data_for_ids(pka_pdb_ids)

        print(f"Number of structures: {len(pka_ligands)}")

        generate_and_save(pka_ligands, data, "pka_compare_synth.csv")

    if KINASES:
        print("==== Kinase =====")

        kinase_query = rcsb.CompositeQuery(
            [
                query_by_title,
                query_by_experimental_method,
                query_by_ligand_mw,
                query_by_non_polymer,
            ],
            "and",
        )
        kinase_pdb_ids = rcsb.search(kinase_query)

        print(f"Number of matches: {len(kinase_pdb_ids)}")

        kinase_ligands = get_data_for_ids(kinase_pdb_ids)

        print(f"Number of structures: {len(kinase_ligands)}")

        generate_and_save(kinase_ligands, data, "kinase_compare_synth.csv")
