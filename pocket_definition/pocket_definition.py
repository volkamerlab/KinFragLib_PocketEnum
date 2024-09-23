import logging

logging.getLogger("numexpr").setLevel(logging.ERROR)

from bravado.client import SwaggerClient
from Bio.Blast.Applications import NcbiblastpCommandline
from opencadd.structure.core import Structure
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess
import argparse
import os, sys

from mmcif.io.BinaryCifReader import BinaryCifReader

if __name__ == "__main__":

    # mapping of three letter code to one letter code of amino acids
    AMINO_ACIDS = {
        "CYS": "C",
        "ASP": "D",
        "SER": "S",
        "GLN": "Q",
        "LYS": "K",
        "ILE": "I",
        "PRO": "P",
        "THR": "T",
        "PHE": "F",
        "ASN": "N",
        "GLY": "G",
        "HIS": "H",
        "LEU": "L",
        "ARG": "R",
        "TRP": "W",
        "ALA": "A",
        "VAL": "V",
        "GLU": "E",
        "TYR": "Y",
        "MET": "M",
        "SEP": "?",
        "TPO": "?",
    }

    # parse command line arguments
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        description="Determines KLIFS binding pocket residues and returns the mapping of residue id to KLIFS numbering scheme",
    )
    parser.add_argument(
        "-f", "--family", help="Name of the kinase family e.g. PKA", required=True
    )
    parser.add_argument(
        "-n", "--name", help="Name of the kinase e.g. Prkaca", required=True
    )
    parser.add_argument(
        "-p",
        "--pdb_structure",
        help="PDB ID of the structure of interest e.g. 3AMB",
        required=True,
    )

    args = parser.parse_args()

    kinase_family = args.family
    kinase_name = args.name
    pdb_code = args.pdb_structure

    # get all KLIFS entries for the kinase of interest

    KLIFS_API_DEFINITIONS = "https://klifs.net/swagger/swagger.json"
    KLIFS_CLIENT = SwaggerClient.from_url(
        KLIFS_API_DEFINITIONS, config={"validate_responses": False}
    )

    kinases = (
        KLIFS_CLIENT.Information.get_kinase_names(kinase_family=kinase_family)
        .response()
        .result
    )
    kinase_klifs_ids = [
        kinase.kinase_ID for kinase in kinases if kinase.name == kinase_name
    ]
    structures = (
        KLIFS_CLIENT.Structures.get_structures_list(kinase_ID=kinase_klifs_ids)
        .response()
        .result
    )

    # write all KLIFS kinase structures to a file for the MSA calculation
    with open("pocket_sequences.fasta", "wt") as fasta_file:
        for i, structure in enumerate(structures):
            fasta_file.write(
                ">"
                + structure.pdb
                + str(structure.structure_ID)
                + str(i)
                + "\n"
                + structure.pocket
                + "\n"
            )

    # perform MSA with clustalo
    output_text = subprocess.run(
        [
            "clustalo",
            "-i",
            "pocket_sequences.fasta",
            "-o",
            "alignment.clustal",
            "--outfmt=clu",
            "--dealign",
            "--force",
        ],
        capture_output=True,
    )
    os.remove("pocket_sequences.fasta")

    # calc simple consensus sequence with cons EMBOSS
    try:
        output_text = subprocess.run(
            ["cons", "-sequence", "alignment.clustal", "-outseq", "consensus.cons"],
            capture_output=True,
        )
    except FileNotFoundError:
        output_text = subprocess.run(
            ["em_cons", "-sequence", "alignment.clustal", "-outseq", "consensus.cons"],
            capture_output=True,
        )
    os.remove("alignment.clustal")

    # parse sequence
    with open("consensus.cons", "rt") as cons_file:
        _temp = [line.strip() for line in cons_file if not line.startswith(">")]
    consensus_sequence = "".join(_temp)

    os.remove("consensus.cons")

    # prepare structure of interest

    # read pdbx from RCSB
    bcr = BinaryCifReader(storeStringsAsBytes=False)
    bcr_results = bcr.deserialize(f"https://models.rcsb.org/{pdb_code}.bcif.gz")

    if len(bcr_results) > 1:
        logging.warning(
            f"Obtained {len(bcr_results)} for {pdb_code} on RCSB, selecting the first result"
        )
    elif not len(bcr_results):
        logging.error(f"No result was found for {pdb_code} on RCSB")
        exit(1)

    # select first result by default (normally it should be len(bcr_results) == 1)
    protein_pdbx = bcr_results[0]

    # get ids of all polymer entities
    sequences = {pol[0]: [] for pol in protein_pdbx.getObj("entity_poly")}
    residue_ids = {pol[0]: [] for pol in protein_pdbx.getObj("entity_poly")}

    # retrieve sequence of all polymer entities
    for residue in protein_pdbx.getObj("entity_poly_seq"):
        polymer_id = residue[0]
        residue_name = residue[2]
        residue_id = residue[1]
        sequences[polymer_id].append(AMINO_ACIDS.get(residue_name) or "?")
        residue_ids[polymer_id].append(residue_id)

    logging.info(
        f"Found {len(sequences)} polymer entities, proceeding on every entity separatly"
    )

    # a pdbx file can comprise several polymer entities
    #   -> perform mapping on each entity with a seq_length >= number of binding pocket residues

    for sequence_id, sequence in sequences.items():
        # skip sequences if  seq_length < number of binding pocket residues
        if len(sequence) < len(consensus_sequence):
            logging.info(
                f"Skipping entity {sequence_id} since number of residues ({len(sequence)}) < {len(consensus_sequence)} (# KLIFS binding pocket residues))"
            )
            continue

        # create two sequence files
        seq1 = SeqRecord(Seq(consensus_sequence), id="seq1")
        seq2 = SeqRecord(Seq("".join(sequence)), id="seq2")
        SeqIO.write(seq1, "seq1.fasta", "fasta")
        SeqIO.write(seq2, "seq2.fasta", "fasta")

        # run BLAST and parse the output as XML
        output = NcbiblastpCommandline(
            query="seq1.fasta", subject="seq2.fasta", outfmt=5
        )()[0]
        blast_result_record = NCBIXML.read(StringIO(output))

        os.remove("seq1.fasta")
        os.remove("seq2.fasta")

        # print some information on the result
        hsp = next(
            hsp
            for aligment in blast_result_record.alignments
            for hsp in aligment.hsps
            if hsp.query_start == 1 and hsp.query_end == 85
        )

        # parse alignment
        hit_from = hsp.sbjct_start
        hit_to = hsp.sbjct_end
        query = hsp.query

        hit = hsp.sbjct

        # get KLIFS residue numbering mapping
        mapping = {}
        assert len(query) == len(residue_ids[sequence_id][hit_from - 1 : hit_to])
        idx = 1
        for q, id in zip(query, residue_ids[sequence_id][hit_from - 1 : hit_to]):
            if q != "-":
                mapping[idx] = id
                idx += 1

        # print full mapping
        print("Mapping (KLIFS numbering: residue ID): ")

        for klifs_idx in mapping:
            print(f"\t{klifs_idx}: {mapping[klifs_idx]}")

        # print subpocket and hinge residues
        print(f"\nHinge (46, 48):\t\t[{mapping[46]}, {mapping[48]}]")
        print(
            f"AP (15, 46, 51, 75):\t{[mapping[15], mapping[46], mapping[51], mapping[75]]}",
        )
        print(
            f"FP (10, 51, 72, 81):\t{[mapping[10], mapping[51], mapping[72], mapping[81]]}",
        )
        print(f"GA (17, 45, 81):\t{[mapping[17], mapping[45], mapping[81]]}")
        print(f"SE (51):\t\t{[mapping[51]]}")
        print(
            f"B1 (28, 38, 43, 81):\t{[mapping[28], mapping[38], mapping[43], mapping[81]]}",
        )
        print(
            f"B2 (18, 24, 70, 83):\t{[mapping[18], mapping[24], mapping[70], mapping[83]]}",
        )
