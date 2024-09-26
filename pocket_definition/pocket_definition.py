import argparse
import logging
import os
import subprocess
import sys
from collections import defaultdict
from io import StringIO

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from bravado.client import SwaggerClient

if __name__ == "__main__":

    # mapping of three letter code to one letter code of amino acids
    AMINO_ACIDS = defaultdict(
        lambda: "?",
        {
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
        },
    )

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
        "--pdb",
        help="Path to PDB structure of interest",
        required=True,
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

    kinase_family = args.family
    kinase_name = args.name
    path_pdb = args.pdb

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
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(kinase_name, path_pdb)

    for chain in structure.get_chains():
        # get protein residues (disregarding waters, ...)
        residues = list(filter(lambda res: res.id[0] == " ", chain.get_residues()))
        sequence = [AMINO_ACIDS[res.get_resname()] for res in residues]

        logging.info(f"Proceed chain {chain.get_id()}")
        # skip sequences if  seq_length < number of binding pocket residues
        if len(sequence) < len(consensus_sequence):
            logging.info(
                f"Skipping entity {chain.get_id()} since number of residues ({len(sequence)}) < {len(consensus_sequence)} (# KLIFS binding pocket residues))"
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

        # get result (highest alignmnent score + covering the whole pocket)
        hsp = next(
            (
                hsp
                for aligment in blast_result_record.alignments
                for hsp in aligment.hsps
                # ensure that query (pocket consensus sequence) is fully covered
            ),
            None,
        )

        if not hsp:
            logging.warning(
                "Could not find an alignment for this sequence covering the whole binding pocket"
            )
            continue

        # parse alignment
        hit_from = hsp.sbjct_start - 1  # 0-starting idx
        hit_to = hsp.sbjct_end
        query = hsp.query
        hit = hsp.sbjct

        logging.info(f"Alignment found (score: {hsp.score}):\n\t{query}\n\t{hit}")

        # get KLIFS residue numbering mapping

        # fill start gaps (not included in alignment)
        mapping = [
            residues[hit_from - hsp.query_start + i].get_id()[1]
            for i in range(1, hsp.query_start)
        ]

        for q, h in zip(query, hit):
            # we are only interested in query (binding pocket) residues
            # thus, we skip gaps ('-') in query sequence
            if q != "-":
                # determine corresponding residue id of given pdb numbering
                if h != "-":
                    # match or missmatch -> add corrsponding residue id and proceed in sequence
                    # otherwise add '?' and do not proceed
                    mapping.append(residues[hit_from].get_id()[1])
                    hit_from += 1
                else:
                    # deletion
                    mapping.append("?")
            else:
                hit_from += 1

        # fill open ends with missmathces
        mapping += [
            residues[hit_to - hsp.query_end + i].get_id()[1]
            for i in range(hsp.query_end, 85)
        ]

        # print full mapping
        print(
            "\n\033[1mMapping \033[0m\t| KLIFS numbering \t|\n\t\t| Residue ID \t\t|\n"
        )

        for i in range(6):
            for klifs_idx in range(i * 15, min((i + 1) * 15, len(mapping))):
                print(f"| {klifs_idx + 1}", end="\t")
            print("|")
            for pdb_idx in mapping[i * 15 : min((i + 1) * 15, len(mapping))]:
                print(f"| {pdb_idx}", end="\t")
            print("|\n")

        # print subpocket and hinge residues
        print("\033[1mSubpocket residues:\033[0m")
        print(f"Hinge (46, 48):\t\t[{mapping[46 + 1]}, {mapping[48 + 1]}]")
        print(
            f"AP (15, 46, 51, 75):\t{[mapping[15 + 1], mapping[46 + 1], mapping[51 + 1], mapping[75 + 1]]}",
        )
        print(
            f"FP (10, 51, 72, 81):\t{[mapping[10 + 1], mapping[51 + 1], mapping[72 + 1], mapping[81 + 1]]}",
        )
        print(
            f"GA (17, 45, 81):\t{[mapping[17 + 1], mapping[45 + 1], mapping[81 + 1]]}"
        )
        print(f"SE (51):\t\t{[mapping[51 + 1]]}")
        print(
            f"B1 (28, 38, 43, 81):\t{[mapping[28 + 1], mapping[38 + 1], mapping[43 + 1], mapping[81 + 1]]}",
        )
        print(
            f"B2 (18, 24, 70, 83):\t{[mapping[18 + 1], mapping[24 + 1], mapping[70 + 1], mapping[83 + 1]]}\n",
        )
