#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
from typing import Any, Dict, List
import numpy as np

# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
    :Parameters:
        path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
    Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(
        description=__doc__, usage="{0} -h".format(sys.argv[0])
    )
    parser.add_argument(
        "-i",
        "-amplicon_file",
        dest="amplicon_file",
        type=isfile,
        required=True,
        help="Amplicon is a compressed fasta file (.fasta.gz)",
    )
    parser.add_argument(
        "-s",
        "-minseqlen",
        dest="minseqlen",
        type=int,
        default=400,
        help="Minimum sequence length for dereplication",
    )
    parser.add_argument(
        "-m",
        "-mincount",
        dest="mincount",
        type=int,
        default=10,
        help="Minimum count for dereplication",
    )
    parser.add_argument(
        "-c",
        "-chunk_size",
        dest="chunk_size",
        type=int,
        default=100,
        help="Chunk size for dereplication",
    )
    parser.add_argument(
        "-k",
        "-kmer_size",
        dest="kmer_size",
        type=int,
        default=8,
        help="kmer size for dereplication",
    )
    parser.add_argument(
        "-o",
        "-output_file",
        dest="output_file",
        type=str,
        default="OTU.fasta",
        help="Output file",
    )
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    """Reads fasta file."""
    sequence = None
    with gzip.open(amplicon_file, "rt") as file:
        for line in file:
            if line[0] != ">":
                sequence += line[:-1]
            else:
                if isinstance(sequence, str):
                    if len(sequence) >= minseqlen:
                        yield sequence
                sequence = ""
    yield sequence


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    output = {}
    for sequence in read_fasta(amplicon_file, minseqlen):
        if sequence not in output.keys():
            output[sequence] = 1
        else:
            output[sequence] += 1

    output = {k: v for k, v in sorted(output.items(), key=lambda item: item[1])}
    # print(len(output.keys()))
    for key, value in list(output.items())[::-1]:
        if value >= mincount:
            yield [key, value]


def get_chunks(sequence, chunk_size):
    if len(sequence) // chunk_size < 4:
        raise ValueError
    else:
        return [
            sequence[i : (i // chunk_size + 1) * chunk_size]
            for i in range(0, len(sequence) - len(sequence) % chunk_size, chunk_size)
        ]


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))


def cut_kmer(sequence, kmer_size):
    for start in range(len(sequence)):
        kmer = sequence[start : start + kmer_size]
        if len(kmer) == kmer_size:
            yield kmer
        else:
            break


def get_identity(alignment_list):
    return (
        sum(
            [
                1
                for i in range(len(alignment_list[0]))
                if alignment_list[0][i] == alignment_list[1][i] != "-"
            ]
        )
        / len(alignment_list[0])
    ) * 100


def compute_id_matrix(chunk_chim, parents):
    """TODO"""
    identity_percentage_matrix = np.zeros((len(chunk_chim), len(parents)))
    for chk_index, chunk in enumerate(chunk_chim):
        for i, parent in enumerate(parents):
            alignment_list = nw.global_align(
                chunk,
                parent["chunks"][chk_index],
                gap_open=-1,
                gap_extend=-1,
                matrix=os.path.abspath(
                    os.path.join(os.path.dirname(__file__), "MATCH")
                ),
            )
            identity_percentage_matrix[chk_index, i] = round(
                get_identity(alignment_list), 2
            )
    return identity_percentage_matrix


def detect_chimera(identity_percentage_matrix):
    """TODO"""
    std_mean = float(np.mean(np.std(identity_percentage_matrix, axis=1)))
    # print(std_mean)
    if std_mean > 5:
        for p_id in range(identity_percentage_matrix.shape[1]):
            sim_to_p = identity_percentage_matrix[:, p_id]
            if len(sim_to_p[sim_to_p != sim_to_p[0]]) > 1:
                return True

    return False


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    references: List[Dict[str, Any]] = []

    for s, c in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        # TODO: 1) Vous diviserez chaque séquence candidate en 4 segments
        # de longueur L=chunk_size.
        chunks = get_chunks(sequence=s, chunk_size=chunk_size)

        if len(references) < 2:
            references.append({"seq": s, "count": c, "chunks": chunks})
            yield (references[-1]["seq"], references[-1]["count"])
        else:
            # TODO: 2) Pour chaque segment, vous identifierez 8 séquences cibles
            # au maximum présentant le nombre le plus grand de k-mer similaires
            # avec notre séquence candidate (au minimum 1 kmer avec notre
            # séquence candidate)
            target_sequencies_dict: Dict[Any, Any] = {}
            for chk_index, chunk in enumerate(chunks):
                target_sequencies_dict[chk_index] = {}
                kmers = [kmer for kmer in cut_kmer(chunk, kmer_size)]
                # print(kmers)
                # For each reference sequence
                for i, target in enumerate(references):
                    target_kmers = [kmer for kmer in cut_kmer(target["chunks"][chk_index], kmer_size)]
                    # print(i, target_kmers)
                    nb_common_kmers = len(common(kmers, target_kmers))
                    if nb_common_kmers > 1:
                        target_sequencies_dict[chk_index][i] = nb_common_kmers

                target_sequencies_dict[chk_index] = {
                    k: v
                    for i, (k, v) in enumerate(
                        sorted(
                            target_sequencies_dict[chk_index].items(),
                            key=lambda item: item[1],
                            reverse=True,
                        )
                    )
                    if i < 8
                }
            # print(target_sequencies_dict)
            # break
            # TODO: 3) Parmi ces 8 séquences cibles identifiées pour chaque
            # segment, vous chercherez si 2 séquences cibles sont présentes
            # dans chacune de ces listes (séquences cibles). Nous appellerons
            # ces séquences: séqences parentes.
            parents = [
                references[index]
                for index in list(
                    set.intersection(
                        *[
                            set(target_sequencies_dict[chk_index].keys())
                            for chk_index in target_sequencies_dict.keys()
                        ]
                    )
                )[:2]
            ]
            if len(parents) == 2:
                # TODO: 4) Vous calculerez le pourcentage d'identité entre les
                # segments des séquences parentes et de la séquence candidate.
                identity_percentage_matrix = compute_id_matrix(chunks, parents)

                # TODO: 5) Si l'écart type moyen des pourcentages d'identité est
                # supérieur à 5 et que 2 segments minimum de notre séquence
                # montrent une similarité différente à un des deux parents, nous
                # identifierons cette séquence comme chimérique.
                if detect_chimera(identity_percentage_matrix):
                    # print("valid check")
                    references.append({"seq": s, "count": c, "chunks": chunks})
                    yield (references[-1]["seq"], references[-1]["count"])
            # else:
            #     references.append({"seq": s, "count": c, "chunks": chunks})
            #     yield (references[-1]["seq"], references[-1]["count"])

def abundance_greedy_clustering(
    amplicon_file, minseqlen, mincount, chunk_size, kmer_size
):
    output = []
    not_chimeric = [
        seq
        for seq in chimera_removal(
            amplicon_file, minseqlen, mincount, chunk_size, kmer_size
        )
    ]

    for index, sequence in enumerate(not_chimeric):
        # Get sequence with higher abundance
        abund_sequences = not_chimeric[:index]
        if len(abund_sequences) < 1:
            output += [sequence]
        else:
            valid = True
            for abund_seq in abund_sequences:
                alignment_list = nw.global_align(
                    sequence[0],
                    abund_seq[0],
                    matrix=os.path.abspath(
                        os.path.join(os.path.dirname(__file__), "MATCH")
                    ),
                )
                similarity = get_identity(alignment_list)

                if similarity > 97:
                    valid = False
                    break
            if valid:
                output += [sequence]

    return output


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i : i + width] for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
    with open(output_file, "w") as file:
        for i, (s, c) in enumerate(OTU_list):
            file.write(f">OTU_{i+1} occurence:{c}\n")
            file.write(f"{fill(s)}\n")


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    otu_list = abundance_greedy_clustering(
        args.amplicon_file,
        args.minseqlen,
        args.mincount,
        args.chunk_size,
        args.kmer_size,
    )

    write_OTU(otu_list, args.output_file)


if __name__ == "__main__":
    main()