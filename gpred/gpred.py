"Anaïs DElASSUS"

import argparse
import sys
import os
import csv
import re
import textwrap
from re import Pattern
from pathlib import Path
from typing import List, Union


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__,
                                     usage=f"{sys.argv[0]} -h")

    parser.add_argument('-i', dest='genome_file', type=isfile, required=True,
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int,
                        default=50,
                        help="Minimum gene length to consider (default 50).")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int,
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif \
                            (default 16).")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes - shine box not \
                            included (default 40).")
    parser.add_argument('-p', dest='predicted_genes_file', type=Path,
                        default=Path("predict_genes.csv"),
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=Path,
                        default=Path("genes.fna"),
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file: Path) -> str:
    """Extract genome sequence from fasta files.

    :param fasta_file: (Path) Path to the fasta file.
    :return: (str) Sequence from the genome. 
    """
    seq_lines = []
    with fasta_file.open("r") as fasta:
        for line in fasta:
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_lines.append(line.strip())
    sequence = "".join(seq_lines).upper()
    return sequence


def find_start(start_regex: Pattern, sequence: str, start: int,
               stop: int) -> Union[int, None]:
    """Find next start codon before a end position.

    :param start_regexp: A regex object that identifies a start codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :param stop: (int) Stop position of the research
    :return: (int) If exist, position of the start codon. Otherwise None. 
    """
    m = start_regex.search(sequence, start, stop)
    if m:
        return m.start()
    return None


def find_stop(stop_regex: Pattern, sequence: str,
              start: int) -> Union[int, None]:
    """Find next stop codon that should be in the same reading phase as the 
    start.

    :param stop_regexp: A regex object that identifies a stop codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :return: (int) If exist, position of the stop codon. Otherwise None. 
    """
    for match in stop_regex.finditer(sequence, start):
        stop_idx = match.start()
        if (stop_idx - start) % 3 == 0:
            return stop_idx
    return None


def has_shine_dalgarno(shine_regex: Pattern, sequence: str, start: int,
                       max_shine_dalgarno_distance: int) -> bool:
    """Find a shine dalgarno motif before the start codon

    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Position of the start in the genome
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine 
    dalgarno to the start position
    :return: (boolean) true -> has a shine dalgarno upstream to the gene, 
    false -> no
    """
    region_end = start - 6
    region_start = start - max_shine_dalgarno_distance
    if region_start < 0:
        return False
    if shine_regex.search(sequence, region_start, region_end):
        return True
    return False


def predict_genes(sequence: str,
                  start_regex: Pattern,
                  stop_regex: Pattern,
                  shine_regex: Pattern,
                  min_gene_len: int, max_shine_dalgarno_distance: int,
                  min_gap: int) -> List:
    """Predict most probable genes

    :param sequence: (str) Sequence from the genome.
    :param start_regexp: A regex object that identifies a start codon.
    :param stop_regexp: A regex object that identifies a stop codon.
    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param min_gene_len: (int) Minimum gene length.
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine 
    dalgarno to the start position.
    :param min_gap: (int) Minimum distance between two genes.
    :return: (list) List of [start, stop] position of each predicted genes.
    """
    predicted_genes = []
    current_pos = 0
    genome_length = len(sequence)
    while current_pos <= genome_length - min_gap:
        start_codon = find_start(
            start_regex, sequence, current_pos, genome_length)
        if start_codon is None:
            break
        stop_codon = find_stop(stop_regex, sequence, start_codon)
        if stop_codon is None:
            current_pos = start_codon + 1
            continue
        length_gene = stop_codon + 3 - start_codon
        if length_gene >= min_gene_len and \
            has_shine_dalgarno(shine_regex, sequence, start_codon + 1,
                               max_shine_dalgarno_distance):
            predicted_genes.append([start_codon + 1, stop_codon + 3])
            current_pos = stop_codon + 3 + min_gap
        else:
            current_pos = start_codon + 1
    return predicted_genes


def write_genes_pos(predicted_genes_file: Path,
                    probable_genes: List[List[int]]) -> None:
    """Write list of gene positions.

    :param predicted_genes_file: (Path) Output file of gene positions.
    :param probable_genes: List of [start, stop] position of each predicted 
    genes.
    """
    with predicted_genes_file.open("w") as genes_predits:
        predict_genes_writer = csv.writer(genes_predits, delimiter=",")
        predict_genes_writer.writerow(["Start", "Stop"])
        predict_genes_writer.writerows(probable_genes)


def write_genes(fasta_file: Path, sequence: str,
                probable_genes: List[List[int]], sequence_rc: str,
                probable_genes_comp: List[List[int]]):
    """Write gene sequence in fasta format

    :param fasta_file: (Path) Output fasta file.
    :param sequence: (str) Sequence of genome file in 5'->3'.
    :param probable_genes: (list) List of [start, stop] position of each 
    predicted genes in 5'->3'.
    :param sequence_rc: (str) Sequence of genome file in 3' -> 5'.
    :param probable_genes_comp: (list)List of [start, stop] position of each 
    predicted genes in 3' -> 5'.
    """
    with open(fasta_file, "w", encoding="utf-8") as fasta:
        for idx, gene_pos in enumerate(probable_genes, start=1):
            fasta.write(
                f">gene_{idx}{os.linesep}"
                f"{textwrap.fill(sequence[gene_pos[0]-1:gene_pos[1]])}"
                f"{os.linesep}"
            )
        for idx, gene_pos in enumerate(probable_genes_comp, start=len(probable_genes)+1):
            fasta.write(
                f">gene_{idx}{os.linesep}"
                f"{textwrap.fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])}"
                f"{os.linesep}"
            )


def reverse_complement(sequence: str) -> str:
    """Get the reverse complement

    :param sequence: (str) DNA Sequence.
    :return: (str) Reverse complemented sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in sequence[::-1]])


# ==============================================================
# Main program
# ==============================================================
def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    args = get_arguments()

    # Read genome sequence
    sequence = read_fasta(args.genome_file)
    seq_len = len(sequence)
    if seq_len == 0:
        sys.exit(f"Error: genome sequence in {args.genome_file} is empty.")

    # Compile regexes
    start_regex = re.compile(r'AT[TG]|[ATCG]TG')
    stop_regex = re.compile(r'TA[GA]|TGA')
    shine_regex = re.compile(r'A?G?GAGG|GGAG|GG.{1}GG')

    # Predict genes on forward strand (5' -> 3')
    probable_genes = predict_genes(
        sequence=sequence,
        start_regex=start_regex,
        stop_regex=stop_regex,
        shine_regex=shine_regex,
        min_gene_len=args.min_gene_len,
        max_shine_dalgarno_distance=args.max_shine_dalgarno_distance,
        min_gap=args.min_gap
    )

    # Predict genes on reverse-complement strand (3' -> 5')
    sequence_rc = reverse_complement(sequence)
    probable_genes_rc = predict_genes(
        sequence=sequence_rc,
        start_regex=start_regex,
        stop_regex=stop_regex,
        shine_regex=shine_regex,
        min_gene_len=args.min_gene_len,
        max_shine_dalgarno_distance=args.max_shine_dalgarno_distance,
        min_gap=args.min_gap
    )

    # Convert positions from reverse-complement coordinates to forward (5'->3')
    probable_genes_comp: List[List[int]] = []
    for start_rc, end_rc in probable_genes_rc:
        start_fwd = seq_len - end_rc + 1
        end_fwd = seq_len - start_rc + 1
        # ensure start <= end (sanity) and clamp to genome bounds
        start_fwd = max(1, min(start_fwd, seq_len))
        end_fwd = max(1, min(end_fwd, seq_len))
        if start_fwd <= end_fwd:
            probable_genes_comp.append([start_fwd, end_fwd])

    # Combine forward + converted reverse lists and sort by increasing start
    # position
    combined_sorted = sorted(
        probable_genes + probable_genes_comp, key=lambda x: x[0])

    # Write predicted gene positions (CSV) and FASTA sequences
    write_genes_pos(args.predicted_genes_file, combined_sorted)
    write_genes(args.fasta_file, sequence, probable_genes,
                sequence_rc, probable_genes_comp)

    # Print summary
    print(f"Genome length: {seq_len} bp")
    print(f"Predicted genes (forward strand): {len(probable_genes)}")
    print(f"Predicted genes (reverse strand, converted): "
          f"{len(probable_genes_comp)}")
    print(f"Total predicted genes (combined): {len(combined_sorted)}")
    print(f"Positions written to: {args.predicted_genes_file}")
    print(f"FASTA written to: {args.fasta_file}")


if __name__ == '__main__':
    main()
