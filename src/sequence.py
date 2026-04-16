"""
FASTA ingestion and six-frame translation using pyfastx and pyhmmer.

Reads input nucleotide FASTA via pyfastx, performs six-frame translation
using pyhmmer's DigitalSequence.translate(), writes translated amino acid
sequences as an indexed FASTA output.
"""

import re

import pyfastx
import pyhmmer.easel as easel


DNA_ALPHABET = easel.Alphabet.dna()
AMINO_ALPHABET = easel.Alphabet.amino()

_COMP_TABLE = str.maketrans("ACGT", "TGCA")
_INVALID_DNA = re.compile(r"[^ACGTN]")

FRAME_SUFFIXES = [
    ("fwd", 0),
    ("fwd", 1),
    ("fwd", 2),
    ("rev", 0),
    ("rev", 1),
    ("rev", 2),
]


def clean_seq(seq):
    """Replace any non-ACGTN characters with N."""
    return _INVALID_DNA.sub("N", seq)


def revcomp(seq):
    """Reverse complement a DNA string."""
    return seq[::-1].translate(_COMP_TABLE)


def open_input(fasta_path):
    """
    Open a FASTA file with pyfastx indexing.

    Returns the pyfastx.Fasta object. The .fxi index is created
    automatically if it doesn't already exist.
    """
    return pyfastx.Fasta(fasta_path, build_index=True)


def get_sequence_metadata(fasta_path):
    """
    Retrieve sequence metadata from the pyfastx SQLite index.

    Returns a dict of {seq_name: seq_length} without re-parsing the FASTA.
    """
    fa = open_input(fasta_path)
    return {rec.name: len(rec) for rec in fa}


def _trim_to_codon(length, frame):
    """Return the number of bases to use so (length - frame) is a multiple of 3."""
    usable = length - frame
    return frame + (usable - usable % 3)


def six_frame_translate_seq(name, sequence):
    """
    Six-frame translate a single nucleotide sequence.

    Args:
        name: sequence identifier (str)
        sequence: nucleotide sequence (str, uppercase)

    Yields:
        (frame_name, aa_sequence_str) tuples for each of the 6 frames.
        frame_name format: "{name}|{direction}{frame+1}" e.g. "seq1|fwd1"
    """
    seq_len = len(sequence)
    rc_seq = revcomp(sequence)

    for direction, frame in FRAME_SUFFIXES:
        trim_len = _trim_to_codon(seq_len, frame)
        if trim_len - frame < 3:
            continue

        source = sequence if direction == "fwd" else rc_seq
        subseq = source[frame:trim_len]

        frame_name = f"{name}|{direction}{frame + 1}"

        ts = easel.TextSequence(
            name=frame_name.encode("ascii"), sequence=subseq
        )
        aa = ts.digitize(DNA_ALPHABET).translate()
        aa_text = aa.textize().sequence

        yield frame_name, aa_text


def load_sequences_dict(fasta_path):
    """
    Load all sequences from a FASTA file into a dict.

    For bulk subsequence extraction -- string slicing from a dict is
    nanoseconds vs ~14ms/call for pyfastx random access.

    Args:
        fasta_path: path to FASTA file

    Returns:
        dict of {seq_name: sequence_str}
    """
    fa = open_input(fasta_path)
    return {rec.name: str(rec.seq) for rec in fa}


def translate_fasta(input_fasta, output_fasta):
    """
    Six-frame translate all sequences in a FASTA file.

    Reads input via pyfastx, translates each sequence in all 6 frames,
    writes the amino acid sequences to output_fasta, and builds a
    pyfastx index on the output.

    Args:
        input_fasta: path to input nucleotide FASTA
        output_fasta: path for output amino acid FASTA

    Returns:
        dict of {original_seq_name: seq_length} from the input
    """
    fa = open_input(input_fasta)
    nucl_lengths = {}

    with open(output_fasta, "w") as fout:
        for rec in fa:
            name = rec.name
            sequence = clean_seq(str(rec.seq).upper())
            nucl_lengths[name] = len(sequence)

            for frame_name, aa_seq in six_frame_translate_seq(name, sequence):
                fout.write(f">{frame_name}\n{aa_seq}\n")

    # Build pyfastx index on the output
    pyfastx.Fasta(output_fasta, build_index=True)

    return nucl_lengths


def parse_frame_suffix(target_name):
    """
    Parse frame info from a translated sequence name.

    Name format: "{seq_id}|{fwd|rev}{1|2|3}"

    Returns:
        (seq_id, strand, frame) where strand is '+'/'-' and frame is 0-indexed.
        Returns (target_name, '.', '.') if not a translated name.
    """
    parts = target_name.rsplit("|", 1)
    if len(parts) != 2:
        return target_name, ".", "."

    seq_id, suffix = parts
    if suffix.startswith("fwd"):
        strand = "+"
    elif suffix.startswith("rev"):
        strand = "-"
    else:
        return target_name, ".", "."

    frame = int(suffix[-1]) - 1  # 1-indexed in name -> 0-indexed
    return seq_id, strand, frame


def aa_to_nucl_coords(env_from, env_to, strand, frame, nucl_length):
    """
    Convert amino acid envelope coordinates to nucleotide coordinates.

    Args:
        env_from: 1-based AA start position
        env_to: 1-based AA end position
        strand: '+' or '-'
        frame: 0-indexed reading frame (0, 1, 2)
        nucl_length: length of the original nucleotide sequence

    Returns:
        (nuc_start, nuc_end) 1-based nucleotide coordinates
    """
    if strand == "+":
        nuc_start = (env_from - 1) * 3 + frame + 1
        nuc_end = env_to * 3 + frame
    elif strand == "-":
        nuc_start = nucl_length - (env_to * 3 + frame) + 1
        nuc_end = nucl_length - ((env_from - 1) * 3 + frame)
    else:
        return env_from, env_to

    return nuc_start, nuc_end
