"""
HMM database loading and optimized profile management.

Loads HMM files, auto-detects alphabet from the ALPH field,
builds OptimizedProfiles once at startup, and provides the dict-based
access needed for ad-hoc model subset searches in pass 2.
"""

import pyhmmer.easel as easel
import pyhmmer.plan7 as plan7


AMINO_ALPHABET = easel.Alphabet.amino()
DNA_ALPHABET = easel.Alphabet.dna()
RNA_ALPHABET = easel.Alphabet.rna()

_ALPHABET_MAP = {
    "amino": AMINO_ALPHABET,
    "DNA":   DNA_ALPHABET,
    "RNA":   RNA_ALPHABET,
}


def peek_alphabet(hmm_path):
    """
    Detect the alphabet of an HMM database by reading the first ALPH field.

    Args:
        hmm_path: path to an HMM file

    Returns:
        easel.Alphabet instance (amino, DNA, or RNA)

    Raises:
        ValueError: if no ALPH line found or unrecognized alphabet
    """
    with open(hmm_path) as f:
        for line in f:
            if line.startswith("ALPH"):
                alph_str = line.split()[1].strip()
                if alph_str in _ALPHABET_MAP:
                    return _ALPHABET_MAP[alph_str]
                raise ValueError(
                    f"Unrecognized HMM alphabet '{alph_str}' in {hmm_path}"
                )
    raise ValueError(f"No ALPH line found in {hmm_path}")


def needs_translation(hmm_path):
    """
    Check whether a database requires amino acid input (i.e. translation).

    Returns True for amino acid HMMs, False for DNA/RNA.
    """
    alphabet = peek_alphabet(hmm_path)
    return alphabet == AMINO_ALPHABET


def load_hmms(hmm_path):
    """
    Load all HMM models from a file.

    Args:
        hmm_path: path to an HMM database file (.hmm)

    Returns:
        list of plan7.HMM objects
    """
    hmms = []
    with plan7.HMMFile(hmm_path) as hf:
        for hmm in hf:
            hmms.append(hmm)
    return hmms


def build_optimized_profiles(hmms, alphabet=None):
    """
    Build OptimizedProfiles from a list of HMMs.

    This is the expensive step that should happen once per worker at startup.
    Each OptimizedProfile has the SSV/MSV filter data pre-computed.

    Args:
        hmms: list of plan7.HMM objects
        alphabet: easel.Alphabet (defaults to amino)

    Returns:
        dict of {model_name: OptimizedProfile}
    """
    if alphabet is None:
        alphabet = AMINO_ALPHABET

    bg = plan7.Background(alphabet)
    profiles = {}

    for hmm in hmms:
        profile = plan7.Profile(hmm.M, alphabet)
        profile.configure(hmm, bg)
        opt = profile.to_optimized()
        profiles[hmm.name] = opt

    return profiles


def load_and_optimize(hmm_path, alphabet=None):
    """
    Load HMMs from file and build optimized profiles in one step.

    Convenience function for worker initialization. Auto-detects
    alphabet if not provided.

    Args:
        hmm_path: path to an HMM database file
        alphabet: easel.Alphabet (auto-detected if None)

    Returns:
        tuple of (list of HMMs, dict of {model_name: OptimizedProfile},
                  easel.Alphabet)
    """
    if alphabet is None:
        alphabet = peek_alphabet(hmm_path)
    hmms = load_hmms(hmm_path)
    profiles = build_optimized_profiles(hmms, alphabet)
    return hmms, profiles, alphabet
