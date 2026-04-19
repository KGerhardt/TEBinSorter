"""
tesorter_output.py — Generate TEsorter-compatible output files.

Produces the same output files as TEsorter from our SQLite database:
  - {prefix}.cls.tsv  — TE classifications (already produced by classifier)
  - {prefix}.dom.gff3 — Domain annotations in GFF3 format
  - {prefix}.dom.tsv  — Domain hit summary table
  - {prefix}.dom.faa  — Domain protein sequences in FASTA
  - {prefix}.cls.lib  — RepeatMasker library (classified sequences)
  - {prefix}.cls.pep  — Classified protein domain sequences
"""

import logging
import re

import pyfastx

from sequence import (parse_frame_suffix, aa_to_nucl_coords,
                      load_sequences_dict, open_input, revcomp)
from classifier import parse_clade_rexdb, parse_clade_gydb

log = logging.getLogger(__name__)


def _parse_gene_clade(model_name, db_name):
    """Extract gene and clade from model name."""
    if db_name.startswith("rexdb") or db_name in ("line", "tir"):
        if ":" in model_name:
            gene = model_name.split(":")[1]
            clade = model_name.split("/")[-1].split(":")[0]
        else:
            gene = model_name
            clade = model_name
    elif db_name == "gydb":
        parts = model_name.split("_", 1)
        gene = parts[0]
        clade = parts[1] if len(parts) > 1 else model_name
    elif db_name.startswith("sine"):
        gene = "SINE"
        clade = "SINE"
    else:
        gene = model_name
        clade = model_name
    return gene, clade


def _format_gff_id(seq_id):
    """Make ID safe for GFF3 (replace special chars)."""
    return re.sub(r"[;=|]", "_", seq_id)


def generate_dom_gff3(conn, prefix, db_name, nucl_lengths, seq_type="nucl"):
    """Generate domain annotation GFF3 file.

    Format: chr  TEsorter  CDS  start  end  score  strand  frame  attributes
    """
    out_path = f"{prefix}.dom.gff3"

    rows = conn.execute("""
        SELECT target_name, query_name, dom_score, i_evalue, acc,
               hmm_from, hmm_to, query_len, env_from, env_to,
               base_seq, strand, frame, domain_type
        FROM legacy_hits
        WHERE database = ?
        ORDER BY base_seq, env_from
    """, (db_name,)).fetchall()

    # Apply TEsorter filters
    with open(out_path, "w") as f:
        for row in rows:
            (target, model, score, evalue, acc,
             hmm_from, hmm_to, model_len, env_from, env_to,
             base_seq, strand, frame_num, domain_type) = row

            hmm_cov = 100.0 * (hmm_to - hmm_from + 1) / model_len
            norm_score = round(score / model_len, 2)

            if not (hmm_cov >= 20 and evalue <= 1e-3 and acc >= 0.5
                    and norm_score >= 0.1):
                continue

            gene, clade = _parse_gene_clade(model, db_name)

            # Convert to nucleotide coordinates
            if seq_type == "nucl" and strand in ("+", "-"):
                nuc_start, nuc_end = aa_to_nucl_coords(
                    env_from, env_to, strand, int(frame_num),
                    nucl_lengths.get(base_seq, 0))
                nuc_frame = int(frame_num)
            else:
                nuc_start, nuc_end = env_from, env_to
                nuc_frame = "."

            gid = f"{_format_gff_id(base_seq)}|{model}"
            name = f"{clade}-{gene.split('-')[-1] if '-' in gene else gene}"
            attr = (f"ID={gid};Name={name};gene={domain_type};"
                    f"clade={clade};evalue={evalue};coverage={hmm_cov:.1f};"
                    f"probability={acc}")

            line = [base_seq, "TEsorter", "CDS",
                    str(nuc_start), str(nuc_end),
                    f"{norm_score}", strand, str(nuc_frame), attr]
            f.write("\t".join(line) + "\n")

    n_lines = sum(1 for _ in open(out_path))
    log.info(f"  {out_path}: {n_lines} domain annotations")
    return out_path


def generate_dom_tsv(conn, prefix, db_name):
    """Generate domain hit summary table."""
    out_path = f"{prefix}.dom.tsv"

    with open(out_path, "w") as f:
        f.write("#id\tlength\tevalue\tcoverge\tprobability\tscore\n")

        for row in conn.execute("""
            SELECT target_name, query_name, dom_score, i_evalue, acc,
                   hmm_from, hmm_to, query_len, env_from, env_to,
                   base_seq
            FROM legacy_hits
            WHERE database = ?
            ORDER BY base_seq, env_from
        """, (db_name,)):
            (target, model, score, evalue, acc,
             hmm_from, hmm_to, model_len, env_from, env_to,
             base_seq) = row

            hmm_cov = 100.0 * (hmm_to - hmm_from + 1) / model_len
            norm_score = round(score / model_len, 2)

            if not (hmm_cov >= 20 and evalue <= 1e-3 and acc >= 0.5
                    and norm_score >= 0.1):
                continue

            gid = f"{_format_gff_id(base_seq)}|{model}"
            seq_len = env_to - env_from + 1

            f.write(f"{gid}\t{seq_len}\t{evalue}\t{hmm_cov:.1f}\t{acc}\t{norm_score}\n")

    log.info(f"  {out_path}")
    return out_path


def generate_dom_faa(conn, prefix, db_name, aa_fasta):
    """Generate domain protein sequences FASTA."""
    out_path = f"{prefix}.dom.faa"

    # Load translated sequences
    seqs = load_sequences_dict(aa_fasta)

    with open(out_path, "w") as f:
        for row in conn.execute("""
            SELECT target_name, query_name, dom_score, i_evalue, acc,
                   hmm_from, hmm_to, query_len, env_from, env_to,
                   base_seq, domain_type
            FROM legacy_hits
            WHERE database = ?
            ORDER BY base_seq, env_from
        """, (db_name,)):
            (target, model, score, evalue, acc,
             hmm_from, hmm_to, model_len, env_from, env_to,
             base_seq, domain_type) = row

            hmm_cov = 100.0 * (hmm_to - hmm_from + 1) / model_len
            norm_score = round(score / model_len, 2)

            if not (hmm_cov >= 20 and evalue <= 1e-3 and acc >= 0.5
                    and norm_score >= 0.1):
                continue

            if target not in seqs:
                continue

            subseq = seqs[target][env_from - 1:env_to]
            gene, clade = _parse_gene_clade(model, db_name)
            gid = f"{_format_gff_id(base_seq)}|{model}"
            name = f"{clade}-{gene.split('-')[-1] if '-' in gene else gene}"
            attr = (f"ID={gid};Name={name};gene={domain_type};"
                    f"clade={clade};evalue={evalue};coverage={hmm_cov:.1f};"
                    f"probability={acc}")

            f.write(f">{gid} {attr}\n{subseq}\n")

    log.info(f"  {out_path}")
    return out_path


def generate_cls_lib(input_fasta, prefix, classifications,
                     no_reverse=False):
    """Generate RepeatMasker library FASTA.

    Each sequence gets its classification appended to the ID:
        >seqid#Order/Superfamily/Clade
    Minus-strand sequences are reverse complemented unless no_reverse.
    """
    out_path = f"{prefix}.cls.lib"

    cls_dict = {r["id"]: r for r in classifications}

    # Load all sequences into memory
    fa = pyfastx.Fasta(input_fasta, build_index=True)
    all_seqs = {rec.name: str(rec.seq) for rec in fa}

    with open(out_path, "w") as f:
        for name, seq in all_seqs.items():
            base_id = name.split("#")[0]
            if name in cls_dict:
                cl = cls_dict[name]
                strand = cl.get("strand", "+")

                parts = [cl["order"]]
                if cl["superfamily"] != "unknown":
                    parts.append(cl["superfamily"])
                if cl["clade"] not in ("unknown", "mixture"):
                    parts.append(cl["clade"])
                cls_str = "/".join(parts)

                if not no_reverse and strand == "-":
                    seq = revcomp(seq.upper())
            else:
                cls_str = "Unknown"

            new_id = base_id + "#" + cls_str
            f.write(f">{new_id}\n{seq}\n")

    log.info(f"  {out_path}")
    return out_path


def generate_cls_pep(conn, prefix, db_name, aa_fasta, classifications):
    """Generate classified protein domain sequences.

    Each domain sequence gets the full classification in the ID:
        >seqid#Order/Superfamily/Clade#gene|clade
    """
    out_path = f"{prefix}.cls.pep"

    cls_dict = {r["id"]: r for r in classifications}
    seqs = load_sequences_dict(aa_fasta)

    with open(out_path, "w") as f:
        for row in conn.execute("""
            SELECT target_name, query_name, dom_score, i_evalue, acc,
                   hmm_from, hmm_to, query_len, env_from, env_to,
                   base_seq, domain_type
            FROM legacy_hits
            WHERE database = ?
            ORDER BY base_seq, env_from
        """, (db_name,)):
            (target, model, score, evalue, acc,
             hmm_from, hmm_to, model_len, env_from, env_to,
             base_seq, domain_type) = row

            hmm_cov = 100.0 * (hmm_to - hmm_from + 1) / model_len
            norm_score = round(score / model_len, 2)

            if not (hmm_cov >= 20 and evalue <= 1e-3 and acc >= 0.5
                    and norm_score >= 0.1):
                continue

            if target not in seqs or base_seq not in cls_dict:
                continue

            subseq = seqs[target][env_from - 1:env_to]
            gene, clade = _parse_gene_clade(model, db_name)
            cl = cls_dict[base_seq]

            parts = [cl["order"]]
            if cl["superfamily"] != "unknown":
                parts.append(cl["superfamily"])
            if cl["clade"] not in ("unknown", "mixture"):
                parts.append(cl["clade"])
            cls_str = "/".join(parts)

            raw_id = base_seq.split("#")[0]
            new_id = f"{raw_id}#{cls_str}#{domain_type}|{clade}"

            gid = f"{_format_gff_id(base_seq)}|{model}"
            attr = (f"ID={gid};Name={clade}-{domain_type};"
                    f"gene={domain_type};clade={clade};"
                    f"evalue={evalue};coverage={hmm_cov:.1f};"
                    f"probability={acc}")

            f.write(f">{new_id} {attr}\n{subseq}\n")

    log.info(f"  {out_path}")
    return out_path


def generate_all_outputs(conn, prefix, db_name, input_fasta, aa_fasta,
                         nucl_lengths, classifications,
                         seq_type="nucl", no_reverse=False, no_library=False):
    """Generate all TEsorter-compatible output files."""
    log.info(f"Generating TEsorter-format outputs for {db_name}")

    generate_dom_gff3(conn, prefix, db_name, nucl_lengths, seq_type)
    generate_dom_tsv(conn, prefix, db_name)

    if aa_fasta:
        generate_dom_faa(conn, prefix, db_name, aa_fasta)
        generate_cls_pep(conn, prefix, db_name, aa_fasta, classifications)

    if not no_library:
        generate_cls_lib(input_fasta, prefix, classifications, no_reverse)
