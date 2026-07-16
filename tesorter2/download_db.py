"""
download_db.py — Fetch and install the TESorter2 HMM databases.

The databases (REXdb, GyDB2, LINE, TIR, AnnoSINE, SINE_SO plus their
similarity-graph sidecars) are distributed separately from the code because
they are hundreds of MB. This installs them into the user data directory
(``$XDG_DATA_HOME/tesorter2/database`` by default) where the pipeline finds
them automatically.

Usage:
    tesorter2-download-db                 # default release, default location
    tesorter2-download-db --db-dir DIR    # install elsewhere
    tesorter2-download-db --url URL       # custom tarball

Set the resulting directory as TESORTER2_DB, or pass --db-dir to tesorter2, if
you install to a non-default location.
"""

import argparse
import os
import sys
import tarfile
import tempfile
import urllib.request

from .paths import user_data_dir

# Canonical database bundle. Set to the release/Zenodo tarball URL once the
# databases are published (GitHub release asset or Zenodo DOI). Until then,
# users must pass --url or point --db-dir at a local copy.
# TODO(maintainer): fill in the published tarball URL before release.
DB_RELEASE_URL = None

# Files expected in a complete installation (used to verify after extraction).
EXPECTED = [
    "REXdb_protein_database_viridiplantae_v4.0_plus_metazoa_v3.1.hmm",
    "GyDB2.hmm", "GyDB2.hmm.info",
    "Kapitonov_et_al.GENE.LINE.hmm",
    "Yuan_and_Wessler.PNAS.TIR.hmm",
    "AnnoSINE_core.hmm", "AnnoSINE.hmm", "SINE_SO.hmm",
]


def _report(done, total):
    if total > 0:
        pct = min(100, done * 100 // total)
        sys.stderr.write(f"\r  downloading... {pct:3d}%")
        sys.stderr.flush()


def download(url, dest_tar):
    def hook(block_num, block_size, total_size):
        _report(block_num * block_size, total_size)
    urllib.request.urlretrieve(url, dest_tar, reporthook=hook)
    sys.stderr.write("\n")


def extract(tar_path, db_dir):
    os.makedirs(db_dir, exist_ok=True)
    with tarfile.open(tar_path) as tf:
        # Flatten: drop any leading directory component so files land in db_dir.
        for member in tf.getmembers():
            if not member.isfile():
                continue
            member.name = os.path.basename(member.name)
            tf.extract(member, db_dir)


def verify(db_dir):
    missing = [f for f in EXPECTED if not os.path.isfile(os.path.join(db_dir, f))]
    return missing


def main():
    parser = argparse.ArgumentParser(
        prog="tesorter2-download-db",
        description="Download and install the TESorter2 HMM databases.")
    parser.add_argument(
        "--db-dir", default=None,
        help="Install location [default: $XDG_DATA_HOME/tesorter2/database]")
    parser.add_argument(
        "--url", default=DB_RELEASE_URL,
        help="Tarball URL to fetch [default: published release]")
    parser.add_argument(
        "--keep-tarball", action="store_true",
        help="Do not delete the downloaded tarball after extraction")
    args = parser.parse_args()

    db_dir = os.path.abspath(args.db_dir) if args.db_dir \
        else os.path.join(user_data_dir(), "database")

    if not args.url:
        sys.exit(
            "No database URL configured. Pass --url with the tarball location, "
            "or point tesorter2 at an existing copy via --db-dir / TESORTER2_DB.")

    print(f"Installing databases into: {db_dir}")
    tmp = tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False)
    tmp.close()
    try:
        download(args.url, tmp.name)
        print("Extracting...")
        extract(tmp.name, db_dir)
    finally:
        if not args.keep_tarball and os.path.exists(tmp.name):
            os.remove(tmp.name)

    missing = verify(db_dir)
    if missing:
        sys.exit("Installed, but these expected files are missing: "
                 + ", ".join(missing))
    print(f"Done. {len(EXPECTED)} databases installed.")
    default_dir = os.path.join(user_data_dir(), "database")
    if db_dir != default_dir:
        print(f"Non-default location: export TESORTER2_DB={db_dir} "
              f"(or pass --db-dir {db_dir} to tesorter2).")


if __name__ == "__main__":
    main()
