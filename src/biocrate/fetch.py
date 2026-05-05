import gzip
import os
import re
import shutil
import urllib.request as request
from .constants import SCOP_MIRROR, CIF_MIRROR, RCSB_PDB_MIRROR, UNIPROT_MIRROR, RCSB_FASTA_MIRROR


def fetch_scope(pdb_id):
    """download SCOPe structure pdb_id"""
    fp = request.urlopen(SCOP_MIRROR)
    txt = fp.read().decode("utf-8")
    fp.close()
    scop_version = sorted(set(re.findall(r"pdbstyle\-\d+\.\d+", txt)))[-1]
    url = "/".join([SCOP_MIRROR, scop_version, pdb_id[2:4], pdb_id + ".ent"])
    pdb_file = pdb_id + ".pdb"
    request.urlretrieve(url, pdb_file)
    return pdb_file


def fetch_rcsb_cif(pdb_id):
    """download pdb_id.cif"""
    pdb_id = pdb_id.lower()
    cif_file = pdb_id + ".cif"
    cif_gz_file = cif_file + ".gz"

    request.urlretrieve(CIF_MIRROR + pdb_id + ".cif.gz", cif_gz_file)

    with gzip.open(cif_gz_file, "rb") as fp_gz:
        with open(cif_file, "wb") as fp_cif:
            shutil.copyfileobj(fp_gz, fp_cif)  # type: ignore
    os.unlink(path=cif_gz_file)
    return cif_file


def fetch_rcsb_pdb(pdb_id):
    """download pdb_id.pdb"""
    pdb_id = pdb_id.lower()
    pdb_file = pdb_id + ".pdb"

    request.urlretrieve(RCSB_PDB_MIRROR + pdb_file, pdb_file)
    return pdb_file


def fetch_uniprot_sequence(accession):
    """retrieve fasta sequence for accession"""
    fp = request.urlopen(UNIPROT_MIRROR + accession + ".fasta")
    txt = fp.read().decode("utf-8")
    fp.close()
    return txt


def fetch_rcsb_sequence(pdb_id):
    """retrieve SEQRES sequence for pdb_id"""
    fp = request.urlopen(RCSB_FASTA_MIRROR + pdb_id)
    txt = fp.read().decode("utf-8")
    fp.close()
    return txt


def fetch_rcsb_chain_sequence(pdb_chain):
    """retrieve SEQRES sequence for pdb_chain"""
    pdb_id = pdb_chain[:4].upper()
    chainID = pdb_chain[4:]
    fasta = fetch_rcsb_sequence(pdb_chain[:4])
    txt = "\n".join([">" + s for s in fasta.split(">") if s.startswith(pdb_id + ":" + chainID)])
    return txt
