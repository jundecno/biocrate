import gzip
import json
import os
import re
import shutil
import urllib.request as request
from .constants import RCSB_LIGAND_MIRROR, SCOP_MIRROR, CIF_MIRROR, RCSB_PDB_MIRROR, UNIPROT_MIRROR, RCSB_FASTA_MIRROR
import requests


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


def fetch_rcsb_ligand(lig_id: str):
    """download lig_id.cif"""
    lig_id = lig_id.upper()
    # 获取json数据
    url = RCSB_LIGAND_MIRROR + lig_id
    fp = request.urlopen(url)
    txt = fp.read().decode("utf-8")
    fp.close()
    return json.loads(txt)


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


def fetch_uniprot_info(uniprot_id, seq: bool = False):
    info = {}

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"

    try:
        r = requests.get(url, timeout=20)
        if r.status_code != 200:
            return info

        data = r.json()

        # protein / enzyme name
        protein_desc = data.get("proteinDescription", {})
        rec_name = protein_desc.get("recommendedName", {})
        full_name = rec_name.get("fullName", {})
        info["enzyme_name"] = full_name.get("value", "")

        # fallback: submittedName
        if not info["enzyme_name"]:
            submitted = protein_desc.get("submissionNames", [])
            if submitted:
                info["enzyme_name"] = submitted[0].get("fullName", {}).get("value", "")

        # organism
        organism = data.get("organism", {})
        info["organism"] = organism.get("scientificName", "")

        # taxonomy ID
        info["tax_id"] = organism.get("taxonId", "")

        # sequence
        if seq:
            info["sequence"] = data.get("sequence", {}).get("value", "")

    except Exception as e:
        print(f"Failed to fetch {uniprot_id}: {e}")

    return info


def fetch_uniparc_info(accession: str, seq: bool = False):
    info = {}

    try:
        # 网页搜索框使用的是普通关键词检索
        response = requests.get(
            "https://rest.uniprot.org/uniparc/search",
            params={
                "query": accession,
                "format": "json",
                "size": 10,
            },
            timeout=30,
        )
        response.raise_for_status()

        results = response.json().get("results", [])
        if not results:
            print(f"No UniParc record found for: {accession}")
            return info

        # 一般精确 accession 搜索的第一个结果就是目标条目
        upi = results[0].get("uniParcId", "")
        if not upi:
            print(f"No UniParc ID found for: {accession}")
            return info

        # 使用 UPI 获取完整记录
        detail_response = requests.get(
            f"https://rest.uniprot.org/uniparc/{upi}.json",
            timeout=30,
        )
        detail_response.raise_for_status()

        record = detail_response.json()

        info["protein_id"] = accession
        info["uniparc_id"] = record.get("uniParcId", upi)
        info["enzyme_name"] = ""
        info["organism"] = ""
        info["tax_id"] = ""
        info["source"] = "UniParc"

        # UniParc 详情中通常可以获取序列
        if seq:
            sequence = record.get("sequence", {})

            if isinstance(sequence, dict):
                info["sequence"] = sequence.get("value", "")
            elif isinstance(sequence, str):
                info["sequence"] = sequence
            else:
                info["sequence"] = ""

        return info

    except requests.RequestException as exc:
        print(f"Failed to fetch UniParc {accession}: {exc}")
        return info
