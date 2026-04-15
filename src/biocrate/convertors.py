from .general import *
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO, MMCIFParser, MMCIFIO  # type: ignore
from Bio.Data.PDBData import protein_letters_3to1_extended, protein_letters_3to1  # type: ignore
from openbabel import pybel  # type: ignore
import pysam
import pybedtools


def json2fasta(json_file: str | PathLike, fasta_file: str | PathLike) -> None:
    uid2seq = json_load(json_file)
    res_list = []
    with open(fasta_file, "w", encoding="utf-8") as f:
        for i, (uid, seq) in enumerate(uid2seq.items()):
            res_list.append(f">{uid}\n{seq}\n")
    txt_dump(fasta_file, content="".join(res_list))


def fasta2json(fasta_file: str | PathLike, json_file: str | PathLike) -> None:
    uid2seq = {}
    for uid, seq in fasta_load(fasta_file):
        uid2seq[uid] = seq
    json_dump(json_file, content=uid2seq)


############ Bioinformatics related functions ############
def pdb2mmcif(pdb_file: str | PathLike, mmcif_file: str | PathLike) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(tmp_name(), pdb_file)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(mmcif_file)


def mmcif2pdb(mmcif_file: str | PathLike, pdb_file: str | PathLike) -> None:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(tmp_name(), mmcif_file)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)


def pdb2fasta(pdb_file: str | PathLike, standard=False) -> dict:
    parser = PDBParser(QUIET=True)
    model = parser.get_structure(tmp_name(), pdb_file)[0]  # type: ignore
    res_dict = {}
    for chain in model:
        if standard:
            chain_str = "".join(protein_letters_3to1.get(res.get_resname(), "X") for res in chain)
        else:
            chain_str = "".join(protein_letters_3to1_extended.get(res.get_resname(), "X") for res in chain)
        res_dict[chain.id] = chain_str
    return res_dict


def mmcif2fasta(mmcif_file: str | PathLike, standard=False) -> dict:
    parser = MMCIFParser(QUIET=True)
    model = parser.get_structure(tmp_name(), mmcif_file)[0]  # type: ignore
    res_dict = {}
    for chain in model:
        if standard:
            chain_str = "".join(protein_letters_3to1.get(res.get_resname(), "X") for res in chain)
        else:
            chain_str = "".join(protein_letters_3to1_extended.get(res.get_resname(), "X") for res in chain)
        res_dict[chain.id] = chain_str
    return res_dict


def sdf2mol2(sdf_file: str | PathLike, mol2_file: str | PathLike) -> None:
    output = pybel.Outputfile("mol2", str(mol2_file), overwrite=True)
    for mol in pybel.readfile("sdf", str(sdf_file)):
        output.write(mol)
    output.close()


def mol2sdf(mol2_file: str | PathLike, sdf_file: str | PathLike) -> None:
    output = pybel.Outputfile("sdf", str(sdf_file), overwrite=True)
    for mol in pybel.readfile("mol2", str(mol2_file)):
        output.write(mol)
    output.close()


def fasta2embl(fasta_file: str | PathLike, embl_file: str | PathLike) -> None:
    records = list(SeqIO.parse(str(fasta_file), "fasta"))
    for record in records:
        if not record.annotations.get("molecule_type"):
            record.annotations["molecule_type"] = "DNA"  # 默认为 DNA
    SeqIO.write(records, str(embl_file), "embl")


def fasta2phylip(fasta_file: str | PathLike, phylip_file: str | PathLike) -> None:
    records = SeqIO.parse(str(fasta_file), "fasta")
    SeqIO.write(records, str(phylip_file), "phylip-relaxed")


def fasta2nexuas(fasta_file: str | PathLike, nexus_file: str | PathLike) -> None:
    records = SeqIO.parse(str(fasta_file), "fasta")
    SeqIO.write(records, str(nexus_file), "nexus")


def fasta2holmes(fasta_file: str | PathLike, stockholm_file: str | PathLike) -> None:
    records = SeqIO.parse(str(fasta_file), "fasta")
    SeqIO.write(records, str(stockholm_file), "stockholm")


def sam2bam(sam_file: str | PathLike, bam_file: str | PathLike) -> None:
    with pysam.AlignmentFile(str(sam_file), "r") as inf:
        with pysam.AlignmentFile(str(bam_file), "wb", template=inf) as outf:
            for read in inf:
                outf.write(read)
    pysam.index(str(bam_file))


def bam2sam(bam_file: str | PathLike, sam_file: str | PathLike) -> None:
    with pysam.AlignmentFile(str(bam_file), "rb") as inf:
        with pysam.AlignmentFile(str(sam_file), "w", template=inf) as outf:
            for read in inf:
                outf.write(read)


def psl2bed(psl_file: str | PathLike, bed_file: str | PathLike) -> None:
    pybedtools.BedTool(psl_file).psl_to_bed().saveas(bed_file)  # type: ignore


def gtf2bed(gtf_file: str | PathLike, bed_file: str | PathLike) -> None:
    pybedtools.BedTool(gtf_file).gtf_to_bed().saveas(bed_file)  # type: ignore


def vcf2bed(vcf_file: str | PathLike, bed_file: str | PathLike) -> None:
    pybedtools.BedTool(vcf_file).to_bed().saveas(bed_file)  # type: ignore
