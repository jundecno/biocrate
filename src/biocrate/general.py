import pickle as pkl
import json
import os
from os import PathLike
from pathlib import Path
from typing import Any
import uuid


def pkl_load(file_path: str | PathLike):
    """
    读取文件内容
    """
    with open(file_path, "rb") as file:
        return pkl.load(file)


def pkl_dump(file_path: str | PathLike, *, content: Any):
    """
    写内容到文件
    """
    with open(file_path, "wb") as file:
        pkl.dump(content, file, protocol=pkl.HIGHEST_PROTOCOL)


def json_load(file_path: str | PathLike):
    """
    读取文件内容
    """
    import orjson

    with open(file_path, "rb") as file:
        return orjson.loads(file.read())


def json_dump(file_path: str | PathLike, *, content: Any):
    """
    写内容到文件
    """
    import orjson

    with open(file_path, "wb") as file:
        file.write(orjson.dumps(content, option=orjson.OPT_INDENT_2 | orjson.OPT_APPEND_NEWLINE))


def lmdb_dump(save_path: str | PathLike, *, content: list | dict, size: int = 512):
    import lmdb
    from tqdm import tqdm

    from .constants import MAP_SIZE

    env = lmdb.open(str(save_path), subdir=False, lock=False, readahead=False, meminit=False, max_readers=64, map_size=size * MAP_SIZE)
    try:
        if isinstance(content, list):
            items = enumerate(content)
        elif isinstance(content, dict):
            items = content.items()
        else:
            raise TypeError("content must be a list or dict")

        with env.begin(write=True) as lmdb_txn:
            for key, value in tqdm(items, total=len(content), desc="Writing to LMDB"):
                lmdb_txn.put(str(key).encode("ascii"), pkl.dumps(value, protocol=pkl.HIGHEST_PROTOCOL))
    finally:
        env.close()


def lmdb_load(file_path: str | PathLike, size: int = 512):
    import lmdb

    from .constants import MAP_SIZE

    env = lmdb.open(str(file_path), subdir=False, lock=False, readahead=False, meminit=False, max_readers=64, map_size=size * MAP_SIZE)
    try:
        with env.begin() as lmdb_txn:
            with lmdb_txn.cursor() as cursor:
                for _, value in cursor:
                    yield pkl.loads(value)
    finally:
        env.close()


def txt_load(file_path: str | PathLike):
    """
    读取文件内容
    """
    with open(file_path, "r", encoding="utf-8") as file:
        return file.read()


def txt_append(file_path: str | PathLike, content: str):
    """
    追加内容到文件
    """
    with open(file_path, "a", encoding="utf-8") as file:
        file.write(content)


def txt_dump(file_path: str | PathLike, content: str):
    """
    写内容到文件
    """
    with open(file_path, "w", encoding="utf-8") as file:
        file.write(content)


def lines_load(file_path: str | PathLike):
    """
    读取文件内容
    """
    with open(file_path, "r", encoding="utf-8") as file:
        res = file.readlines()
    return [x.strip() for x in res]


def lines_dump(file_path: str | PathLike, content: list[str]):
    """
    写内容到文件
    """
    with open(file_path, "w", encoding="utf-8") as file:
        file.write("\n".join(content) + "\n")


def fasta_load(file_path: PathLike, as_dict: bool = True):
    """
    Fast FASTA parser (memory efficient, no dependencies)

    Args:
        file_path: fasta path
        as_dict: True -> return dict, False -> generator

    Returns:
        dict or iterator of (id, seq)
    """
    file_path = Path(file_path)

    def _generator():
        seq_id = None
        seq_parts = []

        with file_path.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    # flush previous
                    if seq_id is not None:
                        yield seq_id, "".join(seq_parts)

                    seq_id = line[1:].split()[0]  # only first token
                    seq_parts = []
                else:
                    seq_parts.append(line)

            # last record
            if seq_id is not None:
                yield seq_id, "".join(seq_parts)

    if as_dict:
        return dict(_generator())

    return _generator()


def fasta_dump(file_path: str | PathLike, *, content: list[tuple[str, str]] | dict[str, str]):
    if isinstance(content, dict):
        content = content.items()  # type: ignore
    with open(file_path, "w", encoding="utf-8") as file:
        for uid, seq in content:
            file.write(f">{uid}\n{seq}\n")


def structure_load(file_path):
    from Bio.PDB import PDBParser, MMCIFParser  # type: ignore

    ext = os.path.splitext(file_path)[-1].lower()
    if ext == ".pdb":
        return PDBParser(QUIET=True).get_structure(tmp_name(), file_path)
    elif ext in [".cif", ".mmcif"]:
        return MMCIFParser(QUIET=True).get_structure(tmp_name(), file_path)
    raise ValueError(f"Unsupported structure format: {ext}")


def structure_dump(structure, file_path):
    from Bio.PDB import PDBIO, MMCIFIO  # type: ignore

    ext = os.path.splitext(file_path)[-1].lower()
    if ext == ".pdb":
        io = PDBIO()
        io.set_structure(structure)
        io.save(file_path)
    elif ext in [".cif", ".mmcif"]:
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(file_path)
    else:
        raise ValueError(f"Unsupported structure format: {ext}")


########################################
def make_dir(path: str | PathLike):
    return Path(path).mkdir(parents=True, exist_ok=True)


def base_name(file_path: str | PathLike):
    return Path(file_path).stem


def tmp_name():
    return uuid.uuid4().hex

def is_uniprot_id(s: str) -> bool:
    from .constants import UNIPROT_ACCESSION_PATTERN

    s = s.strip().upper()
    return bool(UNIPROT_ACCESSION_PATTERN.fullmatch(s))

########################################
# 读取jsonl文件， 添加uniprot访问函数
def jsonl_load(file_path: str | PathLike):
    import orjson

    with open(file_path, "rb") as f:
        for line in f:
            line = line.strip()
            if line:
                yield orjson.loads(line)


def jsonl_dump(data, file_path: str | PathLike):
    import orjson

    with open(file_path, "ab") as f:
        f.write(orjson.dumps(data))
        f.write(b"\n")
