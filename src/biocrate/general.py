import pickle as pkl
import json
import lmdb
from tqdm import tqdm
import os
from os import PathLike
from pathlib import Path
from typing import Any
import orjson
from .constants import MAP_SIZE
import uuid
import pyfastx


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
        pkl.dump(content, file)


def json_load(file_path: str | PathLike):
    """
    读取文件内容
    """
    with open(file_path, "rb") as file:
        return orjson.loads(file.read())


def json_dump(file_path: str | PathLike, *, content: Any):
    """
    写内容到文件
    """
    with open(file_path, "w", encoding="utf-8") as file:
        json.dump(content, file, ensure_ascii=False, indent=4)


def lmdb_dump(save_path: str | PathLike, *, out_list: list, size: int = 512):

    env = lmdb.open(str(save_path), subdir=False, lock=False, readahead=False, meminit=False, max_readers=64, map_size=size * MAP_SIZE)

    with env.begin(write=True) as lmdb_txn:
        for i in tqdm(range(len(out_list)), desc="Writing to LMDB"):
            lmdb_txn.put(str(i).encode("ascii"), pkl.dumps(out_list[i]))


def lmdb_load(file_path: str | PathLike, size: int = 512):
    env = lmdb.open(str(file_path), subdir=False, lock=False, readahead=False, meminit=False, max_readers=64, map_size=size * MAP_SIZE)
    with env.begin() as lmdb_txn:
        with lmdb_txn.cursor() as cursor:
            for _, value in cursor:
                yield pkl.loads(value)


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


def fasta_load(file_path: str | PathLike):
    fa = pyfastx.Fasta(file_path, build_index=False)
    for name, seq in fa:
        yield name, str(seq)


def fasta_dump(file_path: str | PathLike, *, content: list[tuple[str, str]] | dict[str, str]):
    res_list = []
    if isinstance(content, dict):
        content = content.items()  # type: ignore
    for uid, seq in content:
        res_list.append(f">{uid}\n{seq}\n")
    txt_dump(file_path, content="".join(res_list))


########################################
def make_dir(path: str | PathLike):
    return Path(path).mkdir(parents=True, exist_ok=True)


def base_name(file_path: str | PathLike):
    return Path(file_path).stem


def tmp_name():
    return uuid.uuid4().hex
