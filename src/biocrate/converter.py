from .general import *


def json2fasta(json_file: str | PathLike, fasta_file: str | PathLike) -> None:
    uid2seq = json_load(json_file)

    with open(fasta_file, "w", encoding="utf-8") as f:
        output = []
        for i, (uid, seq) in enumerate(uid2seq.items()):
            output.append(f">{uid}\n{seq}\n")

            if i % 1000 == 0:
                f.write("".join(output))
                output = []

        if output:
            f.write("".join(output))
            

def fasta2json(fasta_file: str | PathLike, json_file: str | PathLike) -> None:
    uid2seq = {}
    with open(fasta_file, "r", encoding="utf-8") as f:
        uid = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if uid is not None:
                    uid2seq[uid] = "".join(seq)
                uid = line[1:]
                seq = []
            else:
                seq.append(line)

        if uid is not None:
            uid2seq[uid] = "".join(seq)

    json_dump(json_file, content=uid2seq)