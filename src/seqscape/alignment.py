from __future__ import annotations

import concurrent.futures
import io
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple

import numpy as np
from Bio import SeqIO


def resolve_muscle_exe(path_hint: str) -> str:
    if path_hint and path_hint != "muscle" and Path(path_hint).is_file():
        return str(Path(path_hint).resolve())
    repo_root = Path(__file__).resolve().parents[3]
    candidates = [
        Path("/Users/gerritkoorsen/opt/anaconda3/envs/ciderseq-muscle/bin/muscle"),
        repo_root / ".muscle-env" / "bin" / "muscle",
        repo_root.parent / ".muscle-env" / "bin" / "muscle",
    ]
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
    found = shutil.which("muscle")
    if found:
        return found
    return path_hint or "muscle"


def global_identity_pairwisealigner(ref_seq: str, qry_seq: str) -> float:
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5

    aln = aligner.align(ref_seq, qry_seq)[0]
    idx = aln.indices
    n_cols = idx.shape[1]
    if n_cols == 0:
        return 0.0
    matches = 0
    ungapped_cols = 0
    for i in range(n_cols):
        a = idx[0, i]
        b = idx[1, i]
        if a >= 0 and b >= 0:
            ungapped_cols += 1
            if ref_seq[a] == qry_seq[b]:
                matches += 1
    if ungapped_cols == 0:
        return 0.0
    return 100.0 * matches / ungapped_cols


def global_identity_mafft(ref_seq: str, qry_seq: str, mafft_exe: str) -> float:
    with tempfile.TemporaryDirectory(prefix="seqscape_mafft_") as td:
        in_fa = Path(td) / "pair.fasta"
        with in_fa.open("w") as fh:
            fh.write(">ref\n")
            fh.write(ref_seq + "\n")
            fh.write(">qry\n")
            fh.write(qry_seq + "\n")
        cmd = [mafft_exe, "--quiet", "--thread", "1", str(in_fa)]
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"mafft failed: {p.stderr}")
        recs = list(SeqIO.parse(io.StringIO(p.stdout), "fasta"))
        if len(recs) != 2:
            return 0.0
        a = str(recs[0].seq).upper()
        b = str(recs[1].seq).upper()
        n = min(len(a), len(b))
        if n == 0:
            return 0.0
        matches = 0
        ungapped = 0
        for i in range(n):
            if a[i] != "-" and b[i] != "-":
                ungapped += 1
                if a[i] == b[i]:
                    matches += 1
        if ungapped == 0:
            return 0.0
        return 100.0 * matches / ungapped


def global_identity_muscle(ref_seq: str, qry_seq: str, muscle_exe: str) -> float:
    with tempfile.TemporaryDirectory(prefix="seqscape_muscle_") as td:
        in_fa = Path(td) / "pair.fasta"
        out_fa = Path(td) / "pair.aln.fasta"
        with in_fa.open("w") as fh:
            fh.write(">ref\n")
            fh.write(ref_seq + "\n")
            fh.write(">qry\n")
            fh.write(qry_seq + "\n")

        commands = [
            [muscle_exe, "-align", str(in_fa), "-output", str(out_fa)],
            [muscle_exe, "-in", str(in_fa), "-out", str(out_fa)],
        ]
        stderr_text = ""
        for cmd in commands:
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stderr_text = p.stderr or stderr_text
            if p.returncode == 0 and out_fa.is_file():
                recs = list(SeqIO.parse(str(out_fa), "fasta"))
                if len(recs) != 2:
                    return 0.0
                a = str(recs[0].seq).upper()
                b = str(recs[1].seq).upper()
                n = min(len(a), len(b))
                if n == 0:
                    return 0.0
                matches = 0
                ungapped = 0
                for i in range(n):
                    if a[i] != "-" and b[i] != "-":
                        ungapped += 1
                        if a[i] == b[i]:
                            matches += 1
                if ungapped == 0:
                    return 0.0
                return 100.0 * matches / ungapped
        raise RuntimeError(f"muscle failed: {stderr_text}")


def pairwise_identity(a: str, b: str, aligner: str, mafft: str, muscle: str) -> float:
    if aligner == "mafft":
        return global_identity_mafft(a, b, mafft)
    if aligner == "muscle":
        return global_identity_muscle(a, b, muscle)
    return global_identity_pairwisealigner(a, b)


def build_identity_matrix(
    items: List[dict],
    aligner: str,
    mafft: str,
    muscle: str,
    jobs: int,
    progress_every: int,
) -> tuple[list[str], np.ndarray]:
    ids = [item["id"] for item in items]
    seqs = [item["sequence"] for item in items]
    n = len(seqs)
    matrix = np.eye(n, dtype=float) * 100.0
    pairs: List[Tuple[int, int]] = [(i, j) for i in range(n) for j in range(i + 1, n)]

    def work(pair: Tuple[int, int]) -> Tuple[int, int, float]:
        i, j = pair
        ident = pairwise_identity(seqs[i], seqs[j], aligner, mafft, muscle)
        return i, j, ident

    done = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, jobs)) as ex:
        for i, j, ident in ex.map(work, pairs):
            matrix[i, j] = ident
            matrix[j, i] = ident
            done += 1
            if progress_every > 0 and done % progress_every == 0:
                print(f"progress\t{done}/{len(pairs)}", flush=True)
    return ids, matrix


def identity_to_distance(ident: np.ndarray) -> np.ndarray:
    dist = 1.0 - (ident / 100.0)
    np.fill_diagonal(dist, 0.0)
    return dist
