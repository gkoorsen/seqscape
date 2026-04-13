from __future__ import annotations

import csv
from pathlib import Path


def load_cg_map(path: Path) -> dict[str, str]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    needed = {"identified_id", "orig_ids"}
    if not rows:
        return {}
    missing = needed - set(rows[0].keys())
    if missing:
        raise RuntimeError(f"Missing required columns in {path}: {sorted(missing)}")

    out: dict[str, str] = {}
    for row in rows:
        cg_id = row.get("identified_id", "").strip()
        raw = row.get("orig_ids", "").strip()
        if not cg_id or not raw:
            continue
        tokens = [raw]
        for sep in [",", ";"]:
            expanded: list[str] = []
            for token in tokens:
                expanded.extend(part.strip() for part in token.split(sep) if part.strip())
            tokens = expanded
        for token in tokens:
            out[token] = cg_id
    return out
