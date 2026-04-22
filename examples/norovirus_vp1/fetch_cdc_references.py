"""Fetch CDC Calicivirus Typing Tool reference sequences from NCBI.

Accessions sourced from https://calicivirustypingtool.cdc.gov (Chhabra 2019 nomenclature).
Writes a labeled FASTA: >GII.4_Sydney|JX459908.1
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO

Entrez.email = "g.koorsen@googlemail.com"

# CDC Calicivirus Typing Tool reference accessions, scraped 2026-04-18
# Genogroup I
_GI = {
    "GI.1":  ["EU085529", "M87661", "L23828", "AB447406"],
    "GI.2":  ["AF435807", "FJ515294", "L07418"],
    "GI.3":  ["EF547396", "GQ856470", "GQ856471", "GQ856472", "GQ856473",
               "JX846929", "U04469", "KM289169", "LC122713", "AB187514",
               "AF439267", "KJ196292", "AF145709", "AB447414", "AY038598", "EF529738"],
    "GI.4":  ["GQ856475", "AB042808", "LN854563", "AF394960", "AJ313030"],
    "GI.5":  ["AF414406", "AJ277614", "LC101825", "KJ402295", "AB039774", "MH443711"],
    "GI.6":  ["KP027330", "AF093797", "GQ856463", "GQ856464", "AB354289",
               "LC342057", "AJ277615", "AF538678", "AB081723"],
    "GI.7":  ["KU311161", "AJ844469", "MH130046", "AY675555", "JN603251",
               "FJ383816", "JN899243", "AJ277609"],
    "GI.8":  ["GU299761", "AF538679", "KJ196298"],
    "GI.9":  ["KF586507", "JN183159", "HQ637267"],
}

# Genogroup II — GII.4 variants kept as distinct labels
_GII = {
    "GII.1":  ["KJ194507", "JN797508", "AF425767", "U07611", "AY919139"],
    "GII.2":  ["AY682552", "MF405169", "JX846925", "DQ456824", "X81879",
                "AB089882", "DQ366347", "AY682549", "KF730316", "KY865306",
                "AY134748", "AB281090"],
    "GII.3":  ["MH218579", "KJ194500", "AF190817", "JX846924", "AF539439",
                "JF697282", "U22498", "JN565063", "KT239614", "MF140689",
                "AB385626", "AB190457", "U02030"],
    "GII.4":            ["FJ537135", "FJ537134"],
    "GII.4_Allegany":   ["MT028542", "MW019958", "PQ213438", "MW559992"],
    "GII.4_Apeldoorn":  ["AB445395", "HM748973"],
    "GII.4_Asia":       ["AB220921", "AB220922"],
    "GII.4_Bristol":    ["X76716"],
    "GII.4_DenHaag":    ["EF126965", "EU078417", "AB541348"],
    "GII.4_FarmingtonHills": ["AY502023", "AY485642"],
    "GII.4_Henry":      ["FJ411170"],
    "GII.4_HongKong":   ["MN400355"],
    "GII.4_Hunter":     ["AB220926", "DQ078814", "AY883096"],
    "GII.4_NewOrleans": ["AB933767", "GQ845367", "KX353972", "GU445325"],
    "GII.4_Osaka":      ["AB434770", "EU921388"],
    "GII.4_Richmond":   ["EU078406"],
    "GII.4_SanFrancisco": ["PV275048", "OR262344", "OR262322"],
    "GII.4_Sydney":     ["MW521126", "MK762630", "PQ207921", "PQ195892",
                          "MK752934", "PQ195913", "PQ209098", "PQ218566",
                          "PQ196561", "PQ208322", "PQ207934", "PQ214833",
                          "LC153121", "KX354134", "PQ196529", "OL336352",
                          "JX459908", "JX459907"],
    "GII.4_US9596":     ["AB112306", "AF414424", "AJ004864"],
    "GII.4_Wichita":    ["PQ218943", "PQ207944"],
    "GII.4_Yerseke":    ["AB294790", "EF126963"],
    "GII.5":  ["AJ277607", "AB212306", "AF397156", "AY682550", "KJ196277", "KJ196288"],
    "GII.6":  ["KX158281", "KM198534", "MN248516", "JX989075", "PV174245",
                "AF414410", "PV174244", "KC576910", "AB039778", "AB039777",
                "AJ277620", "HM633213", "GU930737", "AB684664"],
    "GII.7":  ["MH218692", "JX846926", "AF414409", "AJ277608", "AB258331", "KJ196295"],
    "GII.8":  ["AF195848", "AB039780"],
    "GII.9":  ["DQ379715", "AY038599"],
    "GII.10": ["AF427118", "AY237415", "AF504671"],
    "GII.11": ["AB126320", "AB074893"],
    "GII.12": ["LC342059", "EU921353", "AB039775", "GQ845370", "AJ277618"],
    "GII.13": ["KR904229", "AY113106", "DQ379714", "EU921354", "KJ196276"],
    "GII.14": ["GQ856465", "AY130761", "GU017903", "KM289171", "GU594162", "KJ196278"],
    "GII.16": ["AY502006", "AY772730", "AY502010"],
    "GII.17": ["KJ156329", "MK789431", "OP805362", "KX061540", "EF529741",
                "OR685283", "OP689689", "AY502009", "KT589391", "DQ438972",
                "LC037415", "AB983218", "KJ196286"],
    "GII.18": ["AY823304", "AY823305"],
    "GII.19": ["AY823306"],
    "GII.20": ["EU275779", "EU424333", "EU373815", "AB542917"],
    "GII.21": ["AY675554", "KJ196284", "JN899245"],
    "GII.22": ["MG495082", "AB233471", "AB083780"],
    "GII.23": ["KT290889", "MG495080", "KR232647"],
    "GII.24": ["MG495084", "KY225989"],
    "GII.25": ["GQ856469", "MG495083", "HM635128"],
    "GII.26": ["MF352142", "KU306738"],
    "GII.27": ["MK733205", "MG495077", "MG495078"],
}

_GIV = {
    "GIV.1": ["JQ970479", "DQ093067", "AF195847", "KC894731", "KT030674",
               "AF414426", "FM865412", "JQ613567", "LC101824", "AF414427"],
}

_GIX = {
    "GIX.1": ["GQ856474", "AB360387", "AY130762", "AF542090", "KJ196290"],
}

REFERENCE_ACCESSIONS: dict[str, list[str]] = {**_GI, **_GII, **_GIV, **_GIX}


def fetch_all(outfasta: Path, report_tsv: Path, batch: int = 200) -> None:
    # Build flat list preserving genotype label
    pairs: list[tuple[str, str]] = []  # (genotype, accession)
    for genotype, accessions in REFERENCE_ACCESSIONS.items():
        for acc in accessions:
            pairs.append((genotype, acc))

    acc_to_genotype = {acc: gt for gt, acc in pairs}
    all_accessions = [acc for _, acc in pairs]

    print(f"Fetching {len(all_accessions)} accessions in batches of {batch}...", file=sys.stderr)
    n_written = 0
    n_failed = 0

    with outfasta.open("w") as fa, report_tsv.open("w") as rep:
        rep.write("accession\tgenotype\tlength\tstatus\n")
        for i in range(0, len(all_accessions), batch):
            chunk = all_accessions[i: i + batch]
            ids = ",".join(chunk)
            try:
                handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text")
                records = list(SeqIO.parse(handle, "fasta"))
                handle.close()
            except Exception as e:
                print(f"  Batch {i}-{i+batch} error: {e}", file=sys.stderr)
                for acc in chunk:
                    rep.write(f"{acc}\t{acc_to_genotype[acc]}\t0\terror\n")
                n_failed += len(chunk)
                time.sleep(2)
                continue

            fetched_ids = {r.id.split(".")[0]: r for r in records}
            for acc in chunk:
                base = acc.split(".")[0]
                rec = fetched_ids.get(base)
                genotype = acc_to_genotype[acc]
                if rec is None:
                    rep.write(f"{acc}\t{genotype}\t0\tnot_found\n")
                    n_failed += 1
                    continue
                # Label: >GII.4_Sydney|JX459908.1
                fa.write(f">{genotype}|{rec.id}\n{rec.seq}\n")
                rep.write(f"{rec.id}\t{genotype}\t{len(rec.seq)}\tok\n")
                n_written += 1

            time.sleep(0.4)
            print(f"  {min(i+batch, len(all_accessions))}/{len(all_accessions)} done", file=sys.stderr)

    print(f"Written={n_written} failed={n_failed}", file=sys.stderr)


def main() -> int:
    ap = argparse.ArgumentParser(description="Fetch CDC calicivirus typing tool reference sequences.")
    ap.add_argument("--outfasta", type=Path, required=True)
    ap.add_argument("--report-tsv", type=Path, required=True)
    args = ap.parse_args()
    args.outfasta.parent.mkdir(parents=True, exist_ok=True)
    fetch_all(args.outfasta, args.report_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
