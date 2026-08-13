"""Regression tests for the OpenRDP wrapper's parsing and consensus logic.

Each test here corresponds to a defect that caused the IS76 positive control to
be reported as zero events. They are written to fail against the pre-fix code.
"""
from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from seqscape.recombination import (  # noqa: E402
    assign_roles,
    build_consensus,
    parse_openrdp_csv,
    polarise_event,
    read_alignment,
    resolve_openrdp,
    sanitise_alignment_for_openrdp,
)

# Verbatim OpenRDP output for the IS76 control (TYLCV/TYLCSV recombinant).
IS76_CSV = """Method,Start,End,Recombinant,Parent1,Parent2,Pvalue
Geneconv,1382,1492,IS76_LN812978,TYLCSV_NC_003828,-,0.00000
Geneconv,5,1359,IS76_LN812978,TYLCV_IL_AM409201,-,0.03171
Geneconv,1382,1492,TYLCV_IL_AM409201,-,-,0.00000
Geneconv,5,1359,TYLCSV_NC_003828,-,-,0.03037
Rdp,1343,1420,IS76_LN812978,TYLCSV_NC_003828,TYLCV_IL_AM409201,36.363636363636346
"""


class TestParsing(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(__file__).parent / "_tmp_openrdp.csv"
        self.tmp.write_text(IS76_CSV)
        self.events = parse_openrdp_csv(self.tmp)

    def tearDown(self):
        self.tmp.unlink(missing_ok=True)

    def test_all_rows_parsed(self):
        self.assertEqual(len(self.events), 5)

    def test_parent1_is_not_dropped(self):
        """The parser matched only 'major'/'minor'; OpenRDP emits Parent1/Parent2."""
        rdp = [e for e in self.events if e["method"] == "rdp"][0]
        self.assertEqual(rdp["parent_a"], "TYLCSV_NC_003828")
        self.assertEqual(rdp["parent_b"], "TYLCV_IL_AM409201")

    def test_dash_placeholder_becomes_empty(self):
        row = [e for e in self.events
               if e["method"] == "geneconv" and e["recombinant"] == "TYLCV_IL_AM409201"][0]
        self.assertEqual(row["parent_a"], "")
        self.assertEqual(row["parent_b"], "")

    def test_out_of_range_pvalue_flagged(self):
        """RDP emitted 36.36, which is not a probability and must be flagged."""
        rdp = [e for e in self.events if e["method"] == "rdp"][0]
        self.assertTrue(rdp["pvalue_out_of_range"])
        geneconv = [e for e in self.events if e["method"] == "geneconv"][0]
        self.assertFalse(geneconv["pvalue_out_of_range"])


class TestOpenRdpResolution(unittest.TestCase):
    def test_console_script_uses_adjacent_python_when_present(self):
        with tempfile.TemporaryDirectory() as td:
            bindir = Path(td)
            exe = bindir / "openrdp"
            py = bindir / "python"
            exe.write_text("#!/usr/bin/env python3\n")
            py.write_text("#!/bin/sh\n")
            exe.chmod(0o755)
            py.chmod(0o755)

            with patch.dict(os.environ, {"PATH": str(bindir)}):
                self.assertEqual(resolve_openrdp(), [str(py), str(exe)])


class TestConsensus(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(__file__).parent / "_tmp_openrdp2.csv"
        self.tmp.write_text(IS76_CSV)
        self.events = parse_openrdp_csv(self.tmp)

    def tearDown(self):
        self.tmp.unlink(missing_ok=True)

    @staticmethod
    def _is76(cons):
        return [c for c in cons if "IS76_LN812978" in c["members"].split(";")]

    def test_is76_recovered_by_geneconv(self):
        """The positive control must be recovered -- by GENECONV alone.

        Grouping on an exact (recombinant, major, minor) tuple split GENECONV's
        one-parent call from RDP's two-parent call and returned zero events, so
        the grouping fix is still required to get here.

        RDP does NOT support this event: its reported value of 36.36 is a
        Bonferroni-type corrected statistic (G * L/N * binom.sf) whose survival
        term is EXACTLY 1.0 -- on this triplet L/N = 2800/77 = 36.3636..., so
        the reported value is the correction factor alone and the binomial test
        contributes no evidence whatsoever. Requiring two methods here would be
        requiring a detection that the data does not contain.
        """
        cons = build_consensus(self.events, min_methods=1, pvalue_max=0.05)
        hits = self._is76(cons)
        self.assertEqual(len(hits), 1, f"IS76 not recovered cleanly: {cons}")
        self.assertIn("geneconv", hits[0]["methods"])
        self.assertEqual(hits[0]["methods"], "geneconv")
        self.assertEqual(
            hits[0]["members"],
            "IS76_LN812978;TYLCSV_NC_003828;TYLCV_IL_AM409201",
        )

    def test_parentless_rows_are_not_consensus_events(self):
        cons = build_consensus(self.events, min_methods=1, pvalue_max=0.05)
        self.assertTrue(cons, "control event should still be retained")
        for row in cons:
            self.assertGreaterEqual(len(row["members"].split(";")), 2)

    def test_paired_geneconv_fragments_resolve_triplet_roles(self):
        events = [
            {"recombinant": "REC", "parent_a": "DON", "parent_b": "",
             "start": 20, "end": 30, "method": "geneconv", "pvalue": 0.0,
             "pvalue_out_of_range": False},
            {"recombinant": "REC", "parent_a": "BACK", "parent_b": "",
             "start": 0, "end": 20, "method": "geneconv", "pvalue": 0.01,
             "pvalue_out_of_range": False},
            {"recombinant": "BACK", "parent_a": "", "parent_b": "",
             "start": 20, "end": 30, "method": "geneconv", "pvalue": 0.0,
             "pvalue_out_of_range": False},
        ]
        cons = build_consensus(events, min_methods=1, pvalue_max=0.05)
        self.assertEqual(len(cons), 1)
        self.assertEqual(cons[0]["members"], "BACK;DON;REC")

        alignment = {
            "BACK": "A" * 40,
            "DON": "C" * 40,
            "REC": "A" * 20 + "C" * 10 + "A" * 10,
        }
        roles = assign_roles(alignment, cons[0]["members"].split(";"),
                             cons[0]["start_min"], cons[0]["end_max"])
        self.assertEqual(roles["recombinant_resolved"], "REC")
        self.assertEqual(roles["backbone_parent"], "BACK")
        self.assertEqual(roles["donor_parent"], "DON")

    def test_out_of_range_pvalue_is_not_significant(self):
        """A Bonferroni-corrected value above 1 means NOT significant.

        RDP's `pvalue = G * (L/N) * binom.sf(...)` (openrdp/rdp.py:178) is not
        capped at 1; on a 61-genome panel 79% of RDP rows exceed 1. Treating
        those as exempt from the threshold -- on the grounds that a value above
        1 "is not a probability" -- inverts the statistic and admits the LEAST
        significant calls as though they were the most significant. They must
        fail the threshold.
        """
        cons = build_consensus(self.events, min_methods=1, pvalue_max=0.05)
        for c in cons:
            self.assertNotIn("rdp", c["methods"],
                             "RDP support at p=36.36 must not pass a 0.05 threshold")

    def test_out_of_range_pvalue_admitted_when_threshold_permits(self):
        """Clamping to 1 must not hard-drop the row: at pvalue_max=1 it counts."""
        cons = build_consensus(self.events, min_methods=1, pvalue_max=1.0)
        self.assertTrue(any("rdp" in c["methods"] for c in cons))

    def test_breakpoint_tolerance_merges_below_overlap_floor(self):
        """GENECONV 1382-1492 and RDP 1343-1420 reciprocally overlap 0.4935.

        A 0.5 reciprocal-overlap rule alone rejects the real event, so the
        breakpoint-tolerance path must carry it.
        """
        cons = build_consensus(self.events, min_methods=2, pvalue_max=1.0,
                               min_overlap=0.99, breakpoint_tolerance=100)
        self.assertEqual(len(self._is76(cons)), 1)

    def test_distant_intervals_do_not_merge(self):
        """Tolerance must not fuse unrelated events."""
        cons = build_consensus(self.events, min_methods=2, pvalue_max=1.0,
                               min_overlap=0.5, breakpoint_tolerance=5)
        self.assertEqual(len(self._is76(cons)), 0)

    def test_grouping_survives_permuted_recombinant_label(self):
        """OpenRDP's Recombinant column depends on WHICH METHODS were run.

        Verbatim: running `-m rdp` alone labels TYLCSV the recombinant, while all
        six methods label IS76 -- same triplet, same breakpoints, roles permuted
        (deterministic across seeds 3/7/42). Grouping must key on the set of
        sequences involved, not on that label.
        """
        permuted = IS76_CSV.replace(
            "Rdp,1343,1420,IS76_LN812978,TYLCSV_NC_003828,TYLCV_IL_AM409201",
            "Rdp,1343,1420,TYLCSV_NC_003828,IS76_LN812978,TYLCV_IL_AM409201")
        self.assertNotEqual(permuted, IS76_CSV)
        tmp = Path(__file__).parent / "_tmp_perm.csv"
        tmp.write_text(permuted)
        try:
            cons = build_consensus(parse_openrdp_csv(tmp), min_methods=2, pvalue_max=1.0)
        finally:
            tmp.unlink(missing_ok=True)
        hits = self._is76(cons)
        self.assertEqual(len(hits), 1, f"permuted label broke grouping: {cons}")
        self.assertEqual(hits[0]["methods"], "geneconv,rdp")

    def test_min_methods_three_rejects(self):
        cons = build_consensus(self.events, min_methods=3, pvalue_max=0.05)
        self.assertEqual(cons, [])


class TestMergeBound(unittest.TestCase):
    """The consensus merge must not chain transitively across a panel.

    Single-linkage grouping with no bound joins A~B and B~C into one group even
    when A and C share no interval. On the 61-genome panel this collapsed
    142,985 raw events into 8 groups, the largest holding 59 of 61 sequences and
    spanning the whole alignment. A recombination event involves a recombinant
    and two parents, so a group of 59 is not an event.

    This failure mode CANNOT appear on the 3-sequence positive control -- it
    needs a dense panel -- which is why these tests build a synthetic chain.
    """

    @staticmethod
    def _chain_csv(n=8, start=1000, end=1200):
        """A~B, B~C, C~D ... all at identical breakpoints.

        Each row shares exactly two sequences with its neighbour, so every
        adjacent pair is merge-eligible. Unbounded single linkage joins the whole
        chain; a triplet bound stops at three sequences per group.
        """
        rows = ["Method,Start,End,Recombinant,Parent1,Parent2,Pvalue"]
        for i in range(n - 2):
            a, b, c = f"S{i}", f"S{i+1}", f"S{i+2}"
            rows.append(f"Geneconv,{start},{end},{a},{b},{c},0.001")
            rows.append(f"Rdp,{start},{end},{a},{b},{c},0.001")
        return "\n".join(rows) + "\n"

    def setUp(self):
        self.tmp = Path(__file__).parent / "_tmp_chain.csv"
        self.tmp.write_text(self._chain_csv())
        self.events = parse_openrdp_csv(self.tmp)

    def tearDown(self):
        self.tmp.unlink(missing_ok=True)

    def test_no_group_exceeds_a_triplet(self):
        cons = build_consensus(self.events, min_methods=1, pvalue_max=0.05)
        biggest = max(len(c["members"].split(";")) for c in cons)
        self.assertLessEqual(
            biggest, 3,
            f"merge chained to {biggest} sequences; an event is a triplet")

    def test_chain_does_not_collapse_to_one_event(self):
        """Six overlapping triplets must not become a single event."""
        cons = build_consensus(self.events, min_methods=1, pvalue_max=0.05)
        self.assertGreater(
            len(cons), 1,
            "the whole chain collapsed into one group -- transitive chaining")

    def test_raising_the_bound_restores_chaining(self):
        """Confirms the bound is what prevents the collapse, not some other gate."""
        cons = build_consensus(self.events, min_methods=1, pvalue_max=0.05,
                               max_members=99)
        biggest = max(len(c["members"].split(";")) for c in cons)
        self.assertGreater(
            biggest, 3,
            "unbounded merge should chain; if it does not, this test no longer "
            "exercises the bound")

    def test_bound_does_not_break_the_control(self):
        """The bound must not cost the positive control.

        The control is a single triplet, so the bound is inactive on it. A change
        that satisfies the panel-scale tests above but loses IS76 is not a fix.
        """
        tmp = Path(__file__).parent / "_tmp_bound_is76.csv"
        tmp.write_text(IS76_CSV)
        try:
            cons = build_consensus(parse_openrdp_csv(tmp), min_methods=1,
                                   pvalue_max=0.05)
            self.assertTrue(
                [c for c in cons if "IS76_LN812978" in c["members"].split(";")])
        finally:
            tmp.unlink(missing_ok=True)


class TestPolarisation(unittest.TestCase):
    """Parent1/Parent2 are positional; the backbone comes from informative sites."""

    def _aln(self):
        # 40 columns: recombinant follows P_BACK everywhere except cols 20-30,
        # where it follows P_DON.
        back = "A" * 20 + "A" * 10 + "A" * 10
        don = "C" * 20 + "C" * 10 + "C" * 10
        rec = "A" * 20 + "C" * 10 + "A" * 10
        return {"REC": rec, "P_BACK": back, "P_DON": don}

    def test_backbone_and_donor_assigned(self):
        out = polarise_event(self._aln(), "REC", ["P_DON", "P_BACK"], 20, 30)
        self.assertEqual(out["backbone_parent"], "P_BACK")
        self.assertEqual(out["donor_parent"], "P_DON")

    def test_order_of_parents_does_not_matter(self):
        a = polarise_event(self._aln(), "REC", ["P_DON", "P_BACK"], 20, 30)
        b = polarise_event(self._aln(), "REC", ["P_BACK", "P_DON"], 20, 30)
        self.assertEqual(a["backbone_parent"], b["backbone_parent"])
        self.assertEqual(a["donor_parent"], b["donor_parent"])

    def test_no_inversion_yields_no_call(self):
        """If the recombinant follows one parent throughout, it is not a recombinant."""
        aln = {"REC": "A" * 40, "P_BACK": "A" * 40, "P_DON": "C" * 40}
        out = polarise_event(aln, "REC", ["P_BACK", "P_DON"], 20, 30)
        self.assertEqual(out["backbone_parent"], "")
        self.assertEqual(out["donor_parent"], "")

    def test_missing_parent_is_safe(self):
        out = polarise_event(self._aln(), "REC", ["P_BACK", ""], 20, 30)
        self.assertEqual(out["backbone_parent"], "")


class TestReadAlignment(unittest.TestCase):
    def test_gaps_are_preserved(self):
        """io_utils.load_fasta_records strips gaps, breaking breakpoint coordinates."""
        p = Path(__file__).parent / "_tmp_aln.fasta"
        p.write_text(">a\nAC--GT\n>b\nACTTGT\n")
        try:
            aln = read_alignment(p)
            self.assertEqual(aln["a"], "AC--GT")
            self.assertEqual(len(aln["a"]), len(aln["b"]))
        finally:
            p.unlink(missing_ok=True)

    def test_multiline_and_case(self):
        p = Path(__file__).parent / "_tmp_aln2.fasta"
        p.write_text(">a desc here\nac\ngt\n")
        try:
            self.assertEqual(read_alignment(p)["a"], "ACGT")
        finally:
            p.unlink(missing_ok=True)


class TestSanitiseForOpenRDP(unittest.TestCase):
    """OpenRDP rejects anything outside A,T,G,C,-,N; real panels contain R/Y/S/W/K."""

    def setUp(self):
        self.out = Path(__file__).parent / "_tmp_sane.fasta"

    def tearDown(self):
        self.out.unlink(missing_ok=True)

    def test_ambiguity_recoded_and_counted(self):
        aln = {"a": "ACGTRYSWK", "b": "ACGTACGTA"}
        rep = sanitise_alignment_for_openrdp(aln, self.out)
        self.assertEqual(rep["n_sites_recoded"], 5)
        self.assertEqual(rep["n_records_affected"], 1)
        written = read_alignment(self.out)
        self.assertEqual(written["a"], "ACGTNNNNN")
        self.assertEqual(written["b"], "ACGTACGTA")

    def test_output_is_openrdp_safe(self):
        aln = {"a": "ACGT-NRYSWKMBDHV"}
        sanitise_alignment_for_openrdp(aln, self.out)
        self.assertEqual(set(read_alignment(self.out)["a"]) - set("ATGCN-"), set())

    def test_column_coordinates_preserved(self):
        """Breakpoints index into columns, so length must not change."""
        aln = {"a": "AC-RT", "b": "ACGTT"}
        sanitise_alignment_for_openrdp(aln, self.out)
        written = read_alignment(self.out)
        self.assertEqual(len(written["a"]), 5)
        self.assertEqual(written["a"], "AC-NT")

    def test_clean_alignment_reports_zero(self):
        rep = sanitise_alignment_for_openrdp({"a": "ACGTN-"}, self.out)
        self.assertEqual(rep["n_sites_recoded"], 0)
        self.assertEqual(rep["n_records_affected"], 0)


if __name__ == "__main__":
    unittest.main()
