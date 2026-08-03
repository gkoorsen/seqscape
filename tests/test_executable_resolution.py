"""Regression tests for aligner executable resolution.

These cover the defect reported in review, where aligner paths were resolved
through hardcoded absolute paths pointing at the author's local machine instead
of falling back to tools on the user's PATH.
"""

from __future__ import annotations

import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from seqscape.alignment import _is_usable_exe, resolve_mafft_exe, resolve_muscle_exe


def _make_executable(directory: Path, name: str) -> Path:
    exe = directory / name
    exe.write_text("#!/bin/sh\nexit 0\n")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return exe


class ExecutableResolutionTests(unittest.TestCase):
    """No machine-specific path may be baked into aligner resolution."""

    def test_no_absolute_paths_hardcoded_in_resolution_source(self) -> None:
        sources = [
            SRC / "seqscape" / "alignment.py",
            SRC / "seqscape" / "cli.py",
        ]
        for path in sources:
            source = path.read_text()
            for forbidden in ("/Users/", "/home/", "/opt/homebrew", "ciderseq-muscle"):
                self.assertNotIn(
                    forbidden,
                    source,
                    f"{path.name} must not hardcode a machine-specific path ({forbidden})",
                )

    def test_explicit_hint_wins(self) -> None:
        for resolver, name in ((resolve_muscle_exe, "muscle"), (resolve_mafft_exe, "mafft")):
            with tempfile.TemporaryDirectory() as tmp:
                exe = _make_executable(Path(tmp), name)
                self.assertEqual(resolver(str(exe)), str(exe.resolve()))

    def test_falls_back_to_path(self) -> None:
        """A bare tool name must resolve via PATH, not a hardcoded location."""
        for resolver, name in ((resolve_muscle_exe, "muscle"), (resolve_mafft_exe, "mafft")):
            with tempfile.TemporaryDirectory() as tmp:
                exe = _make_executable(Path(tmp), name)
                original = os.environ.get("PATH", "")
                os.environ["PATH"] = tmp
                try:
                    self.assertEqual(Path(resolver(name)).resolve(), exe.resolve())
                finally:
                    os.environ["PATH"] = original

    def test_nonexistent_hint_does_not_raise(self) -> None:
        """An unusable hint is skipped, not fatal; a name is returned for the caller."""
        for resolver, name in ((resolve_muscle_exe, "muscle"), (resolve_mafft_exe, "mafft")):
            with tempfile.TemporaryDirectory() as tmp:
                missing = str(Path(tmp) / "definitely" / "not" / "here" / name)
                original = os.environ.get("PATH", "")
                os.environ["PATH"] = tmp
                try:
                    resolved = resolver(missing)
                finally:
                    os.environ["PATH"] = original
                self.assertIsInstance(resolved, str)
                self.assertTrue(resolved)

    def test_unusable_candidates_are_skipped(self) -> None:
        """A directory and a non-executable file are both rejected, without raising."""
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            self.assertFalse(_is_usable_exe(tmpdir))

            not_executable = tmpdir / "muscle"
            not_executable.write_text("#!/bin/sh\nexit 0\n")
            not_executable.chmod(stat.S_IRUSR | stat.S_IWUSR)
            self.assertFalse(_is_usable_exe(not_executable))

            self.assertFalse(_is_usable_exe(tmpdir / "does-not-exist"))


if __name__ == "__main__":
    unittest.main()
