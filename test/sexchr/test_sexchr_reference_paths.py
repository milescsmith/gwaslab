"""
Chr X / chr Y coverage for Sumstats paths that read reference files (FASTA, VCF).

Assets are generated in a temp dir by ``simulate_ref.py`` (small chrX/chrY FASTA +
tabix-indexed VCF + TSV sumstats).
"""
from __future__ import annotations

import importlib.util
import os
import shutil
import sys
import tempfile
import unittest

_HAVE_BCFTOOLS = shutil.which("bcftools") is not None

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab import Sumstats
from gwaslab.bd.bd_common_data import get_chr_list
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from gwaslab.io.io_fasta import load_and_build_fasta_records
from gwaslab.info.g_Log import Log

_HERE = os.path.dirname(os.path.abspath(__file__))
_mod_name = "gwaslab_sexchr_simulate_ref"
_spec = importlib.util.spec_from_file_location(
    _mod_name,
    os.path.join(_HERE, "simulate_ref.py"),
)
_sim = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
sys.modules[_mod_name] = _sim
_spec.loader.exec_module(_sim)
XY_VARIANTS = _sim.XY_VARIANTS
build_xy_reference = _sim.build_xy_reference
bundled_sumstats_path = _sim.bundled_sumstats_path


def _ref_align_digit(status) -> int:
    try:
        v = int(status)
    except (TypeError, ValueError):
        return -1
    return (v // 10) % 10


try:
    import pysam  # noqa: F401
except ImportError:
    pysam = None  # type: ignore


@unittest.skipUnless(pysam is not None, "pysam required for simulated VCF tabix and reference tests")
class TestSexChrReferencePaths(unittest.TestCase):
    """Reference-backed harmonization / QC on simulated chrX and chrY."""

    @classmethod
    def setUpClass(cls):
        cls._tmpdir = tempfile.mkdtemp(prefix="gwaslab_sexchr_")
        cls.paths = build_xy_reference(cls._tmpdir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._tmpdir, ignore_errors=True)

    def _sumstats(self) -> Sumstats:
        return Sumstats(
            self.paths.sumstats_tsv,
            fmt="auto",
            build="38",
            verbose=False,
        )

    def test_bundled_sumstats_matches_simulator(self):
        self.assertTrue(os.path.isfile(bundled_sumstats_path()))
        s1 = Sumstats(self.paths.sumstats_tsv, fmt="auto", build="38", verbose=False)
        s2 = Sumstats(bundled_sumstats_path(), fmt="auto", build="38", verbose=False)
        self.assertEqual(len(s1.data), len(s2.data))
        self.assertListEqual(s1.data["POS"].tolist(), s2.data["POS"].tolist())

    def test_get_chr_list_includes_numeric_xy_for_fasta_filter(self):
        """Regression: FASTA loader maps chrX→23; chromlist must include 23/24/25."""
        s = set(get_chr_list(add_number=True))
        self.assertIn(23, s)
        self.assertIn(24, s)
        self.assertIn(25, s)

    def test_load_and_build_fasta_records_chrX_chrY(self):
        mapper = ChromosomeMapper(species="homo sapiens", build="38", verbose=False)
        mapper.detect_sumstats_format(__import__("pandas").Series([23, 24]))
        chromlist_set = set(get_chr_list(add_number=True))
        chroms_in = {23, 24}
        log = Log()
        log.verbose = False
        _rec, _starts, lens = load_and_build_fasta_records(
            self.paths.fasta_path,
            chromlist_set,
            chroms_in,
            mapper=mapper,
            log=log,
            verbose=False,
        )
        self.assertIn(23, lens)
        self.assertIn(24, lens)

    def test_check_ref_chrX_chrY(self):
        gl = self._sumstats()
        gl.fix_chr(verbose=False)
        gl.fix_pos(verbose=False)
        gl.check_ref(self.paths.fasta_path, verbose=False)
        d6 = gl.data["STATUS"].map(_ref_align_digit)
        self.assertEqual((d6 == 8).sum(), 0, "no variant should be 'not on reference'")
        self.assertTrue((d6 == 0).all() or ((d6 == 0) | (d6 == 3)).all())

    def test_harmonize_ref_seq_only(self):
        gl = self._sumstats()
        gl.harmonize(basic_check=True, ref_seq=self.paths.fasta_path, verbose=False)
        d6 = gl.data["STATUS"].map(_ref_align_digit)
        self.assertEqual((d6 == 8).sum(), 0)

    @unittest.skipUnless(_HAVE_BCFTOOLS, "harmonize ref_rsid_vcf uses bcftools")
    def test_harmonize_ref_seq_rsid_vcf_and_infer(self):
        gl = self._sumstats()
        gl.harmonize(
            basic_check=True,
            ref_seq=self.paths.fasta_path,
            ref_rsid_vcf=self.paths.vcf_gz,
            ref_infer=self.paths.vcf_gz,
            ref_alt_freq="AF",
            threads=1,
            verbose=False,
        )
        self.assertTrue(gl.meta.get("is_harmonised", False))
        self.assertIn("rsID", gl.data.columns)
        d6 = gl.data["STATUS"].map(_ref_align_digit)
        self.assertEqual((d6 == 8).sum(), 0)

    @unittest.skipUnless(_HAVE_BCFTOOLS, "assign_rsid VCF path uses bcftools")
    def test_assign_rsid_chrX_chrY(self):
        gl = self._sumstats()
        gl.basic_check(verbose=False)
        gl.assign_rsid(ref_rsid_vcf=self.paths.vcf_gz, threads=1, verbose=False)
        ids = set(gl.data["rsID"].astype(str))
        for _c, _p, _nea, _ea, _eaf, rs, _v in XY_VARIANTS:
            self.assertIn(rs, ids)

    def test_check_af_chrX_chrY(self):
        gl = self._sumstats()
        gl.harmonize(
            basic_check=True,
            ref_seq=self.paths.fasta_path,
            ref_infer=self.paths.vcf_gz,
            ref_alt_freq="AF",
            threads=1,
            verbose=False,
        )
        gl.check_af(ref_infer=self.paths.vcf_gz, ref_alt_freq="AF", verbose=False)
        self.assertIn("DAF", gl.data.columns)

    def test_infer_af_chrX_chrY(self):
        gl = self._sumstats()
        gl.harmonize(
            basic_check=True,
            ref_seq=self.paths.fasta_path,
            verbose=False,
        )
        gl.infer_af(ref_infer=self.paths.vcf_gz, ref_alt_freq="AF", verbose=False)
        self.assertTrue(gl.data["EAF"].notna().any())

    def test_infer_eaf_from_maf_chrX_chrY(self):
        gl = self._sumstats()
        gl.harmonize(
            basic_check=True,
            ref_seq=self.paths.fasta_path,
            verbose=False,
        )
        gl.infer_eaf_from_maf(ref_infer=self.paths.vcf_gz, ref_alt_freq="AF", verbose=False)
        self.assertTrue(gl.data["EAF"].notna().all())

    def test_infer_af2_annotation_chrX_chrY(self):
        gl = self._sumstats()
        gl.harmonize(
            basic_check=True,
            ref_seq=self.paths.fasta_path,
            verbose=False,
        )
        gl.infer_af2(path=self.paths.vcf_gz, assign_cols="AF", verbose=False)
        self.assertTrue(gl.data["EAF"].notna().any())

    @unittest.skipUnless(_HAVE_BCFTOOLS, "check_af2 VCF lookup uses bcftools")
    def test_check_af2_chrX_chrY_numeric_chr_lookup(self):
        """Lookup TSV has chrX/chrY strings; sumstats use numeric CHR after fix_chr — DAF is computed."""
        gl = self._sumstats()
        gl.basic_check(verbose=False)
        gl.fix_chr(verbose=False)
        gl.harmonize(
            basic_check=True,
            ref_seq=self.paths.fasta_path,
            verbose=False,
        )
        gl.check_af2(path=self.paths.vcf_gz, ref_alt_freq="AF", verbose=False, force=True)
        self.assertIn("DAF", gl.data.columns)


class TestSexChrRsidToChrposPlaceholder(unittest.TestCase):
    """Document reference paths not covered by the minimal XY fixture."""

    def test_rsid_to_chrpos_needs_hdf5_mod10_files(self):
        self.skipTest(
            "rsid_to_chrpos / rsid_to_chrpos2 expect process_vcf_to_hfd5() HDF5 layout "
            "beside the VCF; add a fixture if we need automated XY coverage for that path."
        )

if __name__ == "__main__":
    unittest.main()
