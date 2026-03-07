"""
Tests for all read_* functions exposed in gwaslab.__init__.

Covers:
  - read_ldsc, read_popcorn, read_greml  (io_read_ldsc)
  - read_tabular, read_bim, read_fam, read_psam, read_pvar, read_bgen_sample (io_read_tabular)
  - read_gtf, read_gtf_file  (io_gtf)
  - read_bed  (io_ucsc_bed)
  - read_bedpe  (io_bedpe)
  - read_bigwig, read_bigwig_intervals, read_bigwig_stats, read_bigbed  (io_bigwig_bigbed)
"""

import os
import sys
import gzip
import unittest
import tempfile
import shutil

import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Setup path so that `gwaslab` can be imported from source tree
# ---------------------------------------------------------------------------
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

# ---------------------------------------------------------------------------
# Imports under test (public API via __init__)
# ---------------------------------------------------------------------------
import gwaslab as gl

# Also import directly from submodules for lower-level tests
from gwaslab.io.io_read_ldsc import read_ldsc, read_popcorn, read_greml
from gwaslab.io.io_read_tabular import (
    _read_tabular as read_tabular,
    read_bim,
    read_fam,
    read_psam,
    read_pvar,
    read_bgen_sample,
)
from gwaslab.io.io_ucsc_bed import read_bed
from gwaslab.io.io_bedpe import read_bedpe
from gwaslab.io.io_gtf import read_gtf, read_gtf_file

# Optional bigwig/bigbed (requires pyBigWig)
try:
    from gwaslab.io.io_bigwig_bigbed import (
        read_bigwig,
        read_bigwig_intervals,
        read_bigwig_stats,
        read_bigbed,
    )
    import pyBigWig

    HAS_PYBIGWIG = True
except ImportError:
    HAS_PYBIGWIG = False


# ===========================================================================
# Helpers
# ===========================================================================
REF_DIR = os.path.join(ROOT, "test", "ref")
LDSC_DIR = os.path.join(REF_DIR, "ldsc_logs")

# Reference PLINK files from 1000 Genomes
BIM_FILE = os.path.join(REF_DIR, "1kg_eas_hg19.chr7_126253550_128253550.bim")
FAM_FILE = os.path.join(REF_DIR, "1kg_eas_hg19.chr7_126253550_128253550.fam")
PSAM_FILE = os.path.join(REF_DIR, "1kg_eas_hg19.chr7_126253550_128253550.psam")
PVAR_FILE = os.path.join(REF_DIR, "1kg_eas_hg19.chr7_126253550_128253550.pvar")


# ===========================================================================
# Test read_ldsc / read_popcorn / read_greml
# ===========================================================================
class TestReadLDSC(unittest.TestCase):
    """Tests for read_ldsc (h2 and rg modes)."""

    def test_read_ldsc_available_in_init(self):
        """read_ldsc should be importable from gwaslab top-level."""
        self.assertTrue(hasattr(gl, "read_ldsc"))

    def test_h2_single_file(self):
        path = os.path.join(LDSC_DIR, "simulated_h2_1.log")
        df = read_ldsc([path], mode="h2")
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 1)
        self.assertIn("h2_obs", df.columns)
        self.assertEqual(df.iloc[0]["h2_obs"], "0.1583")

    def test_h2_multiple_files(self):
        paths = [
            os.path.join(LDSC_DIR, "simulated_h2_1.log"),
            os.path.join(LDSC_DIR, "simulated_h2_2.log"),
        ]
        df = read_ldsc(paths, mode="h2")
        self.assertEqual(len(df), 2)

    def test_rg_single_pair(self):
        path = os.path.join(LDSC_DIR, "simulated_rg.log")
        df = read_ldsc([path], mode="rg")
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 1)
        self.assertIn("rg", df.columns)
        self.assertAlmostEqual(float(df.iloc[0]["rg"]), 0.1601, places=4)

    def test_rg_multiple_pairs(self):
        path = os.path.join(LDSC_DIR, "simulated_rg_multiple.log")
        df = read_ldsc([path], mode="rg")
        self.assertEqual(len(df), 3)

    def test_empty_filelist(self):
        df = read_ldsc([], mode="h2")
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 0)


class TestReadPopcorn(unittest.TestCase):
    """Tests for read_popcorn."""

    def test_available_in_init(self):
        self.assertTrue(hasattr(gl, "read_popcorn"))

    def test_empty_filelist(self):
        df = read_popcorn([])
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 0)
        self.assertIn("pg", df.columns)
        self.assertIn("sfile1", df.columns)


class TestReadGREML(unittest.TestCase):
    """Tests for read_greml."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create a simulated GREML .hsq output
        self.greml_file = os.path.join(self.tmpdir, "test.hsq")
        content = (
            "Source\tVariance\tSE\n"
            "V(G)\t0.1234\t0.0100\n"
            "V(e)\t0.8766\t0.0100\n"
            "Vp\t1.0000\t0.0050\n"
            "Sum of V(G)/Vp\t0.1234\t0.0100\n"
            "Pval\t0.00050\n"
            "n\t50000\n"
        )
        with open(self.greml_file, "w") as f:
            f.write(content)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_available_in_init(self):
        self.assertTrue(hasattr(gl, "read_greml"))

    def test_read_single_greml_file(self):
        df = read_greml([self.greml_file])
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["Sum of V(G)/Vp"], "0.1234")
        self.assertEqual(df.iloc[0]["SE"], "0.0100")
        self.assertEqual(df.iloc[0]["Pval"], "0.00050")
        self.assertEqual(df.iloc[0]["n"], "50000")

    def test_empty_filelist(self):
        df = read_greml([])
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 0)
        self.assertIn("Sum of V(G)/Vp", df.columns)


# ===========================================================================
# Test read_bim / read_fam / read_psam / read_pvar / read_bgen_sample
# ===========================================================================
class TestReadBim(unittest.TestCase):
    """Tests for read_bim (PLINK .bim format)."""

    def test_available_in_init(self):
        self.assertTrue(hasattr(gl, "read_bim"))

    @unittest.skipUnless(os.path.exists(BIM_FILE), "bim reference file missing")
    def test_read_bim_returns_dataframe(self):
        df = read_bim(BIM_FILE)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)

    @unittest.skipUnless(os.path.exists(BIM_FILE), "bim reference file missing")
    def test_read_bim_columns(self):
        df = read_bim(BIM_FILE)
        # PLINK .bim has 6 columns: CHR, SNPID, CM, POS, EA, NEA
        self.assertEqual(df.shape[1], 6)

    @unittest.skipUnless(os.path.exists(BIM_FILE), "bim reference file missing")
    def test_read_bim_first_row(self):
        df = read_bim(BIM_FILE)
        # First row should be chromosome 7
        first_row = df.iloc[0]
        self.assertEqual(int(first_row.iloc[0]), 7)  # CHR


class TestReadFam(unittest.TestCase):
    """Tests for read_fam (PLINK .fam format)."""

    def test_available_in_init(self):
        self.assertTrue(hasattr(gl, "read_fam"))

    @unittest.skipUnless(os.path.exists(FAM_FILE), "fam reference file missing")
    def test_read_fam_returns_dataframe(self):
        df = read_fam(FAM_FILE)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)

    @unittest.skipUnless(os.path.exists(FAM_FILE), "fam reference file missing")
    def test_read_fam_columns(self):
        df = read_fam(FAM_FILE)
        # PLINK .fam has 6 columns: FID, IID, FATHER, MOTHER, SEX, PHENO
        self.assertEqual(df.shape[1], 6)


class TestReadPsam(unittest.TestCase):
    """Tests for read_psam (PLINK2 .psam format)."""

    def test_available_in_init(self):
        self.assertTrue(hasattr(gl, "read_psam"))

    @unittest.skipUnless(os.path.exists(PSAM_FILE), "psam reference file missing")
    @unittest.expectedFailure  # formatbook expects SID/PAT/MAT columns not present in test .psam
    def test_read_psam_returns_dataframe(self):
        df = read_psam(PSAM_FILE)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)


class TestReadPvar(unittest.TestCase):
    """Tests for read_pvar (PLINK2 .pvar format)."""

    def test_available_in_init(self):
        self.assertTrue(hasattr(gl, "read_pvar"))

    @unittest.skipUnless(os.path.exists(PVAR_FILE), "pvar reference file missing")
    def test_read_pvar_returns_dataframe(self):
        df = read_pvar(PVAR_FILE)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)


class TestReadBgenSample(unittest.TestCase):
    """Tests for read_bgen_sample."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create a minimal .sample file (BGEN format)
        self.sample_file = os.path.join(self.tmpdir, "test.sample")
        content = (
            "ID_1 ID_2 missing sex\n"
            "0 0 0 D\n"
            "SAMPLE1 SAMPLE1 0 1\n"
            "SAMPLE2 SAMPLE2 0 2\n"
            "SAMPLE3 SAMPLE3 0 1\n"
        )
        with open(self.sample_file, "w") as f:
            f.write(content)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_available_in_init(self):
        self.assertTrue(hasattr(gl, "read_bgen_sample"))

    @unittest.expectedFailure  # 'bgen_sample' format not yet in formatbook
    def test_read_bgen_sample_returns_dataframe(self):
        df = read_bgen_sample(self.sample_file)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)


# ===========================================================================
# Test read_bed (UCSC BED format, NOT PLINK .bed)
# ===========================================================================
class TestReadBed(unittest.TestCase):
    """Tests for read_bed (UCSC BED format)."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def _write_bed(self, filename, lines):
        path = os.path.join(self.tmpdir, filename)
        with open(path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        return path

    def _write_bed_gz(self, filename, lines):
        path = os.path.join(self.tmpdir, filename)
        with gzip.open(path, "wt") as f:
            for line in lines:
                f.write(line + "\n")
        return path

    def test_available_in_init(self):
        self.assertTrue(hasattr(gl, "read_bed"))

    def test_read_bed3(self):
        """Read a minimal BED3 file."""
        path = self._write_bed("test.bed", [
            "chr1\t0\t1000",
            "chr1\t2000\t3000",
            "chr2\t5000\t6000",
        ])
        df = read_bed(path)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 3)
        self.assertListEqual(list(df.columns[:3]), ["chrom", "chromStart", "chromEnd"])

    def test_read_bed6(self):
        """Read a BED6 file with name, score, strand."""
        path = self._write_bed("test6.bed", [
            "chr1\t0\t1000\tgene_A\t100\t+",
            "chr2\t5000\t6000\tgene_B\t200\t-",
        ])
        df = read_bed(path)
        self.assertEqual(len(df), 2)
        self.assertIn("name", df.columns)
        self.assertIn("score", df.columns)
        self.assertIn("strand", df.columns)

    def test_read_bed_gzipped(self):
        """Read a gzipped BED file."""
        path = self._write_bed_gz("test.bed.gz", [
            "chr1\t0\t1000",
            "chr2\t5000\t6000",
        ])
        df = read_bed(path)
        self.assertEqual(len(df), 2)

    def test_read_bed_with_comments(self):
        """Comment lines starting with # should be skipped."""
        path = self._write_bed("commented.bed", [
            "# this is a comment",
            "chr1\t0\t1000",
            "chr2\t5000\t6000",
        ])
        df = read_bed(path)
        self.assertEqual(len(df), 2)

    def test_read_bed_coordinate_types(self):
        """chromStart and chromEnd should be integers."""
        path = self._write_bed("types.bed", [
            "chr1\t100\t200",
        ])
        df = read_bed(path)
        self.assertTrue(pd.api.types.is_integer_dtype(df["chromStart"]))
        self.assertTrue(pd.api.types.is_integer_dtype(df["chromEnd"]))

    def test_read_bed_file_not_found(self):
        """Should raise FileNotFoundError for missing file."""
        with self.assertRaises(FileNotFoundError):
            read_bed("/nonexistent/path/missing.bed")

    def test_read_bed_invalid_columns(self):
        """File with fewer than 3 columns should raise ValueError."""
        path = self._write_bed("invalid.bed", [
            "chr1\t0",
            "chr2\t5000",
        ])
        with self.assertRaises(ValueError):
            read_bed(path)

    def test_read_bed_usecols(self):
        """usecols should select specific columns."""
        path = self._write_bed("usecols.bed", [
            "chr1\t0\t1000\tgene_A\t100\t+",
        ])
        df = read_bed(path, usecols=[0, 1, 2])
        self.assertEqual(df.shape[1], 3)


# ===========================================================================
# Test read_bedpe
# ===========================================================================
class TestReadBedpe(unittest.TestCase):
    """Tests for read_bedpe (BEDPE format)."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def _write_bedpe(self, filename, lines):
        path = os.path.join(self.tmpdir, filename)
        with open(path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        return path

    def _write_bedpe_gz(self, filename, lines):
        path = os.path.join(self.tmpdir, filename)
        with gzip.open(path, "wt") as f:
            for line in lines:
                f.write(line + "\n")
        return path

    def test_available_in_init(self):
        self.assertTrue(hasattr(gl, "read_bedpe"))

    def test_read_bedpe6(self):
        """Read a minimal BEDPE6 file."""
        path = self._write_bedpe("test.bedpe", [
            "chr1\t0\t1000\tchr1\t5000\t6000",
            "chr2\t1000\t2000\tchr2\t8000\t9000",
        ])
        df = read_bedpe(path)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 2)
        self.assertListEqual(
            list(df.columns[:6]),
            ["chrom1", "chromStart1", "chromEnd1", "chrom2", "chromStart2", "chromEnd2"],
        )

    def test_read_bedpe10(self):
        """Read a full BEDPE10 file with optional fields."""
        path = self._write_bedpe("test10.bedpe", [
            "chr1\t0\t1000\tchr1\t5000\t6000\tinteraction1\t500\t+\t-",
            "chr2\t1000\t2000\tchr2\t8000\t9000\tinteraction2\t300\t+\t+",
        ])
        df = read_bedpe(path)
        self.assertEqual(len(df), 2)
        self.assertIn("name", df.columns)
        self.assertIn("score", df.columns)
        self.assertIn("strand1", df.columns)
        self.assertIn("strand2", df.columns)

    def test_read_bedpe_gzipped(self):
        """Read a gzipped BEDPE file."""
        path = self._write_bedpe_gz("test.bedpe.gz", [
            "chr1\t0\t1000\tchr1\t5000\t6000",
        ])
        df = read_bedpe(path)
        self.assertEqual(len(df), 1)

    def test_read_bedpe_with_comments(self):
        """Comment lines should be skipped."""
        path = self._write_bedpe("commented.bedpe", [
            "# header comment",
            "chr1\t0\t1000\tchr1\t5000\t6000",
        ])
        df = read_bedpe(path)
        self.assertEqual(len(df), 1)

    def test_read_bedpe_coordinate_types(self):
        """Coordinate columns should be integers."""
        path = self._write_bedpe("types.bedpe", [
            "chr1\t100\t200\tchr2\t300\t400",
        ])
        df = read_bedpe(path)
        for col in ["chromStart1", "chromEnd1", "chromStart2", "chromEnd2"]:
            self.assertTrue(pd.api.types.is_integer_dtype(df[col]),
                            f"{col} should be integer")

    def test_read_bedpe_file_not_found(self):
        """Should raise FileNotFoundError for missing file."""
        with self.assertRaises(FileNotFoundError):
            read_bedpe("/nonexistent/path/missing.bedpe")

    def test_read_bedpe_invalid_columns(self):
        """File with fewer than 6 columns should raise ValueError."""
        path = self._write_bedpe("invalid.bedpe", [
            "chr1\t0\t1000\tchr1\t5000",
        ])
        with self.assertRaises(ValueError):
            read_bedpe(path)


# ===========================================================================
# Test read_gtf / read_gtf_file
# ===========================================================================
class TestReadGtf(unittest.TestCase):
    """Tests for read_gtf and read_gtf_file."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def _make_gtf(self, filename, lines, compress=False):
        """Helper to create a GTF file."""
        path = os.path.join(self.tmpdir, filename)
        if compress:
            with gzip.open(path, "wt") as f:
                for line in lines:
                    f.write(line + "\n")
        else:
            with open(path, "w") as f:
                for line in lines:
                    f.write(line + "\n")
        return path

    SAMPLE_GTF_LINES = [
        '1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_biotype "transcribed_unprocessed_pseudogene";',
        '1\thavana\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1"; transcript_id "ENST00000456328"; gene_biotype "transcribed_unprocessed_pseudogene";',
        '1\thavana\texon\t11869\t12227\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1"; transcript_id "ENST00000456328"; gene_biotype "transcribed_unprocessed_pseudogene";',
        '2\tensembl\tgene\t38814\t46870\t.\t-\t.\tgene_id "ENSG00000227232"; gene_name "WASH7P"; gene_biotype "unprocessed_pseudogene";',
    ]

    def test_available_in_init(self):
        self.assertTrue(hasattr(gl, "read_gtf"))
        self.assertTrue(hasattr(gl, "read_gtf_file"))

    def test_read_gtf_basic(self):
        """Read a basic GTF file."""
        path = self._make_gtf("test.gtf", self.SAMPLE_GTF_LINES)
        df = read_gtf(path)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 4)
        self.assertIn("seqname", df.columns)
        self.assertIn("start", df.columns)
        self.assertIn("end", df.columns)

    def test_read_gtf_gzipped(self):
        """Read a gzipped GTF file."""
        path = self._make_gtf("test.gtf.gz", self.SAMPLE_GTF_LINES, compress=True)
        df = read_gtf(path)
        self.assertEqual(len(df), 4)

    def test_read_gtf_feature_filter(self):
        """Filter rows by feature type."""
        path = self._make_gtf("test.gtf", self.SAMPLE_GTF_LINES)
        df = read_gtf(path, features={"gene"})
        self.assertEqual(len(df), 2)
        self.assertTrue((df["feature"] == "gene").all())

    def test_read_gtf_chrom_filter(self):
        """Filter rows by chromosome."""
        path = self._make_gtf("test.gtf", self.SAMPLE_GTF_LINES)
        df = read_gtf(path, chrom="1")
        self.assertGreater(len(df), 0)
        # All rows should be chromosome 1
        self.assertTrue((df["seqname"] == "1").all())

    def test_read_gtf_usecols(self):
        """Select specific columns."""
        path = self._make_gtf("test.gtf", self.SAMPLE_GTF_LINES)
        df = read_gtf(path, usecols=["seqname", "start", "end", "gene_name"])
        self.assertIn("seqname", df.columns)
        self.assertIn("gene_name", df.columns)

    def test_read_gtf_attribute_expansion(self):
        """Attributes should be expanded into separate columns."""
        path = self._make_gtf("test.gtf", self.SAMPLE_GTF_LINES)
        df = read_gtf(path, expand_attribute_column=True)
        self.assertIn("gene_id", df.columns)
        self.assertIn("gene_name", df.columns)
        self.assertEqual(df.iloc[0]["gene_name"], "DDX11L1")

    def test_read_gtf_file_not_found(self):
        """Should raise ValueError for missing file."""
        with self.assertRaises(ValueError):
            read_gtf("/nonexistent/path/missing.gtf")

    def test_read_gtf_file_wrapper(self):
        """read_gtf_file should call read_gtf with standard usecols."""
        path = self._make_gtf("test.gtf", self.SAMPLE_GTF_LINES)
        df = read_gtf_file(path)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)
        # Should contain the standard columns
        for col in ["seqname", "start", "end", "strand", "feature", "gene_name"]:
            self.assertIn(col, df.columns)


# ===========================================================================
# Test read_bigwig / read_bigwig_intervals / read_bigwig_stats / read_bigbed
# ===========================================================================
@unittest.skipUnless(HAS_PYBIGWIG, "pyBigWig not installed")
class TestReadBigWig(unittest.TestCase):
    """Tests for read_bigwig, read_bigwig_intervals, read_bigwig_stats."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create a tiny bigWig file for testing
        self.bw_path = os.path.join(self.tmpdir, "test.bw")
        bw = pyBigWig.open(self.bw_path, "w")
        bw.addHeader([("chr1", 10000), ("chr2", 8000)])
        bw.addEntries(
            ["chr1", "chr1", "chr1"],
            [0, 100, 200],
            ends=[100, 200, 300],
            values=[1.0, 2.0, 3.0],
        )
        bw.close()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_available_in_init(self):
        self.assertTrue(hasattr(gl, "read_bigwig"))
        self.assertTrue(hasattr(gl, "read_bigwig_intervals"))
        self.assertTrue(hasattr(gl, "read_bigwig_stats"))

    def test_read_bigwig_header(self):
        """Without chrom/start/end, should return header info."""
        df = read_bigwig(self.bw_path)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertIn("chrom", df.columns)
        self.assertIn("size", df.columns)
        self.assertEqual(len(df), 2)

    def test_read_bigwig_region_dataframe(self):
        """With chrom/start/end, should return intervals as DataFrame."""
        df = read_bigwig(self.bw_path, chrom="chr1", start=0, end=300)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 3)
        self.assertIn("value", df.columns)
        self.assertAlmostEqual(df.iloc[0]["value"], 1.0)

    def test_read_bigwig_region_list(self):
        """With as_dataframe=False, should return a list."""
        values = read_bigwig(self.bw_path, chrom="chr1", start=0, end=10, as_dataframe=False)
        self.assertIsInstance(values, (list, np.ndarray))
        self.assertEqual(len(values), 10)

    def test_read_bigwig_file_not_found(self):
        """Should raise FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            read_bigwig("/nonexistent/test.bw")

    def test_read_bigwig_invalid_range(self):
        """start >= end should raise ValueError."""
        with self.assertRaises(ValueError):
            read_bigwig(self.bw_path, chrom="chr1", start=300, end=100)

    def test_read_bigwig_intervals(self):
        """read_bigwig_intervals should return intervals as DataFrame."""
        df = read_bigwig_intervals(self.bw_path, chrom="chr1", start=0, end=300)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 3)
        self.assertIn("value", df.columns)

    def test_read_bigwig_intervals_empty_region(self):
        """Region with no data should return empty DataFrame."""
        df = read_bigwig_intervals(self.bw_path, chrom="chr2", start=0, end=100)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 0)

    def test_read_bigwig_stats_mean(self):
        """read_bigwig_stats should return statistics."""
        result = read_bigwig_stats(self.bw_path, chrom="chr1", start=0, end=300, stat="mean")
        self.assertIsNotNone(result)

    def test_read_bigwig_stats_max(self):
        result = read_bigwig_stats(self.bw_path, chrom="chr1", start=0, end=300, stat="max")
        self.assertIsNotNone(result)

    def test_read_bigwig_stats_min(self):
        result = read_bigwig_stats(self.bw_path, chrom="chr1", start=0, end=300, stat="min")
        self.assertIsNotNone(result)


@unittest.skipUnless(HAS_PYBIGWIG, "pyBigWig not installed")
class TestReadBigBed(unittest.TestCase):
    """Tests for read_bigbed."""

    def test_available_in_init(self):
        self.assertTrue(hasattr(gl, "read_bigbed"))

    def test_read_bigbed_file_not_found(self):
        """Should raise FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            read_bigbed("/nonexistent/test.bb")


# ===========================================================================
# Test top-level import completeness
# ===========================================================================
class TestInitExportsAllReadFunctions(unittest.TestCase):
    """Verify all read_* functions are accessible via gwaslab top-level."""

    EXPECTED_READ_FUNCTIONS = [
        "read_ldsc",
        "read_popcorn",
        "read_greml",
        "read_tabular",
        "read_bim",
        "read_fam",
        "read_psam",
        "read_pvar",
        "read_bgen_sample",
        "read_gtf",
        "read_gtf_file",
        "read_bigwig",
        "read_bigwig_intervals",
        "read_bigwig_stats",
        "read_bigbed",
        "read_bed",
        "read_bedpe",
    ]

    def test_all_read_functions_in_namespace(self):
        for name in self.EXPECTED_READ_FUNCTIONS:
            with self.subTest(name=name):
                self.assertTrue(
                    hasattr(gl, name),
                    f"gwaslab.{name} should be accessible but is missing",
                )

    def test_all_are_callable(self):
        for name in self.EXPECTED_READ_FUNCTIONS:
            with self.subTest(name=name):
                self.assertTrue(
                    callable(getattr(gl, name)),
                    f"gwaslab.{name} should be callable",
                )


if __name__ == "__main__":
    unittest.main()
