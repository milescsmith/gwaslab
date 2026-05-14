"""
Sex-chromosome tests for reference-style file loading (chain, BED, BEDPE, GTF) and
ChromosomeMapper notation (23/24, X/Y, chrX/chrY, NC_*). NC strings are hard-coded from
``src/gwaslab/data/chromosomes/chromosomes_nc.json`` via :mod:`nc_constants`.
"""
from __future__ import annotations

import gzip
import importlib.util
import os
import shutil
import sys
import tempfile
import unittest

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)


def _load_nc_constants():
    """Load ``nc_constants.py`` from this directory (works with unittest discover)."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nc_constants.py")
    name = "gwaslab_test_sexchr_nc_constants"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_nc = _load_nc_constants()
NC_HUMAN_19_X = _nc.NC_HUMAN_19_X
NC_HUMAN_19_Y = _nc.NC_HUMAN_19_Y
NC_HUMAN_38_X = _nc.NC_HUMAN_38_X
NC_HUMAN_38_Y = _nc.NC_HUMAN_38_Y

from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from gwaslab.info.g_Log import Log
from gwaslab.io.io_bedpe import bedpe_to_sumstats_coordinates, read_bedpe
from gwaslab.io.io_gtf import read_gtf
from gwaslab.io.io_ucsc_bed import bed_to_sumstats_coordinates, read_bed


class TestSexChrIoFormats(unittest.TestCase):
    """read_* and ChromosomeMapper.detect_reference_format on XY-labeled files."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="gwaslab_sexchr_io_")
        self.log = Log()
        self.log.verbose = False

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write_chain(self, name: str, target_chroms: list[str]) -> str:
        path = os.path.join(self.tmp, name)
        lines = []
        for i, tchrom in enumerate(target_chroms):
            # Minimal valid UCSC chain header; sizes are dummy but parseable.
            lines.append(
                f"chain 1000 {tchrom} 1000000 + 0 1000 {tchrom} 1000000 + 0 1000 {i + 1}\n"
            )
            lines.append("1000\n")
            lines.append("\n")
        with open(path, "w", encoding="utf-8") as f:
            f.write("".join(lines))
        return path

    def test_chain_detect_chr_prefix_xy(self):
        path = self._write_chain("xy.chain", ["chrX", "chrY"])
        m = ChromosomeMapper(log=self.log, verbose=False)
        fmt, prefix = m.detect_reference_format(path)
        self.assertEqual(fmt, "chr")
        self.assertEqual(prefix, "chr")

    def test_chain_detect_numeric_xy(self):
        path = self._write_chain("xy_num.chain", ["23", "24"])
        m = ChromosomeMapper(log=self.log, verbose=False)
        fmt, _prefix = m.detect_reference_format(path)
        self.assertEqual(fmt, "numeric")

    def test_chain_detect_string_name_xy(self):
        path = self._write_chain("xy_str.chain", ["X", "Y"])
        m = ChromosomeMapper(log=self.log, verbose=False)
        fmt, _prefix = m.detect_reference_format(path)
        self.assertIn(fmt, ("string", "numeric"))

    def test_chain_gz_detect_chrX(self):
        plain = self._write_chain("one.chain", ["chrX"])
        gz_path = os.path.join(self.tmp, "one.chain.gz")
        with open(plain, "rb") as fin:
            with gzip.open(gz_path, "wb") as fout:
                fout.write(fin.read())
        m = ChromosomeMapper(log=self.log, verbose=False)
        fmt, prefix = m.detect_reference_format(gz_path)
        self.assertEqual(fmt, "chr")
        self.assertEqual(prefix, "chr")

    def test_read_bed_chrX_chrY_bed_to_sumstats(self):
        path = os.path.join(self.tmp, "p.bed")
        with open(path, "w", encoding="utf-8") as f:
            f.write("chrX\t0\t100\n")
            f.write("chrY\t10\t200\n")
        df = read_bed(path, verbose=False, log=self.log)
        self.assertEqual(len(df), 2)
        self.assertListEqual(df["chrom"].tolist(), ["chrX", "chrY"])
        out = bed_to_sumstats_coordinates(df, log=self.log, verbose=False)
        chrs = set(int(x) for x in out["CHR"].tolist())
        self.assertEqual(chrs, {23, 24})

    def test_read_bed_numeric_23_24_bed_to_sumstats(self):
        path = os.path.join(self.tmp, "n.bed")
        with open(path, "w", encoding="utf-8") as f:
            f.write("23\t0\t100\n")
            f.write("24\t5\t50\n")
        df = read_bed(path, verbose=False, log=self.log)
        out = bed_to_sumstats_coordinates(df, log=self.log, verbose=False)
        self.assertEqual(set(int(x) for x in out["CHR"].tolist()), {23, 24})

    def test_read_bed_X_Y_bed_to_sumstats(self):
        path = os.path.join(self.tmp, "xy.bed")
        with open(path, "w", encoding="utf-8") as f:
            f.write("X\t0\t100\n")
            f.write("Y\t0\t100\n")
        df = read_bed(path, verbose=False, log=self.log)
        out = bed_to_sumstats_coordinates(df, log=self.log, verbose=False)
        self.assertEqual(set(int(x) for x in out["CHR"].tolist()), {23, 24})

    def test_read_bed_nc_hg38_with_build_mapper(self):
        ncx, ncy = NC_HUMAN_38_X, NC_HUMAN_38_Y
        path = os.path.join(self.tmp, "nc.bed")
        with open(path, "w", encoding="utf-8") as f:
            f.write(f"{ncx}\t0\t1000\n")
            f.write(f"{ncy}\t0\t1000\n")
        df = read_bed(path, verbose=False, log=self.log)
        mapper = ChromosomeMapper(
            species="homo sapiens", build="38", log=self.log, verbose=False
        )
        mapper.detect_sumstats_format(pd.Series([ncx, ncy]))
        out = bed_to_sumstats_coordinates(df, mapper=mapper, log=self.log, verbose=False)
        self.assertEqual(len(out), 2)
        # number_to_sumstats keeps NC_* when sumstats layer is NC format
        nums = {mapper.sumstats_to_number(x) for x in out["CHR"].tolist()}
        self.assertEqual(nums, {23, 24})

    def test_read_bedpe_xy_interchromosomal(self):
        path = os.path.join(self.tmp, "xy.bedpe")
        with open(path, "w", encoding="utf-8") as f:
            f.write("chrX\t0\t500\tchrY\t0\t500\tloop1\t1\t+\t-\n")
        df = read_bedpe(path, verbose=False, log=self.log)
        self.assertEqual(len(df), 1)
        out = bedpe_to_sumstats_coordinates(df, log=self.log, verbose=False)
        self.assertEqual(int(out["CHR1"].iloc[0]), 23)
        self.assertEqual(int(out["CHR2"].iloc[0]), 24)

    def test_read_bedpe_numeric_23_24(self):
        path = os.path.join(self.tmp, "n.bedpe")
        with open(path, "w", encoding="utf-8") as f:
            f.write("23\t0\t100\t24\t0\t100\t.\t0\n")
        df = read_bedpe(path, verbose=False, log=self.log)
        out = bedpe_to_sumstats_coordinates(df, log=self.log, verbose=False)
        self.assertEqual(int(out["CHR1"].iloc[0]), 23)
        self.assertEqual(int(out["CHR2"].iloc[0]), 24)

    def test_read_bedpe_X_Y_unprefixed(self):
        path = os.path.join(self.tmp, "xyu.bedpe")
        with open(path, "w", encoding="utf-8") as f:
            f.write("X\t0\t50\tY\t0\t50\t.\t0\n")
        df = read_bedpe(path, verbose=False, log=self.log)
        out = bedpe_to_sumstats_coordinates(df, log=self.log, verbose=False)
        self.assertEqual(int(out["CHR1"].iloc[0]), 23)
        self.assertEqual(int(out["CHR2"].iloc[0]), 24)

    def test_read_bedpe_nc_hg38_with_build_mapper(self):
        ncx, ncy = NC_HUMAN_38_X, NC_HUMAN_38_Y
        path = os.path.join(self.tmp, "nc.bedpe")
        with open(path, "w", encoding="utf-8") as f:
            f.write(f"{ncx}\t0\t100\t{ncy}\t0\t100\t.\t0\n")
        df = read_bedpe(path, verbose=False, log=self.log)
        mapper = ChromosomeMapper(
            species="homo sapiens", build="38", log=self.log, verbose=False
        )
        mapper.detect_sumstats_format(pd.Series([ncx, ncy]))
        out = bedpe_to_sumstats_coordinates(df, mapper=mapper, log=self.log, verbose=False)
        nums1 = mapper.sumstats_to_number(out["CHR1"].iloc[0])
        nums2 = mapper.sumstats_to_number(out["CHR2"].iloc[0])
        self.assertEqual(nums1, 23)
        self.assertEqual(nums2, 24)

    def test_read_gtf_chrX_chrY(self):
        path = os.path.join(self.tmp, "t.gtf")
        with open(path, "w", encoding="utf-8") as f:
            f.write(
                'chrX\ttest\tgene\t100\t200\t.\t+\t.\tgene_id "GX1";\n'
                'chrY\ttest\tgene\t300\t400\t.\t+\t.\tgene_id "GY1";\n'
            )
        df = read_gtf(path)
        self.assertEqual(len(df), 2)
        # read_gtf normalizes chr-prefixed sex chromosomes to X / Y
        self.assertListEqual(df["seqname"].tolist(), ["X", "Y"])

    def test_read_gtf_X_Y_and_numeric_xy(self):
        path = os.path.join(self.tmp, "mix.gtf")
        with open(path, "w", encoding="utf-8") as f:
            f.write(
                'X\ttest\tgene\t1\t2\t.\t+\t.\tgene_id "a";\n'
                'Y\ttest\tgene\t1\t2\t.\t+\t.\tgene_id "b";\n'
                '23\ttest\tgene\t1\t2\t.\t+\t.\tgene_id "c";\n'
                '24\ttest\tgene\t1\t2\t.\t+\t.\tgene_id "d";\n'
            )
        df = read_gtf(path)
        self.assertEqual(len(df), 4)
        # Numeric 23/24 are normalized to X/Y like chrX/chrY
        self.assertEqual(df["seqname"].tolist(), ["X", "Y", "X", "Y"])


class TestSexChrNotationMapper(unittest.TestCase):
    """ChromosomeMapper middle-layer consistency for sex chromosomes (hg38)."""

    def setUp(self):
        self.log = Log()
        self.log.verbose = False

    def _mapper38(self) -> ChromosomeMapper:
        return ChromosomeMapper(
            species="homo sapiens", build="38", log=self.log, verbose=False
        )

    def test_numeric_23_24_roundtrip(self):
        m = self._mapper38()
        m.detect_sumstats_format(pd.Series([23, 24]))
        self.assertEqual(m.sumstats_to_number(23), 23)
        self.assertEqual(m.sumstats_to_number(24), 24)
        self.assertEqual(m.number_to_sumstats(23), 23)
        self.assertEqual(m.number_to_sumstats(24), 24)

    def test_string_X_Y_without_prefix(self):
        m = self._mapper38()
        m.detect_sumstats_format(pd.Series(["X", "Y"]))
        self.assertEqual(m.sumstats_to_number("X"), 23)
        self.assertEqual(m.sumstats_to_number("Y"), 24)

    def test_chr_prefix_Chr_mixed_case(self):
        m = self._mapper38()
        m.detect_sumstats_format(pd.Series(["chrX", "ChrY"]))
        self.assertEqual(m.sumstats_to_number("chrX"), 23)
        self.assertEqual(m.sumstats_to_number("ChrY"), 24)
        self.assertEqual(m.sumstats_to_number("CHRX"), 23)

    def test_string_numeric_23_24_in_sumstats(self):
        m = self._mapper38()
        m.detect_sumstats_format(pd.Series(["23", "24"]))
        self.assertEqual(m.sumstats_to_number("23"), 23)
        self.assertEqual(m.sumstats_to_number("24"), 24)

    def test_nc_notation_hg38_X_Y(self):
        ncx, ncy = NC_HUMAN_38_X, NC_HUMAN_38_Y
        self.assertEqual(ncx, "NC_000023.11")
        self.assertEqual(ncy, "NC_000024.10")
        m = self._mapper38()
        m.detect_sumstats_format(pd.Series([ncx, ncy]))
        self.assertEqual(m.sumstats_to_number(ncx), 23)
        self.assertEqual(m.sumstats_to_number(ncy), 24)
        self.assertEqual(m.to_nc(23), ncx)
        self.assertEqual(m.to_nc(24), ncy)

    def test_nc_notation_hg19_X_Y(self):
        """NC accessions from chromosomes_nc.json homo sapiens build 19."""
        m = ChromosomeMapper(
            species="homo sapiens", build="19", log=self.log, verbose=False
        )
        m.detect_sumstats_format(pd.Series([NC_HUMAN_19_X, NC_HUMAN_19_Y]))
        self.assertEqual(NC_HUMAN_19_X, "NC_000023.10")
        self.assertEqual(NC_HUMAN_19_Y, "NC_000024.9")
        self.assertEqual(m.sumstats_to_number(NC_HUMAN_19_X), 23)
        self.assertEqual(m.sumstats_to_number(NC_HUMAN_19_Y), 24)
        self.assertEqual(m.to_nc(23), NC_HUMAN_19_X)
        self.assertEqual(m.to_nc(24), NC_HUMAN_19_Y)

    def test_sumstats_to_reference_fasta_style_chrX(self):
        m = self._mapper38()
        d = tempfile.mkdtemp()
        fasta = os.path.join(d, "r.fa")
        try:
            with open(fasta, "w", encoding="utf-8") as f:
                f.write(">chrX\nACGT\n>chrY\nACGT\n")
            m.detect_sumstats_format(pd.Series([23, 24]))
            self.assertEqual(
                m.sumstats_to_reference(23, reference_file=fasta, as_string=True), "chrX"
            )
            self.assertEqual(
                m.sumstats_to_reference(24, reference_file=fasta, as_string=True), "chrY"
            )
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_detect_reference_vcf_chrX_chrY_and_reference_to_sumstats(self):
        d = tempfile.mkdtemp()
        vcf = os.path.join(d, "r.vcf")
        try:
            with open(vcf, "w", encoding="utf-8") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("##contig=<ID=chrX>\n##contig=<ID=chrY>\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                f.write("chrX\t1\trs1\tA\tG\t.\t.\t.\n")
            m = ChromosomeMapper(build="38", log=self.log, verbose=False)
            fmt, prefix = m.detect_reference_format(vcf)
            self.assertEqual(fmt, "chr")
            self.assertEqual(prefix, "chr")
            m.detect_sumstats_format(pd.Series([23, 24]))
            self.assertEqual(m.reference_to_sumstats("chrX", reference_file=vcf), 23)
            self.assertEqual(m.reference_to_sumstats("chrY", reference_file=vcf), 24)
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
