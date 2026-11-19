import unittest
from types import SimpleNamespace

import numpy as np

from inspect_splice_annotation import describe_transcript
from prepare_four_regions import boundary_window, decode


class BoundaryTests(unittest.TestCase):
    def test_plus_and_minus_half_open_junction(self):
        fasta = {"chr": "AACCGTAAGGTT"}
        plus, info = boundary_window(fasta, "chr", 4, "+", -4, 4)
        minus, minus_info = boundary_window(fasta, "chr", 4, "-", -4, 4)
        self.assertEqual(decode(plus), "AACCGTAA")
        self.assertEqual(decode(minus), "TTACGGTT")
        self.assertEqual(info["requested_start0"], 0)
        self.assertEqual(minus_info["requested_start0"], 0)
        self.assertEqual(decode(plus[4:6]), "GT")

    def test_chromosome_padding_orientation(self):
        fasta = {"chr": "ACGTACGT"}
        plus, _ = boundary_window(fasta, "chr", 1, "+", -4, 4)
        minus, _ = boundary_window(fasta, "chr", 1, "-", -4, 4)
        self.assertEqual(decode(plus), "NNNACGTA")
        self.assertEqual(decode(minus), "TACGTNNN")
        self.assertEqual(len(plus), 8)
        self.assertEqual(len(minus), 8)

    def test_first_intron_plus(self):
        row = SimpleNamespace(transcript="TX", chrom="chr", strand="+",
                              tss_bed0=10, tes_bed0=59)
        exons = {"TX": [("chr", "+", 10, 20, 1), ("chr", "+", 40, 45, 2),
                        ("chr", "+", 50, 60, 3)]}
        result = describe_transcript(row, exons)
        self.assertEqual(result["first_donor_boundary0"], 20)
        self.assertEqual(result["first_acceptor_boundary0"], 40)
        self.assertEqual(result["first_intron_length"], 20)

    def test_first_intron_minus_not_last_intron(self):
        row = SimpleNamespace(transcript="TX", chrom="chr", strand="-",
                              tss_bed0=59, tes_bed0=10)
        exons = {"TX": [("chr", "-", 10, 20, 3), ("chr", "-", 35, 45, 2),
                        ("chr", "-", 50, 60, 1)]}
        result = describe_transcript(row, exons)
        self.assertEqual(result["first_donor_boundary0"], 50)
        self.assertEqual(result["first_acceptor_boundary0"], 45)
        self.assertEqual(result["first_intron_start0"], 45)
        self.assertEqual(result["first_intron_end0"], 50)

    def test_intronless_and_missing_are_distinct(self):
        row = SimpleNamespace(transcript="TX", chrom="chr", strand="+",
                              tss_bed0=10, tes_bed0=19)
        self.assertEqual(describe_transcript(row, {"TX": [("chr", "+", 10, 20, 1)]})
                         ["splice_status"], "intronless")
        self.assertEqual(describe_transcript(row, {})["splice_status"], "missing_transcript_annotation")

    def test_endpoint_mismatch_rejected(self):
        row = SimpleNamespace(transcript="TX", chrom="chr", strand="+",
                              tss_bed0=11, tes_bed0=19)
        self.assertEqual(describe_transcript(row, {"TX": [("chr", "+", 10, 20, 1)]})
                         ["splice_status"], "transcript_endpoint_mismatch")


if __name__ == "__main__":
    unittest.main()
