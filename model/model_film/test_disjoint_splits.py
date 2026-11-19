"""Pure-standard-library tests; no GPU, real datasets or training needed."""
import unittest
from disjoint_splits import assign_folds, coordinate, intersect_coverage, merge_coverage, overlap_components


def window(index, start, end, chrom="chr1", strand="+", region="tss_sequence"):
    return dict(gene_index=index, gene=f"gene{index}", chrom=chrom, strand=strand,
                region=region, start0=start, end0=end, fold=index % 10 + 1)


class Tests(unittest.TestCase):
    def test_touching_is_not_overlap(self):
        groups, edges, _ = overlap_components(2, [window(0, 10, 20), window(1, 20, 30)])
        self.assertEqual(groups, [[0], [1]])
        self.assertEqual(edges, [])

    def test_opposite_strands_and_different_window_types_overlap(self):
        groups, _, _ = overlap_components(2, [window(0, 10, 20), window(1, 19, 30, strand="-", region="first_acceptor_sequence")])
        self.assertEqual(groups, [[0, 1]])

    def test_different_chromosomes(self):
        groups, _, _ = overlap_components(2, [window(0, 10, 20), window(1, 10, 20, chrom="chr2")])
        self.assertEqual(groups, [[0], [1]])

    def test_transitive(self):
        groups, edges, _ = overlap_components(3, [window(0, 10, 20), window(1, 19, 30), window(2, 29, 40)])
        self.assertEqual(groups, [[0, 1, 2]])
        self.assertEqual(len(edges), 2)

    def test_distant_windows_of_one_gene_connect_groups(self):
        groups, _, _ = overlap_components(3, [window(0, 10, 20), window(0, 100, 120), window(1, 19, 30), window(2, 110, 130)])
        self.assertEqual(groups, [[0, 1, 2]])

    def test_single_gene_kept(self):
        groups, edges, summary = overlap_components(1, [window(0, 10, 20), window(0, 12, 22)])
        self.assertEqual(groups, [[0]])
        self.assertEqual(edges, [])
        self.assertEqual(summary["same_gene_overlapping_window_pairs"], 1)

    def test_balancing_deterministic(self):
        components = [[0, 1, 2]] + [[i] for i in range(3, 103)]
        first = assign_folds(components, ["chr1"] * 103, 20260916)
        self.assertEqual(first, assign_folds(components, ["chr1"] * 103, 20260916))
        self.assertEqual(len(set(first[0][:3])), 1)
        counts = [first[0].count(fold) for fold in range(1, 11)]
        self.assertLessEqual(max(counts) - min(counts), 1)

    def test_independent_coverage_detection(self):
        intervals = [window(0, 10, 20), window(0, 12, 22), window(1, 19, 30)]
        self.assertEqual(intersect_coverage(merge_coverage(intervals, [1]), merge_coverage(intervals, [2])), (1, 3))
        self.assertEqual(intersect_coverage({"chr1": [(10, 20)]}, {"chr1": [(20, 30)]}), (0, 0))

    def test_invalid_coordinate(self):
        for value in ("-1", "inf", "1.5"):
            with self.assertRaises(ValueError):
                coordinate(value)
        self.assertIsNone(coordinate(""))


if __name__ == "__main__":
    unittest.main()
