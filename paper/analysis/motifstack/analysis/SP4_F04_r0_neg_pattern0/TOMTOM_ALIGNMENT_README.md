Current authoritative palette: user explicitly specified A#007F01, C#0200F3, G#F5B029, T#FF0003. Applied exactly to both query and database logos. Supersedes all previous sampled palettes.

Latest palette update: colors sampled from the user R-SAGE-net reference screenshot5fabeee8-37b1-4a71-9751-20cbad9b80d3: A#037F04 green, C#0303F3 blue, G#FCA705 orange-gold, T#FC0102 red. Both query and database logos use the same palette. All matrix values, alignment, labels, fonts and canvas size are unchanged.

Current aspect-ratio update: the upper discovered-pattern logo is1.836in high, increased from0.68in (2.7times) with the same6.96in plot width. The user-authorized canvas is8x5.4in. Native information content and y limits0–2.05bits are unchanged; only the physical height and surrounding layout are adjusted. Tick rules, two-line/one-line annotation layout, Arial12/14 and the axis-free database logo are retained.

Current tick layout: upper pattern displays original positions24–32 and bits ticks0,1,2; lower database logo has no x/y ticks, tick labels or axis titles. Both logo coordinate spines remain hidden. The two-line pattern and one-line database annotation are retained.

Current layout update: query metadata and CWM statistics are exactly two lines; all TOMTOM hit metadata is one line. All annotation text remains Arial12; axis titles Arial14. The canvas is8x5.4in, applying the user's prior explicit instruction that dimensions may be chosen to fit; the width allows retaining every field without reducing font size. The previous 5x4 PDF is replaced.

Text update: reference-style annotation now states Fold04/native pattern ID, 166 seqlets, shown core24–32 and consensus CCCCGCCCC; CWM sum importance score and per-position maximum/minimum calculated only from the native nine-position core CWM. These attribution statistics are distinct from the bits used by the two PPM logos. The p-val/species/TFBS/matrix syntax follows the user reference, using italic H. sapiens and bold MA0685.2; no reference-paper numeric values or plant names were copied.

# SP4 native TOMTOM alignment figure

PDF: ../../figures/SP4_MA0685_2_F04_TOMTOM_ALIGNMENT.pdf

Style: final-figure-small, user-authorized 8x5.4 inch canvas, PDF only, Arial12 tick/annotation text and Arial14 axis titles, pdf.fonttype42; no tight bounding box, no local result downloads. Colors sampled from the user reference: A green, C blue, G yellow-gold, T red.

Query: Fold04 r0_neg_patterns_pattern_0, 50 positions, 166 seqlets. Displayed region is the exact 9-position TOMTOM overlap, original pattern positions 24–32. Target: SP4 JASPAR MA0685.2, full 9 positions, reverse complement as the original TOMTOM Orientation minus. Original query offset -23 maps target position1 to query position24. Both aligned core consensuses are CCCCGCCCC. The target reverse complement is an alignment display operation, not a sign flip of attribution scores.

Native full-query TOMTOM row: p=1.7206e-05, E=0.00215075, q=0.00212697, overlap9, optimal offset-23, orientation-. Values were read from the frozen original result, not recomputed on the displayed crop. The figure shows rounded p/q values. The pattern also has SP1 and SP2 matches; SP4 is the user-selected match and is not claimed to be unique or the top hit.

The two logos use the actual motif probability matrices recorded in the original TOMTOM XML; each column sums to1 after correction of six-significant-digit XML rounding. Letter height is probability times DNA information content (2-entropy), units bits, without small-sample correction. These are motif-comparison logos and neither signed CWM nor a single KLHL23 seqlet attribution trace. Thus the information axes are nonnegative. Original data and the previously generated attribution figures were preserved. Exact matrices, source SHA256s, native row, figure dimensions/fonts and layout QA are in TOMTOM_ALIGNMENT_QA.json.

Orientation/offset documentation: https://meme-suite.org/meme/doc/tomtom-output-format.html
