# SP4-like motif attribution overview

Figure: ../../figures/SP4_MA0685_2_F04_r0_neg_pattern0_importance.pdf

Top: original 600 bp KLHL23 TSS input, five-seed mean actual DeepLIFT contribution per observed base. Lower: the exact 50 bp source seqlet (input [280,330), TSS-relative [-20,+29]); each observed base letter height retains its signed attribution in Y_PC1 units. Grey marks the 9 bp SP4 sequence match, input [296,305), TSS-relative [-4,+4]. The forward transcript-oriented sequence at the hit is CCCCGCCCC; reverse-complement GGGGCGGGG matches MA0685.2 on FIMO strand minus. No sign flip or absolute value was applied.

Representative selected by smallest FIMO p among negative-score SP4 sites, then seqlet index/start/strand, not by maximum attribution magnitude. KLHL23 was selected from pattern F04__r0_neg_patterns_pattern_0 (166 seqlets). This is one illustrative instance, not the average pattern CWM or a population effect estimate. All 166 seqlet one-hot sequences and stored CWM scores were checked against original arrays.

The pattern matches SP1, SP4, and SP2 (original TOMTOM SP4 q=0.00212697; the best original hit is SP1). The requested SP4 label denotes sequence similarity and does not establish SP4 occupancy. The site FIMO p=1.5e-6 is a sequence-match score, not an attribution significance test.

Original model target: signed predicted Y_PC1, inverse transformed; actual FiLM regulation fixed. Method: custom DeepLIFT Rescale with explicit GroupNorm decomposition and symmetric product rule. Twenty dinucleotide-shuffled references per gene, each of five original seeds, no new model inference or training. Scores are model attributions, not physical binding energies and not motif information bits.

Selection/source SHA256: selection.json. Plot numerical/font/page QA: figure_QA.json. Source data and circular motif figures were not changed. PDF only; 8x5 inches (user-authorized size override); Arial 12 pt labels/ticks and 14 pt axis titles; editable PDF text.
