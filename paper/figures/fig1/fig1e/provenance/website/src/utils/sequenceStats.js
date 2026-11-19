// Statistics use only called A/T/G/C bases and adjacent called pairs.
export function sequenceStats(sequence, binSize = 100) {
  const dna = sequence.toUpperCase(), counts = { A: 0, T: 0, G: 0, C: 0 }
  let validPairs = 0, cpg = 0
  for (let i = 0; i < dna.length; i++) {
    if (dna[i] in counts) counts[dna[i]]++
    if (i && dna[i - 1] in counts && dna[i] in counts) {
      validPairs++
      if (dna[i - 1] === 'C' && dna[i] === 'G') cpg++
    }
  }
  const called = Object.values(counts).reduce((sum, value) => sum + value, 0)
  const bins = []
  for (let start = 0; start < dna.length; start += binSize) {
    const fragment = dna.slice(start, start + binSize)
    const valid = [...fragment].filter(base => base in counts)
    const gc = valid.filter(base => base === 'G' || base === 'C').length
    bins.push({ start, end: start + fragment.length - 1, called: valid.length, gc: valid.length ? 100 * gc / valid.length : null })
  }
  return { counts, called, ambiguous: dna.length - called, cpg, validPairs, cpgPercent: validPairs ? 100 * cpg / validPairs : null, bins }
}
