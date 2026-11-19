// Sequence coordinates are one-based inclusive; signal bins are BED0 half-open.
export function baseCoordinate(data, index) {
  return data.display_strand === '-' ? data.genomic_end - index : data.genomic_start + index
}
export function baseOffset(data, index) {
  if (data.endpoint) return (baseCoordinate(data, index) - data.endpoint.anchor) * (data.gene_strand === '-' ? -1 : 1) || 0
  return data.offset + (data.display_strand === data.gene_strand ? index : data.length - index - 1)
}
export function alignBaseSignal(data, response, profile, strand) {
  if (data.length !== data.sequence.length || data.genomic_end - data.genomic_start + 1 !== data.length
    || response.profile !== profile.id || response.strand !== strand || response.assembly !== 'hg38'
    || response.bin_width !== 1 || response.unit !== profile.unit || response.normalization !== 'none'
    || !Array.isArray(response.features) || response.features.length !== data.length) throw Error('Signal and reference sequence do not match at single-base resolution.')
  const values = response.features.map((f, i) => {
    if (f.chr !== data.chrom || f.start !== data.genomic_start - 1 + i || f.end !== f.start + 1
      || !Number.isFinite(f.value) || f.value < 0) throw Error('Incomplete or invalid single-base signal. Please retry.')
    return f.value
  })
  if (data.display_strand === '-') values.reverse()
  return [...data.sequence].map((base, index) => ({ base, index, coordinate: baseCoordinate(data, index), relative: baseOffset(data, index), value: values[index] }))
}
export function signalCeiling(values) {
  const maximum = Math.max(0, ...values)
  if (!maximum) return 1
  const power = 10 ** Math.floor(Math.log10(maximum)), fraction = maximum / power
  return (fraction <= 1 ? 1 : fraction <= 2 ? 2 : fraction <= 5 ? 5 : 10) * power
}
// Glyph paths span exactly y=0..100, so drawn heights are proportional to signal.
export const baseGlyphs = {
  A: 'M0 100L34 0H66L100 100H72L64 76H36L28 100ZM44 54H56L50 32Z',
  T: 'M0 0H100V24H64V100H36V24H0Z',
  C: 'M100 18L78 38C70 26 38 24 30 40C20 60 34 78 50 78C64 78 72 70 78 62L100 82C86 98 70 100 50 100C15 100 0 80 0 50S15 0 50 0C70 0 86 2 100 18Z',
  G: 'M100 18L78 38C70 26 38 24 30 40C20 60 34 78 50 78C61 78 70 74 76 68V60H53V40H100V82C86 98 70 100 50 100C15 100 0 80 0 50S15 0 50 0C70 0 86 2 100 18Z',
}
