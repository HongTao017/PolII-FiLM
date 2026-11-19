// For a single feature, the positive PCA loading is +1 and the score is its Z-score.
export function singleGenePiPca(rows) {
  const observed = rows.filter(row => Number.isFinite(row.log2_pi))
  if (observed.length < 2) return { status: 'unavailable', reason: 'At least two measured PI values are required.' }
  const mean = observed.reduce((sum, row) => sum + row.log2_pi, 0) / observed.length
  const variance = observed.reduce((sum, row) => sum + (row.log2_pi - mean) ** 2, 0) / observed.length
  if (variance === 0) return { status: 'unavailable', reason: 'All measured PI values are identical; no variable principal component exists.' }
  const std = Math.sqrt(variance)
  return { status: 'complete', mean, std, feature_count: 1, component_count: 1, explained_variance_ratio: [1], points: observed.map(row => ({ ...row, pc1: (row.log2_pi - mean) / std, pc2: 0 })) }
}
