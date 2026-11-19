// Linear interpolated quartiles; Tukey whiskers remain within observed values.
export function summarizePi(rows) {
  const values = rows.map(row => row.log2_pi).filter(Number.isFinite).sort((a, b) => a - b)
  if (!values.length) return null
  const quantile = probability => {
    const position = (values.length - 1) * probability
    const low = Math.floor(position), high = Math.ceil(position)
    return values[low] + (values[high] - values[low]) * (position - low)
  }
  const q1 = quantile(.25), median = quantile(.5), q3 = quantile(.75), iqr = q3 - q1
  return { count: values.length, q1, median, q3, low: values.find(value => value >= q1 - 1.5 * iqr), high: values.findLast(value => value <= q3 + 1.5 * iqr) }
}

export function groupPi(rows, key) {
  const grouped = new Map()
  for (const row of rows) {
    if (!Number.isFinite(row.log2_pi)) continue
    const group = row[key] || 'Not recorded'
    if (!grouped.has(group)) grouped.set(group, [])
    grouped.get(group).push(row)
  }
  return [...grouped].sort(([a], [b]) => a.localeCompare(b)).map(([label, profiles]) => ({ label, profiles, ...summarizePi(profiles) }))
}
