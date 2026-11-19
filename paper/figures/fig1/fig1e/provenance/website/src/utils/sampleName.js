// Human-readable labels only. Keep catalog IDs for joins, routes and API calls.
export function sampleName(profile) {
  const alias = String(profile?.alias || '')
  const start = alias.search(/(?:^|_)GSE\d+(?=_|$)/)
  if (start >= 0) return alias.slice(start).replace(/^_/, '')
  return alias || profile?.id || ''
}

export const sampleNameColumn = { label: 'sample_name', value: sampleName }
