import { sampleNameColumn } from './sampleName.js'
export const assays = ['GRO-seq', 'PRO-seq', 'mNET-seq']
export const assayColors = { 'GRO-seq': '#3e77ab', 'PRO-seq': '#e69a9a', 'mNET-seq': '#d69a41' }
export function assayClass(assay) { return assay === 'GRO-seq' ? 'gro' : assay === 'PRO-seq' ? 'pro' : 'mnet' }
export function filterRecords(rows, filters) {
  const query = (filters.q || '').trim().toLowerCase()
  return rows.filter(row => (!filters.assay || row.assay === filters.assay || row.assays?.includes(filters.assay))
    && (!filters.cell || row.cell === filters.cell || row.cells?.includes(filters.cell))
    && (!filters.study || row.gse === filters.study)
    && (!filters.antibody || row.antibody === filters.antibody)
    && (!query || [row.id, row.alias, row.gse, row.cell, row.author, row.antibody, ...(row.cells || []), ...(row.assays || [])].filter(Boolean).join(' ').toLowerCase().includes(query)))
}
export function exportTsv(name, columns, rows) {
  const clean = value => String(value ?? '').replace(/[\t\r\n]/g, ' ')
  const text = '\ufeff' + [columns.map(c => clean(c.label)).join('\t'), ...rows.map(row => columns.map(c => clean(c.value ? c.value(row) : row[c.key])).join('\t'))].join('\n')
  downloadBlob(name, text, 'text/tab-separated-values;charset=utf-8')
}
export function downloadBlob(name, text, type = 'application/json;charset=utf-8') {
  const link = document.createElement('a')
  const url = URL.createObjectURL(new Blob([text], { type }))
  link.href = url; link.download = name; link.click()
  setTimeout(() => URL.revokeObjectURL(url), 1000)
}
export const sampleColumns = [
  { key: 'id', label: 'profile_id' }, sampleNameColumn, { key: 'alias', label: 'original_alias' }, { key: 'assay', label: 'assay' },
  { key: 'gse', label: 'GSE' }, { key: 'cell', label: 'raw_cell_label' }, { key: 'antibody', label: 'antibody_from_alias' },
  { key: 'source_file', label: 'source_file' }, { key: 'status', label: 'curation_status' },
]
