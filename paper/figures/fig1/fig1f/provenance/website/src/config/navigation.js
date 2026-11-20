import { Layers3, Dna, ChartColumn, Download, CircleHelp, Search } from '@lucide/vue'

// One registry keeps the home shortcuts and main navigation aligned.
export const features = [
  { path: '/genes', key: 'genes', label: 'Genes', icon: Dna, description: 'Inspect measured pausing indices across genes, transcripts and cellular contexts.' },
  { path: '/search', key: 'search', label: 'Search', icon: Search, description: 'Find genes, RefSeq transcripts, GEO studies, profiles and cell labels in one search.' },
  { path: '/samples', key: 'browse', label: 'Browse', icon: Layers3, description: 'Explore PRO-seq, GRO-seq and mNET-seq profiles by assay, cell label and antibody subtype.' },
  { path: '/statistics', key: 'statistics', label: 'Statistics', icon: ChartColumn, description: 'Explore the composition of collected studies, assays and raw cell labels.' },
  { path: '/download', key: 'download', label: 'Download', icon: Download, description: 'Download metadata, provenance manifests and filtered results for further analysis.' },
  { path: '/help', key: 'about', label: 'About', icon: CircleHelp, description: 'Learn about database scope, counting units, pausing metrics and data provenance.' },
]
export const mainNavigation = [{ path: '/', label: 'Home' }, ...['genes', 'browse', 'statistics', 'search', 'download', 'about'].map(key => features.find(link => link.key === key))]
export function isNavigationActive(currentPath, path) { return path === '/samples' && currentPath.startsWith('/studies/') ? true : path === '/' ? currentPath === '/' : currentPath === path || currentPath.startsWith(path + '/') }
