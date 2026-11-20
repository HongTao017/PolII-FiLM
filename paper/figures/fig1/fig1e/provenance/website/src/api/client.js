const config = window.NASCENT_CONFIG || {}
export const mode = config.dataMode === 'api' ? 'api' : 'catalog'
export const apiBase = (config.apiBase || '/api/v1').replace(/\/$/, '')
const cache = new Map()

export async function request(url, options = {}) {
  const { allowNotFound, ...fetchOptions } = options
  const controller = new AbortController()
  const timeout = setTimeout(() => controller.abort(), 20000)
  try {
    const response = await fetch(url, { ...fetchOptions, signal: controller.signal })
    if (response.status === 404 && allowNotFound) return null
    if (!response.ok) throw new Error(`Request failed (${response.status})`)
    return await response.json()
  } catch (error) {
    throw new Error(error.name === 'AbortError' ? 'Request timed out. Please try again.' : error.message)
  } finally { clearTimeout(timeout) }
}
function local(name) {
  const key = `./data/${name}.json`
  if (!cache.has(key)) cache.set(key, request(key).catch(error => { cache.delete(key); throw error }))
  return cache.get(key)
}
// Backend responses share the local JSON contract; see docs/API.md.
export function getCatalog() { return mode === 'api' ? request(`${apiBase}/catalog`) : local('catalog') }
export function getProfilePiAnalysis() {
  return mode === 'api' ? request(`${apiBase}/analysis/profile-pi`) : local('analysis/profile-pi')
}
export function getPiSpearmanHeatmap() {
  return mode === 'api' ? request(`${apiBase}/analysis/pi-spearman`) : local('analysis/pi-spearman')
}
export async function getProfilePiValues(version, id) {
  const result = await (mode === 'api'
    ? request(`${apiBase}/analysis/profile-pi/${encodeURIComponent(id)}?version=${encodeURIComponent(version)}`)
    : local(`analysis/profiles-${version}/${id}`))
  if(result.version !== version || result.profile_id !== id || !Array.isArray(result.values)
    || result.values.length !== result.summary?.observed || result.values.some(v=>!Number.isFinite(v))) {
    throw new Error('Profile PI data do not match this analysis release. Reload the page.')
  }
  return result
}
export function validatePiAnalysis(result, catalog) {
  if(result.schema_version !== 1 || result.feature_count !== catalog.summary.transcripts
    || JSON.stringify(result.profile_order) !== JSON.stringify(catalog.metric.profile_order)
    || Object.entries(catalog.source_hashes).some(([name,hash])=>result.source_hashes?.[name]!==hash)) {
    throw new Error('PI analysis and profile catalog versions differ. Reload the page after the data update.')
  }
}
export function getGeneIndex() { return mode === 'api' ? request(`${apiBase}/genes/index`) : local('gene-index') }
export function getGeneCatalog() { return mode === 'api' ? request(`${apiBase}/genes/catalog`) : local('gene-catalog') }
export function getProfilePca() {
  return request(mode === 'api' ? `${apiBase}/analysis/pi-pca` : './GENE/pi-pca.json', { allowNotFound: true })
}
export function getGenePiTsne(symbol, transcript) {
  const base = mode === 'api' ? apiBase : '/api/v1'
  return request(`${base}/analysis/gene-pi-tsne?${new URLSearchParams({ symbol, transcript })}`)
}
export async function getGene(symbol) {
  if (mode === 'api') return request(`${apiBase}/genes/${encodeURIComponent(symbol)}`)
  const index = await getGeneIndex()
  const entry = index.find(gene => gene.symbol.toLowerCase() === symbol.toLowerCase())
  if (!entry) return null
  const chunk = await local(`genes/chunk-${entry.bucket}`)
  return chunk[entry.symbol] || null
}

export function getLibraryQc() { return mode === 'api' ? request(`${apiBase}/analysis/library-qc`) : local('analysis/library-qc') }
export async function getProfileProvenance(id) {
  if (mode === 'api') return request(`${apiBase}/profiles/${encodeURIComponent(id)}/provenance`)
  const result = await local('profile-provenance')
  if (result.schema_version !== 1 || result.profiles?.length !== result.profile_count) throw Error('Invalid profile metadata. Please retry.')
  return { checked_on: result.checked_on, record: result.profiles.find(p=>p.profile_id===id), sources: result.sources }
}

export async function getGeneFeatures(symbol) {
  let result
  if (mode === 'api') {
    result = await request(`${apiBase}/genes/${encodeURIComponent(symbol)}/features`)
  } else {
    const [manifest, index] = await Promise.all([local('gene-features'), getGeneIndex()])
    const entry = index.find(gene => gene.symbol.toLowerCase() === symbol.toLowerCase())
    if (!entry) return { manifest, gene: null }
    const chunk = await local(`gene-features/${manifest.version}/chunk-${entry.bucket}`)
    if (chunk.version !== manifest.version) throw new Error('Feature data versions differ. Reload the page.')
    result = { manifest, gene: chunk.genes[entry.symbol] || null }
  }
  if (result.manifest?.schema_version !== 1 || result.manifest.features?.length !== result.manifest.feature_count
    || (result.gene && (!Array.isArray(result.gene.values) || result.gene.values.length !== result.manifest.feature_count
      || result.gene.values.some(value => value !== null && !Number.isFinite(value))))) {
    throw new Error('Invalid gene feature data. Retry after the data update.')
  }
  return result
}
