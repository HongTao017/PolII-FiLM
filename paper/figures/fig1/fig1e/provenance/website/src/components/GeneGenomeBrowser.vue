<script setup>
import { assayColors } from '../utils/data.js'
import { sampleName } from '../utils/sampleName.js'
import { computed, nextTick, onBeforeUnmount, ref, shallowRef, watch } from 'vue'
import { geneRegions } from '../utils/geneRegions.js'
const props = defineProps({ symbol: { type: String, required: true }, transcript: { type: String, default: '' }, profiles: { type: Array, required: true }, defaultProfileIds: { type: Array, default: () => [] }, compact: { type: Boolean, default: false } })
const manifest = shallowRef(null), host = ref(null), loading = ref(false), applying = ref(false), error = ref(''), selected = ref([]), strand = ref('+'), ready = ref(false)
let browser, igv, generation = 0, trackRevision = 0, disposed = false, pendingTracks = false
const colors = { 'GRO-seq': '#4979ad', 'PRO-seq': assayColors['PRO-seq'], 'mNET-seq': '#d09a42' }
const allowed = computed(() => new Set(props.profiles.map(p => p.id)))
const available = computed(() => (manifest.value?.profiles || []).filter(p => allowed.value.has(p.id)))
const selectedCell = ref(''), pickerOpen = ref(false)
const cellOptions = computed(() => [...new Set(available.value.map(p => p.cell))].sort((a, b) => a.localeCompare(b)).map(cell => ({ cell, count: available.value.filter(p => p.cell === cell).length })))
const cellProfiles = computed(() => available.value.filter(p => !selectedCell.value || p.cell === selectedCell.value))
const activeProfiles = computed(() => available.value.filter(p => selected.value.includes(p.id)))
const selectionKey = computed(() => available.value.map(p => p.id).join('|'))
const activeSummary = ref([])
const activeStrand = ref('+')
const emptySelection = computed(() => ready.value && !applying.value && activeSummary.value.length === 0)
const api = '/api/v1/genome/'
async function get(action, params, timeout = 50000) {
  const controller = new AbortController(), timer = setTimeout(() => controller.abort(), timeout)
  try {
    const response = await fetch(api + action + '?' + new URLSearchParams(params), { signal: controller.signal })
    const value = await response.json()
    if (!response.ok) throw new Error(value.error || `Genome request failed (${response.status}).`)
    return value
  } catch (e) { throw new Error(e.name === 'AbortError' ? 'Signal loading timed out. Zoom in and retry.' : e.message) }
  finally { clearTimeout(timer) }
}
function defaults(profiles = available.value) {
  const preferred = props.defaultProfileIds.filter(id => profiles.some(p => p.id === id)).slice(0, 6)
  if (preferred.length) return preferred
  const assays = new Set()
  return profiles.filter(p => { if (assays.has(p.assay)) return false; assays.add(p.assay); return true }).slice(0, 3).map(p => p.id)
}
function selectCell() {
  selected.value = defaults(cellProfiles.value)
  pickerOpen.value = true
  applyTracks()
}
function destroy() {
  ++trackRevision; pendingTracks = false; ready.value = false
  if (browser && igv) igv.removeBrowser(browser)
  browser = null; activeSummary.value = []
}
function featureSource(profile, token, signalStrand, revision = null) {
  const current = () => token === generation && !disposed && (revision === null || revision === trackRevision)
  return {
    async getHeader() { return {} },
    async getFeatures({ chr, start, end, bpPerPixel }) {
      const chromosomes = manifest.value?.chromosomes
      if (!chromosomes || !chromosomes[chr]) return []
      start = Math.max(0, Math.floor(start)); end = Math.min(chromosomes[chr], Math.ceil(end))
      if (end <= start) return []
      if (end - start > 1000000) {
        if (current()) error.value = 'Zoom in to a region of at most 1 Mb to load annotations and signal.'
        return []
      }
      try {
        const value = await get(profile ? 'signal' : 'annotation', { chr, start, end, ...(profile ? { profile: profile.id, strand: signalStrand, bins: Math.min(2048, Math.max(100, Math.ceil((end - start) / Math.max(1, bpPerPixel || 1)))) } : {}) })
        if (!current()) return []
        const features = value.features
        if (!profile) {
          // IGV custom sources own their transcript row packing.
          const rowEnds = []
          for (const f of features) {
            let row = rowEnds.findIndex(right => right + (bpPerPixel || 1) * 70 < f.start)
            if (row === -1) row = rowEnds.length
            rowEnds[row] = f.end; f.row = row
          }
        }
        return features
      } catch (e) {
        if (!current()) return []
        error.value = `${profile ? sampleName(profile) + ': ' : ''}${e.message}`
        throw e
      }
    }
  }
}
async function applyTracks() {
  ++trackRevision; pendingTracks = true; activeSummary.value = []
  if (!browser || applying.value) return
  const target = browser, token = generation
  applying.value = true; error.value = ''
  try {
    while (pendingTracks && target === browser && token === generation) {
      pendingTracks = false
      error.value = ''
      const revision = trackRevision, profiles = [...activeProfiles.value].slice(0, 6), signalStrand = strand.value
      for (const view of [...target.trackViews]) if (view.track.id?.startsWith('signal-')) target.removeTrack(view.track)
      for (const p of profiles) {
        if (target !== browser || token !== generation) return
        if (revision !== trackRevision) break
        await target.loadTrack({ id: 'signal-' + p.id, name: `${sampleName(p)} · ${p.assay} · ${p.author}`, type: 'wig', height: props.compact ? 56 : 100, color: colors[p.assay] || '#4979ad', autoscale: true, min: 0, graphType: 'bar', visibilityWindow: 300000, supportsWholeGenome: false, featureSource: featureSource(p, token, signalStrand, revision) })
      }
      if (target === browser && token === generation && revision === trackRevision) { activeSummary.value = profiles; activeStrand.value = signalStrand }
    }
  } catch (e) { if (token === generation) error.value = e.message }
  finally { if (token === generation) applying.value = false }
}
async function load() {
  const token = ++generation
  destroy(); manifest.value = null; error.value = ''; loading.value = true; applying.value = false
  try {
    const value = await get('manifest', { symbol: props.symbol, transcript: props.transcript }, 20000)
    if (token !== generation || disposed) return
    manifest.value = value; strand.value = value.strand; selectedCell.value = ''; pickerOpen.value = false; selected.value = defaults()
    await nextTick()
    igv ||= (await import('igv/dist/igv.esm.js')).default
    if (token !== generation || disposed) return
    const reference = Object.fromEntries(Object.entries(value.reference).map(([key, val]) => [key, key.endsWith('URL') ? new URL(val, location.origin).href : val]))
    const created = await igv.createBrowser(host.value, { reference, loadDefaultGenomes: false, locus: value.locus, showNavigation: true, showChromosomeWidget: true, showIdeogram: Boolean(reference.cytobandURL), showSequence: false, showCursorGuide: true, showCenterGuide: false, tracks: [{ id: 'mane-annotation', name: value.annotation, type: 'annotation', displayMode: 'EXPANDED', height: props.compact ? 70 : 170, maxHeight: props.compact ? 110 : 260, color: '#315c9b', visibilityWindow: 300000, supportsWholeGenome: false, featureSource: featureSource(null, token) }] })
    if (token !== generation || disposed) { igv.removeBrowser(created); return }
    browser = created; ready.value = true
    await applyTracks()
  } catch (e) { if (token === generation) error.value = e.message }
  finally { if (token === generation) loading.value = false }
}
async function resetLocus() {
  if (!browser || !manifest.value) return
  error.value = ''
  try { await browser.search(manifest.value.locus) } catch (e) { error.value = e.message }
}
async function showEndpoint(region){
  if(!browser||!manifest.value)return
  error.value=''
  try{const endpoint=manifest.value.endpoint_regions?.[region];if(!endpoint)throw new Error('The selected transcript has no available annotation for this region.');await browser.search(endpoint.locus)}
  catch(e){error.value=e.message}
}
watch(() => [props.symbol, props.transcript], load, { immediate: true })
watch(selectionKey, () => {
  if (!ready.value) return
  if (!cellOptions.value.some(option => option.cell === selectedCell.value)) selectedCell.value = ''
  selected.value = selected.value.filter(id => cellProfiles.value.some(p => p.id === id))
  if (!selected.value.length) selected.value = defaults(cellProfiles.value)
  // Coalesce changes into the latest filter state, without publishing an older result.
  applyTracks()
})
onBeforeUnmount(() => { disposed = true; ++generation; destroy() })
</script>

<template>
  <section class="gene-genome-browser" :class="{'genome-preview':compact}" aria-label="Gene genome browser">
    <header class="genome-heading">
      <div><h2>{{ symbol }} · {{compact ? 'Signal preview' : 'Genome browser'}}</h2><p>{{compact ? 'Example gene structure and sequencing signals.' : 'Explore gene structure and strand-specific sequencing signal.'}}</p></div>
      <RouterLink v-if="compact" class="button secondary" :to="'/genes/'+encodeURIComponent(symbol)">Open full browser</RouterLink><button v-else-if="ready" class="button secondary" type="button" @click="resetLocus">Reset gene view</button>
    </header>
    <component :is="compact ? 'details' : 'div'" v-if="manifest" class="genome-controls"><summary v-if="compact">Preview controls</summary>
      <p class="genome-reference">{{ manifest.reference.name }} · {{ manifest.annotation }} · {{ manifest.transcript }} · Gene strand {{ manifest.strand }}</p>
      <p v-if="!manifest.exact_transcript" class="genome-reference">Requested {{ transcript }} is not present in MANE v1.5. The view uses {{ manifest.transcript }} for this gene.</p>
      <p v-if="manifest.initial_view_truncated" class="genome-reference">This long gene opens around its transcription start site. Pan to explore the full locus.</p>
      <div v-if="!compact" class="genome-endpoints" role="group" aria-label="Transcript endpoint views"><button v-for="r in geneRegions" :key="r.id" type="button" class="button secondary" :disabled="!ready || loading || applying || !manifest.endpoint_regions?.[r.id]" @click="showEndpoint(r.id)">{{r.label}}</button></div>
      <p v-if="!compact" class="genome-reference">Offsets follow the MANE transcript's 5′ → 3′ direction. Splice windows use the first intron in transcription order. ± windows include both endpoints; TES upstream covers −1000 to −1.</p>
      <p v-if="!compact && !manifest.first_intron" class="genome-reference">This MANE transcript has no annotated first intron. Splice-site views are unavailable.</p>
      <div class="genome-actions"><label>Cell line <select v-model="selectedCell" aria-label="Genome cell line" :disabled="!ready || loading" @change="selectCell"><option value="">All cell lines</option><option v-for="option in cellOptions" :key="option.cell" :value="option.cell">{{ option.cell }} ({{ option.count }})</option></select></label><label>Signal strand <select v-model="strand" aria-label="Genome signal strand"><option value="+">+ strand</option><option value="-">− strand</option></select></label><button type="button" class="button primary" :disabled="!ready || applying" @click="applyTracks">{{ applying ? 'Loading tracks…' : 'Apply tracks' }}</button><span>{{ selected.length }} / 6 selected · {{ available.length }} matching signal profiles</span></div>
      <button v-if="compact && ready" class="button secondary preview-reset" type="button" @click="resetLocus">Reset gene view</button><details class="genome-picker" :open="pickerOpen" @toggle="pickerOpen = $event.target.open"><summary>Select signal tracks</summary><p>Initially one profile per available assay. Choose up to six profiles{{ selectedCell ? " for " + selectedCell : " from the current source filters" }}, then apply tracks.</p><div class="genome-options"><label v-for="p in cellProfiles" :key="p.id"><input v-model="selected" type="checkbox" :value="p.id" :disabled="applying || (selected.length >= 6 && !selected.includes(p.id))"/><span><strong>{{ p.assay }} · {{ p.cell }} · {{ p.author }}</strong><small>{{ sampleName(p) }}<span v-if="p.antibody"> · {{ p.antibody }}</span></small></span></label><p v-if="!available.length">No signal profiles match the current filters. Reset the source filters above.</p></div></details>
    </component>
    <div v-if="loading" role="status" class="genome-loading" aria-busy="true"><span/>Loading genome annotation and signal tracks…</div>
    <div v-if="error" role="alert" class="genome-error"><p>{{ error }}</p><button type="button" class="button secondary" @click="load">Retry genome browser</button></div>
    <p v-if="emptySelection && !loading" class="genome-empty">No signal tracks selected. Open “Select signal tracks”, choose profiles and apply.</p>
    <div v-show="manifest" class="genome-scroll"><div ref="host" class="genome-host"/></div>
    <div v-if="activeSummary.length && !compact" class="genome-track-key"><div v-for="p in activeSummary" :key="p.id"><span :style="{ backgroundColor: colors[p.assay] }"/><p><strong>{{ sampleName(p) }} · {{ p.assay }} · {{ p.cell }} · {{ p.author }}</strong><small>{{ p.gse }}<template v-if="p.antibody"> · {{ p.antibody }}</template> · {{ p.unit }} · {{ activeStrand }} strand</small></p></div></div>
    <div v-if="compact && activeSummary.length" class="genome-preview-key"><span v-for="p in activeSummary" :key="p.id" :data-profile="p.id"><i :style="{background:colors[p.assay]}"/>{{p.assay}} · {{p.cell}}</span></div>
    <footer v-if="manifest && !compact" class="genome-note">GRO/PRO-seq: primary reads with MAPQ ≥ 20, using the verified library direction and read selection. mNET-seq: weighted 3′ ends from processed BED6. Signals are raw, binned means with independent track scales; peak heights across assays are not directly comparable. Zoom in to view regions below 300 kb.</footer>
  </section>
</template>

<style scoped>
.gene-genome-browser{margin-top:28px;border:1px solid #dce5ef;border-radius:10px;background:#fff;min-width:0;overflow:hidden;color:#345477}.genome-heading{padding:24px;display:flex;justify-content:space-between;align-items:center;gap:16px;border-bottom:1px solid #e4ebf3}.genome-heading h2{font-size:23px;font-weight:600;line-height:1.4;text-wrap:balance;margin:0 0 6px}.genome-heading p{font-size:14px;color:#6a8099;margin:0;text-wrap:pretty}.genome-controls{padding:18px 24px}.genome-endpoints{display:flex;gap:8px;flex-wrap:wrap;margin:0 0 12px}.genome-endpoints .button{font-size:12px}.genome-reference{font-size:13px;color:#647b95;line-height:1.8;margin:0 0 12px;text-wrap:pretty}.genome-actions{display:flex;gap:14px;align-items:center;flex-wrap:wrap}.genome-actions>label{font-size:14px;display:flex;gap:10px;align-items:center}.genome-actions select{border:1px solid #d4dfeb;border-radius:6px;padding:8px;background:white;color:#345477;font:inherit}.genome-actions>span{font-size:13px;color:#647b95;font-variant-numeric:tabular-nums}.genome-picker{margin-top:16px}.genome-picker summary{font-size:14px;font-weight:500;cursor:pointer;padding:9px 0}.genome-picker>p{font-size:13px;color:#647b95}.genome-options{max-height:250px;overflow:auto;margin-top:10px;border:1px solid #e4ebf3;border-radius:7px}.genome-options>label{display:flex;align-items:center;gap:12px;padding:11px 14px;border-bottom:1px solid #eef2f7;font-size:13px;cursor:pointer;min-height:48px}.genome-options input{accent-color:#4979ad;width:16px;height:16px;flex-shrink:0}.genome-options strong{font-weight:500}.genome-options small{display:block;color:#73869a;font-size:12px;overflow-wrap:anywhere;margin-top:4px}.genome-options>p{padding:16px;font-size:14px}.genome-scroll{overflow-x:auto;max-width:100%;padding:0 12px 12px}.genome-host{min-width:650px;min-height:280px;font-family:Arial,sans-serif;color:#333}.genome-loading{margin:12px 24px;padding:16px;font-size:14px;color:#647b95;background:#f5f8fc;border-radius:6px}.genome-loading>span{display:block;height:90px;background:#eaf0f7;margin-bottom:12px;border-radius:4px}.genome-error{margin:12px 24px;padding:16px;border:1px solid #ead0d0;border-radius:7px;color:#8d4242;font-size:14px;overflow-wrap:anywhere}.genome-error p{margin:0 0 12px}.genome-empty{font-size:14px;margin:0 24px 16px;color:#647b95}.genome-track-key{padding:12px 24px;border-top:1px solid #e4ebf3;display:grid;gap:10px}.genome-track-key>div{display:flex;gap:10px;align-items:flex-start}.genome-track-key>div>span{width:8px;height:8px;border-radius:50%;flex-shrink:0;margin-top:7px}.genome-track-key p{margin:0;font-size:13px;line-height:1.7}.genome-track-key strong{font-weight:500}.genome-track-key small{display:block;font-size:12px;color:#647b95;overflow-wrap:anywhere}.genome-note{padding:16px 24px;background:#f7f9fc;border-top:1px solid #e4ebf3;font-size:12px;line-height:1.9;color:#647b95;text-wrap:pretty}.gene-genome-browser button:disabled{opacity:.55;cursor:wait}.gene-genome-browser :focus-visible{outline:2px solid #4979ad;outline-offset:3px}@media(max-width:760px){.genome-heading{align-items:flex-start;flex-direction:column;padding:20px}.genome-heading h2{font-size:21px}.genome-controls{padding:16px 20px}.genome-actions{gap:12px}.genome-actions>span{flex-basis:100%}.genome-note,.genome-track-key{padding:16px 20px}.genome-error{margin:12px 20px}}
.gene-genome-browser>.genome-heading{background:#204f78;border-bottom-color:#204f78}.gene-genome-browser>.genome-heading h2{color:#fff}.gene-genome-browser>.genome-heading p{color:#dce8f4}
.genome-preview .genome-heading{padding:18px 24px}.genome-preview .genome-heading h2{font-size:20px}.genome-preview .genome-heading p{font-size:13px}
.genome-preview .genome-controls{padding:12px 24px}.genome-preview .genome-controls>summary{cursor:pointer;font-size:13px;font-weight:500}.genome-preview .genome-controls[open]>summary{margin-bottom:14px}
.genome-preview .genome-host{min-height:0}.genome-preview .genome-scroll{padding-bottom:6px}.preview-reset{margin-top:12px;font-size:12px}
.genome-preview-key{display:flex;gap:18px;flex-wrap:wrap;padding:10px 24px;border-top:1px solid #e4ebf3;font-size:12px;color:#345477}.genome-preview-key span{display:inline-flex;align-items:center;gap:7px}.genome-preview-key i{width:8px;height:8px;border-radius:50%;flex-shrink:0}
@media(max-width:760px){.genome-preview .genome-heading{padding:16px 20px}.genome-preview .genome-heading h2{font-size:18px}.genome-preview-key{padding:10px 20px;gap:12px}}
</style>
