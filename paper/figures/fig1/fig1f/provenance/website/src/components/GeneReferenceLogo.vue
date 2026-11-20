<script setup>
import { computed, onBeforeUnmount, ref, watch } from 'vue'
import { sampleName } from '../utils/sampleName.js'
import { downloadBlob } from '../utils/data.js'
import { alignBaseSignal, baseOffset, signalCeiling, baseGlyphs } from '../utils/baseSignal.js'
const props = defineProps({ data: { type: Object, required: true }, selection: { type: Object, default: () => ({}) } })
const emit = defineEmits(['update:selection'])
const start = ref(0), position = ref(1), message = ref(''), count = ref(41), host = ref(null), active = ref(null)
const manifest = ref(null), assay = ref(props.selection.assay || 'PRO-seq'), profileId = ref(props.selection.profileId || ''), cell = ref(props.selection.cell || ''), strandMode = ref(props.selection.strandMode || 'gene'), records = ref(null)
const sourceLoading = ref(false), signalLoading = ref(false), sourceError = ref(''), signalError = ref('')
let observer, manifestController, signalController, manifestRevision = 0, signalRevision = 0
const colors = { A: '#27815e', T: '#bc5260', G: '#997322', C: '#3978ad' }
const availableAssays = computed(() => ['PRO-seq', 'GRO-seq', 'mNET-seq'].filter(a => manifest.value?.profiles.some(p => p.assay === a)))
const profiles = computed(() => (manifest.value?.profiles || []).filter(p => p.assay === assay.value))
const preferred = { 'PRO-seq': 'NP-663C656923', 'GRO-seq': 'NP-410E006606', 'mNET-seq': 'NP-3E3F5FDE87' }
const cells = computed(() => [...new Set(profiles.value.map(p => p.cell))].sort())
const candidates = computed(() => profiles.value.filter(p => !cell.value || p.cell === cell.value))
const profile = computed(() => profiles.value.find(p => p.id === profileId.value))
const signalStrand = computed(() => strandMode.value === 'gene' ? props.data.gene_strand : props.data.gene_strand === '+' ? '-' : '+')
const busy = computed(() => sourceLoading.value || signalLoading.value)
const width = computed(() => 86 + count.value * 18 + 15), baseline = 120, plotHeight = 100
const ceiling = computed(() => signalCeiling((records.value || []).map(r => r.value)))
const letters = computed(() => records.value?.slice(start.value, start.value + count.value) || [])
const fragment = computed(() => props.data.sequence.slice(start.value, start.value + count.value))
const signed = value => value > 0 ? '+' + value.toLocaleString('en-US') : value.toLocaleString('en-US')
const label = i => props.data.endpoint ? signed(baseOffset(props.data, i)) : String(baseOffset(props.data, i) + 1)
const span = computed(() => `${label(start.value)} to ${label(start.value + fragment.value.length - 1)} ${props.data.endpoint ? 'relative to ' + props.data.endpoint.anchor_name : 'gene bp'}`)
const height = value => value / ceiling.value * plotHeight
const format = value => value.toLocaleString('en-US', { maximumFractionDigits: 3 })
const zeroRegion = computed(() => records.value && records.value.every(r => r.value === 0))
async function request(action, params, controller) {
  const timeout = setTimeout(() => controller.abort(), 50000)
  try { const r = await fetch('/api/v1/genome/' + action + '?' + new URLSearchParams(params), { signal: controller.signal }); const value = await r.json(); if (!r.ok) throw Error(value.error || 'Sequencing signal is unavailable.'); return value }
  finally { clearTimeout(timeout) }
}
async function loadSources() {
  manifestController?.abort(); signalController?.abort()
  const token = ++manifestRevision; ++signalRevision
  const controller = manifestController = new AbortController()
  manifest.value = null; records.value = null; sourceError.value = ''; signalError.value = ''; sourceLoading.value = true; signalLoading.value = false
  try {
    const value = await request('manifest', { symbol: props.data.symbol, transcript: props.data.transcript }, controller)
    if (value.transcript !== props.data.transcript || value.reference?.id !== 'hg38' || value.reference.name !== props.data.assembly) throw Error('Signal reference and displayed MANE sequence differ.')
    if (token !== manifestRevision) return
    manifest.value = value
    if (!availableAssays.value.includes(assay.value)) assay.value = availableAssays.value[0] || ''
    if (!cells.value.includes(cell.value)) cell.value = ''
    if (!candidates.value.some(p => p.id === profileId.value)) profileId.value = candidates.value.find(p => p.id === preferred[assay.value])?.id || candidates.value[0]?.id || ''
  } catch (e) { if (token === manifestRevision) sourceError.value = e.name === 'AbortError' ? 'Signal source loading timed out. Please retry.' : e.message }
  finally { if (token === manifestRevision) sourceLoading.value = false }
}
async function loadSignal() {
  signalController?.abort(); const token = ++signalRevision
  records.value = null; signalError.value = ''; signalLoading.value = false; active.value = null
  const p = profile.value, data = props.data, strand = signalStrand.value
  if (!manifest.value || !p) return
  const controller = signalController = new AbortController(); signalLoading.value = true
  try {
    const value = await request('signal', { profile: p.id, chr: data.chrom, start: data.genomic_start - 1, end: data.genomic_end, bins: data.length, strand }, controller)
    const aligned = alignBaseSignal(data, value, p, strand)
    if (token === signalRevision) records.value = aligned
  } catch (e) { if (token === signalRevision) signalError.value = e.name === 'AbortError' ? 'Sequencing signal loading timed out. Please retry.' : e.message }
  finally { if (token === signalRevision) signalLoading.value = false }
}
watch(() => [props.data.symbol, props.data.transcript], loadSources, { immediate: true })
watch([() => props.data, manifest, profileId, signalStrand], loadSignal)
watch([assay, profileId, cell, strandMode], () => emit('update:selection', { assay: assay.value, profileId: profileId.value, cell: cell.value, strandMode: strandMode.value }))
watch(candidates, () => { if (manifest.value && !candidates.value.some(p => p.id === profileId.value)) profileId.value = candidates.value.find(p => p.id === preferred[assay.value])?.id || candidates.value[0]?.id || '' })
watch(assay, () => { if (!cells.value.includes(cell.value)) cell.value = '' })
function move(value) { start.value = Math.max(0, Math.min(props.data.length - 1, value)); position.value = start.value + 1; message.value = ''; active.value = null }
function jump() { if (!Number.isInteger(position.value) || position.value < 1 || position.value > props.data.length) { message.value = `Enter a window base from 1 to ${props.data.length}.`; return } move(position.value - 1) }
watch([() => props.data, count], () => {
  const anchor = props.data.endpoint?.anchor, anchorIndex = anchor == null ? 0 : props.data.display_strand === '-' ? props.data.genomic_end - anchor : anchor - props.data.genomic_start
  move(Math.max(0, Math.min(props.data.length - count.value, anchorIndex - Math.floor(count.value / 2))))
}, { immediate: true })
watch(host, element => { observer?.disconnect(); if (!element) return; observer = new ResizeObserver(entries => { count.value = entries[0].contentRect.width < 420 ? 21 : 41 }); observer.observe(element) })
function download() {
  const header = 'profile_id\tsample_name\tassay\ttranscript\tchromosome\tgenomic_position_1based\tbase\tregion_offset_bp\tsignal_strand\tsignal_value\tunit'
  const rows = letters.value.map(r => [profile.value.id, sampleName(profile.value), profile.value.assay, props.data.transcript, props.data.chrom, r.coordinate, r.base, r.relative, signalStrand.value, r.value, profile.value.unit].join('\t'))
  downloadBlob(`${props.data.symbol}_${sampleName(profile.value)}_ATCG_signal.tsv`, [header, ...rows].join('\n'), 'text/tab-separated-values;charset=utf-8')
}
onBeforeUnmount(() => { ++manifestRevision; ++signalRevision; manifestController?.abort(); signalController?.abort(); observer?.disconnect() })
</script>
<template>
 <section ref="host" class="gene-reference-logo" aria-label="Signal-weighted ATCG">
  <header><h4>Signal-weighted ATCG</h4><p>MANE reference bases · letter height = selected sequencing signal</p></header>
  <div class="logo-signal-controls"><label>Sequencing assay<select v-model="assay" aria-label="ATCG sequencing assay" :disabled="sourceLoading || !availableAssays.length"><option v-for="a in availableAssays" :key="a">{{ a }}</option></select></label><label>Cell line<select v-model="cell" aria-label="ATCG cell line" :disabled="sourceLoading || !profiles.length"><option value="">All cell lines</option><option v-for="c in cells" :key="c">{{ c }}</option></select></label><label class="logo-profile">Sample<select v-model="profileId" aria-label="ATCG signal sample" :disabled="sourceLoading || !candidates.length"><option v-for="p in candidates" :key="p.id" :value="p.id">{{ sampleName(p) }} · {{ p.cell }}</option></select></label><label>Signal strand<select v-model="strandMode" aria-label="ATCG signal strand"><option value="gene">Gene strand ({{ data.gene_strand }})</option><option value="opposite">Opposite strand ({{ data.gene_strand === '+' ? '−' : '+' }})</option></select></label></div>
  <p class="logo-range">{{ span }} · {{ data.transcript }}</p>
  <p v-if="busy" role="status" class="logo-loading">Loading single-base sequencing signal…</p>
  <div v-else-if="sourceError || signalError" class="logo-error" role="alert"><p>{{ sourceError || signalError }}</p><button class="button secondary" @click="sourceError ? loadSources() : loadSignal()">Retry sequencing signal</button></div>
  <p v-else-if="!profiles.length" class="logo-loading">No registered sequencing signal profiles are available for this reference.</p>
  <template v-else-if="records">
   <p class="logo-source">{{ sampleName(profile) }} · {{ profile.unit }} · {{ signalStrand }} strand · 1 bp resolution</p>
   <div class="logo-plot" @mouseleave="active = null">
    <svg class="reference-logo-svg" :viewBox="`0 0 ${width} 178`" role="img" :aria-label="`Signal-weighted reference bases: ${fragment}. ${span}. ${sampleName(profile)}. Signal strand ${signalStrand}.`" :data-sequence="fragment" :data-start="start" :data-profile="profileId" :data-strand="signalStrand" :data-ceiling="ceiling">
     <title>Reference bases scaled by measured sequencing signal</title><desc>The height of each letter is proportional to the raw strand-specific signal at that base. Zero signal has zero height. The small baseline row retains base identity.</desc>
     <g v-for="fraction in [0, .5, 1]" :key="fraction"><line x1="86" :y1="baseline - fraction * plotHeight" :x2="width - 15" :y2="baseline - fraction * plotHeight" stroke="#e5edf5"/><text x="77" :y="baseline - fraction * plotHeight + 4" text-anchor="end" class="logo-tick">{{ format(ceiling * fraction) }}</text></g>
     <text transform="translate(19 72) rotate(-90)" text-anchor="middle" class="logo-tick">{{ profile.assay === 'mNET-seq' ? 'Weighted 3′ ends' : 'Aligned-read depth' }}</text>
     <g v-for="(letter, i) in letters" :key="letter.index" class="signal-base" :data-coordinate="letter.coordinate" :data-base="letter.base" :data-signal="letter.value" :data-height="height(letter.value)" tabindex="0" role="img" :aria-label="`${letter.base}, ${data.chrom}:${letter.coordinate}, ${profile.assay} signal ${letter.value}`" @mouseenter="active = letter" @focus="active = letter" @blur="active = null">
      <title>{{ data.chrom }}:{{ letter.coordinate }} · {{ letter.base }} · {{ profile.assay }} signal {{ letter.value }}</title>
      <rect :x="86 + i * 18" y="16" width="18" height="108" fill="transparent"/>
      <path v-if="letter.value > 0 && baseGlyphs[letter.base]" :d="baseGlyphs[letter.base]" fill-rule="evenodd" :fill="colors[letter.base]" :transform="`translate(${88 + i * 18} ${baseline - height(letter.value)}) scale(.14 ${height(letter.value) / 100})`"/>
      <text :x="95 + i * 18" y="139" text-anchor="middle" :fill="colors[letter.base] || '#64748b'" class="logo-base-label">{{ letter.base }}</text>
      <text v-if="i % 10 === 0 || i === letters.length - 1" :x="95 + i * 18" y="161" text-anchor="middle" class="logo-tick">{{ label(letter.index) }}</text>
     </g>
    </svg>
    <div v-if="active" class="logo-signal-tip" role="status"><b>{{ active.base }} · {{ data.chrom }}:{{ active.coordinate.toLocaleString() }}</b><span>{{ profile.assay }} signal: {{ format(active.value) }} · {{ signalStrand }} strand</span></div>
   </div>
   <p v-if="zeroRegion" class="logo-zero">Zero sequencing signal across this selected region for the chosen sample and strand.</p>
   <div class="logo-navigation"><button class="button secondary" aria-label="Previous logo window" :disabled="start === 0" @click="move(start - count)">Previous</button><form @submit.prevent="jump"><label>Window bp <input v-model.number="position" type="number" min="1" :max="data.length" step="1" aria-label="Logo start at window base"/></label><button class="button secondary" aria-label="Go to logo window">Go</button></form><button class="button secondary" aria-label="Next logo window" :disabled="start + count >= data.length" @click="move(start + count)">Next</button></div>
   <p v-if="message" role="status">{{ message }}</p>
   <div class="logo-signal-foot"><p>Raw signal without normalization. The vertical scale uses the entire selected region and stays fixed while paging. The signal selection is independent of the PI filters.</p><button class="button secondary" @click="download">Download displayed signal</button></div>
  </template>
 </section>
</template>
<style scoped>
.gene-reference-logo{padding:16px 20px;border-top:1px solid #e4eaf1;min-width:0}.gene-reference-logo header{background:#204f78;margin:-16px -20px 16px;padding:16px 20px}.gene-reference-logo h4{margin:0;color:white;font-size:17px;font-weight:600;text-wrap:balance}.gene-reference-logo p{font-size:12px;line-height:1.7;color:#64748b;text-wrap:pretty;margin:5px 0}.gene-reference-logo header p{color:#dce8f4}.logo-range,.logo-source{font-variant-numeric:tabular-nums;overflow-wrap:anywhere}.logo-signal-controls{display:flex;gap:12px;flex-wrap:wrap;margin-bottom:12px}.logo-signal-controls label{display:flex;flex-direction:column;gap:6px;font-size:12px;color:#345776;min-width:0;flex:1 1 160px}.logo-signal-controls .logo-profile{flex:3 1 320px}.logo-signal-controls select{width:100%;min-height:36px;padding:6px 8px;border:1px solid #d5dfeb;border-radius:6px;background:white;color:#345776;font:inherit}.logo-plot{position:relative}.reference-logo-svg{display:block;width:100%;height:auto;max-height:240px;margin:14px 0}.logo-tick{font-family:inherit;font-size:12px;fill:#64748b;font-variant-numeric:tabular-nums}.logo-base-label{font-family:Consolas,ui-monospace,monospace;font-size:12px;font-weight:600}.signal-base{cursor:crosshair}.signal-base:focus-visible{outline:none}.signal-base:focus-visible rect,.signal-base:hover rect{stroke:#a6bfd8;stroke-width:1;fill:#edf3f930}.logo-signal-tip{position:absolute;right:12px;top:0;padding:8px 12px;background:white;border:1px solid #d5dfeb;border-radius:5px;box-shadow:0 3px 10px #0f172a20;font-size:12px;line-height:1.7;color:#345776;pointer-events:none;max-width:90%;font-variant-numeric:tabular-nums}.logo-signal-tip span{display:block}.logo-navigation{display:flex;gap:8px;align-items:center;justify-content:space-between;flex-wrap:wrap}.logo-navigation form,.logo-navigation label{display:flex;align-items:center;gap:6px}.logo-navigation label{font-size:12px;color:#64748b}.logo-navigation input{width:65px;border:1px solid #d5dfeb;border-radius:6px;padding:6px;min-height:34px;color:#345776;font:inherit}.gene-reference-logo .button{font-size:12px;padding:6px 10px;min-height:34px}.logo-signal-foot{display:flex;gap:16px;flex-wrap:wrap;align-items:center;justify-content:space-between;margin-top:12px}.logo-signal-foot p{flex:1 1 400px}.logo-loading,.logo-error{padding:14px;background:#f5f8fc;border-radius:6px}.logo-error{border:1px solid #ead0d0}.gene-reference-logo button:disabled{opacity:.4;cursor:default}.gene-reference-logo :focus-visible{outline:2px solid #417dc6;outline-offset:2px}@media(max-width:600px){.logo-navigation form{order:3;flex-basis:100%;justify-content:center}.logo-tick{font-size:14px}.logo-signal-controls .logo-profile{flex-basis:100%}.logo-signal-tip{font-size:11px}.logo-signal-foot .button{width:100%}}
</style>
