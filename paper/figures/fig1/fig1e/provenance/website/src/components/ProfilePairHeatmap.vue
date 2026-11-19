<script setup>
import { computed, onBeforeUnmount, ref, watch } from 'vue'
import { getPiSpearmanHeatmap, validatePiAnalysis } from '../api/client.js'
import { downloadBlob } from '../utils/data.js'
import { sampleName } from '../utils/sampleName.js'
import LoadState from './LoadState.vue'

const props = defineProps({ profile: Object, catalog: Object, comparisonId: { type: String, default: '' } })
const emit = defineEmits(['update:comparisonId'])
const analysis = ref(null), loading = ref(true), error = ref(''), hover = ref(null), selected = ref([0, 0])
let generation = 0
async function load() {
  const token = ++generation
  loading.value = true; error.value = ''; analysis.value = null
  try {
    const result = await getPiSpearmanHeatmap()
    validatePiAnalysis(result, props.catalog)
    const size = result.profile_order.length
    if (result.method?.correlation !== 'Spearman' || !Number.isInteger(result.method.minimum_overlap)
      || result.spearman?.length !== size || result.overlaps?.length !== size
      || result.spearman.some(row => !Array.isArray(row) || row.length !== size || row.some(r => r !== null && (!Number.isFinite(r) || Math.abs(r) > 1)))
      || result.overlaps.some(row => !Array.isArray(row) || row.length !== size || row.some(n => !Number.isInteger(n) || n < 0))) {
      throw Error('Correlation data do not match this catalog. Reload the page.')
    }
    if (token === generation) analysis.value = result
  } catch (e) { if (token === generation) error.value = e.message }
  finally { if (token === generation) loading.value = false }
}
watch(() => props.catalog, load, { immediate: true })
onBeforeUnmount(() => { ++generation })
const comparison = computed(() => props.catalog.profiles.find(p => p.id === props.comparisonId))
const rows = computed(() => [props.profile, ...(comparison.value && comparison.value.id !== props.profile.id ? [comparison.value] : [])])
const index = computed(() => new Map((analysis.value?.profile_order || []).map((id, i) => [id, i])))
function pairAt(row, column) {
  const a = rows.value[row], b = rows.value[column], i = index.value.get(a.id), j = index.value.get(b.id)
  return { a, b, r: analysis.value?.spearman[i]?.[j] ?? null, n: analysis.value?.overlaps[i]?.[j] ?? 0 }
}
watch([() => props.profile.id, () => props.comparisonId], () => { hover.value = null; selected.value = [0, rows.value.length - 1] }, { immediate: true })
const pair = computed(() => pairAt(...(hover.value || selected.value)))
const number = value => value === null ? 'Unavailable' : value.toFixed(3)
function color(r) {
  if (r === null) return '#e4e8ed'
  const t = Math.min(1, Math.abs(r)), light = [239, 244, 249], dark = [83, 116, 148]
  return `rgb(${light.map((v, i) => Math.round(v + (dark[i] - v) * t)).join(',')})`
}
function ink(r) {
  if (r === null) return '#344b61'
  const [red, green, blue] = color(r).match(/\d+/g).map(v => {
    const c = Number(v) / 255
    return c <= .04045 ? c / 12.92 : ((c + .055) / 1.055) ** 2.4
  })
  const luminance = .2126 * red + .7152 * green + .0722 * blue
  if (1.05 / (luminance + .05) >= 4.5) return '#fff'
  return (luminance + .05) / .082 >= 4.5 ? '#18344c' : '#000'
}
const legend = Array.from({ length: 40 }, (_, i) => color(i / 39))
function label(row, column) {
  const p = pairAt(row, column)
  return `${sampleName(p.a)} versus ${sampleName(p.b)}: signed Spearman rho ${p.r ?? 'unavailable'}, ${p.n} shared records`
}
function download() {
  const fields = ['profile_a', 'sample_name_a', 'profile_b', 'sample_name_b', 'signed_spearman_rho', 'abs_spearman_rho', 'shared_records']
  const lines = rows.value.flatMap((a, i) => rows.value.map((b, j) => {
    const p = pairAt(i, j)
    return [a.id, sampleName(a), b.id, sampleName(b), p.r ?? '', p.r === null ? '' : Math.abs(p.r), p.n]
      .map(v => String(v).replace(/[\t\r\n]/g, ' ')).join('\t')
  }))
  downloadBlob(`${sampleName(props.profile)}_PI_Spearman_comparison.tsv`, [fields.join('\t'), ...lines].join('\n'), 'text/tab-separated-values;charset=utf-8')
}
</script>

<template>
  <section id="profile-correlations" class="panel latest-panel profile-pair-heatmap" aria-label="Profile sample correlation heatmap">
    <div class="panel-head"><div><h2>Sample correlation heatmap</h2><p>Compare shared gene / transcript PI values with another sample.</p></div><button class="button secondary" :disabled="!analysis" @click="download">Download comparison</button></div>
    <div class="pair-heatmap-body">
      <div class="pi-comparison-controls"><label>Comparison sample<select :value="comparisonId" aria-label="Heatmap comparison sample" @change="emit('update:comparisonId', $event.target.value)"><option value="">{{ sampleName(profile) }} · Self-comparison</option><option v-for="p in catalog.profiles.filter(p => p.id !== profile.id)" :key="p.id" :value="p.id">{{ sampleName(p) }} · {{ p.assay }} · {{ p.cell }}</option></select></label><button v-if="comparisonId" class="button secondary" @click="emit('update:comparisonId', '')">Use current sample only</button></div>
      <LoadState :loading="loading" :error="error" @retry="load">
        <template v-if="analysis">
          <p class="pi-method-note pair-mode" role="status">{{ rows.length === 1 ? 'Self-comparison: no other sample selected. This is not an independent replicate.' : 'Both axes show the current sample and the selected comparison sample. The PI distribution and cumulative plots use the same selection.' }}</p>
          <div class="pair-heatmap-layout">
            <div>
              <div class="pair-matrix-scroll" role="region" tabindex="0" aria-label="Scrollable profile correlation matrix">
                <table class="pair-matrix" aria-label="Sample correlation matrix">
                  <caption>Color and cell values: |Spearman ρ|</caption>
                  <thead><tr><th scope="col">Row / column sample</th><th v-for="p in rows" :key="p.id" scope="col"><RouterLink :to="'/samples/' + p.id">{{ sampleName(p) }}</RouterLink><small>{{ p.assay }} · {{ p.cell }}</small></th></tr></thead>
                  <tbody><tr v-for="(a, i) in rows" :key="a.id"><th scope="row"><RouterLink :to="'/samples/' + a.id">{{ sampleName(a) }}</RouterLink><small>{{ a.assay }} · {{ a.cell }}</small></th><td v-for="(b, j) in rows" :key="b.id"><button type="button" :data-row-profile="a.id" :data-column-profile="b.id" :data-rho="pairAt(i, j).r" :data-shared="pairAt(i, j).n" :style="{ background: color(pairAt(i, j).r), color: ink(pairAt(i, j).r) }" :aria-label="label(i, j)" :aria-pressed="selected[0] === i && selected[1] === j" @mouseenter="hover = [i, j]" @mouseleave="hover = null" @focus="selected = [i, j]; hover = null" @click="selected = [i, j]; hover = null"><span>{{ number(pairAt(i, j).r === null ? null : Math.abs(pairAt(i, j).r)) }}</span><small>{{ pairAt(i, j).n.toLocaleString() }} shared</small></button></td></tr></tbody>
                </table>
              </div>
              <div class="pair-color-key" aria-label="Absolute Spearman correlation scale"><span>|Spearman ρ|</span><div><div class="pair-color-strip"><i v-for="(c, i) in legend" :key="i" :style="{ background: c }"/></div><div class="pair-color-ticks"><span v-for="v in ['0', '0.25', '0.50', '0.75', '1']" :key="v">{{ v }}</span></div></div><span><i class="pair-unavailable"/>Unavailable</span></div>
            </div>
            <aside class="pair-inspection" aria-label="Profile correlation pair details">
              <h3>{{ pair.a.id === pair.b.id ? 'Self-comparison' : 'Sample comparison' }}</h3>
              <dl><div><dt>Signed Spearman ρ</dt><dd data-pair-rho>{{ number(pair.r) }}</dd></div><div><dt>Shared gene / transcript records</dt><dd data-pair-shared>{{ pair.n.toLocaleString() }}</dd></div></dl>
              <p><b>Row:</b> {{ sampleName(pair.a) }}<br><b>Column:</b> {{ sampleName(pair.b) }}</p>
              <p v-if="pair.r === null">Correlation is unavailable because there are fewer than {{ analysis.method.minimum_overlap }} shared records or a constant PI vector.</p>
              <p v-else-if="pair.a.id === pair.b.id">The same measured PI values are compared with themselves.</p>
              <p v-else>Color shows correlation magnitude; the signed value preserves positive or inverse relationships.</p>
            </aside>
          </div>
          <p class="pi-method-note">Original log₂ PI · Spearman ranks on records measured in both samples · minimum {{ analysis.method.minimum_overlap }} shared records. Missing values are excluded, true zeros retained. Hover, focus or select a cell to inspect it.</p>
        </template>
      </LoadState>
      <RouterLink class="table-link" to="/statistics#sample-correlations">Open the full sample correlation heatmap →</RouterLink>
    </div>
  </section>
</template>

<style scoped>
.profile-pair-heatmap{min-width:0;font-variant-numeric:tabular-nums}.profile-pair-heatmap .panel-head{gap:16px;flex-wrap:wrap}.profile-pair-heatmap h2,.pair-inspection h3{text-wrap:balance}.pair-heatmap-body{padding:22px}.pair-mode{margin:0 0 18px}.pair-heatmap-layout{display:grid;grid-template-columns:minmax(0,1.5fr) minmax(230px,1fr);gap:24px;align-items:start}.pair-matrix-scroll{max-width:100%;overflow:auto}.pair-matrix{width:100%;table-layout:fixed;border-collapse:separate;border-spacing:4px;font-size:12px}.pair-matrix caption{text-align:left;margin-bottom:10px;color:#627b9b}.pair-matrix th,.pair-matrix td{padding:0;border:0;background:transparent;white-space:normal;vertical-align:middle}.pair-matrix th{color:#244e78;overflow-wrap:anywhere;font-weight:500;line-height:1.7;text-align:center}.pair-matrix th:first-child{width:36%;text-align:left;padding-right:10px}.pair-matrix thead th{padding-bottom:12px}.pair-matrix th small{display:block;color:#627b9b;font-size:11px;margin-top:4px}.pair-matrix td button{display:flex;flex-direction:column;align-items:center;justify-content:center;width:100%;min-height:100px;gap:8px;border:2px solid transparent;border-radius:5px;font:inherit;cursor:pointer;line-height:1.5}.pair-matrix td button span{font-size:21px;font-weight:600}.pair-matrix td button small{font-size:11px}.pair-matrix td button:hover{box-shadow:inset 0 0 0 2px currentColor}.pair-matrix td button[aria-pressed=true]{border-color:#18344c}.profile-pair-heatmap :focus-visible{outline:2px solid #417dc6;outline-offset:2px}.pair-color-key{display:flex;gap:12px;align-items:center;flex-wrap:wrap;font-size:11px;color:#627b9b;margin-top:16px}.pair-color-key>div{width:200px;max-width:100%}.pair-color-strip{display:flex;height:10px}.pair-color-strip i{flex:1}.pair-color-ticks{display:flex;justify-content:space-between;margin-top:5px}.pair-unavailable{display:inline-block;width:11px;height:11px;background:#e4e8ed;margin-right:5px;vertical-align:middle}.pair-inspection{background:#f5f8fc;border:1px solid #dce5ef;border-radius:8px;padding:18px;overflow-wrap:anywhere}.pair-inspection h3{font-size:15px;margin:0 0 18px;color:#244e78}.pair-inspection dl{margin:0;display:grid;gap:18px}.pair-inspection dt{font-size:12px;color:#627b9b}.pair-inspection dd{font-size:24px;margin:5px 0 0;color:#244e78}.pair-inspection p{font-size:12px;color:#627b9b;line-height:1.8;text-wrap:pretty}.pair-heatmap-body>.pi-method-note{margin-top:18px}.pair-heatmap-body>.table-link{display:inline-block;margin-top:12px;font-size:12px}@media(max-width:900px){.pair-heatmap-layout{grid-template-columns:1fr}.pair-inspection dl{grid-template-columns:1fr 1fr}}@media(max-width:600px){.pair-heatmap-body{padding:14px}.pair-matrix td button{min-height:84px}.pair-matrix td button span{font-size:18px}.pair-matrix td button small{font-size:10px}.pair-matrix th{font-size:11px}.pair-matrix th small{font-size:10px}.pair-inspection{padding:14px}.pair-inspection dd{font-size:21px}}
</style>
