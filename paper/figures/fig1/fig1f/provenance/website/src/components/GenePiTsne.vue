<script setup>
import { sampleName, sampleNameColumn } from '../utils/sampleName.js'
import { computed, ref, watch, onBeforeUnmount } from 'vue'
import { Download } from '@lucide/vue'
import { getGenePiTsne } from '../api/client.js'
import { exportTsv } from '../utils/data.js'
import PcaScatter from './PcaScatter.vue'
const props = defineProps({ catalog: Object, symbol: String, transcript: String, rows: Array, allRows: Array })
const analysis = ref(null), loading = ref(false), error = ref(''), activeId = ref('')
let generation = 0, timer, wake
function cancel() { generation++; clearTimeout(timer); wake?.(); wake = null }
function validate(data) {
  if (data.algorithm !== 't-SNE' || data.scope !== 'single_gene_pi' || data.metric !== 'log2_pi' || data.feature_count !== 1 || data.symbol !== props.symbol || data.transcript !== props.transcript || data.catalog_release !== props.catalog.release) throw new Error('t-SNE input does not match this gene and PI catalog.')
  if (Object.keys(props.catalog.source_hashes).some(key => data.source_hashes?.[key] !== props.catalog.source_hashes[key])) throw new Error('t-SNE sources differ from the current PI matrices.')
  if (data.status === 'unavailable') return data
  if (data.status !== 'complete' || !Array.isArray(data.points) || data.points.length !== props.allRows.length || data.profile_count !== data.points.length || !data.parameters || !data.preprocessing?.input) throw new Error('t-SNE result is incomplete.')
  const values = new Map(props.allRows.map(row => [row.id, row.log2_pi])), unique = new Set()
  for (const point of data.points) {
    if (!values.has(point.id) || unique.has(point.id) || values.get(point.id) !== point.log2_pi || !Number.isFinite(point.x) || !Number.isFinite(point.y)) throw new Error('t-SNE coordinates or original PI values are invalid.')
    unique.add(point.id)
  }
  return data
}
async function load() {
  cancel()
  const id = generation
  analysis.value = null; error.value = ''; activeId.value = ''; loading.value = true
  try {
    for (let attempt = 0; attempt < 60; attempt++) {
      const data = await getGenePiTsne(props.symbol, props.transcript)
      if (id !== generation) return
      if (data.status === 'computing') {
        if (data.symbol !== props.symbol || data.transcript !== props.transcript) throw new Error('Analysis job does not match this transcript.')
        await new Promise(resolve => { wake = resolve; timer = setTimeout(() => { wake = null; resolve() }, 1500) })
        if (id !== generation) return
        continue
      }
      analysis.value = validate(data)
      return
    }
    throw new Error('t-SNE is taking longer than expected. Please retry shortly.')
  } catch (failure) { if (id === generation) error.value = failure.message }
  finally { if (id === generation) loading.value = false }
}
watch(() => [props.symbol, props.transcript, props.catalog], load, { immediate: true })
watch(() => props.rows, () => { activeId.value = '' })
onBeforeUnmount(cancel)
const points = computed(() => {
  const metadata = new Map(props.catalog.profiles.map(row => [row.id, row]))
  return (analysis.value?.points || []).map(point => ({ ...metadata.get(point.id), ...point, pc1: point.x, pc2: point.y }))
})
const visibleIds = computed(() => props.rows.map(row => row.id))
const active = computed(() => points.value.find(point => point.id === activeId.value))
function download() {
  exportTsv(`${props.symbol}_${props.transcript}_PI_tSNE.tsv`, [
    { key: 'id', label: 'profile_id' }, sampleNameColumn, { key: 'alias', label: 'source_alias' }, { key: 'assay', label: 'assay' },
    { key: 'cell', label: 'raw_cell_label' }, { key: 'antibody', label: 'antibody' }, { key: 'gse', label: 'study' },
    { key: 'log2_pi', label: 'log2_PI' }, { key: 'x', label: 'tSNE_1' }, { key: 'y', label: 'tSNE_2' },
  ], points.value.filter(point => visibleIds.value.includes(point.id)))
}
</script>
<template>
  <section class="gene-pi-tsne" aria-labelledby="gene-tsne-title" :aria-busy="loading">
    <div class="tsne-heading"><div><h2 id="gene-tsne-title">{{ symbol }} · PI t-SNE</h2><p>{{ transcript }} · Independent embedding of this gene's measured PI values.</p></div><button v-if="analysis?.status === 'complete'" class="button secondary" type="button" :disabled="!rows.length" @click="download"><Download :size="16"/>Download t-SNE</button></div>
    <div v-if="error" class="error-state" role="alert"><strong>Unable to load PI t-SNE</strong><p>{{ error }}</p><button class="button secondary" type="button" @click="load">Retry t-SNE</button></div>
    <div v-else-if="loading" class="tsne-pending" role="status"><strong>Computing this gene's PI t-SNE…</strong><p>Using observed log₂ PI values. The result is cached for future visits.</p></div>
    <div v-else-if="analysis?.status === 'unavailable'" class="tsne-pending"><strong>Insufficient PI variation for t-SNE</strong><p>{{ analysis.reason }}</p></div>
    <template v-else-if="analysis?.status === 'complete'">
      <p class="tsne-context">{{ analysis.profile_count }} measured profiles · 1 PI feature · {{ rows.length }} shown · Perplexity {{ analysis.parameters.perplexity }} · Seed {{ analysis.parameters.random_state }}</p>
      <p class="tsne-input-note">Input: this gene's log₂ PI only. Colors identify assay or cell labels; t-SNE does not assign clusters. Use the PI distributions below to check the original differences.</p>
      <p v-if="!rows.length" class="tsne-no-matches">No measured profiles match these filters. Reset the filters to show points.</p>
      <div class="tsne-grid"><PcaScatter :points="points" :visible-ids="visibleIds" algorithm="t-SNE" group-by="assay" :active-id="activeId" @inspect="activeId = $event"/><PcaScatter :points="points" :visible-ids="visibleIds" algorithm="t-SNE" group-by="cell" :active-id="activeId" @inspect="activeId = $event"/></div>
      <div class="tsne-inspector"><template v-if="active"><strong>{{ sampleName(active) }}</strong><div><span>{{ active.assay }}</span><span>{{ active.cell }}</span><span>{{ active.antibody || 'No antibody subtype' }}</span><span>log₂ PI: {{ active.log2_pi.toFixed(4) }}</span><RouterLink :to="`/samples/${active.id}`">View profile ↗</RouterLink></div></template><template v-else><strong>Compare PI neighbors</strong><p>Hover or focus a sample to highlight it in both views. Select a legend group to highlight an assay or cell line.</p></template></div>
      <details class="tsne-methods"><summary>t-SNE methods and interpretation</summary><dl><dt>Input</dt><dd>{{ analysis.preprocessing.input }}</dd><dt>Scaling</dt><dd>{{ analysis.preprocessing.scaling }}</dd><dt>Missing values</dt><dd>{{ analysis.preprocessing.missing_values }}</dd><dt>Initialization</dt><dd>{{ analysis.preprocessing.initialization }}</dd><dt>Settings</dt><dd>Exact Euclidean t-SNE · {{ analysis.parameters.max_iter }} maximum iterations · Learning rate: {{ analysis.parameters.learning_rate }} · Early exaggeration: {{ analysis.parameters.early_exaggeration }}</dd><dt>Interpretation</dt><dd>{{ analysis.preprocessing.interpretation }}</dd><dt>Study effects</dt><dd>{{ analysis.preprocessing.batch_correction }}</dd><dt>Neighborhood check</dt><dd>Trustworthiness ({{ analysis.quality.trustworthiness_neighbors }} neighbors): {{ analysis.quality.trustworthiness.toFixed(3) }} · KL divergence: {{ analysis.quality.kl_divergence.toFixed(4) }}. These describe the embedding, not biological cluster validity.</dd></dl><p>Both charts use the same coordinates fitted from all measured profiles for this transcript. Filters only change visible samples; they do not refit the embedding.</p></details>
    </template>
  </section>
</template>
<style scoped>
.gene-pi-tsne{margin:26px 0 36px;min-width:0}.tsne-heading{display:flex;align-items:center;justify-content:space-between;gap:16px;margin-bottom:14px}.tsne-heading h2{font-size:24px;font-weight:600;color:#204f78;line-height:1.4;margin:0;text-wrap:balance}.tsne-heading p{font-size:13px;color:#64748b;line-height:1.7;margin:7px 0 0;text-wrap:pretty}.tsne-context{font-size:13px;color:#64748b;font-variant-numeric:tabular-nums;margin:0 0 12px}.tsne-input-note{font-size:13px;color:#475569;background:#f0f6fc;border-left:3px solid #7398bb;padding:12px 16px;margin:0 0 18px;line-height:1.8;text-wrap:pretty}.tsne-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:22px;min-width:0}.tsne-pending{padding:30px 24px;border:1px solid #e4eaf1;border-radius:10px;background:#f8fafc;min-height:150px}.tsne-pending strong{font-size:16px;font-weight:500;color:#334155}.tsne-pending p,.tsne-no-matches{font-size:13px;color:#64748b;line-height:1.8;text-wrap:pretty}.tsne-inspector{border:1px solid #e4eaf1;border-radius:10px;background:#f8fafc;min-height:104px;padding:18px 22px;margin-top:16px;overflow-wrap:anywhere}.tsne-inspector strong{font-size:14px;font-weight:500;color:#334155}.tsne-inspector>div{display:flex;flex-wrap:wrap;gap:8px 16px;margin-top:8px;font-size:13px;color:#64748b;font-variant-numeric:tabular-nums}.tsne-inspector a{color:#3e77ab}.tsne-inspector p{font-size:13px;color:#64748b;line-height:1.7;margin:8px 0 0;text-wrap:pretty}.tsne-methods{font-size:13px;color:#64748b;line-height:1.8;margin-top:12px}.tsne-methods summary{cursor:pointer;color:#3e77ab;padding:8px 0}.tsne-methods dl{display:grid;grid-template-columns:150px minmax(0,1fr);gap:10px 16px}.tsne-methods dt{font-weight:500;color:#334155}.tsne-methods dd{margin:0;text-wrap:pretty}.tsne-methods p{text-wrap:pretty}
@media(max-width:860px){.tsne-grid{grid-template-columns:1fr}.tsne-heading{align-items:flex-start;flex-wrap:wrap}.tsne-heading h2{font-size:22px}.tsne-methods dl{grid-template-columns:1fr;gap:4px}.tsne-methods dd{margin-bottom:8px}}
</style>
