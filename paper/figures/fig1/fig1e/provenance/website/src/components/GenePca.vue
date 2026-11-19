<script setup>
import { sampleName, sampleNameColumn } from '../utils/sampleName.js'
import { computed, ref, watch } from 'vue'
import { Download, ScatterChart } from '@lucide/vue'
import { getProfilePca } from '../api/client.js'
import { exportTsv } from '../utils/data.js'
import PcaScatter from './PcaScatter.vue'
const props = defineProps({ catalog: { type: Object, required: true } })
const analysis = ref(null), loading = ref(false), error = ref(''), activeId = ref('')
let requestId = 0
function validate(data) {
  if (data.scope !== 'all_gene_pi' || data.metric !== 'log2_pi' || data.catalog_release !== props.catalog.release) throw new Error('PCA data does not match the current PI catalog.')
  if (!Array.isArray(data.points)) throw new Error('PCA profile coordinates are missing.')
  if (JSON.stringify(data.profile_order) !== JSON.stringify(props.catalog.metric.profile_order) || data.points.length !== data.profile_order.length) throw new Error('PCA profile order does not match the catalog.')
  if (Object.keys(props.catalog.source_hashes).some(key => data.source_hashes?.[key] !== props.catalog.source_hashes[key])) throw new Error('PCA sources differ from the current PI matrices.')
  if (!Number.isInteger(data.feature_count) || data.feature_count < 2 || !Array.isArray(data.points) || data.points.length < 3) throw new Error('PCA requires at least three profiles and two PI features.')
  if (!Array.isArray(data.explained_variance_ratio) || data.explained_variance_ratio.length !== 2 || data.explained_variance_ratio.some(value => !Number.isFinite(value) || value <= 0 || value > 1) || data.explained_variance_ratio.reduce((a,b)=>a+b,0) > 1.000001) throw new Error('PCA variance information is invalid.')
  const ids = new Set(props.catalog.profiles.map(profile => profile.id)), unique = new Set()
  for (const point of data.points) {
    if (!ids.has(point.id) || unique.has(point.id) || !Number.isFinite(point.pc1) || !Number.isFinite(point.pc2)) throw new Error('PCA coordinates or profile identifiers are invalid.')
    unique.add(point.id)
  }
  if (!data.preprocessing?.input || !data.preprocessing?.scaling || !data.preprocessing?.missing_values) throw new Error('PCA method information is incomplete.')
  return data
}
async function load() {
  const id = ++requestId
  analysis.value = null; error.value = ''; activeId.value = ''
  if (!props.catalog) { loading.value = false; return }
  loading.value = true
  try { const data = await getProfilePca(); if (id === requestId) analysis.value = data ? validate(data) : null }
  catch (failure) { if (id === requestId) error.value = failure.message }
  finally { if (id === requestId) loading.value = false }
}
watch(() => props.catalog, load, { immediate: true })
const points = computed(() => {
  const metadata = new Map(props.catalog.profiles.map(profile => [profile.id, profile]))
  return (analysis.value?.points || []).map(point => ({ ...metadata.get(point.id), ...point }))
})
const active = computed(() => points.value.find(point => point.id === activeId.value))
function download() {
  exportTsv('profile_PI_PCA.tsv', [
    { key: 'id', label: 'profile_id' }, sampleNameColumn, { key: 'alias', label: 'source_alias' },
    { key: 'assay', label: 'assay' }, { key: 'cell', label: 'raw_cell_label' },
    { key: 'gse', label: 'study' }, { key: 'pc1', label: 'PC1' }, { key: 'pc2', label: 'PC2' },
  ], points.value)
}
</script>
<template>
  <section class="gene-pca" aria-labelledby="gene-pca-title">
    <div class="gene-pca-heading"><div><h2 id="gene-pca-title">Profile PCA</h2><p>Explore sample clustering across the all-gene PI matrix.</p></div><button v-if="analysis" class="button secondary" type="button" @click="download"><Download :size="16"/>Download PCA</button></div>
    <div v-if="error" class="error-state" role="alert"><strong>Unable to load profile PCA</strong><p>{{ error }}</p><button class="button secondary" type="button" @click="load">Retry PCA</button></div>
    <template v-else-if="analysis">
      <p class="pca-context">{{ points.length }} profiles · {{ analysis.feature_count.toLocaleString() }} gene / transcript features · Same PCA coordinates, two groupings.</p>
      <div class="gene-pca-grid"><PcaScatter :points="points" group-by="assay" :variance="analysis.explained_variance_ratio" :active-id="activeId" @inspect="activeId = $event"/><PcaScatter :points="points" group-by="cell" :variance="analysis.explained_variance_ratio" :active-id="activeId" @inspect="activeId = $event"/></div>
      <div class="pca-inspector"><template v-if="active"><strong>{{ sampleName(active) }}</strong><div><span>{{ active.assay }}</span><span>{{ active.cell }}</span><span>{{ active.gse || 'Study unavailable' }}</span><RouterLink :to="`/samples/${active.id}`">View profile ↗</RouterLink></div></template><template v-else><strong>Explore individual profiles</strong><p>Hover or focus a point to compare it in both plots. Select a legend group to highlight it.</p></template></div>
      <details class="pca-methods"><summary>PCA methods and PI source</summary><dl><dt>Input</dt><dd>{{ analysis.preprocessing.input }}</dd><dt>Feature coverage</dt><dd>{{ analysis.preprocessing.coverage }}</dd><dt>Transformation</dt><dd>{{ analysis.preprocessing.transformation }}</dd><dt>Feature scaling</dt><dd>{{ analysis.preprocessing.scaling }}</dd><dt>Missing data</dt><dd>{{ analysis.preprocessing.missing_values }}</dd><dt>Batch adjustment</dt><dd>{{ analysis.preprocessing.batch_correction }}</dd></dl><p>Colors show existing labels, not inferred clusters. Assay, cell line, study and antibody effects may overlap.</p></details>
    </template>
    <div v-else class="gene-pca-grid" :aria-busy="loading"><section v-for="group in ['sequencing assay', 'cell line']" :key="group" class="panel pca-pending-card"><div class="panel-head"><div><h3>PCA by {{ group }}</h3><p>All-gene PI matrix</p></div></div><div v-if="loading" class="pca-skeleton" role="status"><span/>Loading PI profile PCA…</div><div v-else class="pca-pending"><ScatterChart :size="36" aria-hidden="true"/><h4>PI PCA is not available yet</h4><p>The all-gene PI matrix PCA has not been connected yet.</p><button type="button" class="button secondary" @click="load">Reload PCA</button></div></section></div>
  </section>
</template>
<style scoped>
.gene-pca{margin:28px 0 32px;min-width:0}.gene-pca-heading{display:flex;align-items:center;justify-content:space-between;gap:16px;margin-bottom:14px}.gene-pca-heading h2{font-size:24px;font-weight:600;line-height:1.35;color:#204f78;text-wrap:balance;margin:0}.gene-pca-heading p{font-size:13px;color:#64748b;line-height:1.7;margin:7px 0 0;text-wrap:pretty}.pca-context{font-size:13px;color:#64748b;margin:0 0 16px;font-variant-numeric:tabular-nums}.gene-pca-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:22px;min-width:0}
.pca-inspector{padding:18px 22px;margin-top:16px;min-height:100px;border:1px solid #e4eaf1;border-radius:10px;background:#f8fafc;overflow-wrap:anywhere}.pca-inspector strong{font-size:14px;font-weight:500;color:#334155}.pca-inspector p,.pca-inspector>div{font-size:13px;color:#64748b;margin:8px 0 0;line-height:1.65}.pca-inspector>div{display:flex;flex-wrap:wrap;gap:8px 18px;font-variant-numeric:tabular-nums}.pca-inspector a{color:#3e77ab}.pca-methods{font-size:13px;color:#64748b;margin-top:16px;line-height:1.7}.pca-methods summary{cursor:pointer;color:#3e77ab;padding:8px 0}.pca-methods dl{display:grid;grid-template-columns:140px minmax(0,1fr);gap:8px 14px;margin:12px 0;overflow-wrap:anywhere}.pca-methods dt{font-weight:500;color:#334155}.pca-methods dd{margin:0}.pca-methods p{text-wrap:pretty}
.pca-pending-card{min-width:0;overflow:hidden;margin:0}.pca-pending-card .panel-head{padding:22px 24px}.pca-pending-card h3{font-size:19px;font-weight:600;color:#204f78;margin:0;text-wrap:balance}.pca-pending-card .panel-head p{font-size:13px;margin-top:6px}.pca-pending{min-height:290px;display:flex;flex-direction:column;align-items:center;justify-content:center;text-align:center;padding:28px;color:#7b9bb9}.pca-pending h4{font-size:16px;font-weight:500;color:#334155;margin:16px 0 8px}.pca-pending p{font-size:13px;line-height:1.8;max-width:340px;margin:0 0 16px;color:#64748b;text-wrap:pretty}.pca-pending a{font-size:13px;line-height:1.5}.pca-skeleton{min-height:290px;padding:28px;display:flex;flex-direction:column;align-items:center;justify-content:center;gap:18px;color:#64748b;font-size:13px}.pca-skeleton span{display:block;width:100%;height:168px;border:1px solid #e4eaf1;border-radius:8px;background:#f8fafc}
@media(max-width:860px){.gene-pca-grid{grid-template-columns:1fr}.gene-pca-heading{align-items:flex-start;flex-wrap:wrap}.gene-pca-heading h2{font-size:22px}.pca-methods dl{grid-template-columns:1fr;gap:4px}.pca-methods dd{margin-bottom:8px}.pca-pending{min-height:250px}}
</style>
