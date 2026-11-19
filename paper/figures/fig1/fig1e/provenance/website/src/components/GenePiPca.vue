<script setup>
import { sampleName, sampleNameColumn } from '../utils/sampleName.js'
import { computed, ref, watch } from 'vue'
import { Download } from '@lucide/vue'
import { singleGenePiPca } from '../utils/gene-pca.js'
import { exportTsv } from '../utils/data.js'
import PcaScatter from './PcaScatter.vue'
const props = defineProps({ symbol: String, transcript: String, rows: Array, allRows: Array })
const activeId = ref('')
watch(() => [props.symbol, props.transcript, props.rows], () => { activeId.value = '' })
const analysis = computed(() => singleGenePiPca(props.allRows))
const points = computed(() => analysis.value.points || [])
const visibleIds = computed(() => props.rows.map(row => row.id))
const active = computed(() => points.value.find(point => point.id === activeId.value))
function download() {
  exportTsv(`${props.symbol}_${props.transcript}_PI_PCA.tsv`, [
    { key: 'id', label: 'profile_id' }, sampleNameColumn, { key: 'alias', label: 'source_alias' }, { key: 'assay', label: 'assay' },
    { key: 'cell', label: 'raw_cell_label' }, { key: 'antibody', label: 'antibody' }, { key: 'gse', label: 'study' },
    { key: 'log2_pi', label: 'log2_PI' }, { key: 'pc1', label: 'PC1' },
  ], points.value.filter(point => visibleIds.value.includes(point.id)))
}
</script>
<template>
  <section class="gene-pi-pca" aria-labelledby="gene-pi-pca-title">
    <div class="pi-pca-heading"><div><h2 id="gene-pi-pca-title">{{ symbol }} · PI PCA</h2><p>{{ transcript }} · PCA of this gene's observed PI values across profiles.</p></div><button v-if="analysis.status === 'complete'" class="button secondary" type="button" :disabled="!rows.length" @click="download"><Download :size="16"/>Download gene PCA</button></div>
    <div v-if="analysis.status === 'unavailable'" class="pi-pca-pending"><strong>PCA is unavailable for this transcript</strong><p>{{ analysis.reason }}</p></div>
    <template v-else>
      <p class="pi-pca-context">{{ points.length }} measured profiles · 1 PI feature · 1 principal component · {{ rows.length }} shown</p>
      <p class="pi-pca-input-note">This gene has one PI feature: PC1 explains 100% of its variance. Points lie on a display baseline because there is no second principal component. Colors identify existing assay or cell labels.</p>
      <p v-if="!rows.length" class="pi-pca-no-matches">No measured profiles match these filters. Reset the filters to show points.</p>
      <div class="pi-pca-grid"><PcaScatter :points="points" :visible-ids="visibleIds" :variance="[1]" one-dimensional group-by="assay" :active-id="activeId" @inspect="activeId = $event"/><PcaScatter :points="points" :visible-ids="visibleIds" :variance="[1]" one-dimensional group-by="cell" :active-id="activeId" @inspect="activeId = $event"/></div>
      <div class="pi-pca-inspector"><template v-if="active"><strong>{{ sampleName(active) }}</strong><div><span>{{ active.assay }}</span><span>{{ active.cell }}</span><span>{{ active.antibody || 'No antibody subtype' }}</span><span>log₂ PI: {{ active.log2_pi.toFixed(4) }}</span><span>PC1: {{ active.pc1.toFixed(4) }}</span><RouterLink :to="`/samples/${active.id}`">View profile ↗</RouterLink></div></template><template v-else><strong>Compare PI profiles</strong><p>Hover or focus a sample to highlight it in both views. Select a legend group to highlight an assay or cell line.</p></template></div>
      <details class="pi-pca-methods"><summary>Gene PCA methods</summary><dl><dt>Input</dt><dd>Only {{ symbol }} / {{ transcript }}: one original log₂ PI value per measured profile.</dd><dt>Missing values</dt><dd>Excluded, never imputed or replaced with zero. Original PI values are preserved.</dd><dt>Scaling and PCA</dt><dd>Population Z-score across all measured profiles. With one feature, the covariance matrix is 1 × 1; the positive loading is +1 and PC1 equals the standardized PI. Mean: {{ analysis.mean.toFixed(6) }}; population SD: {{ analysis.std.toFixed(6) }}.</dd><dt>Display</dt><dd>The vertical position is zero for every point, solely a display baseline. There is no PC2 or artificial vertical jitter. The export contains the original PI and PC1 only.</dd><dt>Filters</dt><dd>Both views use the same scores and scale. Filters change visible profiles without recalculating PCA.</dd><dt>Interpretation</dt><dd>PCA does not assign clusters. This one-feature view shows PI similarity; assay, cell line, study and antibody effects can overlap. No batch correction is applied. The all-gene matrix PCA remains on the gene search page.</dd></dl></details>
    </template>
  </section>
</template>
<style scoped>
.gene-pi-pca{margin:26px 0 36px;min-width:0}.pi-pca-heading{display:flex;align-items:center;justify-content:space-between;gap:16px;margin-bottom:14px}.pi-pca-heading h2{font-size:24px;font-weight:600;color:#204f78;line-height:1.4;margin:0;text-wrap:balance}.pi-pca-heading p{font-size:13px;color:#64748b;line-height:1.7;margin:7px 0 0;text-wrap:pretty}.pi-pca-context{font-size:13px;color:#64748b;font-variant-numeric:tabular-nums;margin:0 0 12px}.pi-pca-input-note{font-size:13px;color:#475569;background:#f0f6fc;border-left:3px solid #7398bb;padding:12px 16px;margin:0 0 18px;line-height:1.8;text-wrap:pretty}.pi-pca-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:22px;min-width:0}.pi-pca-pending{padding:30px 24px;border:1px solid #e4eaf1;border-radius:10px;background:#f8fafc;min-height:150px}.pi-pca-pending strong{font-size:16px;font-weight:500;color:#334155}.pi-pca-pending p,.pi-pca-no-matches{font-size:13px;color:#64748b;line-height:1.8;text-wrap:pretty}.pi-pca-inspector{border:1px solid #e4eaf1;border-radius:10px;background:#f8fafc;min-height:104px;padding:18px 22px;margin-top:16px;overflow-wrap:anywhere}.pi-pca-inspector strong{font-size:14px;font-weight:500;color:#334155}.pi-pca-inspector>div{display:flex;flex-wrap:wrap;gap:8px 16px;margin-top:8px;font-size:13px;color:#64748b;font-variant-numeric:tabular-nums}.pi-pca-inspector a{color:#3e77ab}.pi-pca-inspector p{font-size:13px;color:#64748b;line-height:1.7;margin:8px 0 0;text-wrap:pretty}.pi-pca-methods{font-size:13px;color:#64748b;line-height:1.8;margin-top:12px}.pi-pca-methods summary{cursor:pointer;color:#3e77ab;padding:8px 0}.pi-pca-methods dl{display:grid;grid-template-columns:150px minmax(0,1fr);gap:10px 16px}.pi-pca-methods dt{font-weight:500;color:#334155}.pi-pca-methods dd{margin:0;text-wrap:pretty}
@media(max-width:860px){.pi-pca-grid{grid-template-columns:1fr}.pi-pca-heading{align-items:flex-start;flex-wrap:wrap}.pi-pca-heading h2{font-size:22px}.pi-pca-methods dl{grid-template-columns:1fr;gap:4px}.pi-pca-methods dd{margin-bottom:8px}}
</style>
