<script setup>
import { sampleName } from '../utils/sampleName.js'
import { computed, ref, watch } from 'vue'
import { assays, assayColors } from '../utils/data.js'
import GenePiPlot from './GenePiPlot.vue'
import EmptyState from './EmptyState.vue'
const props = defineProps({ rows: Array, allRows: Array, symbol: String, transcript: String, profileCount: Number })
defineEmits(['reset'])
const activeId = ref('')
watch(() => [props.symbol, props.transcript, props.rows], () => { activeId.value = '' })
const measured = computed(() => props.rows.filter(row => Number.isFinite(row.log2_pi)))
const range = computed(() => {
  const values = props.allRows.map(row => row.log2_pi).filter(Number.isFinite)
  if (!values.length) return [-1, 1]
  const low = Math.min(...values), high = Math.max(...values), padding = Math.max((high - low) * .08, .25)
  return [low - padding, high + padding]
})
const active = computed(() => measured.value.find(row => row.id === activeId.value))
</script>
<template>
  <section id="gene-results" class="gene-pi-groups" aria-labelledby="gene-pi-title">
    <div class="gene-pi-heading"><h2 id="gene-pi-title">{{ symbol }} · PI across assays, cell lines and authors</h2><p>{{ transcript }} · {{ allRows.length }} of {{ profileCount }} profiles measured · {{ rows.length }} shown with current source filters.</p></div>
    <div v-if="measured.length" class="gene-pi-legend" aria-label="Sequencing assay colors"><span v-for="assay in assays" :key="assay"><i :style="{ background: assayColors[assay] }" aria-hidden="true"/>{{ assay }}</span></div>
    <div class="gene-pi-layout">
      <template v-if="measured.length">
        <div class="panel gene-pi-combined" role="group" aria-label="PI by sequencing assay and author">
          <GenePiPlot :rows="measured" :range="range" group-by="assay" :active-id="activeId" @inspect="activeId = $event"/>
          <GenePiPlot class="gene-pi-authors" :rows="measured" :range="range" group-by="author" :active-id="activeId" @inspect="activeId = $event"/>
        </div>
        <div class="gene-pi-cell-column"><GenePiPlot class="gene-pi-cell-panel" :rows="measured" :range="range" group-by="cell" :active-id="activeId" @inspect="activeId = $event"/></div>
      </template>
      <EmptyState v-else title="No measured PI for these filters" description="Choose another transcript or reset the source filters."><button type="button" class="button secondary" @click="$emit('reset')">Reset filters</button></EmptyState>
    </div>
    <template v-if="measured.length">
      <div class="gene-pi-inspector"><template v-if="active"><strong>{{ sampleName(active) }}</strong><div><span>{{ active.assay }}</span><span>{{ active.cell }}</span><span>{{ active.author || 'Author not recorded' }}</span><span>{{ active.antibody || 'No antibody subtype' }}</span><span>{{ active.gse || 'Study unavailable' }}</span><span>log₂ PI: {{ active.log2_pi.toFixed(4) }}</span><RouterLink :to="`/samples/${active.id}`">View profile ↗</RouterLink></div></template><template v-else><strong>Compare this gene across profiles</strong><p>Hover or focus a point to highlight the same profile in all three views. Scroll within the cell-line and author plots to explore all groups.</p></template></div>
      <details class="gene-pi-methods"><summary>How to read these distributions</summary><p>Each point is an original measured log₂ PI for {{ symbol }} / {{ transcript }}. Boxes show the middle 50% of values, the center line is the median, and whiskers extend to the most extreme observations within 1.5 × IQR. Single-profile groups show a point without a distribution estimate. All three charts share the same numeric scale. Vertical point offsets only separate overlapping profiles.</p><p>Missing values are excluded, never filled or replaced with zero. These views show one-dimensional PI distributions; no PCA or cluster assignments are calculated for an individual gene. Author groups use the original source labels and may combine studies, assays or cell lines. Laboratory affiliations have not been verified. Assay, cell line, study and antibody effects may overlap.</p></details>
    </template>
  </section>
</template>
<style scoped>
.gene-pi-groups{margin:24px 0 32px;min-width:0;scroll-margin-top:110px}.gene-pi-heading h2{font-size:23px;font-weight:600;color:#204f78;line-height:1.4;margin:0;text-wrap:balance}.gene-pi-heading p{font-size:13px;color:#64748b;margin:8px 0 14px;line-height:1.7;font-variant-numeric:tabular-nums;text-wrap:pretty}.gene-pi-legend{display:flex;flex-wrap:wrap;gap:10px 22px;margin:0 0 16px;color:#475569;font-size:13px}.gene-pi-legend>span{display:inline-flex;align-items:center;gap:7px}.gene-pi-legend i{display:inline-block;width:10px;height:10px;border-radius:50%}.gene-pi-layout{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:22px;min-width:0}.gene-pi-inspector{border:1px solid #e4eaf1;border-radius:10px;background:#f8fafc;min-height:104px;padding:18px 22px;margin-top:16px;overflow-wrap:anywhere}.gene-pi-inspector strong{font-size:14px;font-weight:500;color:#334155}.gene-pi-inspector>div{display:flex;flex-wrap:wrap;gap:8px 16px;margin-top:8px;color:#64748b;font-size:13px;font-variant-numeric:tabular-nums}.gene-pi-inspector a{color:#3e77ab}.gene-pi-inspector p{font-size:13px;color:#64748b;line-height:1.7;margin:8px 0 0;text-wrap:pretty}.gene-pi-methods{font-size:13px;color:#64748b;line-height:1.8;margin-top:12px}.gene-pi-methods summary{cursor:pointer;color:#3e77ab;padding:8px 0}.gene-pi-methods p{text-wrap:pretty}
.gene-pi-layout :deep(.panel-head){padding:18px 20px;min-height:0}
.gene-pi-combined{min-width:0;display:flex;flex-direction:column}
.gene-pi-combined :deep(.gene-pi-card + .gene-pi-card){border-top:1px solid #e4eaf1!important}.gene-pi-combined :deep(.gene-pi-rows){min-height:0}.gene-pi-combined :deep(.gene-pi-authors .gene-pi-rows){max-height:280px}
.gene-pi-cell-column{position:relative;min-width:0;min-height:0}.gene-pi-layout :deep(.gene-pi-cell-panel){position:absolute;inset:0;height:100%;min-height:0}.gene-pi-layout :deep(.gene-pi-cell-panel .gene-pi-rows){flex:1;height:0;min-height:260px;max-height:none}
@media(max-width:860px){.gene-pi-cell-column{position:static}.gene-pi-layout :deep(.gene-pi-cell-panel){position:static;height:auto}.gene-pi-layout :deep(.gene-pi-cell-panel .gene-pi-rows){flex:none;height:auto;max-height:480px}}
@media(max-width:860px){.gene-pi-layout{grid-template-columns:1fr}.gene-pi-heading h2{font-size:21px}.gene-pi-inspector{padding:16px}}
</style>
