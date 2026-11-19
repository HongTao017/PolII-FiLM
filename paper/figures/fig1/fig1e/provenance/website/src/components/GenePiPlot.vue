<script setup>
import { sampleName } from '../utils/sampleName.js'
import { computed } from 'vue'
import { assayColors } from '../utils/data.js'
import { groupPi } from '../utils/gene-pi.js'
const props = defineProps({ rows: Array, groupBy: String, range: Array, activeId: String })
const emit = defineEmits(['inspect'])
const groups = computed(() => {
  const result = groupPi(props.rows, props.groupBy)
  if (props.groupBy === 'cell') {
    const byMedian = (a, b) => b.median - a.median || a.label.localeCompare(b.label)
    const byCount = [...result].sort((a, b) => b.count - a.count || a.label.localeCompare(b.label))
    // Match the preview's eight most represented cells; retain every other group below.
    return [...byCount.slice(0, 8).sort(byMedian), ...byCount.slice(8).sort(byMedian)]
  }
  if (props.groupBy === 'assay') {
    const order = ['GRO-seq', 'PRO-seq', 'mNET-seq']
    const rank = label => { const index = order.indexOf(label); return index < 0 ? order.length : index }
    return result.sort((a, b) => rank(a.label) - rank(b.label) || a.label.localeCompare(b.label))
  }
  return result
})
const x = value => 3 + (value - props.range[0]) / (props.range[1] - props.range[0]) * 94
function jitter(id) { return [...id].reduce((hash, letter) => (hash * 31 + letter.charCodeAt(0)) >>> 0, 0) % 13 - 6 }
const labels = {
  assay: { title: 'PI by sequencing assay', description: 'mNET-seq profiles are combined across antibody subtypes.', column: 'Sequencing assay', count: 'assay groups' },
  cell: { title: 'PI by cell line', description: 'Original cell labels · points colored by assay.', column: 'Cell line', count: 'cell labels' },
  author: { title: 'PI by author / research source', description: 'Source author labels · points colored by sequencing assay.', column: 'Author / source', count: 'author labels' },
}
const label = computed(() => labels[props.groupBy])
const title = computed(() => label.value.title)
const active = computed(() => props.rows.find(row => row.id === props.activeId))
const ticks = computed(() => Array.from({ length: 5 }, (_, index) => props.range[0] + index * (props.range[1] - props.range[0]) / 4))
</script>
<template>
  <section class="panel gene-pi-card" :aria-label="title" :data-group-by="groupBy">
    <div class="panel-head"><div><h3>{{ title }}</h3><p>{{ label.description }}</p></div></div>
    <div class="gene-pi-columns" aria-hidden="true"><span>{{ label.column }}</span><span>Observed log₂ PI</span><span>n</span><span>Median</span></div>
    <div class="gene-pi-rows" :class="{ 'gene-pi-cells': groupBy !== 'assay' }" :aria-label="`${title} groups`" tabindex="0">
      <div v-for="group in groups" :key="group.label" class="gene-pi-row" :data-group="group.label" :data-count="group.count" :data-median="group.median">
        <span class="gene-pi-label" :title="group.label">{{ group.label }}</span>
        <div class="gene-pi-track">
          <span v-for="tick in ticks" :key="tick" class="gene-pi-grid" :style="{ left: `${x(tick)}%` }" aria-hidden="true"/>
          <span v-if="range[0] <= 0 && range[1] >= 0" class="gene-pi-zero" :style="{ left: `${x(0)}%` }" aria-hidden="true"/><svg v-if="group.count >= 5" viewBox="0 0 100 50" preserveAspectRatio="none" class="gene-pi-box" aria-hidden="true"><line :x1="x(group.low)" :x2="x(group.high)" y1="12" y2="12"/><line :x1="x(group.low)" :x2="x(group.low)" y1="7" y2="17"/><line :x1="x(group.high)" :x2="x(group.high)" y1="7" y2="17"/><rect :x="x(group.q1)" y="6" :width="Math.max(.15, x(group.q3)-x(group.q1))" height="12"/><line :x1="x(group.median)" :x2="x(group.median)" y1="6" y2="18" class="gene-pi-median"/></svg>
          <a v-for="row in group.profiles" :key="row.id" :href="`#/samples/${row.id}`" class="gene-pi-point" :class="{ 'is-active': activeId === row.id, 'is-dimmed': activeId && activeId !== row.id }" :style="{ left: `${x(row.log2_pi)}%`, transform: `translate(-50%, ${jitter(row.id)}px)`, background: assayColors[row.assay] }" :aria-label="`${sampleName(row)}, ${row.assay}, ${row.cell}, ${row.author || 'Author not recorded'}, ${row.gse || 'Study not recorded'}, log2 PI ${row.log2_pi.toFixed(4)}. View profile.`" :data-profile="row.id" :data-value="row.log2_pi" @mouseenter="emit('inspect', row.id)" @mouseleave="emit('inspect', '')" @focus="emit('inspect', row.id)" @blur="emit('inspect', '')"/>
        </div>
        <span class="gene-pi-count">{{ group.count }}</span><span class="gene-pi-value" :title="`Median log₂ PI: ${group.median}`">{{ group.median.toFixed(2) }}</span>
      </div>
      <div class="gene-pi-row gene-pi-axis-row" aria-hidden="true"><span/><div class="gene-pi-axis"><span v-for="tick in ticks" :key="tick" :style="{ left: `${x(tick)}%` }">{{ tick.toFixed(1) }}</span></div><span/><span/></div>
    </div>
    <p class="gene-pi-axis-title">log₂(Pause Index)</p>
    <p class="gene-pi-card-foot">{{ groups.length }} {{ label.count }} · {{ rows.length }} measured profiles<span v-if="active"> · Highlight: {{ active.cell }}</span></p>
  </section>
</template>
<style scoped>
.gene-pi-card{min-width:0;overflow:hidden;margin:0;display:flex;flex-direction:column}.gene-pi-card .panel-head{padding:22px 24px;min-height:106px}.gene-pi-card h3{font-size:19px;font-weight:600;line-height:1.4;color:#204f78;margin:0;text-wrap:balance}.gene-pi-card .panel-head p{font-size:13px;line-height:1.6;margin-top:6px;text-wrap:pretty}
.gene-pi-columns,.gene-pi-row{display:grid;grid-template-columns:106px minmax(0,1fr) 26px 52px;gap:8px;align-items:center}.gene-pi-columns{padding:16px 20px 10px;color:#64748b;font-size:12px}.gene-pi-columns>span:nth-last-child(-n+2){text-align:right}.gene-pi-rows{padding:0 20px;max-height:480px;min-height:260px;overflow-y:auto;scrollbar-width:thin;scrollbar-gutter:stable}.gene-pi-cells{min-height:260px}.gene-pi-row{min-height:56px}.gene-pi-row[data-group]:nth-child(odd){background:#f6f9fd}.gene-pi-row[data-group]{padding:2px 0}.gene-pi-label{font-size:12px;line-height:1.5;color:#475569;overflow-wrap:anywhere}.gene-pi-track{position:relative;height:50px;min-width:0}.gene-pi-zero{position:absolute;top:0;height:100%;border-left:1px dashed #91a9be;pointer-events:none}.gene-pi-value{font-size:12px;text-align:right;color:#64748b;font-variant-numeric:tabular-nums}.gene-pi-grid{position:absolute;height:100%;border-left:1px solid #edf1f6;pointer-events:none}.gene-pi-box{position:absolute;inset:0;width:100%;height:100%;pointer-events:none;overflow:visible}.gene-pi-box line{stroke:#829ab0;stroke-width:1;vector-effect:non-scaling-stroke}.gene-pi-box rect{fill:#e8eff7;fill-opacity:.8;stroke:#9eb2c7;stroke-width:1;vector-effect:non-scaling-stroke}.gene-pi-box .gene-pi-median{stroke:#345776;stroke-width:2}.gene-pi-count{font-size:13px;color:#64748b;text-align:right;font-variant-numeric:tabular-nums}
.gene-pi-point{position:absolute;top:31px;width:9px;height:9px;border:1px solid #fff;border-radius:50%;opacity:.85;transition:opacity 150ms ease-out;z-index:1}.gene-pi-point.is-active{outline:2px solid #172e47;outline-offset:2px;opacity:1;z-index:2}.gene-pi-point.is-dimmed{opacity:.18}.gene-pi-point:focus-visible{outline:2px solid #172e47;outline-offset:3px;z-index:2}.gene-pi-axis-row{position:sticky;bottom:0;min-height:30px;background:#fff;z-index:3}.gene-pi-axis{position:relative;height:26px;border-top:1px solid #b8c9d8}.gene-pi-axis>span{position:absolute;top:5px;transform:translateX(-50%);font-size:12px;color:#64748b;font-variant-numeric:tabular-nums}.gene-pi-axis-title{text-align:center;font-size:14px;color:#334155;margin:8px 20px 16px}.gene-pi-card-foot{font-size:12px;color:#64748b;padding:13px 20px;background:#f8fafc;border-top:1px solid #e4eaf1;margin:auto 0 0;text-wrap:pretty}
@media(max-width:760px){.gene-pi-card .panel-head{padding:20px;min-height:0}.gene-pi-card h3{font-size:18px}.gene-pi-columns,.gene-pi-row{grid-template-columns:80px minmax(0,1fr) 22px 46px;gap:6px}.gene-pi-columns{padding-inline:14px}.gene-pi-rows{padding-inline:14px}.gene-pi-label{font-size:12px}.gene-pi-value{font-size:11px}.gene-pi-columns{font-size:11px}.gene-pi-axis>span{font-size:11px}.gene-pi-point{width:8px;height:8px}.gene-pi-card-foot{padding-inline:14px}}
@media(prefers-reduced-motion:reduce){.gene-pi-point{transition:none}}
.gene-pi-card>.panel-head{background:#204f78;border-bottom-color:#204f78}.gene-pi-card>.panel-head h3{color:#fff}.gene-pi-card>.panel-head p{color:#dce8f4}
</style>
