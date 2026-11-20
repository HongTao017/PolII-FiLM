<script setup>
import { sampleName } from '../utils/sampleName.js'
import { computed, ref, watch } from 'vue'
import { assayColors } from '../utils/data.js'

const props = defineProps({
  points: { type: Array, required: true },
  groupBy: { type: String, required: true },
  variance: { type: Array, default: () => [] },
  algorithm: { type: String, default: 'PCA' },
  visibleIds: { type: Array, default: null },
  activeId: { type: String, default: '' },
  oneDimensional: { type: Boolean, default: false },
})
const emit = defineEmits(['inspect'])
const selectedGroup = ref('')
const visiblePoints = computed(() => props.visibleIds === null ? props.points : props.points.filter(point => props.visibleIds.includes(point.id)))
const groups = computed(() => [...new Set(visiblePoints.value.map(point => point[props.groupBy]))].sort())
const colorGroups = computed(() => [...new Set(props.points.map(point => point[props.groupBy]))].sort())
watch(groups, value => { if (!value.includes(selectedGroup.value)) selectedGroup.value = '' })
const colors = ['#3e77ab', '#d69a41', '#218879', '#b45968', '#6877ad', '#91653c', '#348da2', '#657c3e', '#c46e46', '#7d658f', '#456180']
const color = group => props.groupBy === 'assay' ? assayColors[group] : colors[colorGroups.value.indexOf(group) % colors.length]
const shape = group => props.groupBy === 'assay' ? 0 : Math.floor(colorGroups.value.indexOf(group) / colors.length)
function domain(key) {
  const values = props.points.map(point => point[key])
  const low = Math.min(...values), high = Math.max(...values)
  const padding = (high - low || 1) * .09
  return [low - padding, high + padding]
}
const xDomain = computed(() => domain('pc1'))
const yDomain = computed(() => domain('pc2'))
const x = value => 66 + (value - xDomain.value[0]) / (xDomain.value[1] - xDomain.value[0]) * 372
const y = value => 288 - (value - yDomain.value[0]) / (yDomain.value[1] - yDomain.value[0]) * 262
const ticks = range => Array.from({ length: 5 }, (_, index) => range[0] + (range[1] - range[0]) * index / 4)
const tickLabel = value => Math.abs(value) >= 100 ? value.toFixed(0) : value.toFixed(1)
const title = computed(() => `${props.algorithm} by ${props.groupBy === 'assay' ? 'sequencing assay' : 'cell line'}`)
const xLabel = computed(() => props.algorithm === 't-SNE' ? 't-SNE 1' : `PC1 (${(props.variance[0] * 100).toFixed(1)}%)`)
const yLabel = computed(() => props.oneDimensional ? 'Display baseline (no PC2)' : props.algorithm === 't-SNE' ? 't-SNE 2' : `PC2 (${(props.variance[1] * 100).toFixed(1)}%)`)
const description = computed(() => props.groupBy === 'assay' ? 'Colored by GRO-seq, PRO-seq and mNET-seq.' : 'Colored by the original cell labels in the catalog.')
const dimmed = point => Boolean((selectedGroup.value && point[props.groupBy] !== selectedGroup.value) || (props.activeId && point.id !== props.activeId))
const activePoint = computed(() => visiblePoints.value.find(point => point.id === props.activeId))
const pointLabel = point => `${sampleName(point)}, ${point.assay}, ${point.cell}, ${xLabel.value} ${point.pc1.toFixed(2)}${props.oneDimensional ? '' : `, ${yLabel.value} ${point.pc2.toFixed(2)}`}${Number.isFinite(point.log2_pi) ? `, log2 PI ${point.log2_pi.toFixed(4)}` : ''}. View profile.`
</script>

<template>
  <section class="panel embedding-card" :class="algorithm === 't-SNE' ? 'tsne-card' : 'pca-card'" :aria-label="title">
    <div class="panel-head"><div><h3>{{ title }}</h3><p>{{ description }}</p></div></div>
    <div class="pca-plot-wrap">
      <svg class="pca-scatter" viewBox="0 0 464 350" role="group" :aria-label="`${title}. Each point is one profile. Both plots use the same coordinates.`">
        <g aria-hidden="true">
          <g v-for="tick in ticks(xDomain)" :key="`x-${tick}`">
            <line :x1="x(tick)" :x2="x(tick)" y1="26" y2="288" class="pca-gridline" />
            <text :x="x(tick)" y="311" text-anchor="middle">{{ tickLabel(tick) }}</text>
          </g>
          <g v-for="tick in oneDimensional ? [0] : ticks(yDomain)" :key="`y-${tick}`">
            <line x1="66" x2="438" :y1="y(tick)" :y2="y(tick)" class="pca-gridline" />
            <text x="56" :y="y(tick) + 5" text-anchor="end">{{ tickLabel(tick) }}</text>
          </g>
          <line v-if="xDomain[0] < 0 && xDomain[1] > 0" :x1="x(0)" :x2="x(0)" y1="26" y2="288" class="pca-zero" />
          <line v-if="yDomain[0] < 0 && yDomain[1] > 0" x1="66" x2="438" :y1="y(0)" :y2="y(0)" class="pca-zero" />
          <path d="M66 26 V288 H438" class="pca-axis" />
          <text x="252" y="342" text-anchor="middle" class="pca-axis-title">{{ xLabel }}</text>
          <text x="18" y="157" transform="rotate(-90 18 157)" text-anchor="middle" class="pca-axis-title">{{ yLabel }}</text>
        </g>
        <a v-for="point in visiblePoints" :key="point.id" :href="`#/samples/${point.id}`" class="embedding-point" :class="{ 'pca-point': algorithm === 'PCA', 'tsne-point': algorithm === 't-SNE', 'is-dimmed': dimmed(point), 'is-active': point.id === activeId }" :aria-label="pointLabel(point)" :data-profile="point.id" :data-value="point.log2_pi" :data-x="point.pc1" :data-y="point.pc2" :data-group="point[groupBy]" @mouseenter="emit('inspect', point.id)" @mouseleave="emit('inspect', '')" @focus="emit('inspect', point.id)" @blur="emit('inspect', '')">
          <title>{{ sampleName(point) }} · {{ point.assay }} · {{ point.cell }}</title>
          <circle v-if="shape(point[groupBy]) === 0" :cx="x(point.pc1)" :cy="y(point.pc2)" r="5.5" :fill="color(point[groupBy])" />
          <rect v-else-if="shape(point[groupBy]) === 1" :x="x(point.pc1)-5" :y="y(point.pc2)-5" width="10" height="10" rx="1" :fill="color(point[groupBy])" />
          <path v-else :d="`M${x(point.pc1)} ${y(point.pc2)-6.5} l6 11 h-12 Z`" :fill="color(point[groupBy])" />
        </a>
        <circle v-if="activePoint" :cx="x(activePoint.pc1)" :cy="y(activePoint.pc2)" r="9" fill="none" stroke="#172e47" stroke-width="1.5" pointer-events="none" aria-hidden="true"/>
      </svg>
    </div>
    <div class="pca-legend-heading"><span>{{ groupBy === 'assay' ? 'Sequencing assay' : 'Cell line' }} · {{ groups.length }} groups</span><button v-if="selectedGroup" type="button" @click="selectedGroup = ''">Reset highlight</button></div>
    <div class="pca-legend" :class="{ 'pca-legend-cells': groupBy === 'cell' }" :aria-label="`${title} legend`">
      <button v-for="group in groups" :key="group" type="button" :aria-pressed="selectedGroup === group" :class="{ 'is-muted': selectedGroup && selectedGroup !== group }" @click="selectedGroup = selectedGroup === group ? '' : group">
        <svg viewBox="0 0 16 16" aria-hidden="true"><circle v-if="shape(group) === 0" cx="8" cy="8" r="5" :fill="color(group)"/><rect v-else-if="shape(group) === 1" x="3" y="3" width="10" height="10" rx="1" :fill="color(group)"/><path v-else d="M8 2 L14 13 H2 Z" :fill="color(group)"/></svg>
        <span>{{ group }}</span><small>{{ visiblePoints.filter(point => point[groupBy] === group).length }}</small>
      </button>
    </div>
  </section>
</template>

<style scoped>
.embedding-card{min-width:0;margin:0;overflow:hidden;display:flex;flex-direction:column}
.embedding-card .panel-head{padding:22px 24px;min-height:104px;align-items:flex-start}
.embedding-card h3{font-size:19px;line-height:1.4;font-weight:600;color:#204f78;text-wrap:balance;margin:0}
.embedding-card .panel-head p{font-size:13px;line-height:1.6;margin-top:6px;text-wrap:pretty}
.pca-plot-wrap{padding:10px 16px 0;min-width:0}
.pca-scatter{display:block;width:100%;height:auto;max-height:420px;font-family:inherit;overflow:visible}
.pca-scatter text{fill:#64748b;font-size:14px;font-variant-numeric:tabular-nums}
.pca-scatter .pca-axis-title{fill:#334155;font-size:15px;font-weight:500}
.pca-gridline{stroke:#edf1f6;stroke-width:1}
.pca-zero{stroke:#b8c9d8;stroke-dasharray:4 4}
.pca-axis{fill:none;stroke:#a9bacb;stroke-width:1}
.embedding-point{opacity:.85;cursor:pointer;transition:opacity 150ms ease-out;outline:none}
.embedding-point>*{stroke:#fff;stroke-width:.9}
.embedding-point.is-dimmed{opacity:.15}
.embedding-point.is-active,.embedding-point:focus-visible{opacity:1}
.embedding-point.is-active>*,.embedding-point:focus-visible>*{stroke:#172e47;stroke-width:2.5}
.pca-legend-heading{display:flex;justify-content:space-between;align-items:center;gap:8px;padding:12px 24px 6px;min-height:40px;font-size:12px;color:#64748b}
.pca-legend-heading button{font-size:12px;color:#3e77ab;background:transparent;padding:4px 0;white-space:nowrap}
.pca-legend{display:flex;align-content:flex-start;flex-wrap:wrap;gap:6px;padding:0 24px 20px;min-height:164px}
.pca-legend-cells{max-height:164px;overflow-y:auto;scrollbar-width:thin}
.pca-legend button{display:inline-flex;align-items:center;gap:6px;min-height:32px;max-width:100%;padding:5px 8px;border:1px solid #e4eaf1;border-radius:6px;background:#fff;color:#334155;font-size:12px;line-height:1.4}
.pca-legend button[aria-pressed=true]{border-color:#3e77ab;background:#f0f6fc}
.pca-legend button.is-muted{opacity:.45}
.pca-legend svg{width:16px;height:16px;flex-shrink:0}
.pca-legend small{font-size:11px;color:#64748b;font-variant-numeric:tabular-nums}
@media(max-width:760px){.embedding-card .panel-head{padding:20px;min-height:0}.embedding-card h3{font-size:18px}.pca-plot-wrap{padding:8px 8px 0}.pca-legend-heading{padding-inline:20px}.pca-legend{padding-inline:20px;min-height:0;max-height:164px}.pca-legend button{min-height:36px}}
@media(prefers-reduced-motion:reduce){.embedding-point{transition:none}}
</style>
