<script setup>
import { computed, ref } from 'vue'
import { sequenceStats } from '../utils/sequenceStats.js'
import SequenceRatioPie from './SequenceRatioPie.vue'
import GcRegionProfile from './GcRegionProfile.vue'
import { useChartReveal } from '../composables/chartReveal.js'
const props = defineProps({ data: { type: Object, required: true } })
const stats = computed(() => sequenceStats(props.data.sequence))
const { host, visible } = useChartReveal()
const gcPercent = computed(() => stats.value.called ? 100 * (stats.value.counts.G + stats.value.counts.C) / stats.value.called : null)
const activeBase = ref('')
const colors = { A: '#12a89b', T: '#df7387', G: '#dfa32f', C: '#417dc6' }
const pct = value => value === null ? 'Not available' : value.toFixed(2) + '%'

</script>
<template>
  <section class="region-stats" aria-label="Selected sequence region statistics" :data-region="data.region" :data-sequence="data.sequence" :data-strand="data.display_strand">
    <h4>Selected-region sequence statistics</h4>
    <p class="stats-scope">{{ data.length.toLocaleString('en-US') }} displayed bases · {{ data.region === 'gene' ? 'Current 500 bp window (or remaining bases)' : 'Entire selected region' }} · {{ data.display_strand }} strand</p>
    <div :key="data.sequence" ref="host" class="region-pies" :class="{'is-revealed':visible}">
      <div class="base-chart"><h5>ATGC composition</h5><div v-for="(color, base) in colors" :key="base" class="base-stat" :data-base="base" :data-count="stats.counts[base]" :class="{'is-active':activeBase===base}" tabindex="0" @mouseenter="activeBase=base" @mouseleave="activeBase=''" @focus="activeBase=base" @blur="activeBase=''"><span>{{ base }}</span><svg viewBox="0 0 100 10" preserveAspectRatio="none" role="img" :aria-label="`${base}: ${stats.counts[base]} bases, ${pct(stats.called ? 100 * stats.counts[base] / stats.called : null)}`"><rect width="100" height="10" fill="#edf2f8"/><rect :width="stats.called ? 100 * stats.counts[base] / stats.called : 0" height="10" :fill="color"/></svg><strong>{{ pct(stats.called ? 100 * stats.counts[base] / stats.called : null) }}<small>{{ stats.counts[base].toLocaleString('en-US') }} bases</small></strong></div><div v-if="activeBase" class="base-tooltip" role="tooltip"><strong>{{activeBase}}</strong> · {{stats.counts[activeBase].toLocaleString('en-US')}} bases · {{pct(stats.called ? 100*stats.counts[activeBase]/stats.called : null)}}</div><p>{{ stats.called.toLocaleString('en-US') }} called bases · {{ stats.ambiguous.toLocaleString('en-US') }} N / ambiguous excluded</p></div>
      <SequenceRatioPie label="GC content" :numerator="stats.counts.G+stats.counts.C" :denominator="stats.called" unit="bases" :percent="gcPercent" primary="G + C" remainder="A + T" color="#12a89b" :animated="visible"><template #value><span data-stat="gc-percent">{{gcPercent === null ? 'Not available' : gcPercent.toFixed(1) + '%'}}</span></template>{{(stats.counts.G+stats.counts.C).toLocaleString('en-US')}} G/C bases / {{stats.called.toLocaleString('en-US')}} called bases.</SequenceRatioPie>
      <SequenceRatioPie label="CpG dinucleotides" :numerator="stats.cpg" :denominator="stats.validPairs" unit="pairs" :percent="stats.cpgPercent" :precision="2" primary="CG pairs" remainder="Other valid pairs" color="#df9b32" :animated="visible"><template #value><span data-stat="cpg-percent">{{pct(stats.cpgPercent)}}</span></template><strong data-stat="cpg-count">{{stats.cpg.toLocaleString('en-US')}}</strong> CG pairs / {{stats.validPairs.toLocaleString('en-US')}} adjacent A/T/G/C pairs.</SequenceRatioPie>
    </div>
    <p class="pie-methods">GC uses called bases; CpG uses adjacent called pairs, counted with overlap. N and ambiguous bases are excluded from each denominator.</p>
    <GcRegionProfile :data="data" :bins="stats.bins"/>
  </section>
</template>
<style scoped>
.region-pies{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:24px;margin-top:16px}.region-stats .pie-methods{margin:4px 0 16px}.base-chart{min-width:0;position:relative}.base-stat{border-radius:4px}.base-stat.is-active{background:#eef3f9}.base-stat:focus-visible{outline:2px solid #417dc6;outline-offset:2px}.base-stat svg rect:last-child{transform-origin:left center;transition:transform 180ms ease-out}.base-stat.is-active svg rect:last-child{transform:scaleY(1.3)}.base-tooltip{position:absolute;top:30px;left:20px;z-index:4;background:#fff;border:1px solid #d5dfeb;border-radius:5px;box-shadow:0 3px 10px #0f172a20;padding:8px 12px;pointer-events:none;font-size:12px;color:#345776}.is-revealed .base-stat svg{transform-origin:left center;animation:composition-reveal 600ms ease-out both}@keyframes composition-reveal{from{opacity:0;transform:scaleX(0)}to{opacity:1;transform:scaleX(1)}}@media(prefers-reduced-motion:reduce){.is-revealed .base-stat svg{animation:none}.base-stat svg rect:last-child{transition:none}}@media(max-width:860px){.region-pies{grid-template-columns:1fr;gap:14px}}
.region-stats{padding:16px 20px;border-top:1px solid #e4eaf1;min-width:0}.region-stats h4{font-size:16px;color:#23415f;margin:0 0 6px;text-wrap:balance}.region-stats h5{font-size:13px;color:#345776;margin:0 0 12px;font-weight:600}.region-stats p{font-size:11px;color:#64748b;line-height:1.7;margin:6px 0;text-wrap:pretty}.stats-grid{display:grid;grid-template-columns:1fr;gap:24px;margin:16px 0}.base-stat{display:grid;grid-template-columns:18px minmax(0,1fr) 90px;gap:10px;align-items:center;margin:10px 0;font-size:12px;color:#345776}.base-stat svg{width:100%;height:10px}.base-stat strong{text-align:right;font-variant-numeric:tabular-nums;font-weight:500}.base-stat small{display:block;font-size:10px;color:#64748b;font-weight:400}@media(max-width:600px){.stats-grid{grid-template-columns:1fr;gap:16px}}
</style>
