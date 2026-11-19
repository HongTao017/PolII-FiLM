<script setup>
import { computed, ref, watch } from 'vue'
import { useChartReveal } from '../composables/chartReveal.js'
const props = defineProps({ data: Object, bins: Array })
const { host, visible } = useChartReveal()
const activeIndex = ref(null)
watch(() => props.data, () => { activeIndex.value = null })
const position = index => {
  const v = props.data, coordinate = v.display_strand === '-' ? v.genomic_end - index : v.genomic_start + index
  if (!v.endpoint) return coordinate.toLocaleString('en-US')
  const offset = (coordinate - v.endpoint.anchor) * (v.gene_strand === '-' ? -1 : 1)
  return offset > 0 ? '+' + offset : String(offset)
}
const points = computed(() => props.bins.map(bin => ({ ...bin, x: 44 + 670 * (bin.start + bin.end + 1) / (2 * props.data.sequence.length), y: bin.gc === null ? null : 146 - bin.gc * 1.2 })))
const active = computed(() => activeIndex.value === null ? null : points.value[activeIndex.value])
const label = bin => `${props.data.endpoint?.anchor_name || props.data.chrom} ${position(bin.start)} to ${position(bin.end)} · ${bin.called} called bases · GC ${bin.gc === null ? 'Not available' : bin.gc.toFixed(2) + '%'}`
const motion = i => ({ '--gc-delay': `${i * 600 / points.value.length}ms`, '--gc-duration': `${600 / points.value.length + 80}ms` })
</script>
<template>
  <div :key="data.sequence" ref="host" class="gc-profile" :class="{'is-revealed':visible}">
    <h5>GC across the selected region</h5><p>Non-overlapping 100 bp bins · final bin may be shorter · GC / called bases in each bin</p>
    <svg viewBox="0 0 740 182" role="img" aria-label="GC content across the selected sequence, 100 base bins" @mouseleave="activeIndex=null">
      <g v-for="tick in [0,50,100]" :key="tick"><line x1="44" :y1="146-tick*1.2" x2="714" :y2="146-tick*1.2" stroke="#e2e9f2"/><text x="36" :y="150-tick*1.2" text-anchor="end">{{tick}}%</text></g>
      <g v-for="(bin,i) in points" :key="bin.start" :data-bin="i" :data-gc="bin.gc ?? ''" :style="motion(i)">
        <line v-if="i && bin.y !== null && points[i-1].y !== null" class="gc-segment" :style="{transformOrigin:`${points[i-1].x}px ${points[i-1].y}px`}" :x1="points[i-1].x" :y1="points[i-1].y" :x2="bin.x" :y2="bin.y" stroke="#417dc6" stroke-width="2"/>
        <g v-if="bin.y !== null" class="gc-dot"><circle class="gc-dot-visible" :class="{'is-active':activeIndex===i}" :style="{transformOrigin:`${bin.x}px ${bin.y}px`}" :cx="bin.x" :cy="bin.y" r="3" fill="#417dc6"/></g>
        <circle v-if="bin.y !== null" class="gc-hit" :cx="bin.x" :cy="bin.y" r="9" fill="transparent" tabindex="0" role="img" :aria-label="label(bin)" @mouseenter="activeIndex=i" @mouseleave="activeIndex=null" @focus="activeIndex=i" @blur="activeIndex=null" @click="activeIndex=i"/>
      </g>
      <text x="44" y="169">{{ position(0) }}</text><text x="714" y="169" text-anchor="end">{{ position(data.length-1) }}</text>
    </svg>
    <div v-if="active" class="gc-tooltip" role="tooltip"><strong>GC {{active.gc.toFixed(2)}}%</strong><span>{{data.endpoint?.anchor_name || data.chrom}} {{position(active.start)}} to {{position(active.end)}}</span><small>{{active.called}} called bases</small></div>
    <p>{{ data.endpoint ? 'Offsets relative to ' + data.endpoint.anchor_name : 'Genomic coordinates' }} · displayed sequence direction, 5′ → 3′</p>
  </div>
</template>
<style scoped>
.gc-profile{border-top:1px solid #e4eaf1;padding-top:14px;position:relative;min-width:0}.gc-profile h5{font-size:13px;color:#345776;margin:0 0 12px;font-weight:600}.gc-profile p{font-size:11px;color:#64748b;line-height:1.7;margin:6px 0;text-wrap:pretty}.gc-profile svg{width:100%;height:auto;display:block;max-height:190px}.gc-profile text{font-size:12px;fill:#64748b;font-variant-numeric:tabular-nums}.gc-segment{pointer-events:none}.gc-segment,.gc-dot{opacity:0}.is-revealed .gc-segment{animation:gc-line-reveal var(--gc-duration) ease-out var(--gc-delay) both}.is-revealed .gc-dot{animation:gc-point-reveal 180ms ease-out var(--gc-delay) both}.gc-dot-visible{transition:transform 180ms ease-out;pointer-events:none}.gc-dot-visible.is-active{transform:scale(1.8);stroke:#fff;stroke-width:1}.gc-hit{cursor:pointer}.gc-hit:focus-visible{outline:none;stroke:#23415f;stroke-width:1}.gc-tooltip{position:absolute;top:55px;right:12px;max-width:calc(100% - 24px);z-index:4;padding:8px 12px;border:1px solid #d5dfeb;background:white;border-radius:5px;box-shadow:0 3px 10px #0f172a20;pointer-events:none;font-size:12px;line-height:1.7;color:#345776}.gc-tooltip span,.gc-tooltip small{display:block;font-variant-numeric:tabular-nums}.gc-tooltip small{font-size:10px;color:#64748b}@keyframes gc-line-reveal{from{opacity:0;transform:scale(0)}to{opacity:1;transform:scale(1)}}@keyframes gc-point-reveal{from{opacity:0}to{opacity:1}}@media(prefers-reduced-motion:reduce){.gc-segment,.gc-dot{opacity:1}.is-revealed .gc-segment,.is-revealed .gc-dot{animation:none}.gc-dot-visible{transition:none;transform:none!important}}@media(max-width:600px){.gc-profile text{font-size:17px}}
</style>
