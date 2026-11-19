<script setup>
import { computed } from 'vue'
const props = defineProps({ analysis: Object, primary: Object, comparison: Object, cumulative: Boolean })
const left=64, right=560, top=20, bottom=254
const edges=computed(()=>props.analysis.bin_edges)
const x=v=>left+(v-edges.value[0])/(edges.value.at(-1)-edges.value[0])*(right-left)
const ceiling=computed(()=>props.cumulative ? 1 : Math.ceil(100*Math.max(.01,...[props.primary,props.comparison].filter(Boolean).flatMap(p=>p.summary.histogram.map(n=>n/(p.summary.observed||1)))))/100)
const y=v=>bottom-v/ceiling.value*(bottom-top)
const ticks=computed(()=>Array.from({length:5},(_,i)=>edges.value[0]+i*(edges.value.at(-1)-edges.value[0])/4))
function step(profile){
  const values=profile?.values||[]
  if(!values.length)return ''
  let path=`M ${left} ${bottom}`, i=0
  while(i<values.length){let j=i+1;while(j<values.length&&values[j]===values[i])j++;path+=` H ${x(values[i]).toFixed(2)} V ${y(j/values.length).toFixed(2)}`;i=j}
  return path+` H ${right}`
}
function histogram(profile){
  if(!profile?.summary.observed)return ''
  let path=`M ${left} ${bottom}`
  profile.summary.histogram.forEach((n,i)=>{path+=` V ${y(n/profile.summary.observed).toFixed(2)} H ${x(edges.value[i+1]).toFixed(2)}`})
  return path+` V ${bottom}`
}
</script>
<template>
  <svg class="pi-distribution-svg" viewBox="0 0 580 316" role="img" :aria-label="cumulative ? 'Empirical cumulative distribution of observed log2 PI values' : 'Histogram of observed log2 PI values'">
    <title>{{ cumulative ? 'PI cumulative distribution' : 'Single-profile PI distribution' }}</title>
    <desc>Each observation is a measured gene / transcript record. Missing values are excluded. Download the observations for exact values.</desc>
    <g v-for="i in [0,1,2,3,4]" :key="i"><line :x1="left" :x2="right" :y1="y(ceiling*i/4)" :y2="y(ceiling*i/4)" stroke="#e9eff6"/><text :x="left-10" :y="y(ceiling*i/4)+4" text-anchor="end">{{ (ceiling*i/4*100).toFixed(cumulative?0:1) }}%</text></g>
    <path :d="`M${left} ${top}V${bottom}H${right}`" fill="none" stroke="#9cb9d3"/>
    <g v-for="tick in ticks" :key="tick"><text :x="x(tick)" :y="bottom+23" text-anchor="middle">{{ tick.toFixed(1) }}</text></g>
    <path v-if="!cumulative" :d="histogram(primary)" fill="#3e77ab" fill-opacity=".22" stroke="#3e77ab" stroke-width="2"/>
    <path v-else :d="step(primary)" fill="none" stroke="#3e77ab" stroke-width="2.4"/>
    <path v-if="comparison" :d="cumulative ? step(comparison) : histogram(comparison)" fill="none" stroke="#d69a41" stroke-width="2" stroke-dasharray="6 3"/>
    <text x="310" y="307" text-anchor="middle" class="pi-axis-title">log₂(Pause Index)</text>
    <text transform="translate(16 142) rotate(-90)" text-anchor="middle" class="pi-axis-title">{{ cumulative?'Cumulative proportion':'Proportion of records' }}</text>
  </svg>
</template>
