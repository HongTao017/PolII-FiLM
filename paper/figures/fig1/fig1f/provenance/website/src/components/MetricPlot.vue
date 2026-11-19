<script setup>
import { sampleName } from '../utils/sampleName.js'
import { computed } from 'vue'
import { assayColors } from '../utils/data.js'
const props=defineProps({ rows:Array })
const groups=computed(()=>[...new Set(props.rows.map(p=>p.assay==='mNET-seq'?`${p.assay} / ${p.antibody||'Unspecified'}`:p.assay))])
const range=computed(()=>{const values=props.rows.map(p=>p.log2_pi);return [Math.floor(Math.min(...values,0)-.5),Math.ceil(Math.max(...values,0)+.5)]})
const x=value=>230+(value-range.value[0])/(range.value[1]-range.value[0])*620
const ticks=computed(()=>Array.from({length:6},(_,i)=>range.value[0]+i*(range.value[1]-range.value[0])/5))
const height=computed(()=>groups.value.length*42+65)
function y(row,index){const group=row.assay==='mNET-seq'?`${row.assay} / ${row.antibody||'Unspecified'}`:row.assay;return groups.value.indexOf(group)*42+30+(index%5-2)*3}
</script>
<template><div class="metric-plot"><svg :viewBox="`0 0 930 ${height}`" role="img" aria-label="Measured log2 pausing indices by assay and antibody subtype"><g v-for="tick in ticks" :key="tick"><line :x1="x(tick)" :x2="x(tick)" y1="12" :y2="height-40" stroke="#e3e5e9"/><text :x="x(tick)" :y="height-20" text-anchor="middle" fill="#7b8086" font-size="12">{{ tick.toFixed(1) }}</text></g><line :x1="x(0)" :x2="x(0)" y1="12" :y2="height-40" stroke="#9aa1aa" stroke-dasharray="4 4"/><text v-for="(group,index) in groups" :key="group" x="12" :y="index*42+35" fill="#484f59" font-size="12">{{ group }}</text><circle v-for="(row,index) in rows" :key="row.id" :cx="x(row.log2_pi)" :cy="y(row,index)" r="4.3" :fill="assayColors[row.assay]" fill-opacity=".8"><title>{{ sampleName(row) }} · {{ row.log2_pi.toFixed(4) }}</title></circle></svg><p class="small">Each point represents a measured profile. The x-axis shows log2(Pause Index). Assays and antibodies are stratified; no cross-assay differential test is performed.</p></div></template>
