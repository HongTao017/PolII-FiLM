<script setup>
import {computed,ref} from 'vue'
const props=defineProps({catalog:{type:Object,required:true}})
const emit=defineEmits(['select'])
const total=computed(()=>props.catalog.profiles.length)
const palette=['#739bc6','#96c979','#f3cb64','#e57d78','#79bed6','#58a789','#a988c7','#e493c5','#7190ce','#e3aa6e','#8dc4b0','#bca66c','#b892d5','#e3a0a6','#9ac8e0','#a8bd74','#6eafbd','#ca9dc3','#afbacf','#ecb982','#8bc6a2','#c5aecc','#98b0cd','#d5c878','#89b9a9','#d7a1a0','#9cadcf','#b3d494','#deb488','#a7cbd3','#c5a3bb']
const page=ref(1),pageSize=10,hovered=ref('')
const cells=computed(()=>{
 const counts=new Map()
 for(const p of props.catalog.profiles){const name=p.cell||'Not recorded';counts.set(name,(counts.get(name)||0)+1)}
 let offset=0
 return [...counts].sort((a,b)=>b[1]-a[1]||a[0].localeCompare(b[0])).map(([name,count],i)=>{
  const start=offset,totalCount=total.value;offset+=count
  const startAngle=-Math.PI/2+start/totalCount*Math.PI*2,endAngle=-Math.PI/2+offset/totalCount*Math.PI*2
  const point=angle=>[150+118*Math.cos(angle),150+118*Math.sin(angle)].join(' ')
  const d=count===totalCount?'M 150 32 A 118 118 0 1 1 150 268 A 118 118 0 1 1 150 32 Z':'M 150 150 L '+point(startAngle)+' A 118 118 0 '+(count/totalCount>.5?1:0)+' 1 '+point(endAngle)+' Z'
  return {name,count,percentage:count/totalCount*100,color:palette[i%palette.length],d}
 })
})
const pages=computed(()=>Math.max(1,Math.ceil(cells.value.length/pageSize)))
const visible=computed(()=>cells.value.slice((page.value-1)*pageSize,page.value*pageSize))
const active=computed(()=>cells.value.find(c=>c.name===hovered.value))
function changePage(delta){page.value=Math.min(pages.value,Math.max(1,page.value+delta));hovered.value=''}
const label=c=>'Browse '+c.name+': '+c.count+' profiles ('+c.percentage.toFixed(1)+'%)'
</script>
<template>
 <section class="panel cell-composition" aria-label="Cell line composition">
  <div class="panel-head"><div><h2>Cell line composition</h2><p>Share of profile records by original cell label</p></div></div>
  <div v-if="total" class="cell-pie-layout">
   <figure class="cell-pie-figure">
    <svg class="cell-pie-svg" viewBox="0 0 300 300" role="img" :aria-label="'Cell line proportions: '+total+' profiles across '+cells.length+' original cell labels. Use the legend buttons to browse each cell line.'">
     <title>Cell line proportions across {{total}} profile records</title>
     <path v-for="item in cells" :key="item.name" class="cell-pie-slice" :d="item.d" :fill="item.color" :data-cell="item.name" :data-count="item.count" :data-share="item.percentage" :opacity="hovered && hovered!==item.name ? .48 : 1" stroke="#fff" stroke-width=".8" @pointerenter="hovered=item.name" @pointerleave="hovered=''" @click="emit('select',{cell:item.name})"><title>{{item.name}}: {{item.count}} profiles · {{item.percentage.toFixed(1)}}%</title></path>
    </svg>
    <div v-if="active" class="cell-pie-tooltip" role="status"><strong>{{active.name}}</strong><span>{{active.count}} profiles · {{active.percentage.toFixed(1)}}%</span></div>
    <figcaption>{{total}} profiles · {{cells.length}} cell labels</figcaption>
   </figure>
   <div class="cell-pie-legend">
    <ul class="cell-pie-legend-list">
     <li v-for="item in visible" :key="item.name">
      <button type="button" :aria-label="label(item)" :title="label(item)" :data-cell="item.name" :class="{'cell-legend-active':hovered===item.name}" @pointerenter="hovered=item.name" @pointerleave="hovered=''" @focus="hovered=item.name" @blur="hovered=''" @click="emit('select',{cell:item.name})"><i :style="{background:item.color}" aria-hidden="true"/><span class="cell-legend-name">{{item.name}}</span><span class="cell-legend-value">{{item.count}} <small>{{item.percentage.toFixed(1)}}%</small></span></button>
     </li>
    </ul>
    <div class="cell-legend-pagination"><button type="button" aria-label="Previous cell legend page" :disabled="page===1" @click="changePage(-1)">‹</button><span aria-live="polite">{{page}} / {{pages}}</span><button type="button" aria-label="Next cell legend page" :disabled="page===pages" @click="changePage(1)">›</button></div>
   </div>
  </div>
  <p v-else class="cell-pie-empty">No cell records are available in this catalog.</p>
  <div class="panel-foot">Select a slice or legend item to browse profiles →</div>
 </section>
</template>
<style scoped>
.cell-composition{display:flex;flex-direction:column;min-width:0;overflow:hidden}
.cell-pie-layout{display:grid;grid-template-columns:minmax(170px,1fr) minmax(190px,1.05fr);gap:12px;align-items:center;padding:20px 20px 16px;flex:1}
.cell-pie-figure{margin:0;min-width:0;position:relative;text-align:center}
.cell-pie-svg{display:block;width:100%;max-width:300px;aspect-ratio:1;margin:0 auto;overflow:visible}
.cell-pie-slice{cursor:pointer}
.cell-pie-figure figcaption{font-size:12px;color:#647b95;line-height:1.7;font-variant-numeric:tabular-nums}
.cell-pie-tooltip{position:absolute;left:50%;top:0;transform:translateX(-50%);display:flex;flex-direction:column;gap:3px;max-width:100%;width:max-content;padding:10px 13px;border:1px solid #98bedf;border-radius:5px;background:white;box-shadow:0 2px 5px #102e4f18;color:#244e78;pointer-events:none;font-size:12px;text-align:left;font-variant-numeric:tabular-nums}
.cell-pie-tooltip strong{font-size:13px;font-weight:600;overflow-wrap:anywhere}
.cell-pie-legend{min-width:0}
.cell-pie-legend-list{list-style:none;padding:0;margin:0;display:grid;grid-template-columns:1fr;gap:2px;min-height:300px;align-content:start}
.cell-pie-legend-list li{min-width:0}
.cell-pie-legend-list button{display:grid;grid-template-columns:10px minmax(0,1fr) auto;gap:8px;align-items:center;width:100%;min-height:29px;padding:4px 5px;border:1px solid transparent;border-radius:4px;background:white;color:#345477;text-align:left;font:inherit;font-size:13px;cursor:pointer}
.cell-pie-legend-list button:hover,.cell-pie-legend-list .cell-legend-active{border-color:#d5e3ef;background:#f4f8fc}
.cell-pie-legend-list i{width:10px;height:10px;border-radius:2px}
.cell-legend-name{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.cell-legend-value{font-size:12px;font-variant-numeric:tabular-nums;white-space:nowrap}
.cell-legend-value small{font-size:11px;color:#647b95;margin-left:5px}
.cell-legend-pagination{display:flex;align-items:center;justify-content:center;gap:10px;margin-top:8px;color:#647b95;font-size:12px;font-variant-numeric:tabular-nums}
.cell-legend-pagination button{background:white;border:1px solid #dce5ef;color:#244e78;border-radius:4px;width:27px;height:27px;font-size:20px;cursor:pointer}
.cell-legend-pagination button:disabled{opacity:.4;cursor:default}
.cell-composition button:focus-visible{outline:2px solid #3e77ab;outline-offset:2px}
.cell-pie-empty{padding:24px;color:#647b95;font-size:14px}
.cell-composition .panel-foot{margin-top:auto}
@media(max-width:520px){.cell-pie-layout{grid-template-columns:1fr;padding:16px;gap:18px}.cell-pie-svg{max-width:250px}.cell-pie-legend-list{grid-template-columns:repeat(2,minmax(0,1fr));min-height:150px;gap:5px 10px}.cell-pie-legend-list button{font-size:12px;gap:5px}.cell-legend-value{font-size:11px}.cell-legend-value small{font-size:10px;margin-left:2px}}
</style>

