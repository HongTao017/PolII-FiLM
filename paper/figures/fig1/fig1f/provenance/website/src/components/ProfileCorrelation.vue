<script setup>
import { sampleName } from '../utils/sampleName.js'
import {ref,computed,watch,nextTick,onMounted} from 'vue'
import {getPiSpearmanHeatmap,validatePiAnalysis} from '../api/client.js'
import {assays,assayColors,downloadBlob} from '../utils/data.js'
import LoadState from './LoadState.vue'
import EmptyState from './EmptyState.vue'
const props=defineProps({catalog:Object})
const analysis=ref(null),loading=ref(true),error=ref(''),assay=ref(''),cell=ref(''),study=ref(''),sort=ref('cluster'),zoom=ref('100'),canvas=ref(null),first=ref(''),second=ref(''),hover=ref(null)
const coordinates={x:190,y:260,side:1120,size:1600,labelGap:30}
const assayPalette={'GRO-seq':'#3659d9','PRO-seq':assayColors['PRO-seq'],'mNET-seq':'#a44ac5'}
const palette=['#4477aa','#ee7733','#228833','#cc3377','#66ccee','#aa3377','#999933','#004488','#ddaa33','#882255','#44aa99','#6699cc','#bb5566','#117733','#aa4499','#88ccee','#332288','#cc6677','#999955','#661100','#44bb99','#77aadd','#ee8866','#eedd88','#ffaabb','#99ddff','#44aa77','#bbbbbb','#7788aa','#ddcc77','#885588']
const cells=computed(()=>[...new Set(props.catalog.profiles.map(p=>p.cell))].sort())
const studies=computed(()=>[...new Set(props.catalog.profiles.map(p=>p.gse))].sort((a,b)=>a.localeCompare(b,undefined,{numeric:true})))
const index=computed(()=>new Map((analysis.value?.profile_order||[]).map((id,i)=>[id,i])))
const profiles=computed(()=>new Map(props.catalog.profiles.map(p=>[p.id,p])))
const filtered=computed(()=>(analysis.value?.profile_order||[]).map(id=>profiles.value.get(id)).filter(p=>(!assay.value||p.assay===assay.value)&&(!cell.value||p.cell===cell.value)&&(!study.value||p.gse===study.value)))
const tree=computed(()=>analysis.value?.clusters[filtered.value.map(p=>index.value.get(p.id)).join(',')])
const clustered=computed(()=>sort.value==='cluster'&&tree.value?.available)
const rows=computed(()=>clustered.value?tree.value.order.map(i=>profiles.value.get(analysis.value.profile_order[i])):sort.value==='cluster'?filtered.value:filtered.value.slice().sort((a,b)=>String(a[sort.value]).localeCompare(String(b[sort.value]),undefined,{numeric:true})||a.id.localeCompare(b.id)))
const pair=computed(()=>{
 const [a,b]=hover.value||[first.value,second.value],i=index.value.get(a),j=index.value.get(b)
 if(i===undefined||j===undefined||!analysis.value)return null
 return {a:profiles.value.get(a),b:profiles.value.get(b),r:analysis.value.spearman[i][j],n:analysis.value.overlaps[i][j]}
})
const command='python reproduce-spearman.py --base-url '+window.location.origin+' --data ./pi-input --output ./pi-spearman.json'
const cellColor=name=>palette[cells.value.indexOf(name)%palette.length]
async function load(){
 loading.value=true;error.value=''
 try{const result=await getPiSpearmanHeatmap();validatePiAnalysis(result,props.catalog);if(result.method.correlation!=='Spearman'||!result.spearman||!result.clusters)throw Error('Spearman heatmap data are unavailable.');analysis.value=result}catch(e){error.value=e.message}finally{loading.value=false}
}
load()
onMounted(()=>{if(window.location.hash.endsWith('#sample-correlations'))document.getElementById('sample-correlations')?.scrollIntoView()})
// Reference rainbow scale encodes magnitude. Signed coefficients remain inspectable.
function color(r){
 if(r===null)return '#e4e8ed'
 const t=Math.abs(r),c=o=>Math.round(255*Math.max(0,Math.min(1,1.5-Math.abs(4*t-o))))
 return 'rgb('+[c(3),c(2),c(1)].join(',')+')'
}
function draw(mark=true){
 const node=canvas.value,a=analysis.value,list=rows.value
 if(!node||!a||list.length<2)return
 const ctx=node.getContext('2d'),{x,y,side,size,labelGap}=coordinates,unit=side/list.length
 ctx.clearRect(0,0,size,size);ctx.fillStyle='#fff';ctx.fillRect(0,0,size,size)
 ctx.textBaseline='middle';ctx.textAlign='left';ctx.fillStyle='#244e78';ctx.font='bold 21px Arial'
 ctx.fillText('Sample PI correlations · '+list.length+' profiles',35,25)
 ctx.font='15px Arial';ctx.fillText(clustered.value?'Average linkage · distance = 1 − |Spearman ρ|':'Metadata order · no clustering tree',35,53)
 ctx.font='12px Arial'
 cells.value.forEach((name,i)=>{
  const lx=35+(i%8)*126,ly=83+Math.floor(i/8)*23
  ctx.fillStyle=cellColor(name);ctx.fillRect(lx,ly-5,9,9);ctx.fillStyle='#344b61';ctx.fillText(name,lx+14,ly)
 })
 assays.forEach((name,i)=>{ctx.fillStyle=assayPalette[name];ctx.fillRect(550+i*145,165,12,12);ctx.fillStyle='#344b61';ctx.fillText(name,568+i*145,171)})
 const bx=1130,by=74,bw=370
 ctx.font='16px Arial';ctx.fillStyle='#244e78';ctx.fillText('|Spearman correlation|',bx,49)
 for(let k=0;k<bw;k++){ctx.fillStyle=color(k/(bw-1));ctx.fillRect(bx+k,by,1,18)}
 ctx.font='13px Arial';ctx.textAlign='center'
 ;[0,.25,.5,.75,1].forEach(v=>{ctx.fillStyle='#344b61';ctx.fillText(v.toFixed(2),bx+v*bw,by+34)})
 ctx.textAlign='left';ctx.fillStyle='#526e8d'
 ctx.fillText(a.feature_count.toLocaleString('en-US')+' gene / transcript records',bx,133)
 ctx.fillText('Gray: unavailable · labels: author / source',bx,158)
 list.forEach((p,i)=>list.forEach((q,j)=>{ctx.fillStyle=color(a.spearman[index.value.get(p.id)][index.value.get(q.id)]);ctx.fillRect(x+j*unit,y+i*unit,unit+.2,unit+.2)}))
 if(clustered.value){
  const positions=new Map(list.map((p,i)=>[index.value.get(p.id),y+(i+.5)*unit]))
  const points=tree.value.indices.map(i=>({x:x-36,y:positions.get(i)}))
  const maxHeight=Math.max(...tree.value.linkage.map(r=>r[2]),1e-12)
  ctx.strokeStyle='#a3acb6';ctx.lineWidth=1
  tree.value.linkage.forEach(([left,right,height])=>{
   const l=points[left],r=points[right],parent={x:x-36-height/maxHeight*125,y:(l.y+r.y)/2}
   ctx.beginPath();ctx.moveTo(l.x,l.y);ctx.lineTo(parent.x,l.y);ctx.lineTo(parent.x,r.y);ctx.lineTo(r.x,r.y);ctx.stroke();points.push(parent)
  })
  ctx.font='12px Arial';ctx.fillStyle='#526e8d';ctx.textAlign='center';ctx.fillText('1 − |ρ|',x-100,y-46)
  ctx.fillText(maxHeight.toFixed(2),x-161,y-26);ctx.fillText('0',x-36,y-26)
 }
 ctx.font=Math.min(16,unit*.87)+'px Arial'
 list.forEach((p,i)=>{
  const o=i*unit,label=(p.author||'Unknown')+' · '+String(index.value.get(p.id)+1).padStart(3,'0')
  ctx.fillStyle=cellColor(p.cell);ctx.fillRect(x-30,y+o,12,unit);ctx.fillRect(x+o,y-30,unit,12)
  ctx.fillStyle=assayPalette[p.assay];ctx.fillRect(x-15,y+o,12,unit);ctx.fillRect(x+o,y-15,unit,12)
  ctx.fillStyle='#344b61';ctx.textAlign='left';ctx.fillText(label,x+side+labelGap,y+(i+.5)*unit,230)
  ctx.save();ctx.translate(x+(i+.5)*unit,y+side+labelGap);ctx.rotate(Math.PI/2);ctx.fillText(label,0,0,size-y-side-labelGap-50);ctx.restore()
 })
 ctx.strokeStyle='#8e9daa';ctx.lineWidth=1;ctx.strokeRect(x,y,side,side)
 if(mark&&pair.value){
  const i=list.findIndex(p=>p.id===pair.value.a.id),j=list.findIndex(p=>p.id===pair.value.b.id)
  if(i>=0&&j>=0){ctx.strokeStyle='#fff';ctx.lineWidth=2;ctx.strokeRect(x+j*unit,y+i*unit,unit,unit)}
 }
 ctx.font='14px Arial';ctx.fillStyle='#526e8d';ctx.textAlign='left'
 ctx.fillText('Shared observations per pair; missing PI excluded. Absolute correlation does not distinguish positive from inverse relationships.',35,1577)
}
watch([rows,analysis],async()=>{
 hover.value=null
 if(!rows.value.some(p=>p.id===first.value))first.value=rows.value[0]?.id||''
 if(!rows.value.some(p=>p.id===second.value))second.value=rows.value[1]?.id||rows.value[0]?.id||''
 await nextTick();draw()
},{immediate:true})
watch([pair,canvas],()=>nextTick(draw))
function pointer(e,lock=false){
 const rect=canvas.value.getBoundingClientRect(),{x,y,side,size}=coordinates
 const px=(e.clientX-rect.left)/rect.width*size,py=(e.clientY-rect.top)/rect.height*size
 if(px<x||px>=x+side||py<y||py>=y+side){hover.value=null;return}
 const ids=[rows.value[Math.floor((py-y)/side*rows.value.length)].id,rows.value[Math.floor((px-x)/side*rows.value.length)].id]
 if(lock){first.value=ids[0];second.value=ids[1];hover.value=null}else if(hover.value?.[0]!==ids[0]||hover.value?.[1]!==ids[1])hover.value=ids
}
function exportPairs(){
 const fields=['profile_a','sample_name_a','alias_a','assay_a','cell_a','GSE_a','author_a','profile_b','sample_name_b','alias_b','assay_b','cell_b','GSE_b','author_b','abs_spearman_rho','signed_spearman_rho','shared_records'],lines=[]
 rows.value.forEach((a,i)=>rows.value.slice(i).forEach(b=>{
  const x=index.value.get(a.id),y=index.value.get(b.id),r=analysis.value.spearman[x][y]
  lines.push([a.id,sampleName(a),a.alias,a.assay,a.cell,a.gse,a.author,b.id,sampleName(b),b.alias,b.assay,b.cell,b.gse,b.author,r===null?'':Math.abs(r),r??'',analysis.value.overlaps[x][y]].map(v=>String(v??'').replace(/[\t\r\n]/g,' ')).join('\t'))
 }))
 downloadBlob('NascentAtlas_PI_Spearman_correlations.tsv',[fields.join('\t'),...lines].join('\n'),'text/tab-separated-values;charset=utf-8')
}
function exportImage(){
 draw(false);const a=document.createElement('a');a.download='NascentAtlas_PI_Spearman_clustermap.png';a.href=canvas.value.toDataURL('image/png');a.click();draw()
}
function exportMethod(){
 const a=analysis.value
 downloadBlob('NascentAtlas_Spearman_method.json',JSON.stringify({method:a.method,source_hashes:a.source_hashes,input_sha256:a.input_sha256,feature_count:a.feature_count,filters:{assay:assay.value,cell:cell.value,gse:study.value},display_order:rows.value.map(p=>p.id),clustering:clustered.value?tree.value:null},null,2),'application/json')
}
</script>
<template>
 <section id="sample-correlations" class="panel latest-panel clustered-pi-heatmap" aria-label="Sample PI correlations">
  <div class="panel-head"><div><h2>Sample correlation heatmap</h2><p>Hierarchical clustering of shared PI records · absolute Spearman correlation.</p></div><button class="button secondary" :disabled="!analysis||!rows.length" @click="exportPairs">Download pairwise results</button></div>
  <div class="pi-heatmap-body"><LoadState :loading="loading" :error="error" @retry="load"><template v-if="analysis">
   <div class="pi-heatmap-controls">
    <label>Sequencing assay<select v-model="assay" aria-label="Heatmap assay"><option value="">All assays</option><option v-for="a in assays" :key="a">{{a}}</option></select></label>
    <label>Cell line<select v-model="cell" aria-label="Heatmap cell line"><option value="">All cell labels</option><option v-for="c in cells" :key="c">{{c}}</option></select></label>
    <label>GEO study<select v-model="study" aria-label="Heatmap GEO study"><option value="">All studies</option><option v-for="g in studies" :key="g">{{g}}</option></select></label>
    <label>Profile order<select v-model="sort" aria-label="Heatmap sort order"><option value="cluster">Hierarchical clustering</option><option value="assay">Assay</option><option value="cell">Cell line</option><option value="gse">GSE study</option><option value="author">Author / source</option></select></label>
   </div>
   <p class="pi-method-note">{{rows.length}} profiles · {{analysis.feature_count.toLocaleString()}} gene / transcript records. {{clustered?'Average linkage on 1 − |ρ|; both axes share the clustered order.':sort==='cluster'?(tree?.reason||'No profiles match these filters.'):'Metadata order selected; the clustering tree is hidden.'}} Outer strips: cell line. Inner strips: sequencing assay.</p>
   <details class="pi-calculation-details">
    <summary>How this plot is calculated</summary>
    <ol>
     <li><b>Input:</b> {{analysis.feature_count.toLocaleString()}} gene / transcript records × {{analysis.profile_order.length}} profiles, using the original imported log₂ PI matrix. Each position on either axis is a profile.</li>
     <li><b>Shared observations:</b> for each pair, retain only records measured in both profiles. Require at least {{analysis.method.minimum_overlap}} shared records. Missing values are excluded, never replaced with zero; constant vectors have no correlation.</li>
     <li><b>Spearman ρ:</b> rank the shared values separately for each profile, assign average ranks to ties, then calculate the Pearson correlation of those ranks.</li>
     <li><b>Color:</b> display |ρ| from 0 (blue) to 1 (red). A negative and a positive correlation of the same magnitude have the same color. Inspect a pair or download the results to see its signed ρ.</li>
     <li><b>Tree:</b> use distance 1 − |ρ| and average-linkage hierarchical clustering (SciPy; optimal ordering off). Filtering selects a tree computed for that exact subset. Unavailable distances are never imputed. Metadata sorting hides the tree.</li>
    </ol>
    <p class="pi-method-note">No additional normalization, PCA imputation or batch correction. Across all profiles, pairwise overlap ranges from {{analysis.summary.minimum_shared_records.toLocaleString()}} to {{analysis.summary.maximum_shared_records.toLocaleString()}} records. Assay, cell and author labels describe the existing catalog; correlations cannot establish which factor caused a difference.</p>
    <p class="pi-method-note">Implementation: NumPy {{analysis.method.numpy_version}} · SciPy {{analysis.method.scipy_version}}. <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html" target="_blank" rel="noopener">Spearman reference</a> · <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html" target="_blank" rel="noopener">Average-linkage reference</a></p>
    <details><summary>Input files and reproducibility</summary><pre>{{JSON.stringify(analysis.source_hashes,null,2)}}</pre><p class="pi-method-note">Input SHA-256: <code>{{analysis.input_sha256}}</code></p><p>Install NumPy and SciPy, download the script, then run:</p><pre>pip install numpy=={{analysis.method.numpy_version}} scipy=={{analysis.method.scipy_version}}
{{command}}</pre><div class="pi-plot-actions"><a class="button secondary" href="./data/analysis/reproduce-spearman.py" download>Download calculation script</a><button class="button secondary" @click="exportMethod">Download method &amp; tree</button></div></details>
   </details>
   <template v-if="rows.length>=2">
    <div class="pi-plot-actions"><label>Plot size<select v-model="zoom" aria-label="Heatmap zoom"><option value="100">Fit to panel</option><option value="150">150%</option><option value="200">200%</option></select></label><button class="button secondary" @click="exportImage">Download heatmap PNG</button><span class="pi-method-note">Hover to inspect · click to keep a pair · enlarge to read labels</span></div>
    <div class="pi-clustermap-scroll" tabindex="0" role="region" aria-label="Scrollable sample clustermap"><canvas ref="canvas" class="pi-heatmap-canvas" :style="{width:zoom+'%'}" width="1600" height="1600" role="img" :aria-label="'Absolute Spearman correlation heatmap for '+rows.length+' profiles, with '+(clustered?'an average-linkage dendrogram, ':'')+'cell and assay annotation strips and author labels. Use the row and column selectors below to inspect any pair with a keyboard.'" @mousemove="pointer($event)" @mouseleave="hover=null" @click="pointer($event,true)"/></div>
    <aside v-if="pair" class="pi-pair-details" aria-label="Selected correlation pair"><div><h3>Pair details</h3><strong>{{pair.r===null?'Unavailable':Math.abs(pair.r).toFixed(3)}}</strong><p>|Spearman ρ| · <span class="pi-signed-rho">Signed ρ: {{pair.r===null?'Unavailable':pair.r.toFixed(3)}}</span><br>{{pair.n.toLocaleString()}} shared records</p><p v-if="pair.r===null" class="pi-method-note">Insufficient overlap or a constant PI vector.</p></div><dl><dt>Row profile</dt><dd><RouterLink :to="'/samples/'+pair.a.id">{{ sampleName(pair.a) }}</RouterLink><br>{{pair.a.assay}} · {{pair.a.cell}}<br>{{pair.a.gse}} · {{pair.a.author||'Source unknown'}}</dd></dl><dl><dt>Column profile</dt><dd><RouterLink :to="'/samples/'+pair.b.id">{{ sampleName(pair.b) }}</RouterLink><br>{{pair.b.assay}} · {{pair.b.cell}}<br>{{pair.b.gse}} · {{pair.b.author||'Source unknown'}}</dd></dl></aside>
    <div class="pi-heatmap-controls pair-controls"><label>Row profile<select v-model="first" aria-label="Correlation row profile" @change="hover=null"><option v-for="(p,i) in rows" :key="p.id" :value="p.id">{{i+1}} · {{ sampleName(p) }} · {{p.assay}} · {{p.cell}}</option></select></label><label>Column profile<select v-model="second" aria-label="Correlation column profile" @change="hover=null"><option v-for="(p,i) in rows" :key="p.id" :value="p.id">{{i+1}} · {{ sampleName(p) }} · {{p.assay}} · {{p.cell}}</option></select></label></div>
   </template>
   <EmptyState v-else title="Choose at least two profiles" description="Broaden the assay, cell or study filters to compare profiles."/>
  </template></LoadState></div>
 </section>
</template>
<style scoped>
.pi-calculation-details{border:1px solid #dce5ef;border-radius:8px;padding:14px 18px;margin:18px 0;color:#344f6a;font-size:14px;line-height:1.8}
.pi-calculation-details summary{cursor:pointer;color:#244e78;font-weight:600}
.pi-calculation-details li{margin:10px 0;text-wrap:pretty}
.pi-calculation-details pre{font:12px/1.7 monospace;white-space:pre-wrap;overflow-wrap:anywhere;background:#f5f8fc;padding:14px;border-radius:6px}
.pi-calculation-details code{overflow-wrap:anywhere}
.pi-plot-actions{display:flex;align-items:center;gap:14px;flex-wrap:wrap;margin:16px 0}
.pi-plot-actions label{display:flex;align-items:center;gap:10px;font-size:13px;color:#244e78}
.pi-plot-actions select{font:inherit;color:inherit;background:white;border:1px solid #cddae9;border-radius:6px;padding:8px}
.pi-clustermap-scroll{max-width:100%;overflow:auto;border:1px solid #e3eaf2;border-radius:6px}
.pi-clustermap-scroll:focus-visible{outline:2px solid #3659d9;outline-offset:3px}
.pi-heatmap-canvas{max-width:none;min-width:780px;touch-action:auto}
.pi-pair-details{display:grid;grid-template-columns:minmax(180px,.8fr) minmax(0,1fr) minmax(0,1fr);gap:24px;margin-top:20px;font-variant-numeric:tabular-nums}
.pi-pair-details h3{margin:0 0 10px;font-size:16px}.pi-pair-details p{font-size:13px;line-height:1.8}
.pi-pair-details dl{display:block}.pi-pair-details dt{margin-top:0}
@media(max-width:700px){.pi-pair-details{grid-template-columns:1fr;gap:16px}.pi-calculation-details{padding:12px}.pi-calculation-details ol{padding-left:20px}.pi-plot-actions .button{font-size:12px}}
</style>

