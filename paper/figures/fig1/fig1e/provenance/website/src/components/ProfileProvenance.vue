<script setup>
import { computed, onBeforeUnmount, ref, watch } from 'vue'
import { Download, ExternalLink } from '@lucide/vue'
import { getProfileProvenance } from '../api/client.js'
import { downloadBlob } from '../utils/data.js'
const props=defineProps({profile:Object})
const result=ref(null),loading=ref(false),error=ref('')
let revision=0
async function load(){const id=++revision;result.value=null;loading.value=true;error.value='';try{const data=await getProfileProvenance(props.profile.id);if(!data.record||data.record.alias!==props.profile.alias||data.record.gse!==props.profile.gse)throw Error('Profile metadata and catalog differ. Reload the page.');if(id===revision)result.value=data}catch(e){if(id===revision)error.value=e.message}finally{if(id===revision)loading.value=false}}
watch(()=>props.profile.id,load,{immediate:true});onBeforeUnmount(()=>{++revision})
const record=computed(()=>result.value?.record)
const unique=key=>[...new Set(record.value?.runs.map(r=>r[key]).filter(Boolean)||[])]
const gsms=computed(()=>unique('gsm')),titles=computed(()=>unique('sample_title')),platforms=computed(()=>unique('instrument')),organisms=computed(()=>unique('organism')),projects=computed(()=>unique('bioproject')),biosamples=computed(()=>unique('biosample'))
const layouts=computed(()=>unique('layout').map(v=>v==='SINGLE'?'Single-end':v==='PAIRED'?'Paired-end':v))
const lengths=computed(()=>[...new Set(record.value?.runs.flatMap(r=>r.raw_read_lengths)||[])])
const subseries=computed(()=>record.value?.runs.filter(r=>r.geo_series!==props.profile.gse)||[])
function download(){downloadBlob(`${props.profile.alias}_provenance.json`,JSON.stringify({checked_on:result.value.checked_on,record:record.value,sources:result.value.sources},null,2))}
</script>
<template>
 <section class="panel profile-provenance" aria-label="Profile provenance">
  <div class="panel-head"><div><h2>Provenance and experimental labels</h2><p v-if="record">GEO / ENA identifiers checked · {{record.runs.length}} sequencing run{{record.runs.length===1?'':'s'}}</p></div><button v-if="record" class="button secondary" @click="download"><Download :size="15"/>Download metadata</button></div>
  <p v-if="loading" class="provenance-message" role="status">Loading verified run metadata…</p>
  <div v-else-if="error" class="provenance-message" role="alert"><p>{{error}}</p><button class="button secondary" @click="load">Retry metadata</button></div>
  <template v-else-if="record">
   <dl class="detail-list wide">
    <dt>Original alias</dt><dd class="mono break-all">{{profile.alias}}</dd>
    <dt>GEO study</dt><dd><RouterLink :to="`/studies/${profile.gse}`" class="table-link">{{profile.gse}}</RouterLink> <a :href="`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${profile.gse}`" target="_blank" rel="noopener" class="inline-icon" aria-label="Open original GEO record"><ExternalLink :size="14"/></a><small v-if="subseries.length">Run metadata uses subseries {{[...new Set(subseries.map(r=>r.geo_series))].join(', ')}}; the GEO parent/child relationship was verified.</small></dd>
    <dt>Public sample title</dt><dd><div v-for="title in titles" :key="title">{{title}}</div></dd>
    <dt>Assay / cell label</dt><dd>{{profile.assay}} / {{profile.cell}}<small>Catalog labels retained from the source matrix.</small></dd>
    <template v-if="record.geo_review"><dt>Verified cell / condition</dt><dd>{{record.geo_review.cell_line}} · {{record.geo_review.culture_condition}}</dd><dt>Research source</dt><dd>{{record.geo_review.contact}} · {{record.geo_review.laboratory}}<small>{{record.geo_review.organization}}</small></dd></template>
    <dt>Source author label</dt><dd>{{profile.author || 'Not recorded'}}</dd>
    <dt>Organism</dt><dd>{{organisms.join(' / ') || 'Not recorded'}}</dd>
    <dt>Antibody subtype</dt><dd v-if="profile.antibody">{{profile.antibody}}<small>Pol II CTD subtype parsed from the source alias; antibody reagent not independently reviewed.</small></dd><dd v-else-if="profile.assay!=='mNET-seq'">Not applicable to Pol II CTD subtype grouping<small>This field describes mNET-seq Pol II states, rather than all antibodies used during library preparation.</small></dd><dd v-else>Not recorded in the source alias</dd>
    <dt>Source matrix</dt><dd class="mono break-all">{{profile.source_file}}</dd>
    <dt>GEO sample (GSM)</dt><dd class="accession-links"><a v-for="gsm in gsms" :key="gsm" :href="`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gsm}`" target="_blank" rel="noopener">{{gsm}} ↗</a><span v-if="!gsms.length">Not verified</span></dd>
    <dt>Sequencing runs (SRR)</dt><dd class="accession-links"><a v-for="r in record.runs" :key="r.run" :href="`https://www.ncbi.nlm.nih.gov/sra/${r.run}`" target="_blank" rel="noopener">{{r.run}} ↗</a></dd>
    <dt>Sequencing platform</dt><dd>{{platforms.join(' / ') || 'Not recorded'}}</dd>
    <dt>Library layout</dt><dd>{{layouts.join(' / ') || 'Not recorded'}}<small v-if="lengths.length">Raw read length: {{lengths.join(' / ')}} bp · existing FastQC reports</small></dd>
    <dt>BioProject / BioSample</dt><dd class="accession-links"><a v-for="id in projects" :key="id" :href="`https://www.ncbi.nlm.nih.gov/bioproject/${id}`" target="_blank" rel="noopener">{{id}} ↗</a><a v-for="id in biosamples" :key="id" :href="`https://www.ncbi.nlm.nih.gov/biosample/${id}`" target="_blank" rel="noopener">{{id}} ↗</a></dd>
    <dt>Reference assembly / annotation</dt><dd v-if="record.track_reference">{{record.track_reference.reference_name}} · {{record.track_reference.annotation_name}}<small>Verified reference for the connected genome tracks and sequence display.</small></dd><dd v-else>Track reference not available</dd>
    <template v-if="record.geo_review"><dt>GEO original processing</dt><dd>{{record.geo_review.original_processing_assembly}}<small>{{record.geo_review.original_processing_note}}</small></dd></template>
    <dt>Metric</dt><dd>log₂(Pause Index) · Existing matrix values<small>The displayed PI values are imported without recalculation.</small></dd>
   </dl>
   <details class="metadata-sources"><summary>Metadata sources and verification scope</summary><p>Checked {{result.checked_on}}. Local sample-to-run mappings come from the existing QC metadata; GSE, organism, platform and library layout are cross-checked against official ENA run records. GSM identifiers come from ENA experiment aliases or verified GEO sample-to-SRA experiment relations. The current track reference comes from the registered genome manifest. Source matrix values and labels are unchanged.</p><ul><li v-for="r in record.runs" :key="r.run"><a :href="`https://www.ebi.ac.uk/ena/browser/view/${r.run}`" target="_blank" rel="noopener">ENA {{r.run}}</a> · {{r.gsm}} · {{r.geo_series}} · {{r.instrument}} · {{r.layout}}<small>GSM mapping: {{r.gsm_source}}</small></li></ul><p v-if="record.geo_review"><a :href="record.geo_review.source_url" target="_blank" rel="noopener">Reviewed GEO sample {{record.geo_review.gsm}}</a> · Cell line, culture condition, research source and original processing reference.</p><p>PI calculation windows, normalization and original matrix annotation still require their own processing records. Technical runs may be merged into a single catalog profile.</p></details>
  </template>
 </section>
</template>
<style scoped>
.profile-provenance{min-width:0}.profile-provenance .panel-head{gap:16px;flex-wrap:wrap}.profile-provenance .panel-head p{font-size:12px}.profile-provenance .button{font-size:12px}.profile-provenance dd{overflow-wrap:anywhere}.profile-provenance dd small{display:block;font-size:11px;line-height:1.7;color:#64748b;margin-top:5px;text-wrap:pretty}.accession-links{display:flex;flex-wrap:wrap;gap:8px 14px;align-items:baseline}.accession-links a,.metadata-sources a{font-size:12px;color:#3e77ab}.metadata-sources{padding:16px 24px;border-top:1px solid #e4eaf1;font-size:12px;color:#64748b;line-height:1.8}.metadata-sources summary{cursor:pointer;color:#345776}.metadata-sources p{text-wrap:pretty}.metadata-sources ul{padding-left:18px}.provenance-message{padding:20px 24px;font-size:13px;color:#64748b}.profile-provenance :focus-visible{outline:2px solid #417dc6;outline-offset:2px}@media(max-width:760px){.profile-provenance .detail-list.wide{grid-template-columns:1fr;gap:6px}.profile-provenance dd{margin-bottom:10px}.metadata-sources{padding:16px 20px}}
</style>
<style scoped>
.metadata-sources li small{display:block;color:#64748b;font-size:11px}
</style>
