<script setup>
import { sampleName } from '../utils/sampleName.js'
import { ref, computed, watch, onBeforeUnmount } from 'vue'
import { getProfilePiAnalysis, getProfilePiValues, validatePiAnalysis } from '../api/client.js'
import { downloadBlob } from '../utils/data.js'
import LoadState from './LoadState.vue'
import EmptyState from './EmptyState.vue'
import PiDistributionPlot from './PiDistributionPlot.vue'
import ProfilePairHeatmap from './ProfilePairHeatmap.vue'
const props=defineProps({ profile:Object, catalog:Object })
const analysis=ref(null),primary=ref(null),comparison=ref(null),compareId=ref(''),loading=ref(true),error=ref(''),compareError=ref(''),compareLoading=ref(false)
let generation=0,comparisonGeneration=0
async function load(){
  const token=++generation;++comparisonGeneration;loading.value=true;error.value='';analysis.value=null;primary.value=null;comparison.value=null;compareId.value='';compareError.value='';compareLoading.value=false
  try{const a=await getProfilePiAnalysis();validatePiAnalysis(a,props.catalog);const p=await getProfilePiValues(a.version,props.profile.id);if(token!==generation)return;analysis.value=a;primary.value=p}
  catch(e){if(token===generation)error.value=e.message}
  finally{if(token===generation)loading.value=false}
}
watch(()=>props.profile.id,load,{immediate:true})
async function loadComparison(){
  const token=++comparisonGeneration;comparison.value=null;compareError.value='';compareLoading.value=false
  if(!compareId.value||!analysis.value)return
  compareLoading.value=true
  try{const result=await getProfilePiValues(analysis.value.version,compareId.value);if(token===comparisonGeneration)comparison.value=result}
  catch(e){if(token===comparisonGeneration)compareError.value=e.message}
  finally{if(token===comparisonGeneration)compareLoading.value=false}
}
watch([compareId,analysis],loadComparison)
onBeforeUnmount(()=>{++generation;++comparisonGeneration})
const compareProfile=computed(()=>props.catalog.profiles.find(p=>p.id===compareId.value))
const summary=computed(()=>primary.value?.summary)
const format=v=>v==null?'Not available':v.toLocaleString(undefined,{maximumFractionDigits:2})
function download(){downloadBlob(`${sampleName(props.profile)}_observed_log2_PI.tsv`,'profile_id\tsample_name\tlog2_PI\n'+primary.value.values.map(v=>`${props.profile.id}\t${sampleName(props.profile)}\t${v}`).join('\n'),'text/tab-separated-values;charset=utf-8')}
</script>
<template>
  <section class="profile-pi-analysis latest-panel" aria-label="Profile PI analysis">
    <div class="pi-section-heading"><div><h2>PI across gene / transcript records</h2><p>One profile, all measured records in the imported PI matrix.</p></div><button class="button secondary" :disabled="!primary?.values.length" @click="download">Download observed PI</button></div>
    <LoadState :loading="loading" :error="error" @retry="load"><template v-if="summary">
      <dl class="pi-summary-strip"><div><dt>Measured records</dt><dd>{{ format(summary.observed) }}</dd></div><div><dt>Missing records</dt><dd>{{ format(summary.missing) }}</dd></div><div><dt>Median log₂ PI</dt><dd>{{ format(summary.median) }}</dd></div><div><dt>Interquartile range</dt><dd>{{ format(summary.q25) }} to {{ format(summary.q75) }}</dd></div></dl>
      <div class="pi-comparison-controls"><label>Compare with another profile<select v-model="compareId" aria-label="Compare PI profile"><option value="">Current profile only · heatmap self-comparison</option><option v-for="p in catalog.profiles.filter(p=>p.id!==profile.id)" :key="p.id" :value="p.id">{{p.assay}} · {{p.cell}} · {{ sampleName(p) }}</option></select></label><span v-if="compareLoading" role="status">Loading comparison…</span><span v-if="compareError" role="alert">{{compareError}} <button class="reset-button" @click="loadComparison">Retry comparison</button></span></div>
      <div class="pi-curve-legend"><span><i class="pi-legend-primary"/>{{ sampleName(profile) }} · {{format(summary.observed)}} measured</span><span v-if="comparison"><i class="pi-legend-comparison"/>{{sampleName(compareProfile)}} · {{format(comparison.summary.observed)}} measured</span></div>
      <div v-if="summary.observed" class="detail-grid"><section class="panel"><div class="panel-head"><div><h3>Single-profile PI distribution</h3><p>60 common bins · proportion of measured records</p></div></div><PiDistributionPlot :analysis="analysis" :primary="primary" :comparison="comparison"/></section><section class="panel"><div class="panel-head"><div><h3>PI cumulative distribution</h3><p>Fraction of measured records with log₂ PI ≤ x</p></div></div><PiDistributionPlot :analysis="analysis" :primary="primary" :comparison="comparison" cumulative/></section></div>
      <EmptyState v-else title="No observed PI values" description="This profile has no finite PI values in the current matrix." to="/samples"/>
      <p class="pi-method-note">Original log₂ PI values; no additional log transform or missing-value imputation. True zero values are retained ({{format(summary.zeros)}} records). Each curve uses its own observed-record denominator. The ECDF is exact, including ties. Imported PI calculation protocols remain pending review.</p>
    </template></LoadState>
    <ProfilePairHeatmap :profile="profile" :catalog="catalog" v-model:comparison-id="compareId"/>
  </section>
</template>
