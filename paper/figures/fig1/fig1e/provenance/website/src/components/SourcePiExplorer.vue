<script setup>
import { computed, onBeforeUnmount, ref, watch } from 'vue'
import { Search } from '@lucide/vue'
import { getGene, getGeneIndex } from '../api/client.js'
import ExperimentalSources from './ExperimentalSources.vue'
import GenePiGroups from './GenePiGroups.vue'
import GeneGenomeBrowser from './GeneGenomeBrowser.vue'
import GeneFeatures from './GeneFeatures.vue'
import GeneSequence from './GeneSequence.vue'
const props = defineProps({ catalog: { type: Object, required: true }, rows: { type: Array, required: true }, allRows: { type: Array, default: () => [] }, symbol: { type: String, default: '' }, transcript: { type: String, default: '' } })
defineEmits(['reset'])
const sourceTable = ref(null), matchingRows = ref(props.rows)
const query = ref('SAMD11'), loadedGene = ref(null), selectedTranscript = ref(''), loading = ref(false), error = ref('')
let requestId = 0
const isDetail = computed(() => Boolean(props.symbol))
async function loadDistribution() {
  if (isDetail.value) return
  const value = query.value.trim()
  const id = ++requestId
  loadedGene.value = null; selectedTranscript.value = ''; error.value = ''
  if (!value) { error.value = 'Enter a gene symbol or RefSeq transcript.'; loading.value = false; return }
  loading.value = true
  try {
    const index = await getGeneIndex(), key = value.toLowerCase()
    const entry = index.find(row => row.symbol.toLowerCase() === key || row.transcripts.some(t => t.toLowerCase() === key))
    if (!entry) throw new Error('Gene or transcript not found. Try SAMD11, NOC2L or TP53.')
    const gene = await getGene(entry.symbol)
    if (!gene?.transcripts?.length) throw new Error('No PI data is available for this gene.')
    if (id === requestId) { loadedGene.value = gene; selectedTranscript.value = gene.transcripts.find(t => t.id.toLowerCase() === key)?.id || gene.transcripts[0].id }
  } catch (failure) { if (id === requestId) error.value = failure.message }
  finally { if (id === requestId) loading.value = false }
}
watch(isDetail, detail => { if (!detail) loadDistribution() }, { immediate: true })
onBeforeUnmount(() => { ++requestId })
const symbol = computed(() => props.symbol || loadedGene.value?.symbol || '')
const transcript = computed(() => props.transcript || selectedTranscript.value)
const current = computed(() => loadedGene.value?.transcripts.find(row => row.id === selectedTranscript.value))
const allMetrics = computed(() => {
  if (isDetail.value) return props.allRows
  if (!current.value) return []
  const byId = new Map(props.catalog.profiles.map(profile => [profile.id, profile]))
  return props.catalog.metric.profile_order.map((id, i) => ({ ...byId.get(id), log2_pi: current.value.log2_pi[i] })).filter(row => row.id && Number.isFinite(row.log2_pi))
})
const measured = computed(() => {
  const ids = new Set(matchingRows.value.map(row => row.id))
  return allMetrics.value.filter(row => ids.has(row.id))
})
function reset() { sourceTable.value?.reset() }
</script>
<template>
  <div class="source-pi-explorer">
    <ExperimentalSources ref="sourceTable" :rows="rows" :symbol="symbol && isDetail ? symbol : ''" :transcript="isDetail ? transcript : ''" @filtered="matchingRows = $event" @reset="$emit('reset')"/>
    <div class="source-distribution-controls">
      <form v-if="!isDetail" class="distribution-gene-form" @submit.prevent="loadDistribution"><label for="distribution-gene">Gene for PI distributions</label><div class="distribution-gene-input"><Search :size="17" aria-hidden="true"/><input id="distribution-gene" v-model="query" aria-label="Gene for PI distributions" placeholder="Gene or RefSeq transcript"/><button type="submit" class="button primary">Show PI</button></div></form>
      <label v-if="!isDetail && loadedGene" class="distribution-transcript"><span>Transcript</span><select v-model="selectedTranscript" aria-label="PI distribution transcript"><option v-for="t in loadedGene.transcripts" :key="t.id" :value="t.id">{{ t.id }}</option></select></label>
      <p v-if="!loading && !error && symbol" class="source-distribution-status" role="status">{{ measured.length }} of {{ matchingRows.length }} filtered profiles have measured PI<span v-if="matchingRows.length > measured.length"> · {{ matchingRows.length - measured.length }} missing values excluded</span>. All three charts use all matching profiles across table pages.</p>
    </div>
    <GeneGenomeBrowser v-if="symbol && !loading && !error" :symbol="symbol" :transcript="transcript" :profiles="matchingRows"/>
    <div v-if="loading" class="distribution-loading" role="status" aria-busy="true"><div/><p>Loading gene PI distributions…</p></div>
    <div v-else-if="error" class="error-state" role="alert"><strong>Unable to load PI distributions</strong><p>{{ error }}</p><button class="button secondary" type="button" @click="loadDistribution">Retry PI</button></div>
    <GenePiGroups v-else-if="symbol" :symbol="symbol" :transcript="transcript" :rows="measured" :all-rows="allMetrics" :profile-count="catalog.profiles.length" @reset="reset"/>
    <template v-if="symbol && !loading && !error">
      <GeneFeatures v-if="isDetail" :symbol="symbol" :transcript="transcript"/>
      <GeneSequence class="source-sequence-row" :symbol="symbol" :transcript="transcript"/>
    </template>
  </div>
</template>
<style scoped>
.source-sequence-row{margin-top:24px}.source-sequence-row :deep(.sequence-error),.source-sequence-row :deep(.sequence-loading){min-height:0}
.source-pi-explorer{min-width:0}.source-distribution-controls{display:flex;align-items:flex-end;flex-wrap:wrap;gap:14px 20px;margin:8px 0 0}.distribution-gene-form{flex:1 1 330px;max-width:570px}.distribution-gene-form>label,.distribution-transcript>span{display:block;font-size:13px;font-weight:500;color:#536b83;margin-bottom:8px}.distribution-gene-input{border:1px solid #d5e0eb;background:white;display:flex;align-items:center;gap:10px;padding:5px 6px 5px 13px;border-radius:7px;min-height:50px}.distribution-gene-input>svg{flex-shrink:0;color:#7b92a9}.distribution-gene-input>input{flex:1;min-width:0;font-size:14px;color:#334e68;height:36px}.distribution-gene-input>.button{min-height:38px;font-size:13px}.distribution-transcript select{min-height:50px;border:1px solid #d5e0eb;background:white;border-radius:7px;color:#334e68;font-size:13px;padding:0 12px;max-width:100%}.source-distribution-status{flex-basis:100%;font-size:13px;line-height:1.8;color:#64748b;font-variant-numeric:tabular-nums;text-wrap:pretty}.distribution-loading{border:1px solid #e4eaf1;border-radius:10px;padding:24px;margin:20px 0;color:#64748b;font-size:13px}.distribution-loading>div{height:150px;border-radius:6px;background:#f1f5f9;margin-bottom:14px}.source-pi-explorer :deep(.gene-pi-groups){margin-top:18px}@media(max-width:760px){.distribution-gene-form{flex-basis:100%}.distribution-transcript{width:100%}.distribution-transcript select{width:100%}}
</style>
