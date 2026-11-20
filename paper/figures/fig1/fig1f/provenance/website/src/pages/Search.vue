<script setup>
import { computed, ref, watch } from 'vue'
import { useRoute, useRouter } from 'vue-router'
import { Search, ArrowUpRight } from '@lucide/vue'
import { useCatalog } from '../composables/catalog.js'
import { getGeneIndex } from '../api/client.js'
import { filterRecords } from '../utils/data.js'
import PageTitle from '../components/PageTitle.vue'
import LoadState from '../components/LoadState.vue'
import StudyTable from '../components/StudyTable.vue'
import SampleTable from '../components/SampleTable.vue'
import EmptyState from '../components/EmptyState.vue'
const {catalog,loading,error,reload}=useCatalog(),route=useRoute(),router=useRouter()
const input=ref(route.query.q||''),geneIndex=ref([]),geneError=ref('')
const query=computed(()=>String(route.query.q||'').trim())
watch(()=>route.query.q,value=>input.value=value||'')
getGeneIndex().then(v=>geneIndex.value=v).catch(e=>geneError.value=e.message)
const studies=computed(()=>query.value?filterRecords(catalog.value?.studies||[],{q:query.value}):[])
const profiles=computed(()=>query.value?filterRecords(catalog.value?.profiles||[],{q:query.value}):[])
const genes=computed(()=>query.value?geneIndex.value.filter(g=>[g.symbol,...g.transcripts].some(v=>v.toLowerCase().includes(query.value.toLowerCase()))).slice(0,24):[])
function submit(){router.push({path:'/search',query:input.value.trim()?{q:input.value.trim()}:{}})}
</script>
<template><PageTitle eyebrow="UNIFIED SEARCH" title="One search, connected data" description="Search genes, transcripts, GSE accessions, original profile aliases or cell labels."/><form class="hero-search search-page-form" @submit.prevent="submit"><Search :size="21"/><input v-model="input" aria-label="Unified search query" placeholder="Try TP53, NM_000546, GSE81662 or K562"/><button type="submit">Search <ArrowUpRight :size="17"/></button></form><LoadState :loading="loading" :error="error" @retry="reload"><template v-if="catalog"><EmptyState v-if="!query" title="Start with a keyword" description="Find pausing metrics by gene, trace a study by accession or explore coverage by cell."/><template v-else><p class="search-result-label">Results for “{{ query }}”</p><div v-if="geneError" class="notice">Unable to load the gene index: {{ geneError }}</div><section v-if="genes.length" class="panel search-section"><div class="panel-head"><h2>Genes / transcripts</h2><span class="small">Showing up to 24 matches</span></div><div class="gene-results"><RouterLink v-for="gene in genes" :key="gene.symbol" :to="`/genes/${encodeURIComponent(gene.symbol)}`"><span><strong>{{ gene.symbol }}</strong><small class="mono">{{ gene.transcripts.join(' / ') }}</small></span><ArrowUpRight :size="18"/></RouterLink></div></section><section v-if="studies.length" class="panel search-section"><div class="panel-head"><h2>Studies · {{ studies.length }}</h2><RouterLink :to="{path:'/samples',query:{q:query,view:'study'}}" class="text-link">View all →</RouterLink></div><StudyTable :rows="studies.slice(0,5)"/></section><section v-if="profiles.length" class="panel search-section"><div class="panel-head"><h2>Profiles · {{ profiles.length }}</h2><RouterLink :to="{path:'/samples',query:{q:query}}" class="text-link">View all →</RouterLink></div><SampleTable :rows="profiles.slice(0,5)"/></section><EmptyState v-if="!studies.length&&!profiles.length&&!genes.length&&!geneError" title="No matching records" description="GSM/SRR mappings are not yet available in this 105-profile catalog. Search by GSE, cell label or gene name."/></template></template></LoadState></template>
