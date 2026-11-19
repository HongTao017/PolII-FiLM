<script setup>
import { sampleName } from '../utils/sampleName.js'
import { computed, ref, watch } from 'vue'
import { useRoute, useRouter } from 'vue-router'
import { Search, Download, RotateCcw } from '@lucide/vue'
import { useCatalog } from '../composables/catalog.js'
import { assays, filterRecords, exportTsv, sampleColumns } from '../utils/data.js'
import ActiveFilters from '../components/ActiveFilters.vue'
import LoadState from '../components/LoadState.vue'
import SampleTable from '../components/SampleTable.vue'
import LibraryQcTable from '../components/LibraryQcTable.vue'
import Pagination from '../components/Pagination.vue'
import EmptyState from '../components/EmptyState.vue'
const { catalog, loading, error, reload }=useCatalog()
const route=useRoute(),router=useRouter(),page=ref(1),pageSize=ref(10),sortKey=ref(''),sortDirection=ref('asc'),cellSearch=ref('')
const parse=v=>[...new Set((Array.isArray(v)?v:[v||'']).flatMap(x=>x.split(',')).filter(Boolean))]
const view=ref(['study','qc'].includes(route.query.view)?route.query.view:'list')
const q=ref(route.query.q||''),selectedAssays=ref(parse(route.query.assay)),selectedCells=ref(parse(route.query.cell)),antibody=ref(route.query.antibody||''),study=ref(route.query.study||'')
watch(()=>route.query,v=>{q.value=v.q||'';selectedAssays.value=parse(v.assay);selectedCells.value=parse(v.cell);antibody.value=v.antibody||'';study.value=v.study||'';view.value=['study','qc'].includes(v.view)?v.view:'list'})
watch([q,selectedAssays,selectedCells,antibody,study,view],()=>{page.value=1;const values={q:q.value,assay:selectedAssays.value.join(','),cell:selectedCells.value.join(','),antibody:antibody.value,study:study.value,view:view.value==='list'?'':view.value};router.replace({query:Object.fromEntries(Object.entries(values).filter(([,v])=>v))})})
watch([pageSize,sortKey,sortDirection],()=>{page.value=1})
const filtered=computed(()=>filterRecords(catalog.value?.profiles||[],{q:q.value,antibody:antibody.value,study:study.value}).filter(p=>(!selectedAssays.value.length||selectedAssays.value.includes(p.assay))&&(!selectedCells.value.length||selectedCells.value.includes(p.cell))))
const sorted=computed(()=>sortKey.value?[...filtered.value].sort((a,b)=>String(sortKey.value==='id'?sampleName(a):a[sortKey.value]||'').localeCompare(String(sortKey.value==='id'?sampleName(b):b[sortKey.value]||''),undefined,{numeric:true})*(sortDirection.value==='asc'?1:-1)):filtered.value)
const rows=computed(()=>sorted.value.slice((page.value-1)*pageSize.value,page.value*pageSize.value))
const studies=computed(()=>[...new Set(catalog.value?.profiles.map(p=>p.gse)||[])].sort((a,b)=>a.localeCompare(b,undefined,{numeric:true})))
const groups=computed(()=>[...new Set(filtered.value.map(p=>p.gse))].sort((a,b)=>a.localeCompare(b,undefined,{numeric:true})).map(gse=>({gse,profiles:filtered.value.filter(p=>p.gse===gse)})))
const shownGroups=computed(()=>groups.value.slice((page.value-1)*pageSize.value,page.value*pageSize.value))
const antibodies=computed(()=>[...new Set(catalog.value?.profiles.map(p=>p.antibody).filter(Boolean))].sort())
const assayCounts=computed(()=>Object.fromEntries(assays.map(a=>[a,catalog.value?.profiles.filter(p=>p.assay===a).length||0])))
const cellOptions=computed(()=>{const counts=new Map();for(const p of catalog.value?.profiles||[])counts.set(p.cell,(counts.get(p.cell)||0)+1);return [...counts].map(([name,count])=>({name,count})).sort((a,b)=>b.count-a.count||a.name.localeCompare(b.name)).filter(c=>c.name.toLowerCase().includes(cellSearch.value.toLowerCase()))})
function reset(){q.value='';selectedAssays.value=[];selectedCells.value=[];antibody.value='';study.value='';cellSearch.value=''}
const activeFilters=computed(()=>[{key:'q',label:'Search',value:q.value},...selectedAssays.value.map(a=>({key:'assay:'+a,label:'Assay',value:a})),...selectedCells.value.map(c=>({key:'cell:'+c,label:'Cell',value:c})),{key:'antibody',label:'Antibody',value:antibody.value},{key:'study',label:'Study',value:study.value}].filter(item=>item.value))
function clearFilter(key){if(key.startsWith('assay:'))selectedAssays.value=selectedAssays.value.filter(a=>a!==key.slice(6));else if(key.startsWith('cell:'))selectedCells.value=selectedCells.value.filter(c=>c!==key.slice(5));else{const filters={q,antibody,study};if(filters[key])filters[key].value=''}}
function sort(key){if(sortKey.value===key)sortDirection.value=sortDirection.value==='asc'?'desc':'asc';else{sortKey.value=key;sortDirection.value='asc'}}
</script>
<template>
 <div class="browse-page">
  <header class="browse-banner"><h1>Browse</h1><p>Human nascent transcription profiles</p></header>
  <div class="container browse-content"><LoadState :loading="loading" :error="error" @retry="reload"><div v-if="catalog" class="browse-layout">
   <aside class="panel browse-sidebar" aria-label="Profile filters">
    <div class="browse-filter-top"><h2>Search &amp; filter</h2><button type="button" class="browse-reset" @click="reset"><RotateCcw :size="13" aria-hidden="true"/>Reset</button></div>
    <label class="browse-search"><span class="sr-only">Search profiles</span><input v-model="q" aria-label="Search profiles" placeholder="Profile, author or GSE"/><Search :size="18" aria-hidden="true"/></label>
    <fieldset class="browse-facet"><legend>Sequencing assay</legend><label v-for="a in assays" :key="a" class="browse-choice"><input v-model="selectedAssays" type="checkbox" :value="a" :aria-label="`Assay ${a}`"/><span>{{a}}</span><small>{{assayCounts[a]}}</small></label></fieldset>
    <fieldset class="browse-facet"><legend>Cell line</legend><input v-model="cellSearch" class="browse-cell-search" aria-label="Find a cell line" placeholder="Find a cell line…"/><div class="browse-cell-options"><label v-for="c in cellOptions" :key="c.name" class="browse-choice"><input v-model="selectedCells" type="checkbox" :value="c.name" :aria-label="`Cell line ${c.name}`"/><span>{{c.name}}</span><small>{{c.count}}</small></label><p v-if="!cellOptions.length" class="browse-no-cells">No cell labels found</p></div></fieldset>
    <label class="browse-select"><span>Antibody subtype</span><select v-model="antibody" aria-label="Filter antibody subtype"><option value="">All subtypes</option><option v-for="a in antibodies" :key="a">{{a}}</option></select></label>
    <label class="browse-select"><span>GEO study</span><select v-model="study" aria-label="Filter GEO study"><option value="">All studies</option><option v-for="g in studies" :key="g">{{g}}</option></select></label>
   </aside>
   <section class="browse-results" aria-label="Profile catalog">
    <div class="browse-results-top"><div class="browse-count"><div><strong>{{filtered.length}}</strong> profiles found</div><span>{{groups.length}} GEO studies</span></div><button class="button secondary" :disabled="!filtered.length" @click="exportTsv('NascentAtlas_profiles_filtered.tsv',sampleColumns,filtered)"><Download :size="15"/>Export filtered records</button></div>
    <ActiveFilters :items="activeFilters" @remove="clearFilter"/>
    <div class="browse-toolbar"><div class="pi-view-switch" role="group" aria-label="Profile browsing mode"><button class="button secondary" :aria-pressed="view==='list'" @click="view='list'">Flat list</button><button class="button secondary" :aria-pressed="view==='study'" @click="view='study'">Group by GSE</button><button class="button secondary" :aria-pressed="view==='qc'" @click="view='qc'">Library QC</button></div><RouterLink class="table-link" to="/statistics#sample-correlations">Sample correlation heatmap →</RouterLink></div>
    <section v-if="view==='list'" class="panel browse-table-panel"><SampleTable v-if="filtered.length" variant="browse" sortable :rows="rows" :sort-key="sortKey" :sort-direction="sortDirection" @sort="sort"/><EmptyState v-else title="No matching profiles" description="Try adjusting the assay, cell or study filters."/><div class="browse-pagination"><label>Rows per page<select v-model.number="pageSize" aria-label="Rows per page"><option :value="10">10</option><option :value="25">25</option><option :value="50">50</option></select></label><Pagination v-model="page" :total="filtered.length" :size="pageSize"/></div></section>
    <LibraryQcTable v-else-if="view==='qc'" id="library-qc" :profiles="filtered"/>
    <div v-else><details v-for="group in shownGroups" :key="group.gse" class="panel pi-study-group browse-study-group" open><summary>{{group.gse}}<span>{{group.profiles.length}} profiles · {{[...new Set(group.profiles.map(p=>p.assay))].join(', ')}} · {{[...new Set(group.profiles.map(p=>p.cell))].join(', ')}}</span></summary><RouterLink class="table-link pi-study-link" :to="`/studies/${group.gse}`">View study provenance →</RouterLink><SampleTable variant="browse" :rows="group.profiles"/></details><EmptyState v-if="!filtered.length" title="No matching studies" description="Try adjusting the filters."/><div class="browse-pagination"><label>Studies per page<select v-model.number="pageSize" aria-label="Studies per page"><option :value="10">10</option><option :value="25">25</option><option :value="50">50</option></select></label><Pagination v-model="page" :total="groups.length" :size="pageSize"/></div></div>
    <p class="browse-count-note">Counts represent profile records in the current catalog.</p>
   </section>
  </div></LoadState></div>
 </div>
</template>
<style scoped src="../browse-layout.css"></style>
