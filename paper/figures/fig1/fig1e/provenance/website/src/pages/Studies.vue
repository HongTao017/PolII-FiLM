<script setup>
import { computed, ref, watch } from 'vue'
import { useRoute, useRouter } from 'vue-router'
import { Search, Download, RotateCcw } from '@lucide/vue'
import { useCatalog } from '../composables/catalog.js'
import { assays, exportTsv } from '../utils/data.js'
import ActiveFilters from '../components/ActiveFilters.vue'
import LoadState from '../components/LoadState.vue'
import StudyTable from '../components/StudyTable.vue'
import Pagination from '../components/Pagination.vue'
import EmptyState from '../components/EmptyState.vue'
const { catalog, loading, error, reload }=useCatalog()
const route=useRoute(),router=useRouter(),page=ref(1),pageSize=ref(10),cellSearch=ref(''),sortKey=ref(''),sortDirection=ref('asc')
const parse=v=>[...new Set((Array.isArray(v)?v:[v||'']).flatMap(x=>x.split(',')).filter(Boolean))]
const q=ref(route.query.q||''),selectedAssays=ref(parse(route.query.assay)),selectedCells=ref(parse(route.query.cell)),author=ref(route.query.author||''),study=ref(route.query.study||'')
watch(()=>route.query,v=>{q.value=v.q||'';selectedAssays.value=parse(v.assay);selectedCells.value=parse(v.cell);author.value=v.author||'';study.value=v.study||''})
watch([q,selectedAssays,selectedCells,author,study],()=>{page.value=1;const values={q:q.value,assay:selectedAssays.value.join(','),cell:selectedCells.value.join(','),author:author.value,study:study.value};router.replace({query:Object.fromEntries(Object.entries(values).filter(([,v])=>v))})})
watch([pageSize,sortKey,sortDirection],()=>{page.value=1})
const studies=computed(()=>catalog.value?.studies||[])
const filtered=computed(()=>{const query=q.value.trim().toLowerCase();return studies.value.filter(s=>(!query||[s.gse,...s.authors,...s.cells,...s.assays].join(' ').toLowerCase().includes(query))&&(!selectedAssays.value.length||s.assays.some(a=>selectedAssays.value.includes(a)))&&(!selectedCells.value.length||s.cells.some(c=>selectedCells.value.includes(c)))&&(!author.value||s.authors.includes(author.value))&&(!study.value||s.gse===study.value))})
const sorted=computed(()=>sortKey.value?[...filtered.value].sort((a,b)=>{const key=sortKey.value,av=Array.isArray(a[key])?a[key].join(' / '):a[key],bv=Array.isArray(b[key])?b[key].join(' / '):b[key];return(typeof av==='number'?av-bv:String(av||'').localeCompare(String(bv||''),undefined,{numeric:true}))*(sortDirection.value==='asc'?1:-1)}):filtered.value)
const rows=computed(()=>sorted.value.slice((page.value-1)*pageSize.value,page.value*pageSize.value))
const profileCount=computed(()=>filtered.value.reduce((total,s)=>total+s.profile_count,0))
const assayCounts=computed(()=>Object.fromEntries(assays.map(a=>[a,studies.value.filter(s=>s.assays.includes(a)).length])))
const cellOptions=computed(()=>{const counts=new Map();for(const s of studies.value)for(const c of new Set(s.cells))counts.set(c,(counts.get(c)||0)+1);return [...counts].map(([name,count])=>({name,count})).sort((a,b)=>b.count-a.count||a.name.localeCompare(b.name)).filter(c=>c.name.toLowerCase().includes(cellSearch.value.toLowerCase()))})
const authors=computed(()=>[...new Set(studies.value.flatMap(s=>s.authors))].sort())
const accessions=computed(()=>studies.value.map(s=>s.gse).sort((a,b)=>a.localeCompare(b,undefined,{numeric:true})))
function reset(){q.value='';selectedAssays.value=[];selectedCells.value=[];author.value='';study.value='';cellSearch.value=''}
const activeFilters=computed(()=>[{key:'q',label:'Search',value:q.value},...selectedAssays.value.map(a=>({key:'assay:'+a,label:'Assay',value:a})),...selectedCells.value.map(c=>({key:'cell:'+c,label:'Cell',value:c})),{key:'author',label:'Author',value:author.value},{key:'study',label:'Study',value:study.value}].filter(item=>item.value))
function clearFilter(key){if(key.startsWith('assay:'))selectedAssays.value=selectedAssays.value.filter(a=>a!==key.slice(6));else if(key.startsWith('cell:'))selectedCells.value=selectedCells.value.filter(c=>c!==key.slice(5));else{const values={q,author,study};if(values[key])values[key].value=''}}
function sort(key){if(sortKey.value===key)sortDirection.value=sortDirection.value==='asc'?'desc':'asc';else{sortKey.value=key;sortDirection.value='asc'}}
function download(){exportTsv('NascentAtlas_studies_filtered.tsv',[{key:'gse',label:'GSE'},{label:'authors',value:s=>s.authors.join(';')},{label:'assays',value:s=>s.assays.join(';')},{label:'raw_cell_labels',value:s=>s.cells.join(';')},{key:'profile_count',label:'profile_records'},{key:'source_url',label:'source_url'},{key:'status',label:'status'}],sorted.value)}
</script>
<template>
 <div class="browse-page studies-browse-page">
  <header class="browse-banner"><h1>Studies</h1><p>Human nascent transcription GEO studies</p></header>
  <div class="container browse-content"><LoadState :loading="loading" :error="error" @retry="reload"><div v-if="catalog" class="browse-layout">
   <aside class="panel browse-sidebar" aria-label="Study filters">
    <div class="browse-filter-top"><h2>Search &amp; filter</h2><button type="button" class="browse-reset" @click="reset"><RotateCcw :size="13" aria-hidden="true"/>Reset</button></div>
    <label class="browse-search"><span class="sr-only">Search studies</span><input v-model="q" aria-label="Search studies" placeholder="GSE, author or cell line"/><Search :size="18" aria-hidden="true"/></label>
    <fieldset class="browse-facet"><legend>Sequencing assay</legend><label v-for="a in assays" :key="a" class="browse-choice"><input v-model="selectedAssays" type="checkbox" :value="a" :aria-label="`Assay ${a}`"/><span>{{a}}</span><small>{{assayCounts[a]}}</small></label></fieldset>
    <fieldset class="browse-facet"><legend>Cell line</legend><input v-model="cellSearch" class="browse-cell-search" aria-label="Find a cell line" placeholder="Find a cell line…"/><div class="browse-cell-options"><label v-for="c in cellOptions" :key="c.name" class="browse-choice"><input v-model="selectedCells" type="checkbox" :value="c.name" :aria-label="`Cell line ${c.name}`"/><span>{{c.name}}</span><small>{{c.count}}</small></label><p v-if="!cellOptions.length" class="browse-no-cells">No cell labels found</p></div></fieldset>
    <label class="browse-select"><span>Author / source</span><select v-model="author" aria-label="Study author"><option value="">All authors</option><option v-for="a in authors" :key="a">{{a}}</option></select></label>
    <label class="browse-select"><span>GEO study</span><select v-model="study" aria-label="Filter GEO study"><option value="">All studies</option><option v-for="g in accessions" :key="g">{{g}}</option></select></label>
   </aside>
   <section class="browse-results" aria-label="Study catalog">
    <div class="browse-results-top"><div class="browse-count"><div><strong>{{filtered.length}}</strong> studies found</div><span>{{profileCount}} catalog profiles</span></div><button class="button secondary" :disabled="!filtered.length" @click="download"><Download :size="15"/>Export filtered records</button></div>
    <ActiveFilters :items="activeFilters" @remove="clearFilter"/>
    <div class="browse-toolbar"><span class="study-list-label">GEO study list</span><RouterLink class="table-link" to="/samples">Browse sequencing profiles →</RouterLink></div>
    <section class="panel browse-table-panel"><StudyTable v-if="filtered.length" variant="browse" sortable :rows="rows" :sort-key="sortKey" :sort-direction="sortDirection" @sort="sort"/><EmptyState v-else title="No matching studies" description="Try adjusting the assay, cell line, author or GSE filters."/><div class="browse-pagination"><label>Rows per page<select v-model.number="pageSize" aria-label="Rows per page"><option :value="10">10</option><option :value="25">25</option><option :value="50">50</option></select></label><Pagination v-model="page" :total="filtered.length" :size="pageSize"/></div></section>
    <p class="browse-count-note">Facet counts represent studies. Profile counts include all catalog records in each matching study.</p>
   </section>
  </div></LoadState></div>
 </div>
</template>
<style scoped src="../browse-layout.css"></style>
<style scoped>.study-list-label{font-size:12px;color:#536b83}</style>
