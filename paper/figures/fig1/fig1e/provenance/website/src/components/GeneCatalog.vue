<script setup>
import { computed, onBeforeUnmount, ref, watch } from 'vue'
import { Download, ArrowDown, ArrowUp, Search, RotateCcw } from '@lucide/vue'
import { useRoute, useRouter } from 'vue-router'
import { getGeneCatalog } from '../api/client.js'
import { assays, exportTsv } from '../utils/data.js'
import AssayBadge from './AssayBadge.vue'
import Pagination from './Pagination.vue'
import EmptyState from './EmptyState.vue'
import ActiveFilters from './ActiveFilters.vue'
const route=useRoute(), router=useRouter()
const parse=v=>[...new Set((Array.isArray(v)?v:[v||'']).flatMap(x=>x.split(',')).filter(Boolean))]
const query=ref(String(route.query.q||'')), selectedAssays=ref(parse(route.query.assay)), strand=ref(String(route.query.strand||'')), minimum=ref([0,10,50,100].includes(Number(route.query.min))?Number(route.query.min):0)
const data = ref(null), loading = ref(false), error = ref(''), chromosome = ref(String(route.query.chrom||'')), page = ref(1), size = ref(10), sortKey = ref('symbol'), direction = ref('asc')
let revision = 0
async function load() {
  const id = ++revision; loading.value = true; error.value = ''
  try {
    const result = await getGeneCatalog()
    if (result.schema_version !== 1 || result.rows?.length !== result.gene_count) throw Error('Invalid gene catalog. Please retry.')
    if (id === revision) data.value = result
  } catch (e) { if (id === revision) error.value = e.message }
  finally { if (id === revision) loading.value = false }
}
load(); onBeforeUnmount(() => { ++revision })
const chromosomes = computed(() => [...new Set(data.value?.rows.map(row => row.chrom).filter(Boolean))].sort((a,b) => a.localeCompare(b, 'en', { numeric: true })))
const filtered = computed(() => {
  const q = query.value.trim().toLowerCase()
  return (data.value?.rows || []).filter(row => (!chromosome.value || row.chrom === chromosome.value) && (!selectedAssays.value.length || row.assays.some(a=>selectedAssays.value.includes(a))) && (!strand.value || row.strand===strand.value) && row.observed>=minimum.value && (!q || [row.symbol, ...row.transcripts].some(v => v.toLowerCase().includes(q))))
})
const sorted = computed(() => [...filtered.value].sort((a,b) => {
  const av = a[sortKey.value], bv = b[sortKey.value]
  if (av === null) return bv === null ? a.symbol.localeCompare(b.symbol) : 1
  if (bv === null) return -1
  const value = typeof av === 'number' ? av - bv : av.localeCompare(bv, 'en', { numeric: true })
  return (direction.value === 'asc' ? value : -value) || a.symbol.localeCompare(b.symbol)
}))
const shown = computed(() => sorted.value.slice((page.value-1)*size.value, page.value*size.value))
watch(()=>route.query,v=>{query.value=String(v.q||'');selectedAssays.value=parse(v.assay);chromosome.value=String(v.chrom||'');strand.value=String(v.strand||'');minimum.value=[0,10,50,100].includes(Number(v.min))?Number(v.min):0})
watch([query,selectedAssays,chromosome,strand,minimum],()=>{page.value=1;const values={q:query.value,assay:selectedAssays.value.join(','),chrom:chromosome.value,strand:strand.value,min:minimum.value?String(minimum.value):''};router.replace({query:Object.fromEntries(Object.entries(values).filter(([,v])=>v))})})
watch([size,sortKey,direction],()=>{page.value=1})
const assayCounts=computed(()=>Object.fromEntries(assays.map(a=>[a,data.value?.rows.filter(r=>r.assays.includes(a)).length||0])))
const activeFilters=computed(()=>[{key:'q',label:'Search',value:query.value},...selectedAssays.value.map(a=>({key:'assay:'+a,label:'Assay',value:a})),{key:'chrom',label:'Chromosome',value:chromosome.value},{key:'strand',label:'Strand',value:strand.value},{key:'min',label:'PI profiles',value:minimum.value?'≥ '+minimum.value:''}].filter(v=>v.value))
function clearFilter(key){if(key.startsWith('assay:'))selectedAssays.value=selectedAssays.value.filter(a=>a!==key.slice(6));else if(key==='min')minimum.value=0;else{const values={q:query,chrom:chromosome,strand};if(values[key])values[key].value=''}}
function sort(key) { if (sortKey.value === key) direction.value = direction.value === 'asc' ? 'desc' : 'asc'; else { sortKey.value = key; direction.value = key === 'symbol' ? 'asc' : 'desc' } }
function reset(){query.value='';selectedAssays.value=[];chromosome.value='';strand.value='';minimum.value=0}
function download() {
  exportTsv('NascentAtlas_genes_filtered.tsv', [{key:'symbol',label:'gene'}, {key:'transcript',label:'summary_transcript'}, {label:'indexed_transcripts',value:r=>r.transcripts.join(';')}, {key:'chrom',label:'source_chromosome'}, {key:'strand',label:'source_strand'}, {key:'observed',label:'measured_PI_profiles'}, {key:'median',label:'median_log2_PI'}, {label:'assays',value:r=>r.assays.join(';')}, {key:'cells',label:'cell_labels_with_PI'}, {key:'studies',label:'GEO_studies_with_PI'}], sorted.value)
}
</script>
<template>
 <div class="browse-page genes-browse-page">
  <header class="browse-banner"><h1>Genes</h1><p>Human genes and RefSeq transcripts</p></header>
  <div class="container browse-content">
   <p v-if="loading" class="gene-list-message" role="status">Loading gene catalog…</p>
   <div v-else-if="error" class="gene-list-message" role="alert"><p>{{error}}</p><button class="button secondary" @click="load">Retry gene catalog</button></div>
   <div v-else-if="data" class="browse-layout">
    <aside class="panel browse-sidebar" aria-label="Gene filters">
     <div class="browse-filter-top"><h2>Search &amp; filter</h2><button type="button" class="browse-reset" @click="reset"><RotateCcw :size="13" aria-hidden="true"/>Reset</button></div>
     <label class="browse-search"><span class="sr-only">Search a gene or RefSeq transcript</span><input v-model="query" aria-label="Search a gene or RefSeq transcript" placeholder="Gene or RefSeq transcript"/><Search :size="18" aria-hidden="true"/></label>
     <fieldset class="browse-facet"><legend>Sequencing assay</legend><label v-for="a in assays" :key="a" class="browse-choice"><input v-model="selectedAssays" type="checkbox" :value="a" :aria-label="`Assay ${a}`"/><span>{{a}}</span><small>{{assayCounts[a].toLocaleString('en-US')}}</small></label></fieldset>
     <label class="browse-select"><span>Chromosome</span><select v-model="chromosome" aria-label="Gene chromosome"><option value="">All chromosomes</option><option v-for="c in chromosomes" :key="c">{{c}}</option></select></label>
     <label class="browse-select"><span>Strand</span><select v-model="strand" aria-label="Gene strand"><option value="">All strands</option><option value="+">+ strand</option><option value="-">− strand</option></select></label>
     <label class="browse-select"><span>Measured PI profiles</span><select v-model.number="minimum" aria-label="Minimum measured PI profiles"><option :value="0">All genes</option><option :value="10">At least 10 profiles</option><option :value="50">At least 50 profiles</option><option :value="100">At least 100 profiles</option></select></label>
    </aside>
    <section class="browse-results gene-catalog" aria-label="Gene catalog">
     <div class="browse-results-top"><div class="browse-count"><div><strong>{{filtered.length.toLocaleString('en-US')}}</strong> genes found</div><span>{{data.gene_count.toLocaleString('en-US')}} genes · {{data.transcript_count.toLocaleString('en-US')}} indexed transcripts</span></div><button class="button secondary" :disabled="!filtered.length" @click="download"><Download :size="15"/>Export filtered records</button></div>
     <ActiveFilters :items="activeFilters" @remove="clearFilter"/>
     <div class="browse-toolbar"><span class="gene-list-label">Gene list · one gene per row</span><RouterLink class="table-link" to="/samples">Browse sequencing profiles →</RouterLink></div>
     <section class="panel browse-table-panel">
      <div v-if="shown.length" class="table-scroll gene-catalog-scroll" tabindex="0" role="region" aria-label="Gene list table, scroll horizontally for all columns"><table><thead><tr><th scope="col" :aria-sort="sortKey==='symbol' ? (direction==='asc'?'ascending':'descending') : 'none'"><button @click="sort('symbol')">Gene<component :is="direction==='asc'?ArrowUp:ArrowDown" v-if="sortKey==='symbol'" :size="12"/></button></th><th scope="col">RefSeq transcript</th><th scope="col">Chromosome</th><th scope="col">Strand</th><th scope="col" :aria-sort="sortKey==='observed' ? (direction==='asc'?'ascending':'descending') : 'none'"><button @click="sort('observed')">Measured PI profiles<component :is="direction==='asc'?ArrowUp:ArrowDown" v-if="sortKey==='observed'" :size="12"/></button></th><th scope="col" :aria-sort="sortKey==='median' ? (direction==='asc'?'ascending':'descending') : 'none'"><button @click="sort('median')">Median log₂ PI<component :is="direction==='asc'?ArrowUp:ArrowDown" v-if="sortKey==='median'" :size="12"/></button></th><th scope="col">Sequencing assays</th><th scope="col">Details</th></tr></thead><tbody><tr v-for="row in shown" :key="row.symbol" :data-gene="row.symbol"><th scope="row"><RouterLink class="table-link" :to="`/genes/${encodeURIComponent(row.symbol)}`">{{row.symbol}}</RouterLink></th><td>{{row.transcript}}<small v-if="row.transcripts.length>1">{{row.transcripts.length}} indexed transcripts</small></td><td>{{row.chrom || 'Not recorded'}}</td><td>{{row.strand || '—'}}</td><td class="gene-number" :data-observed="row.observed">{{row.observed}} / {{data.profile_count}}</td><td class="gene-number" :data-median="row.median ?? ''">{{row.median===null ? 'Not measured' : row.median.toFixed(2)}}</td><td><div class="gene-assays"><AssayBadge v-for="assay in row.assays" :key="assay" :assay="assay"/><span v-if="!row.assays.length">Not measured</span></div></td><td><RouterLink class="table-link" :to="`/genes/${encodeURIComponent(row.symbol)}`" :aria-label="`View ${row.symbol} gene details`">View ↗</RouterLink></td></tr></tbody></table></div>
      <EmptyState v-else title="No matching genes" description="Try adjusting the gene, chromosome, assay or PI coverage filters."><button class="button secondary" @click="reset">Reset filters</button></EmptyState>
      <div class="browse-pagination"><label>Rows per page<select v-model.number="size" aria-label="Gene rows per page"><option :value="10">10</option><option :value="25">25</option><option :value="50">50</option></select></label><Pagination v-model="page" :total="filtered.length" :size="size"/></div>
     </section>
     <p class="browse-count-note">Assay filters select genes with measured PI for those assays. Summaries use all measured profiles for the listed transcript; missing PI values are excluded. Chromosome and strand come from the source matrices.</p>
    </section>
   </div>
  </div>
 </div>
</template>
<style scoped src="../browse-layout.css"></style>
<style scoped>
.genes-browse-page .browse-select{align-content:start}.gene-catalog{min-width:0}.gene-list-label{font-size:12px;color:#536b83}.gene-catalog table{min-width:920px;width:100%;font-size:12px}.gene-catalog td,.gene-catalog th{padding:16px 14px}.gene-catalog thead th{background:#3e77ab;color:white;font-weight:600}.gene-catalog thead button{display:inline-flex;align-items:center;gap:6px;text-align:left;font:inherit;background:none;color:inherit;padding:0;min-height:24px}.gene-catalog th:first-child{position:sticky;left:0;min-width:100px;z-index:1}.gene-catalog thead th:first-child{background:#3e77ab;z-index:2}.gene-catalog tbody th{font-weight:500;text-align:left;border-bottom:0}.gene-catalog tbody td{border-bottom:0;color:#52677e}.gene-catalog tbody tr:nth-child(odd),.gene-catalog tbody tr:nth-child(odd) th{background:#f5f7f9}.gene-catalog tbody tr:nth-child(even),.gene-catalog tbody tr:nth-child(even) th{background:white}.gene-catalog tbody tr:hover td,.gene-catalog tbody tr:hover th{background:#eaf2f9}.gene-catalog .table-link{color:#3e77ab;font-size:12px}.gene-number{font-variant-numeric:tabular-nums;white-space:nowrap}.gene-assays{display:flex;gap:5px;flex-wrap:wrap;max-width:160px}.gene-catalog td small{display:block;font-size:10px;color:#64748b}.gene-catalog :deep(.assay-badge){white-space:nowrap}.gene-list-message{padding:20px;color:#64748b;font-size:13px}.gene-catalog :focus-visible{outline:2px solid #417dc6;outline-offset:2px}@media(max-width:760px){.gene-catalog th:first-child{min-width:95px}.gene-catalog td,.gene-catalog th{padding:14px 12px}}
</style>
