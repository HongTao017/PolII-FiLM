<script setup>
import { sampleNameColumn } from '../utils/sampleName.js'
import { computed,ref,watch } from 'vue'
import { useRoute,useRouter } from 'vue-router'
import { Search,Download,Dna } from '@lucide/vue'
import { useCatalog } from '../composables/catalog.js'
import { getGene,getGeneIndex } from '../api/client.js'
import { assays,filterRecords,exportTsv } from '../utils/data.js'
import PageTitle from '../components/PageTitle.vue'
import ActiveFilters from '../components/ActiveFilters.vue'
import LoadState from '../components/LoadState.vue'
import EmptyState from '../components/EmptyState.vue'
import SourcePiExplorer from '../components/SourcePiExplorer.vue'
import GeneCatalog from '../components/GeneCatalog.vue'
const {catalog,loading:catalogLoading,error:catalogError,reload}=useCatalog(),route=useRoute(),router=useRouter()
const gene=ref(null),loading=ref(false),error=ref(''),input=ref(''),transcript=ref('')
const assay=ref(route.query.assay||''),cell=ref(route.query.cell||''),antibody=ref('')
let requestId=0
async function load(symbol){const id=++requestId;gene.value=null;error.value='';transcript.value='';input.value=symbol||'';if(!symbol){loading.value=false;return}loading.value=true;try{const value=await getGene(symbol);if(id===requestId){gene.value=value;transcript.value=value?.transcripts[0]?.id||''}}catch(e){if(id===requestId)error.value=e.message}finally{if(id===requestId)loading.value=false}}
watch(()=>route.params.symbol,load,{immediate:true})
watch(()=>route.query,v=>{assay.value=v.assay||'';cell.value=v.cell||''})
async function search(value=input.value){if(!value.trim())return;try{const index=await getGeneIndex();const q=value.trim().toLowerCase();const entry=index.find(g=>g.symbol.toLowerCase()===q||g.transcripts.some(t=>t.toLowerCase()===q));router.push({path:`/genes/${encodeURIComponent(entry?.symbol||value.trim())}`,query:{...(assay.value?{assay:assay.value}:{}),...(cell.value?{cell:cell.value}:{})}})}catch(e){error.value=e.message}}
const current=computed(()=>gene.value?.transcripts.find(t=>t.id===transcript.value))
const metrics=computed(()=>current.value&&catalog.value?catalog.value.metric.profile_order.map((id,index)=>{const p=catalog.value.profiles.find(profile=>profile.id===id);return p?{...p,log2_pi:current.value.log2_pi[index]}:null}).filter(p=>p&&p.log2_pi!==null):[])
const filtered=computed(()=>filterRecords(metrics.value,{assay:assay.value,cell:cell.value,antibody:antibody.value}))
const antibodies=computed(()=>[...new Set(metrics.value.map(p=>p.antibody).filter(Boolean))].sort())
const coordinate=computed(()=>current.value?.coordinates||[])
function download(){exportTsv(`${gene.value.symbol}_${transcript.value}_log2PI.tsv`,[{key:'id',label:'profile_id'},sampleNameColumn,{key:'alias',label:'source_alias'},{key:'assay',label:'assay'},{key:'cell',label:'raw_cell_label'},{key:'author',label:'author_source_label'},{key:'gse',label:'GSE'},{key:'antibody',label:'antibody_from_alias'},{key:'log2_pi',label:'log2_pause_index'}],filtered.value)}
const activeFilters=computed(()=>[{key:'assay',label:'Assay',value:assay.value},{key:'cell',label:'Cell',value:cell.value},{key:'antibody',label:'Antibody',value:antibody.value}].filter(item=>item.value))
function clearFilter(key){const filters={assay,cell,antibody};if(filters[key])filters[key].value=''}
function resetFilters(){assay.value='';cell.value='';antibody.value=''}
</script>
<template><GeneCatalog v-if="!route.params.symbol"/><template v-else><PageTitle eyebrow="GENE EXPLORER" :title="gene?.symbol||'Gene catalog'" :description="route.params.symbol ? 'Explore measured pausing metrics across assays, cellular contexts and transcripts.' : 'Browse genes and RefSeq transcripts, then select a gene to explore its data.'"><button v-if="gene" class="button secondary" :disabled="!filtered.length" @click="download"><Download :size="16"/>Download gene results</button></PageTitle><form class="hero-search search-page-form" @submit.prevent="search()"><Search :size="21"/><input v-model="input" aria-label="Search a gene or RefSeq transcript" placeholder="Enter a gene or RefSeq transcript, e.g. NOC2L"/><button type="submit">View gene <Dna :size="17"/></button></form><LoadState :loading="loading||catalogLoading" :error="error||catalogError" @retry="()=>{reload();load(route.params.symbol)}"><template v-if="gene&&catalog"><div class="notice">Values come from the existing log2 PI matrices. Assembly, annotation and calculation protocols await verification; cross-assay views are exploratory.</div><div class="filter-bar gene-filters"><label><span>Transcript</span><select v-model="transcript" aria-label="Select transcript"><option v-for="t in gene.transcripts" :key="t.id" :value="t.id">{{ t.id }}</option></select></label><label><span>Assay</span><select v-model="assay" aria-label="Filter assay"><option value="">All assays</option><option v-for="a in assays" :key="a">{{ a }}</option></select></label><label><span>Cell</span><select v-model="cell" aria-label="Filter cell label"><option value="">All cell labels</option><option v-for="c in catalog.cells" :key="c">{{ c }}</option></select></label><label><span>Antibody</span><select v-model="antibody" aria-label="Filter antibody subtype"><option value="">All subtypes</option><option v-for="a in antibodies" :key="a">{{ a }}</option></select></label></div><ActiveFilters :items="activeFilters" @remove="clearFilter"/><SourcePiExplorer :catalog="catalog" :rows="filtered" :all-rows="metrics" :symbol="gene.symbol" :transcript="transcript" @reset="resetFilters"/><section class="panel latest-panel"><div class="panel-head"><h2>Coordinates from source matrices</h2><span class="small">BED0 · Reference version unverified</span></div><div v-if="coordinate.length" class="table-scroll"><table><thead><tr><th>Source assay</th><th>Chromosome</th><th>Strand</th><th>TSS (BED0)</th><th>TES (BED0)</th></tr></thead><tbody><tr v-for="c in coordinate" :key="c.source_assay"><td>{{ c.source_assay }}</td><td>{{ c.chrom }}</td><td>{{ c.strand }}</td><td class="mono">{{ c.tss_bed0 }}</td><td class="mono">{{ c.tes_bed0 }}</td></tr></tbody></table></div><EmptyState v-else title="No coordinates in this source matrix" description="Connect transcript annotations to provide coordinates."/></section></template><EmptyState v-else-if="route.params.symbol" title="Gene or transcript not found" description="Try NOC2L, TP53 or an exact RefSeq accession."/></LoadState></template></template>
