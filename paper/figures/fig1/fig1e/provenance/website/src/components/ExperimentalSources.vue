<script setup>
import { sampleName, sampleNameColumn } from '../utils/sampleName.js'
import { computed, ref, watch } from 'vue'
import { Download, Search } from '@lucide/vue'
import AssayBadge from './AssayBadge.vue'
import Pagination from './Pagination.vue'
import EmptyState from './EmptyState.vue'
import { exportTsv, filterRecords } from '../utils/data.js'
const props = defineProps({ rows: { type: Array, required: true }, symbol: { type: String, default: '' }, transcript: { type: String, default: '' }, title: { type: String, default: 'Experimental sources' }, description: { type: String, default: 'Trace each profile to its sequencing assay, cell line and research source.' }, pageSize: { type: Number, default: 10 } })
const emit = defineEmits(['reset', 'filtered'])
const query = ref(''), assay = ref(''), cell = ref(''), author = ref(''), study = ref(''), page = ref(1)
const hasPi = computed(() => Boolean(props.symbol))
const options = key => [...new Set(props.rows.map(row => row[key]).filter(Boolean))].sort((a, b) => a.localeCompare(b))
const authors = computed(() => options('author')), studies = computed(() => options('gse'))
const assays = computed(() => options('assay')), cells = computed(() => options('cell'))
const filtered = computed(() => filterRecords(props.rows, { q: query.value, assay: assay.value, cell: cell.value, study: study.value }).filter(row => !author.value || row.author === author.value))
watch(filtered, rows => emit('filtered', rows), { immediate: true })
defineExpose({ reset })
const shown = computed(() => filtered.value.slice((page.value - 1) * props.pageSize, page.value * props.pageSize))
const studyCount = computed(() => new Set(filtered.value.map(row => row.gse).filter(Boolean)).size)
function clear() { query.value = ''; assay.value = ''; cell.value = ''; author.value = ''; study.value = ''; page.value = 1 }
function reset() { clear(); emit('reset') }
watch([query, assay, cell, author, study, () => props.rows], () => page.value = 1)
watch(() => [props.symbol, props.transcript], clear)
function download() {
  const columns = [
    ...(hasPi.value ? [{ key: 'gene', label: 'gene' }, { key: 'transcript', label: 'transcript' }] : []),
    { key: 'id', label: 'profile_id' }, sampleNameColumn, { key: 'alias', label: 'source_alias' },
    { key: 'assay', label: 'sequencing_assay' }, { key: 'cell', label: 'raw_cell_label' },
    { key: 'author', label: 'author_source_label' }, { key: 'gse', label: 'GSE' },
    { key: 'antibody', label: 'antibody_from_alias' }, { key: 'source_file', label: 'source_matrix' },
    ...(hasPi.value ? [{ key: 'log2_pi', label: 'log2_PI' }] : []),
  ]
  exportTsv(`${hasPi.value ? `${props.symbol}_${props.transcript}_` : ''}experimental_sources.tsv`, columns, filtered.value.map(row => ({ ...row, gene: props.symbol, transcript: props.transcript })))
}
</script>
<template>
  <section class="panel source-table" aria-labelledby="source-table-title">
    <div class="panel-head source-heading"><div><h2 id="source-table-title">{{title}}</h2><p>{{description}}</p></div><button class="button secondary" type="button" :disabled="!filtered.length" @click="download"><Download :size="16"/>Download sources</button></div>
    <div class="source-filters">
      <label v-if="!hasPi"><span>Sequencing assay</span><select v-model="assay" aria-label="Source sequencing assay"><option value="">All assays</option><option v-for="value in assays" :key="value">{{ value }}</option></select></label>
      <label v-if="!hasPi"><span>Cell line</span><select v-model="cell" aria-label="Source cell line"><option value="">All cell labels</option><option v-for="value in cells" :key="value">{{ value }}</option></select></label>
      <label><span>Author / source</span><select v-model="author" aria-label="Source author"><option value="">All authors</option><option v-for="value in authors" :key="value">{{ value }}</option></select></label>
      <label><span>GEO study</span><select v-model="study" aria-label="Source GEO study"><option value="">All studies</option><option v-for="value in studies" :key="value">{{ value }}</option></select></label>
      <label class="source-query"><span>Find a source</span><div><Search :size="16" aria-hidden="true"/><input v-model="query" aria-label="Search experimental sources" placeholder="Author, GSE or profile" type="search"/></div></label>
    </div>
    <p class="source-context" role="status">{{ filtered.length }} profiles · {{ studyCount }} studies<span v-if="hasPi"> · {{ symbol }} / {{ transcript }} · Current gene filters apply</span></p>
    <div v-if="filtered.length" class="table-scroll source-scroll" tabindex="0" role="region" aria-label="Experimental sources table, scroll horizontally for all columns"><table><caption class="source-caption">Author labels come from the current source metadata; laboratory affiliations have not been verified.</caption><thead><tr><th scope="col">Profile</th><th scope="col">Sequencing assay</th><th scope="col">Cell line</th><th scope="col">Author / source label</th><th scope="col">GEO study</th><th scope="col">Antibody</th><th v-if="hasPi" scope="col">log₂ PI</th></tr></thead><tbody><tr v-for="row in shown" :key="row.id" :data-profile="row.id"><td><RouterLink :to="`/samples/${row.id}`" class="table-link">{{ sampleName(row) }}</RouterLink></td><td><AssayBadge :assay="row.assay"/></td><td>{{ row.cell || 'Not recorded' }}</td><td>{{ row.author || 'Not recorded' }}</td><td><template v-if="row.gse"><RouterLink :to="`/studies/${row.gse}`" class="table-link">{{ row.gse }}</RouterLink><a class="geo-link" :href="`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${encodeURIComponent(row.gse)}`" target="_blank" rel="noopener noreferrer" :aria-label="`Open ${row.gse} at NCBI GEO`">NCBI GEO ↗</a></template><span v-else>Not recorded</span></td><td>{{ row.antibody || 'Not recorded' }}</td><td v-if="hasPi" class="number source-pi">{{ Number.isFinite(row.log2_pi) ? row.log2_pi.toFixed(4) : 'Not measured' }}</td></tr></tbody></table></div>
    <EmptyState v-else title="No matching experimental sources" description="Reset the filters to explore the available profiles."><button class="button secondary" type="button" @click="reset">Reset filters</button></EmptyState>
    <Pagination v-if="filtered.length" v-model="page" :total="filtered.length" :size="pageSize"/>
  </section>
</template>
<style scoped>
.source-table{margin:28px 0;min-width:0;overflow:hidden}.source-heading{gap:18px;align-items:center}.source-heading h2{font-size:23px;font-weight:600;color:#204f78;line-height:1.4;text-wrap:balance}.source-heading p{font-size:13px;line-height:1.7;color:#64748b;margin:7px 0 0;text-wrap:pretty}.source-filters{display:flex;flex-wrap:wrap;gap:16px;padding:20px 24px 12px}.source-filters>label{display:flex;flex:1 1 155px;flex-direction:column;gap:7px;min-width:0}.source-filters>label>span{font-size:12px;font-weight:500;color:#536b83}.source-filters select,.source-query>div{height:42px;border:1px solid #d5e0eb;border-radius:6px;background:white;color:#334e68;font-size:13px;width:100%;min-width:0}.source-filters select{padding:0 10px}.source-query>div{display:flex;align-items:center;gap:8px;padding:0 11px}.source-query svg{flex-shrink:0;color:#7b92a9}.source-query input{width:100%;height:100%;min-width:0;color:#334e68;font-size:13px}.source-context{font-size:12px;color:#64748b;padding:0 24px 17px;line-height:1.7;font-variant-numeric:tabular-nums;text-wrap:pretty}.source-scroll{max-width:100%;overflow:auto}.source-scroll table{min-width:1050px;font-size:14px}.source-scroll th{font-size:12px;line-height:1.6;font-weight:500;color:#536b83;white-space:nowrap}.source-scroll td{line-height:1.65;vertical-align:top;padding-top:17px;padding-bottom:17px}.source-scroll .table-link{font-size:13px;color:#3e77ab}.source-alias{font-size:11px;color:#7b8da2;white-space:nowrap;max-width:200px;overflow:hidden;text-overflow:ellipsis;margin-top:4px}.source-scroll .assay-badge{font-size:12px}.source-scroll .geo-link{display:block;color:#6684a1;font-size:11px;margin-top:4px;width:fit-content}.source-scroll .geo-link:hover{text-decoration:underline}.source-pi{color:#204f78;font-weight:500}.source-caption{caption-side:bottom;text-align:left;white-space:normal;padding:15px 24px;color:#738398;font-size:12px;line-height:1.7;background:#f8fafc;text-wrap:pretty}
@media(max-width:760px){.source-heading{flex-wrap:wrap;align-items:flex-start}.source-heading h2{font-size:21px}.source-filters{padding:18px 16px 12px;gap:12px}.source-filters>label{flex-basis:calc(50% - 12px)}.source-filters>.source-query{flex-basis:100%}.source-context{padding-left:16px;padding-right:16px}.source-caption{padding:14px 16px}}
.source-table>.source-heading{background:#204f78;border-bottom-color:#204f78}.source-table>.source-heading h2{color:#fff}.source-table>.source-heading p{color:#dce8f4}
</style>
