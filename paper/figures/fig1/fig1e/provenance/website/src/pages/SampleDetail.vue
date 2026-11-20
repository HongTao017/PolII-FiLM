<script setup>
import { sampleName } from '../utils/sampleName.js'
import { computed } from 'vue'
import { useRoute } from 'vue-router'
import { ArrowLeft, Activity } from '@lucide/vue'
import { useCatalog } from '../composables/catalog.js'
import PageTitle from '../components/PageTitle.vue'
import LoadState from '../components/LoadState.vue'
import EmptyState from '../components/EmptyState.vue'
import AssayBadge from '../components/AssayBadge.vue'
import ProfilePiAnalysis from '../components/ProfilePiAnalysis.vue'
import LibraryQcTable from '../components/LibraryQcTable.vue'
import ProfileProvenance from '../components/ProfileProvenance.vue'
const { catalog,loading,error,reload }=useCatalog(),route=useRoute()
const profile=computed(()=>catalog.value?.profiles.find(p=>p.id===route.params.id))
</script>
<template><LoadState :loading="loading" :error="error" @retry="reload"><template v-if="profile"><PageTitle eyebrow="PROFILE DETAIL" :title="sampleName(profile)" description="A profile record from an existing pausing-index matrix."><AssayBadge :assay="profile.assay"/></PageTitle><RouterLink to="/samples" class="back-link"><ArrowLeft :size="15"/>Back to profile catalog</RouterLink><ProfileProvenance :profile="profile"/><LibraryQcTable id="library-qc" class="latest-panel" :profiles="[profile]" detail/><ProfilePiAnalysis :profile="profile" :catalog="catalog"/><div class="latest-panel"><section class="panel"><div class="panel-head"><h2><Activity :size="17"/>Genome signal tracks</h2></div><div class="pi-heatmap-body"><p class="pi-method-note">Open a gene to explore the connected reference annotation and profile signal tracks.</p><RouterLink class="table-link" :to="{path:'/genes',query:{assay:profile.assay,cell:profile.cell}}">Open genome signal browser →</RouterLink></div></section></div><div class="notice">Select {{ profile.assay }} and {{ profile.cell }} in the gene explorer to view pausing metrics from the existing matrices.<RouterLink :to="{path:'/genes',query:{assay:profile.assay,cell:profile.cell}}">Open gene explorer →</RouterLink></div></template><EmptyState v-else-if="catalog" title="Profile not found" description="Check the record ID or return to the profile catalog." to="/samples"/></LoadState></template>

<style scoped>
:deep(.page-title>div:first-child){min-width:0}:deep(.page-title h1){overflow-wrap:anywhere}
</style>
