<script setup>
import { ref } from 'vue'
import { useRouter } from 'vue-router'
import { ArrowRight, Search, X } from '@lucide/vue'
import { useCatalog } from '../composables/catalog.js'
import LoadState from '../components/LoadState.vue'
import GeneGenomeBrowser from '../components/GeneGenomeBrowser.vue'
import QuickLinks from '../components/QuickLinks.vue'
import CollectionCharts from '../components/CollectionCharts.vue'
const { catalog, loading, error, reload } = useCatalog()
const router = useRouter()
const query = ref('')
const searchInput = ref(null)
const exampleProfileIds = ['NP-410E006606','NP-3E3F5FDE87','NP-663C656923']
function search(q = query.value) { if(q.trim()) router.push({ path:'/search', query:{ q:q.trim() } }) }
function clearSearch() { query.value = ''; searchInput.value?.focus() }
</script>
<template>
  <section class="website-introduction" aria-labelledby="website-introduction-title">
    <div class="container website-introduction-inner">
      <div class="website-introduction-overview">
        <h1 id="website-introduction-title">NascentAtlas</h1>
        <p class="website-introduction-subtitle">Human Nascent Transcription Database</p>
        <p class="website-introduction-description">Explore human PRO-seq, GRO-seq and mNET-seq profiles and promoter-proximal pausing.</p>
        <form class="hero-search reference-search" role="search" aria-label="Search database" @submit.prevent="search()">
          <Search :size="23" class="search-icon" aria-hidden="true"/>
          <input ref="searchInput" v-model="query" aria-label="Search genes, studies or cells" placeholder="Search genes, studies, profiles or cells"/>
          <button v-if="query" type="button" class="search-clear" aria-label="Clear search" @click="clearSearch"><X :size="18" aria-hidden="true"/></button>
          <button type="submit" class="search-submit">Search <ArrowRight :size="18" aria-hidden="true"/></button>
        </form>
        <RouterLink to="/help" class="website-introduction-link">Learn More <ArrowRight :size="17" aria-hidden="true"/></RouterLink>
      </div>
    </div>
  </section>
  <div class="container home-content">
    <LoadState :loading="loading" :error="error" @retry="reload">
      <template v-if="catalog">
        <section class="home-distribution" aria-label="Collection distribution charts">
          <CollectionCharts :catalog="catalog" @select="filters=>router.push({path:'/samples',query:filters})"/>
          <div class="distribution-heading distribution-footer">
            <RouterLink to="/statistics" class="text-link">Explore collection statistics <ArrowRight :size="16" aria-hidden="true"/></RouterLink>
          </div>
        </section>
        <GeneGenomeBrowser class="home-genome-preview" compact symbol="SAMD11" transcript="NM_001385641.1" :profiles="catalog.profiles" :default-profile-ids="exampleProfileIds"/>
      </template>
    </LoadState>
  </div>
  <QuickLinks/>
</template>
<style scoped>
.home-distribution { margin-top: 32px; }
.home-genome-preview { margin: 32px 0 24px; }
.distribution-footer { margin: 24px auto 0; }
@media(max-width:760px){.home-genome-preview{margin:24px 0 32px}}
</style>
