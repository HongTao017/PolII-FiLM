<script setup>
import { ref, watch } from 'vue'
import { useRoute } from 'vue-router'
import { Dna, ArrowUpRight, Menu, X, Database } from '@lucide/vue'
import { mainNavigation, isNavigationActive } from './config/navigation.js'
const route = useRoute()
const menuOpen = ref(false)
const links = mainNavigation
watch(() => route.path, () => { menuOpen.value = false })
</script>
<template><a class="skip-link" href="#main">Skip to main content</a><header class="site-header"><div class="header-inner"><RouterLink to="/" class="brand" @click="menuOpen=false"><span class="brand-mark"><Dna :size="27" /></span><div>Nascent<span>Atlas</span><small>NASCENT TRANSCRIPTION ATLAS</small></div></RouterLink><button class="mobile-menu" aria-label="Toggle navigation" :aria-expanded="menuOpen" @click="menuOpen=!menuOpen"><X v-if="menuOpen" :size="22" /><Menu v-else :size="22" /></button><nav :class="{open:menuOpen}" aria-label="Main navigation"><RouterLink v-for="link in links" :key="link.path" :to="link.path" :class="{active:isNavigationActive(route.path,link.path)}" :aria-current="isNavigationActive(route.path,link.path) ? 'page' : undefined" @click="menuOpen=false">{{ link.label }}</RouterLink></nav></div></header><main id="main" class="main-container" :class="route.path==='/' ? 'home-main' : ['/samples','/genes'].includes(route.path) ? 'browse-main' : 'page-main'"><RouterView /></main><footer class="site-footer"><div class="container footer-top"><div><div class="footer-brand"><Dna :size="20" />NascentAtlas</div><p>Discover, explore and trace public nascent transcription data.</p></div><div class="footer-links"><RouterLink to="/help">Methods and data</RouterLink><RouterLink to="/download">Download catalog</RouterLink><a href="https://www.ncbi.nlm.nih.gov/geo/" target="_blank" rel="noopener">NCBI GEO <ArrowUpRight :size="13" /></a></div></div><div class="container footer-bottom"><span><Database :size="13" /> PRO-seq · GRO-seq · mNET-seq</span><span>Snapshot 2026-10-09 · Front end v0.1.0</span></div></footer></template>
