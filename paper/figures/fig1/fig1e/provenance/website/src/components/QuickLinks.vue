<script setup>
import { ArrowUpRight } from '@lucide/vue'
import { features } from '../config/navigation.js'
</script>
<template>
  <section id="quick-links" class="quick-links-section" aria-labelledby="quick-links-title">
    <div class="container quick-links-inner">
      <header class="quick-links-heading"><h2 id="quick-links-title">Quick Links</h2><p>Explore available data and tools in NascentAtlas.</p></header>
      <nav class="quick-links-grid" aria-label="Quick links">
        <RouterLink v-for="link in features" :key="link.key" :to="link.path" class="quick-link-card" :aria-label="link.label" :aria-describedby="`quick-${link.key}-description`">
          <span class="quick-link-icon"><component :is="link.icon" :size="26" :stroke-width="1.8" aria-hidden="true" /></span>
          <h3>{{ link.label }}</h3><p :id="`quick-${link.key}-description`">{{ link.description }}</p>
          <span class="quick-link-action" aria-hidden="true">Explore <ArrowUpRight :size="15" /></span>
        </RouterLink>
      </nav>
    </div>
  </section>
</template>
