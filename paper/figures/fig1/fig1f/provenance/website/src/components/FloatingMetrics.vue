<script setup>
import { ArrowUpRight } from '@lucide/vue'
defineProps({ cards: Array })
function move(event) {
  if (event.pointerType !== 'mouse' || window.matchMedia('(prefers-reduced-motion: reduce)').matches) return
  const box = event.currentTarget.getBoundingClientRect()
  const x = Math.max(-.5, Math.min(.5, (event.clientX - box.left) / box.width - .5))
  const y = Math.max(-.5, Math.min(.5, (event.clientY - box.top) / box.height - .5))
  event.currentTarget.style.setProperty('--tilt-x', `${-y}deg`)
  event.currentTarget.style.setProperty('--tilt-y', `${x}deg`)
}
function reset(event) { event.currentTarget.style.removeProperty('--tilt-x'); event.currentTarget.style.removeProperty('--tilt-y') }
</script>
<template>
  <div class="metric-grid floating-metrics" aria-label="Collection summary">
    <RouterLink v-for="card in cards" :key="card.label" :to="card.to" class="metric-card floating-metric-card" :aria-label="card.accessible" @pointermove="move" @pointerleave="reset" @blur="reset">
      <div class="metric-top"><span>{{ card.label }}</span><span class="metric-icon"><component :is="card.icon" :size="19" aria-hidden="true" /></span></div>
      <strong>{{ card.value }}</strong>
      <span class="metric-action" aria-hidden="true">{{ card.action }} <ArrowUpRight :size="15" /></span>
    </RouterLink>
  </div>
</template>
