<script setup>
import { computed } from 'vue'
import { assays, assayColors } from '../utils/data.js'
import CellComposition from './CellComposition.vue'
import { useChartReveal } from '../composables/chartReveal.js'
const props = defineProps({ catalog: Object })
const emit = defineEmits(['select'])
const counts = computed(() => assays.map(assay => ({ assay, count: props.catalog.profiles.filter(p => p.assay === assay).length })))
const { host: assayHost, visible: assayVisible } = useChartReveal({ once: true, threshold: 0.12 })
const { host: cellHost, visible: cellVisible } = useChartReveal({ once: true, threshold: 0.12 })
</script>
<template>
  <div class="chart-grid">
    <div ref="assayHost" class="collection-chart-reveal" :class="{ 'is-visible': assayVisible }">
    <section class="panel assay-composition">
      <div class="panel-head">
        <div><h2>Assay composition</h2><p>Each profile is counted as a separate catalog record</p></div>
      </div>
      <div class="assay-bars">
        <button v-for="(item,index) in counts" :key="item.assay" class="bar-row" :style="{'--bar-delay': `${index*70}ms`}" :aria-label="`Filter ${item.assay}`" @click="emit('select',{ assay:item.assay })">
          <span>{{ item.assay }}</span>
          <div class="bar-track"><div :style="{width: `${item.count / catalog.profiles.length * 100}%`, background: assayColors[item.assay]}" /></div>
          <b>{{ item.count }}</b><small>{{ (item.count/catalog.profiles.length*100).toFixed(1) }}%</small>
        </button>
      </div>
      <div class="panel-foot">Select a bar to explore profiles for that assay →</div>
    </section>
    </div>
    <div ref="cellHost" class="collection-chart-reveal" :class="{ 'is-visible': cellVisible }" style="--reveal-delay:120ms">
      <CellComposition :catalog="catalog" @select="emit('select',$event)"/>
    </div>
  </div>
</template>

<style scoped>
.collection-chart-reveal{display:flex;min-width:0;opacity:0}
.collection-chart-reveal>.panel{width:100%}
.collection-chart-reveal.is-visible{animation:collection-fade-in-up 650ms ease-out var(--reveal-delay,0ms) both}
.collection-chart-reveal:not(.is-visible):focus-within{opacity:1;transform:none;animation:none}
@keyframes collection-fade-in-up{from{opacity:0;transform:translateY(20px)}to{opacity:1;transform:translateY(0)}}
@media(prefers-reduced-motion:reduce){.collection-chart-reveal{opacity:1;transform:none}.collection-chart-reveal.is-visible{animation:none}}
.assay-composition{display:flex;flex-direction:column}
.assay-composition>.assay-bars{flex:1;display:grid;grid-template-rows:repeat(3,minmax(0,1fr));align-items:center;gap:16px;padding:20px 24px;min-height:330px}
.assay-composition .bar-row{min-height:64px;padding:12px 0;font-size:13px}
.assay-composition .bar-track{height:32px}
.assay-composition .bar-row b{font-size:15px}
.assay-composition .bar-row small{font-size:12px}
.assay-composition>.panel-foot{margin-top:auto}
@media(max-width:760px){.assay-composition>.assay-bars{min-height:240px;padding:16px 18px;gap:8px}.assay-composition .bar-row{min-height:56px;font-size:12px}.assay-composition .bar-track{height:28px}}
</style>
