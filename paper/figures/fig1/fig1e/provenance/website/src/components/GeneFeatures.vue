<script setup>
import { computed, onBeforeUnmount, ref, watch } from 'vue'
import { Download, Search } from '@lucide/vue'
import { getGeneFeatures } from '../api/client.js'
import { exportTsv } from '../utils/data.js'
const props = defineProps({ symbol: String, transcript: String })
const result = ref(null), loading = ref(false), error = ref(''), category = ref(''), query = ref('')
const view = ref('bars')
const featureHost = ref(null), revealed = ref(new Set()), activeFeature = ref('')
let featureObserver
let revision = 0
async function load() {
  const token = ++revision
  result.value = null; loading.value = true; error.value = ''; category.value = ''; query.value = ''
  try { const value = await getGeneFeatures(props.symbol); if (token === revision) result.value = value }
  catch (e) { if (token === revision) error.value = e.message }
  finally { if (token === revision) loading.value = false }
}
watch(() => [props.symbol, props.transcript], load, { immediate: true })
onBeforeUnmount(() => { ++revision; featureObserver?.disconnect() })
const matched = computed(() => result.value?.gene?.transcript === props.transcript)
const rows = computed(() => matched.value ? result.value.manifest.features.map((feature, index) => {
  const raw = result.value.gene.values[index]
  const value = raw === null ? null : feature.transform === 'length' ? Math.max(0, 10 ** raw - .1) : feature.transform === 'percent' ? raw * 100 : raw
  return { ...feature, raw, value, category: result.value.manifest.groups[feature.group] }
}) : [])
const filtered = computed(() => rows.value.filter(row => (!category.value || row.group === category.value) && `${row.label} ${row.name}`.toLowerCase().includes(query.value.trim().toLowerCase())))
watch([featureHost, filtered, view], () => {
  featureObserver?.disconnect(); revealed.value = new Set(); activeFeature.value = ''
  if (!featureHost.value || view.value !== 'bars') return
  featureObserver = new IntersectionObserver(entries => {
    const next = new Set(revealed.value)
    for (const entry of entries) {
      if (entry.isIntersecting) next.add(entry.target.dataset.feature)
      else next.delete(entry.target.dataset.feature)
    }
    revealed.value = next
  }, { threshold: 0.2 })
  featureHost.value.querySelectorAll('.feature-chart-row').forEach(row => featureObserver.observe(row))
}, { flush: 'post' })
const format = row => row.value === null ? 'Not available' : row.value.toLocaleString('en-US', { minimumFractionDigits: row.precision, maximumFractionDigits: row.precision })
const colors = { B: '#417dc6', Q: '#12a89b', I: '#df9b32', P: '#db738c' }
function scale(unit, values) {
  if (unit === '%') return 100
  const maximum = Math.max(0, ...values.filter(value => value !== null))
  if (!maximum) return 1
  const magnitude = 10 ** Math.floor(Math.log10(maximum))
  return Math.ceil(maximum / magnitude) * magnitude
}
const charts = computed(() => Object.entries(result.value?.manifest.groups || {}).map(([key, label]) => {
  const selected = filtered.value.filter(row => row.group === key)
  const units = [...new Set(selected.map(row => row.unit))].map(unit => {
    const values = selected.filter(row => row.unit === unit)
    return { unit, maximum: scale(unit, values.map(row => row.value)), rows: values }
  })
  return { key, label, units, count: selected.length, color: colors[key] }
}).filter(chart => chart.count))
const width = (row, unit) => row.value === null ? 0 : Math.max(0, Math.min(100, row.value / unit.maximum * 100))
const escapeXml = value => String(value).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&apos;' }[c]))
function downloadImage() {
  const chartWidth = 600, gap = 24, margin = 32
  const heights = charts.value.map(chart => 64 + chart.units.reduce((sum, unit) => sum + (42 + unit.rows.length * 48), 0))
  let y = 114
  const parts = []
  for (let i = 0; i < charts.value.length; i += 2) {
    for (let j = i; j < Math.min(i + 2, charts.value.length); j++) {
      const chart = charts.value[j], x = margin + (j % 2) * (chartWidth + gap)
      parts.push(`<rect x="${x}" y="${y}" width="${chartWidth}" height="${heights[j]}" rx="12" fill="white" stroke="#dce5ef"/><text x="${x + 24}" y="${y + 34}" font-size="20" font-weight="600" fill="${chart.color}">${escapeXml(chart.label)}</text>`)
      let rowY = y + 64
      for (const unit of chart.units) {
        parts.push(`<text x="${x + 24}" y="${rowY}" font-size="13" fill="#64748b">${escapeXml(unit.unit)} · scale 0–${unit.maximum.toLocaleString('en-US')}</text>`)
        rowY += 24
        for (const row of unit.rows) {
          parts.push(`<g><title>${escapeXml(row.name)}</title><text x="${x + 24}" y="${rowY}" font-size="14">${escapeXml(row.label)}</text><text x="${x + chartWidth - 24}" y="${rowY}" text-anchor="end" font-size="14" font-weight="600">${escapeXml(format(row))} ${row.value === null ? '' : escapeXml(row.unit)}</text><rect x="${x + 24}" y="${rowY + 10}" width="${chartWidth - 48}" height="10" rx="3" fill="#eef3f8"/><rect x="${x + 24}" y="${rowY + 10}" width="${(chartWidth - 48) * width(row, unit) / 100}" height="10" rx="3" fill="${chart.color}"/></g>`)
          rowY += 48
        }
        rowY += 18
      }
    }
    y += Math.max(...heights.slice(i, i + 2)) + gap
  }
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="1288" height="${y + 46}" viewBox="0 0 1288 ${y + 46}"><rect width="100%" height="100%" fill="#f4f7fb"/><g font-family="Arial, sans-serif" fill="#23415f"><text x="32" y="44" font-size="26" font-weight="600">${escapeXml(props.symbol)} · Extracted gene features</text><text x="32" y="74" font-size="15">${escapeXml(props.transcript)} · ${filtered.value.length} actual reference features · ${escapeXml(result.value.manifest.provenance.assembly)}</text>${parts.join('')}<text x="32" y="${y + 20}" font-size="13" fill="#64748b">Separate scales by category and unit. Bar lengths show actual values; missing values are unavailable.</text></g></svg>`
  const url = URL.createObjectURL(new Blob([svg], { type: 'image/svg+xml;charset=utf-8' }))
  const link = document.createElement('a'); link.href = url; link.download = `${props.symbol}_${props.transcript}_features.svg`; link.click()
  setTimeout(() => URL.revokeObjectURL(url), 1000)
}
function download() {
  exportTsv(`${props.symbol}_${props.transcript}_features.tsv`, [
    { label: 'gene', value: () => props.symbol }, { label: 'transcript', value: () => props.transcript },
    { key: 'category', label: 'category' }, { key: 'name', label: 'feature' },
    { key: 'raw', label: 'original_feature_value' },
    { label: 'display_value', value: row => row.value === null ? '' : Number(row.value.toFixed(row.precision === 0 ? 0 : 8)) },
    { key: 'unit', label: 'display_unit' }, { key: 'transform', label: 'display_transform' },
  ], filtered.value)
}
</script>

<template>
  <section class="panel gene-features-card" aria-label="Gene actual feature values">
    <header class="panel-head"><div><h2>{{ symbol }} · Extracted gene features</h2><p>Actual feature values for {{ transcript }}.</p></div><div v-if="matched" class="feature-downloads"><button class="button secondary" :disabled="!filtered.length" @click="downloadImage"><Download :size="15"/>Download SVG</button><button class="button secondary" :disabled="!filtered.length" @click="download"><Download :size="15"/>Download features</button></div></header>
    <p v-if="loading" class="feature-message" role="status">Loading extracted features…</p>
    <div v-else-if="error" class="feature-message" role="alert"><p>{{ error }}</p><button class="button secondary" @click="load">Retry features</button></div>
    <p v-else-if="!matched" class="feature-message">No extracted feature record for this transcript<span v-if="result?.gene">. Available record: {{ result.gene.transcript }}</span>.</p>
    <template v-else>
      <div class="feature-controls"><label>Display mode<select v-model="view" aria-label="Feature display mode"><option value="bars">Horizontal bars</option><option value="table">Value table</option></select></label><label>Category<select v-model="category" aria-label="Feature category"><option value="">All categories</option><option v-for="(label, key) in result.manifest.groups" :key="key" :value="key">{{ label }}</option></select></label><label class="feature-search"><span>Find a feature</span><div><Search :size="16" aria-hidden="true"/><input v-model="query" aria-label="Search gene features" placeholder="GC, intron, CpG…"/></div></label><span class="feature-count">{{ filtered.length }} / {{ rows.length }} features</span></div>
      <p v-if="view !== 'table'" class="feature-scale-note">Bar lengths show actual values. Each category and unit has a separate scale. Scroll inside the charts to see more features.</p>
      <div v-if="view !== 'table'" ref="featureHost" class="feature-charts">
        <section v-for="chart in charts" :key="chart.key" class="feature-chart" :style="{'--feature-color':chart.color}" :aria-label="chart.label + ' feature chart'">
          <h3>{{ chart.label }}<span>{{ chart.count }} features</span></h3>
          <div v-for="unit in chart.units" :key="unit.unit" class="feature-unit-chart">
            <div class="feature-scale"><strong>{{ unit.unit }}</strong><span>0–{{ unit.maximum.toLocaleString('en-US') }}</span></div>
            <div class="feature-unit-rows"><div v-for="row in unit.rows" :key="row.name" class="feature-chart-row" :class="{'is-revealed':revealed.has(row.name),'is-active':activeFeature===row.name}" :data-feature="row.name" :data-raw="row.raw ?? ''" :data-value="row.value ?? ''" tabindex="0" role="group" :aria-label="`${row.label}: ${format(row)} ${row.value === null ? '' : row.unit}`" :aria-describedby="activeFeature===row.name ? 'feature-tip-'+row.name : undefined" @mouseenter="activeFeature=row.name" @mouseleave="activeFeature=''" @focus="activeFeature=row.name" @blur="activeFeature=''">
              <div class="feature-row-label"><span>{{ row.label }}</span><strong class="feature-value" :class="{'feature-missing':row.value === null}">{{ format(row) }}<small v-if="row.value !== null">{{ row.unit }}</small></strong></div>
              <svg class="feature-bar" viewBox="0 0 100 10" preserveAspectRatio="none" role="img" :aria-label="`${row.label}: ${format(row)} ${row.value === null ? '' : row.unit}; scale 0 to ${unit.maximum} ${unit.unit}`"><rect width="100" height="10" rx="1" fill="#edf2f8"/><g class="feature-growth"><rect class="feature-fill" :width="width(row, unit)" height="10" rx="1" fill="var(--feature-color)"/></g></svg>
              <div v-if="activeFeature===row.name" :id="'feature-tip-'+row.name" class="feature-tooltip" role="tooltip"><strong>{{row.label}}</strong><span>{{format(row)}} {{row.value===null ? '' : row.unit}}</span><small>{{row.name}} · Original value: {{row.raw ?? 'Not available'}}</small></div>
            </div></div>
          </div>
        </section>
      </div>
      <div v-else class="table-scroll feature-table"><table><thead><tr><th>Feature</th><th>Category</th><th>Value</th><th>Unit</th></tr></thead><tbody><tr v-for="row in filtered" :key="row.name" :data-feature="row.name" :data-raw="row.raw ?? ''" :data-value="row.value ?? ''"><td><strong>{{ row.label }}</strong><small>{{ row.name }}</small></td><td>{{ row.category }}</td><td class="feature-value" :class="{'feature-missing':row.value === null}">{{ format(row) }}</td><td>{{ row.unit }}</td></tr></tbody></table></div>
      <p v-if="!filtered.length" class="feature-message">No matching features. Clear the search or choose another category.</p>
      <footer class="feature-footer">MANE-based reference features · Lengths in bp · Fractions as percentages · Missing values remain unavailable.</footer>
      <details class="feature-methods"><summary>Feature values and sources</summary><p>{{ result.manifest.provenance.annotation }} · {{ result.manifest.provenance.assembly }} · {{ result.gene.transcript }} · {{ result.gene.chrom }} ({{ result.gene.strand }}).</p><p>Length fields stored as log₁₀(length + 0.1) are restored to bp and verified against the original core table. Percentages are the stored fractions × 100. The download preserves the original feature values. These reference features do not change when assay or cell filters change.</p><dl><template v-for="(hash, name) in result.manifest.source_hashes" :key="name"><dt>{{ name }}</dt><dd>{{ hash }}</dd></template></dl></details>
    </template>
  </section>
</template>

<style scoped>
.feature-chart-row{position:relative;border-radius:4px}.feature-chart-row.is-active{background:#eef3f9}.feature-chart-row:focus-visible{outline:2px solid #417dc6;outline-offset:3px}.feature-growth{transform-origin:left center;opacity:0}.is-revealed .feature-growth{animation:feature-bar-reveal 600ms ease-out both}.feature-fill{transform-origin:left center;transition:transform 180ms ease-out}.feature-chart-row.is-active .feature-fill{transform:scaleY(1.25)}.feature-chart-row .feature-bar{overflow:visible}.feature-tooltip{position:absolute;bottom:calc(100% + 4px);right:0;max-width:100%;z-index:4;background:#fff;border:1px solid #d5dfeb;border-radius:5px;box-shadow:0 3px 10px #0f172a20;padding:8px 12px;pointer-events:none;font-size:12px;line-height:1.7;color:#345776;overflow-wrap:anywhere}.feature-tooltip span,.feature-tooltip small{display:block;font-variant-numeric:tabular-nums}.feature-tooltip small{font-size:10px;color:#64748b}@keyframes feature-bar-reveal{from{opacity:0;transform:scaleX(0)}to{opacity:1;transform:scaleX(1)}}@media(prefers-reduced-motion:reduce){.feature-growth{opacity:1}.is-revealed .feature-growth{animation:none}.feature-fill{transition:none;transform:none!important}}
.feature-downloads{display:flex;gap:8px;flex-wrap:wrap}.feature-scale-note{margin:0;padding:0 20px 16px;font-size:12px;color:#64748b;text-wrap:pretty}.feature-charts{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px;padding:0 20px 20px;align-items:start;max-height:640px;overflow:auto;scrollbar-width:thin}.feature-chart{border:1px solid #dce5ef;border-radius:8px;padding:18px;min-width:0}.feature-chart h3{display:flex;align-items:baseline;justify-content:space-between;gap:12px;font-size:16px;color:#23415f;margin:0 0 16px;text-wrap:balance}.feature-chart h3 span{font-size:11px;font-weight:400;color:#64748b;white-space:nowrap}.feature-unit-chart+.feature-unit-chart{margin-top:18px;padding-top:14px;border-top:1px solid #e5ebf2}.feature-scale{display:flex;justify-content:space-between;gap:12px;font-size:11px;color:#64748b;margin-bottom:8px;font-variant-numeric:tabular-nums}.feature-scale strong{font-weight:500}.feature-chart-row{padding:7px 0}.feature-row-label{display:flex;justify-content:space-between;gap:12px;margin-bottom:6px;font-size:12px;color:#345776;align-items:baseline}.feature-row-label>span{min-width:0;text-wrap:pretty}.feature-value{font-variant-numeric:tabular-nums;font-size:12px;white-space:nowrap;font-weight:600}.feature-value small{margin-left:5px;font-weight:400;font-size:10px;color:#64748b}.feature-missing{color:#8392a4;font-weight:400}.feature-bar{display:block;width:100%;height:10px}@media(max-width:860px){.feature-charts{grid-template-columns:1fr;gap:16px}.feature-chart{padding:14px}.feature-downloads{width:100%}}@media(max-width:480px){.feature-charts{padding:0 12px 16px}.feature-controls select{max-width:100%}}
.gene-features-card{margin-top:24px;min-width:0}.gene-features-card .panel-head{padding:18px 20px;gap:14px;flex-wrap:wrap}.gene-features-card .panel-head p{font-size:13px}.gene-features-card .button{font-size:12px}.feature-controls{display:flex;align-items:flex-end;gap:16px;padding:8px 20px 18px;flex-wrap:wrap}.feature-controls label{display:flex;flex-direction:column;gap:8px;font-size:12px;color:#536b83}.feature-controls select,.feature-search>div{border:1px solid #d5dfeb;border-radius:6px;background:#fff;min-height:38px;color:#345776}.feature-controls select{padding:7px 10px;font:inherit}.feature-search>div{display:flex;align-items:center;gap:8px;padding:7px 10px}.feature-search input{min-width:0;font:inherit;color:#345776;flex:1}.feature-count{font-size:12px;color:#64748b;padding-bottom:10px;font-variant-numeric:tabular-nums}.feature-table{max-height:440px;scrollbar-width:thin}.feature-table table{min-width:650px;font-size:13px;white-space:normal}.feature-table th{position:sticky;top:0;background:#f5f8fc;font-size:12px;z-index:1}.feature-table td,.feature-table th{padding:12px 20px}.feature-table td:first-child{width:42%}.feature-table strong{font-weight:500;color:#345776}.feature-table small{display:block;color:#64748b;font-size:10px;overflow-wrap:anywhere;margin-top:3px}.feature-table tbody tr:nth-child(odd){background:#f7f9fc}.feature-table .feature-value{text-align:right;font-weight:600;color:#345776;white-space:nowrap;font-variant-numeric:tabular-nums}.feature-table th:nth-child(3){text-align:right}.feature-table .feature-missing{font-weight:400;color:#8392a4;font-size:12px}.feature-footer{padding:12px 20px;border-top:1px solid #e4eaf1;background:#f8fafc;font-size:12px;line-height:1.7;color:#64748b;text-wrap:pretty}.feature-message{padding:18px 20px;color:#64748b;font-size:13px;line-height:1.8}.feature-message p{margin:0 0 12px}.feature-methods{padding:12px 20px;font-size:12px;color:#64748b;line-height:1.8}.feature-methods summary{cursor:pointer;color:#345776;font-weight:500}.feature-methods p{margin:10px 0;text-wrap:pretty}.feature-methods dl{font-size:11px}.feature-methods dd{margin:0 0 8px;overflow-wrap:anywhere;font-family:ui-monospace,monospace}.gene-features-card :focus-visible{outline:2px solid #4979ad;outline-offset:2px}@media(max-width:760px){.feature-controls{gap:12px}.feature-search{flex:1;min-width:180px}.feature-table{max-height:400px}.feature-table table{min-width:0;table-layout:fixed}.feature-table th:nth-child(2),.feature-table td:nth-child(2){display:none}.feature-table td:first-child,.feature-table th:first-child{width:58%;overflow-wrap:anywhere}.feature-table th:nth-child(3){width:28%}.feature-table th:last-child{width:14%}.feature-table td,.feature-table th{padding:12px 10px}.feature-table .feature-missing{white-space:normal}.feature-table td:last-child{overflow-wrap:anywhere;font-size:11px}.feature-count{flex-basis:100%;padding-bottom:0}}
</style>
