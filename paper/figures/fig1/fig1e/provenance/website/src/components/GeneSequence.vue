<script setup>
import { computed, onBeforeUnmount, ref, watch } from 'vue'
import { geneRegions } from '../utils/geneRegions.js'
import GeneReferenceLogo from './GeneReferenceLogo.vue'
import SequenceRegionStats from './SequenceRegionStats.vue'
const props = defineProps({ symbol: { type: String, required: true }, transcript: { type: String, default: '' } })
const data = ref(null), loading = ref(false), error = ref(''), actionMessage = ref(''), orientation = ref('gene'), position = ref(1), body = ref(null), basesPerLine = ref(40)
const region = ref('tss600'), signalSequenceSelection = ref({})
const signed = value => value > 0 ? `+${number(value)}` : number(value)
let controller, requestId = 0, requestedOffset = 0, observer
const number = value => Number(value).toLocaleString('en-US')
const rows = computed(() => {
  if (!data.value) return []
  const result = [], value = data.value
  for (let index = 0; index < value.sequence.length; index += basesPerLine.value) {
    const bases = [...value.sequence.slice(index, index + basesPerLine.value)].map((base, offset) => {
      const i = index + offset
      const coordinate = value.display_strand === '-' ? value.genomic_end - i : value.genomic_start + i
      return { base, index: i, coordinate, genePosition: value.display_strand === value.gene_strand ? value.offset + i + 1 : value.offset + value.length - i, relative: value.endpoint ? (coordinate-value.endpoint.anchor)*(value.gene_strand==='-'?-1:1) : null }
    })
    const groups = []
    for (let i = 0; i < bases.length; i += 10) groups.push(bases.slice(i, i + 10))
    result.push({ index, coordinate: bases[0].coordinate, groups })
  }
  return result
})
const fasta = computed(() => {
  const v = data.value
  if (!v) return ''
  const interval = v.endpoint ? `${v.endpoint.anchor_name}_offsets=${v.endpoint.relative_start}:${v.endpoint.relative_end}` : `gene_bp=${v.offset+1}-${v.offset+v.length}`
  return `>${v.symbol}|${v.transcript}|${v.assembly}|${v.chrom}:${v.genomic_start}-${v.genomic_end}|strand=${v.display_strand}|reference_genomic_DNA|region=${v.region}|${interval}\n${v.sequence.match(/.{1,60}/g).join('\n')}\n`
})
async function load(offset = requestedOffset) {
  controller?.abort(); controller = new AbortController()
  const id = ++requestId, signal = controller.signal
  requestedOffset = offset; loading.value = true; error.value = ''; actionMessage.value = ''
  const timeout = setTimeout(() => controller?.signal === signal && controller.abort(), 20000)
  try {
    const response = await fetch('/api/v1/genome/sequence?' + new URLSearchParams({ symbol: props.symbol, transcript: props.transcript, offset, length: 500, orientation: orientation.value, region: region.value }), { signal })
    const value = await response.json()
    if (!response.ok) throw new Error(value.error || 'Unable to read the reference sequence.')
    if (id === requestId) { data.value = value; position.value = value.offset + 1 }
  } catch (e) {
    if (id === requestId) { data.value = null; error.value = e.name === 'AbortError' ? 'Sequence loading timed out. Please retry.' : e.message }
  } finally { clearTimeout(timeout); if (id === requestId) loading.value = false }
}
function jump() {
  if (!data.value || loading.value) return
  if (!Number.isInteger(position.value) || position.value < 1 || position.value > data.value.total_length) {
    actionMessage.value = `Enter a gene position from 1 to ${number(data.value.total_length)}.`; return
  }
  load(position.value - 1)
}
async function copy() {
  if (!data.value || loading.value) return
  try {
    if (navigator.clipboard && window.isSecureContext) await navigator.clipboard.writeText(data.value.sequence)
    else {
      const active = document.activeElement, input = document.createElement('textarea')
      input.value = data.value.sequence; input.readOnly = true; input.style.position = 'fixed'; input.style.opacity = '0'; input.setAttribute('aria-label', 'Sequence to copy'); document.body.append(input)
      let copied = false
      try { input.select(); copied = document.execCommand('copy') } finally { input.remove(); active?.focus() }
      if (!copied) throw new Error('Copy is unavailable. Use FASTA to download this window.')
    }
    actionMessage.value = `${number(data.value.length)} bases copied.`
  } catch (e) { actionMessage.value = e.message || 'Copy is unavailable. Use FASTA to download this window.' }
}
function download() {
  if (!data.value || loading.value) return
  const url = URL.createObjectURL(new Blob([fasta.value], { type: 'text/plain;charset=utf-8' })), link = document.createElement('a')
  link.href = url; link.download = `${data.value.symbol}_${data.value.transcript}_${data.value.chrom}_${data.value.genomic_start}-${data.value.genomic_end}_${data.value.display_strand === '-' ? 'minus' : 'plus'}.fasta`
  link.click(); setTimeout(() => URL.revokeObjectURL(url), 1000)
}
watch(() => [props.symbol, props.transcript], () => { data.value = null; position.value = 1; requestedOffset = 0; orientation.value = 'gene'; load(0) }, { immediate: true })
watch(body, element => {
  observer?.disconnect()
  if (!element) return
  observer = new ResizeObserver(entries => { basesPerLine.value = Math.min(60, Math.max(10, Math.floor((entries[0].contentRect.width - 86) / 94) * 10)) })
  observer.observe(element)
})
onBeforeUnmount(() => { ++requestId; controller?.abort(); observer?.disconnect() })
</script>

<template>
  <section class="panel gene-sequence-card" aria-label="Gene DNA sequence">
    <div class="panel-head"><div><h3>{{ symbol }} · Sequence regions</h3><p>TSS · first intron splice sites · TES · MANE reference.</p></div></div>
    <div class="sequence-region"><label>Sequence region<select v-model="region" aria-label="DNA sequence region" @change="load(0)"><option v-for="r in geneRegions" :key="r.id" :value="r.id">{{r.label}}</option><option value="gene">Gene body · 500 bp windows</option></select></label></div>
    <template v-if="data && !loading">
      <div class="sequence-meta"><strong>{{ data.transcript }} · {{ data.chrom }}:{{ number(data.genomic_start) }}–{{ number(data.genomic_end) }} ({{ data.display_strand }})</strong><p v-if="data.endpoint">{{data.endpoint.anchor_name}} offsets {{signed(data.endpoint.relative_start)}} to {{signed(data.endpoint.relative_end)}} · {{number(data.length)}} bases · {{data.assembly}}</p><p v-else>Gene bp {{ number(data.offset + 1) }}–{{ number(data.offset + data.length) }} of {{ number(data.total_length) }} · {{ data.assembly }}</p><p v-if="data.endpoint?.clipped">Window clipped at the chromosome boundary.</p></div>
      <div class="sequence-actions"><label><span class="sr-only">Sequence orientation</span><select v-model="orientation" aria-label="Sequence orientation" @change="load()"><option value="gene">Gene strand (5′ → 3′)</option><option value="reference">Reference strand (+)</option></select></label><button class="button secondary" type="button" aria-label="Copy displayed DNA sequence" @click="copy">Copy</button><button class="button secondary" type="button" aria-label="Download window FASTA" @click="download">FASTA</button></div>
      <SequenceRegionStats :data="data"/>
      <GeneReferenceLogo :data="data" v-model:selection="signalSequenceSelection"/>
      <details class="sequence-raw" :open="data.region==='gene'"><summary>View ATGC sequence · {{number(data.length)}} bases</summary>
      <div class="sequence-composition" aria-label="Displayed sequence composition"><span v-for="base in ['A', 'T', 'G', 'C']" :key="base"><b :class="'base-' + base">{{ base }}</b>{{ number(data.counts[base]) }}</span><span v-if="data.counts.N || data.counts.other">N / other {{ number(data.counts.N + data.counts.other) }}</span></div>
      <div class="sequence-column-labels" aria-hidden="true"><span>Position</span><span>{{ data.display_strand === '-' ? 'Reverse complement · 5′ → 3′' : 'Reference (+) · 5′ → 3′' }}</span></div>
      <div ref="body" class="sequence-body" tabindex="0" role="region" aria-label="Displayed DNA bases" :data-sequence="data.sequence" :data-offset="data.offset" :data-strand="data.display_strand">
        <div v-for="row in rows" :key="row.index" class="sequence-line"><span class="sequence-coordinate" aria-hidden="true">{{ number(row.coordinate) }}</span><span class="sequence-groups"><span v-for="(group, groupIndex) in row.groups" :key="groupIndex" class="sequence-base-group"><span v-for="letter in group" :key="letter.index" class="sequence-base" :class="'base-' + letter.base" :data-coordinate="letter.coordinate" :data-relative="letter.relative" :title="`${data.chrom}:${number(letter.coordinate)} · ${letter.base} · ${data.endpoint ? data.endpoint.anchor_name+' '+signed(letter.relative) : 'Gene bp '+number(letter.genePosition)}`">{{ letter.base }}</span></span></span></div>
      </div>
      <div v-if="data.region==='gene'" class="sequence-navigation"><button class="button secondary" type="button" aria-label="Previous DNA sequence window" :disabled="data.offset === 0" @click="load(Math.max(0, data.offset - 500))">Previous</button><form @submit.prevent="jump"><label>Gene bp <input v-model.number="position" aria-label="Start at gene base" type="number" min="1" :max="data.total_length" step="1" inputmode="numeric"/></label><button type="submit" class="button secondary">Go</button></form><button class="button secondary" type="button" aria-label="Next DNA sequence window" :disabled="data.offset + data.length >= data.total_length" @click="load(data.offset + data.length)">Next</button></div>
      </details>
    </template>
    <div v-if="loading" class="sequence-loading" role="status" aria-busy="true"><span/><span/><span/><span/><p>Loading reference sequence…</p></div>
    <div v-if="error" class="sequence-error" role="alert"><strong>Sequence unavailable</strong><p>{{ error }}</p><button class="button secondary" type="button" @click="load()">Retry sequence</button></div>
    <p v-if="actionMessage && !loading" class="sequence-action-status" role="status">{{ actionMessage }}</p>
    <footer class="sequence-foot"><template v-if="data && !data.exact_transcript">{{ transcript }} is absent from MANE v1.5; this view uses {{ data.transcript }}. </template><template v-if="data?.region?.startsWith('splice')">Splice offset 0 is the first intronic base at the donor and the first exonic base at the acceptor. </template>Reference DNA is independent of profile filters. GC uses A/T/G/C bases in the displayed window.</footer>
  </section>
</template>

<style scoped>
.gene-sequence-card{display:flex;flex-direction:column;min-width:0;margin:0;overflow:hidden}.gene-sequence-card .panel-head{padding:22px 24px;min-height:106px}.gene-sequence-card h3{font-size:19px;font-weight:600;line-height:1.4;color:#204f78;margin:0;text-wrap:balance}.gene-sequence-card .panel-head p{font-size:13px;line-height:1.6;margin-top:6px;text-wrap:pretty}.sequence-raw{border-top:1px solid #e4eaf1;padding-top:10px}.sequence-raw summary{cursor:pointer;margin:0 20px 15px;font-size:13px;font-weight:500;color:#345776}.sequence-region{padding:16px 20px 0}.sequence-region label{display:flex;flex-direction:column;gap:8px;font-size:13px;color:#345776}.sequence-region select{width:100%;border:1px solid #d5dfeb;border-radius:6px;min-height:38px;padding:6px 9px;background:white;color:#345776;font:inherit}.sequence-meta{padding:16px 20px 10px;font-size:12px;color:#64748b;font-variant-numeric:tabular-nums;overflow-wrap:anywhere}.sequence-meta strong{font-size:13px;font-weight:500;color:#345776}.sequence-meta p{margin:5px 0 0;line-height:1.7}.sequence-actions{padding:0 20px 12px;display:flex;gap:8px;flex-wrap:wrap;align-items:center}.sequence-actions label{flex:1;min-width:175px}.sequence-actions select{width:100%;border:1px solid #d5dfeb;border-radius:6px;color:#345776;min-height:36px;padding:6px 9px;background:#fff;font:inherit;font-size:12px}.gene-sequence-card .button{font-size:12px;min-height:36px;padding:7px 12px}.sequence-composition{display:flex;flex-wrap:wrap;gap:8px 16px;font-size:12px;color:#64748b;padding:0 20px 14px;font-variant-numeric:tabular-nums}.sequence-composition>span{display:inline-flex;align-items:center;gap:6px}.sequence-composition b{font:600 14px ui-monospace,SFMono-Regular,Consolas,monospace}.sequence-column-labels{display:grid;grid-template-columns:86px 1fr;color:#6a7d90;font-size:11px;padding:0 20px 8px}.sequence-body{font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-size:14px;line-height:27px;max-height:310px;min-height:240px;overflow:auto;margin:0 20px 14px;scrollbar-width:thin;font-variant-numeric:tabular-nums}.sequence-line{display:grid;grid-template-columns:86px minmax(0,1fr);min-height:27px}.sequence-coordinate{font-size:11px;color:#7b8b9c;user-select:none}.sequence-groups{display:flex;gap:1ch;white-space:nowrap}.sequence-base-group{display:inline-flex}.sequence-base{width:1ch;text-align:center;font-weight:500;color:#64748b}.sequence-base:hover{background:#edf2f7;border-radius:2px}.base-A{color:#27815e}.base-T{color:#bc5260}.base-G{color:#997322}.base-C{color:#3978ad}.sequence-navigation{display:flex;align-items:center;justify-content:space-between;gap:8px;padding:2px 20px 16px;flex-wrap:wrap}.sequence-navigation form{display:flex;gap:6px;align-items:center}.sequence-navigation label{display:flex;align-items:center;gap:6px;font-size:12px;color:#64748b}.sequence-navigation input{width:85px;border:1px solid #d5dfeb;min-height:36px;border-radius:6px;padding:6px 8px;color:#345776;font:inherit;font-variant-numeric:tabular-nums}.gene-sequence-card button:disabled{opacity:.4;cursor:default}.sequence-foot{font-size:12px;line-height:1.7;padding:13px 20px;background:#f8fafc;border-top:1px solid #e4eaf1;color:#64748b;margin-top:auto;text-wrap:pretty}.sequence-action-status{margin:0;padding:0 20px 14px;font-size:12px;line-height:1.7;color:#345776;text-wrap:pretty}.sequence-loading{padding:28px 20px;min-height:410px}.sequence-loading>span{display:block;max-width:100%;height:20px;background:#edf2f7;border-radius:4px;margin:16px 0}.sequence-loading>span:nth-child(even){max-width:85%}.sequence-loading p{font-size:13px;color:#64748b;margin-top:24px}.sequence-error{padding:24px 20px;min-height:350px;color:#8d4242;font-size:13px}.sequence-error p{line-height:1.8;text-wrap:pretty}.gene-sequence-card :focus-visible{outline:2px solid #4979ad;outline-offset:2px}.sr-only{position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}@media(max-width:760px){.gene-sequence-card .panel-head{padding:20px;min-height:0}.gene-sequence-card h3{font-size:18px}.sequence-body{max-height:310px}.sequence-navigation{gap:10px}.sequence-navigation form{order:3;flex-basis:100%;justify-content:center}}
.gene-sequence-card>.panel-head{background:#204f78;border-bottom-color:#204f78}.gene-sequence-card>.panel-head h3{color:#fff}.gene-sequence-card>.panel-head p{color:#dce8f4}
</style>
