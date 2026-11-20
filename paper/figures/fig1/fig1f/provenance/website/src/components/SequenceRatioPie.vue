<script setup>
import { computed, ref } from 'vue'
const props = defineProps({ label: String, percent: { type: Number, default: null }, precision: { type: Number, default: 1 }, primary: String, remainder: String, color: String, animated: Boolean, numerator: Number, denominator: Number, unit: String })
const activeId = ref('')
const display = computed(() => props.percent === null ? 'Not available' : props.percent.toFixed(props.precision) + '%')
function arc(start, fraction) {
  if (fraction >= 1) return 'M80 18 A62 62 0 1 1 80 142 A62 62 0 1 1 80 18 Z'
  const a = start * Math.PI * 2 - Math.PI / 2, b = (start + fraction) * Math.PI * 2 - Math.PI / 2
  return `M80 80 L${80 + 62*Math.cos(a)} ${80 + 62*Math.sin(a)} A62 62 0 ${fraction > .5 ? 1 : 0} 1 ${80 + 62*Math.cos(b)} ${80 + 62*Math.sin(b)} Z`
}
const sectors = computed(() => {
  if (props.percent === null) return []
  const fraction = props.percent / 100
  return [{ id: 'primary', label: props.primary, value: props.percent, count: props.numerator, start: 0, color: props.color }, { id: 'remainder', label: props.remainder, value: 100-props.percent, count: props.denominator-props.numerator, start: fraction, color: '#d7e3f0' }].filter(s => s.value > 0).map(s => {
    const angle = (s.start + s.value / 200) * Math.PI * 2 - Math.PI/2
    return { ...s, path: arc(s.start, s.value / 100), offset: `translate(${Math.cos(angle)*4}px, ${Math.sin(angle)*4}px) scale(1.035)` }
  })
})
const active = computed(() => sectors.value.find(s => s.id === activeId.value))
</script>
<template>
  <figure class="ratio-pie" :class="{'is-revealed':animated}" @mouseleave="activeId=''">
    <figcaption>{{ label }}</figcaption>
    <div class="ratio-pie-content">
      <svg viewBox="0 0 160 160" role="img" :aria-label="`${label}: ${display}. ${primary} versus ${remainder}.`" :data-percent="percent ?? ''">
        <g class="pie-motion"><circle v-if="percent === null" cx="80" cy="80" r="62" fill="#e3ebf4"/><path v-for="sector in sectors" :key="sector.id" :d="sector.path" :fill="sector.color" class="pie-sector" :class="{'is-active':activeId===sector.id,'is-muted':activeId && activeId!==sector.id}" :style="{transform:activeId===sector.id ? sector.offset : undefined}" :data-sector="sector.id" @mouseenter="activeId=sector.id" @mouseleave="activeId=''" @click="activeId=sector.id"><title>{{sector.label}}: {{sector.value.toFixed(precision)}}%</title></path></g>
      </svg>
      <div class="ratio-pie-values"><strong class="ratio-percent"><slot name="value">{{ display }}</slot></strong><div class="pie-legend"><button v-for="sector in sectors" :key="sector.id" type="button" :class="{'legend-active':activeId===sector.id}" :aria-label="`${label}, ${sector.label}: ${sector.value.toFixed(precision)} percent, ${sector.count} ${unit}`" @mouseenter="activeId=sector.id" @mouseleave="activeId=''" @focus="activeId=sector.id" @blur="activeId=''" @click="activeId=sector.id"><i :style="{background:sector.color}"/>{{ sector.label }}</button></div><p><slot/></p></div>
    </div>
    <div v-if="active" class="pie-tooltip" role="tooltip"><strong>{{active.label}}</strong><span>{{active.count.toLocaleString('en-US')}} {{unit}} · {{active.value.toFixed(precision)}}%</span></div>
  </figure>
</template>
<style scoped>
.ratio-pie{margin:0;min-width:0;position:relative}.ratio-pie figcaption{font-size:14px;font-weight:600;color:#345776}.ratio-pie-content{display:flex;align-items:center;gap:14px}.ratio-pie svg{display:block;width:150px;max-width:45%;flex-shrink:0;overflow:visible}.ratio-pie-values{min-width:0}.ratio-percent{font-size:27px;font-weight:600;color:#23415f;font-variant-numeric:tabular-nums}.pie-legend{display:flex;gap:6px 10px;flex-wrap:wrap;font-size:11px;line-height:1.7;color:#64748b;margin:8px 0}.pie-legend button{display:inline-flex;align-items:center;gap:6px;font:inherit;color:inherit;border:0;background:transparent;padding:2px 4px;border-radius:3px;text-align:left;cursor:pointer}.pie-legend button.legend-active{background:#eef3f9;color:#23415f}.pie-legend button:focus-visible{outline:2px solid #417dc6;outline-offset:2px}.pie-legend i{width:9px;height:9px;border-radius:2px;flex-shrink:0}.ratio-pie p{font-size:11px;color:#64748b;line-height:1.7;margin:6px 0;text-wrap:pretty}.pie-motion{transform-origin:80px 80px;opacity:0}.is-revealed .pie-motion{animation:sequence-pie-reveal 600ms ease-out both}.pie-sector{transform-origin:80px 80px;transition:transform 180ms ease-out,opacity 180ms ease-out;stroke:white;stroke-width:1;cursor:pointer}.pie-sector.is-active{stroke-width:2}.pie-sector.is-muted{opacity:.55}.pie-tooltip{position:absolute;top:30px;left:8px;z-index:4;background:#fff;border:1px solid #d5dfeb;border-radius:5px;box-shadow:0 3px 10px #0f172a20;padding:8px 12px;pointer-events:none;font-size:12px;line-height:1.7;color:#345776;max-width:calc(100% - 16px)}.pie-tooltip span{display:block;font-variant-numeric:tabular-nums}@keyframes sequence-pie-reveal{from{opacity:0;transform:scale(0) rotate(-60deg)}to{opacity:1;transform:scale(1) rotate(0)}}@media(prefers-reduced-motion:reduce){.pie-motion{opacity:1;transform:none}.is-revealed .pie-motion{animation:none}.pie-sector{transition:none;transform:none!important}}@media(max-width:480px){.ratio-pie svg{width:130px}.ratio-percent{font-size:24px}.ratio-pie-content{gap:10px}}
</style>
