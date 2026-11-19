<script setup>
import { computed } from 'vue'
import figures from '../config/about-figures.json'
const props = defineProps({ name: String, alt: String, caption: String })
const figure = computed(() => figures[props.name])
</script>
<template><figure v-if="figure" class="about-figure"><a :href="figure.path" target="_blank" rel="noopener" :aria-label="`Open full-size example: ${alt}`"><img :src="figure.path" :alt="alt" :width="figure.width" :height="figure.height" loading="lazy" decoding="async"/></a><figcaption>{{ caption }} <span>Open image to enlarge.</span></figcaption></figure></template>
<style scoped>
.about-figure{margin:22px 0 8px}.about-figure>a{display:block;border:1px solid #dce5ef;border-radius:6px;overflow:hidden;background:#f8fafc}.about-figure img{display:block;width:100%;height:auto;max-height:480px;object-fit:contain}.about-figure figcaption{font-size:12px;color:#64748b;line-height:1.8;margin:10px 0;text-wrap:pretty}.about-figure figcaption span{color:#3e77ab}.about-figure>a:focus-visible{outline:2px solid #417dc6;outline-offset:3px}
</style>
