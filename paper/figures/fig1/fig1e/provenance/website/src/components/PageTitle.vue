<script setup>
defineProps({ eyebrow: { default: 'NAScent transcription atlas' }, title: String, description: String })
</script>
<template><header class="page-title"><div class="container page-title-inner"><div class="page-title-copy"><p class="eyebrow">{{ eyebrow }}</p><h1>{{ title }}</h1><p v-if="description" class="subtitle">{{ description }}</p></div><div v-if="$slots.default" class="page-actions"><slot /></div></div></header></template>
