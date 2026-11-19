<script setup>
import { Database, ArrowRight } from '@lucide/vue'
defineProps({ title: { default: 'No results available' }, description: String, to: String, action: { default: 'Browse catalog' } })
</script>
<template><div class="empty-state"><Database :size="28" /><h3>{{ title }}</h3><p>{{ description }}</p><RouterLink v-if="to" :to="to" class="text-link">{{ action }} <ArrowRight :size="15" /></RouterLink><slot /></div></template>
