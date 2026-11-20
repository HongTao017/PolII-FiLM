<script setup>
import { X } from '@lucide/vue'
defineProps({ items: { type: Array, default: () => [] } })
defineEmits(['remove'])
</script>
<template>
  <TransitionGroup name="filter-chip" tag="div" class="active-filters" aria-label="Active filters">
    <button v-for="item in items" :key="item.key" type="button" class="active-filter-chip" :aria-label="`Remove ${item.label.toLowerCase()} filter ${item.value}`" @click="$emit('remove', item.key)">
      <span>{{ item.label }}</span><b :title="item.value">{{ item.value }}</b><X :size="12" aria-hidden="true" />
    </button>
  </TransitionGroup>
</template>
