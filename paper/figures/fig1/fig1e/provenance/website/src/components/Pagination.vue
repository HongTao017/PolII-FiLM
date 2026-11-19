<script setup>
import { computed } from 'vue'
import { ChevronLeft, ChevronRight } from '@lucide/vue'
const props = defineProps({ modelValue: Number, total: Number, size: { default: 10 } })
defineEmits(['update:modelValue'])
const pages = computed(() => Math.max(1, Math.ceil(props.total / props.size)))
</script>
<template><div class="pagination"><span>{{ total }} records · Page {{ modelValue }} of {{ pages }}</span><div><button aria-label="Previous page" :disabled="modelValue <= 1" @click="$emit('update:modelValue', modelValue - 1)"><ChevronLeft :size="16" /></button><span>{{ modelValue }}</span><button aria-label="Next page" :disabled="modelValue >= pages" @click="$emit('update:modelValue', modelValue + 1)"><ChevronRight :size="16" /></button></div></div></template>
