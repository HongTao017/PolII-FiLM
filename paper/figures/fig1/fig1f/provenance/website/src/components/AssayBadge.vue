<script setup>
import { assayClass } from '../utils/data.js'
defineProps({ assay: String, to: [String, Object] })
</script>
<template>
  <RouterLink v-if="to" :to="to" class="assay-badge interactive-tag" :class="assayClass(assay)" :aria-label="`Browse ${assay} profiles`"><span />{{ assay }}</RouterLink>
  <span v-else class="assay-badge" :class="assayClass(assay)"><span />{{ assay }}</span>
</template>
