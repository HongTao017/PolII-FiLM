<script setup>
defineProps({ loading: Boolean, error: String })
defineEmits(['retry'])
</script>
<template><div v-if="loading" class="loading-state" role="status"><span class="spinner" />Loading the data catalog…</div><div v-else-if="error" class="error-state" role="alert"><strong>Unable to load data</strong><p>{{ error }}</p><button class="button secondary" @click="$emit('retry')">Reload</button></div><slot v-else /></template>
