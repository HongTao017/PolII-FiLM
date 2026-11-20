<script setup>
import PageTitle from '../components/PageTitle.vue'
import EmptyState from '../components/EmptyState.vue'
</script>
<template><PageTitle eyebrow="NAVIGATION" title="Page not found" description="Return to the database to continue exploring."/><EmptyState title="Page not found" description="This address does not match a page on this site." to="/" action="Back to home"/></template>
