import { ref } from 'vue'
import { getCatalog } from '../api/client.js'
const catalog = ref(null)
const loading = ref(false)
const error = ref('')
let inflight
export function useCatalog() {
  async function reload() {
    if (inflight) return inflight
    loading.value = true; error.value = ''
    inflight = getCatalog().then(value => { catalog.value = value }).catch(e => { error.value = e.message }).finally(() => { loading.value = false; inflight = null })
    return inflight
  }
  if (!catalog.value && !inflight) reload()
  return { catalog, loading, error, reload }
}
