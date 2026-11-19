import { onBeforeUnmount, ref, watch } from 'vue'

// Run the short reveal when the chart enters the viewport, including on scroll.
export function useChartReveal({ once = false, threshold = 0.05 } = {}) {
  const host = ref(null), visible = ref(false)
  let observer
  watch(host, element => {
    observer?.disconnect()
    visible.value = false
    if (!element) return
    if (typeof IntersectionObserver === 'undefined') { visible.value = true; return }
    observer = new IntersectionObserver(entries => {
      const entered = entries[0].isIntersecting
      if (once && visible.value) return
      visible.value = entered
      if (once && entered) observer.unobserve(element)
    }, { threshold })
    observer.observe(element)
  }, { flush: 'post' })
  onBeforeUnmount(() => observer?.disconnect())
  return { host, visible }
}
