import { createRouter, createWebHashHistory } from 'vue-router'

const routes = [
  { path: '/', name: 'home', component: () => import('./pages/Home.vue'), meta: { title: 'Overview' } },
  { path: '/search', name: 'search', component: () => import('./pages/Search.vue'), meta: { title: 'Unified search' } },
  { path: '/studies', name: 'studies', redirect: to => { const { author, ...query } = to.query; return { path: '/samples', query: { ...query, ...(!query.q && author ? { q: author } : {}), view: 'study' } } } },
  { path: '/studies/:id', name: 'study', component: () => import('./pages/StudyDetail.vue'), meta: { title: 'Study details' } },
  { path: '/samples', name: 'samples', component: () => import('./pages/Samples.vue'), meta: { title: 'Profile catalog' } },
  { path: '/samples/:id', name: 'sample', component: () => import('./pages/SampleDetail.vue'), meta: { title: 'Profile details' } },
  { path: '/genes', name: 'genes', component: () => import('./pages/Gene.vue'), meta: { title: 'Genes' } },
  { path: '/genes/:symbol', name: 'gene', component: () => import('./pages/Gene.vue'), meta: { title: 'Gene details' } },
  { path: '/statistics', name: 'statistics', component: () => import('./pages/Statistics.vue'), meta: { title: 'Collection statistics' } },
  { path: '/download', name: 'download', component: () => import('./pages/Download.vue'), meta: { title: 'Download center' } },
  { path: '/help', name: 'help', component: () => import('./pages/Help.vue'), meta: { title: 'About and methods' } },
  { path: '/:pathMatch(.*)*', name: 'notfound', component: () => import('./pages/NotFound.vue'), meta: { title: 'Page not found' } },
]
const router = createRouter({ history: createWebHashHistory(), routes, scrollBehavior: (to) => to.hash ? ({ el: to.hash, top: window.innerWidth <= 760 ? 92 : 108, behavior: window.matchMedia('(prefers-reduced-motion: reduce)').matches ? 'auto' : 'smooth' }) : ({ top: 0 }) })
router.afterEach((to) => { document.title = `${to.meta.title} · NascentAtlas` })
export default router
