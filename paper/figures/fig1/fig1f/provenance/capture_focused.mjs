import { chromium } from 'playwright'
import { mkdir,writeFile } from 'node:fs/promises'
const out='G:/序列分析数据/tmp/pdfs/nascentatlas_focus_20261025'
await mkdir(out,{recursive:true})
const browser=await chromium.launch({headless:true,executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',args:['--no-proxy-server']})
const page=await browser.newPage({viewport:{width:1440,height:1100},deviceScaleFactor:4,reducedMotion:'reduce'})
page.setDefaultTimeout(60000)
const origin='http://172.16.12.164:8088',figures=[],errors=[]
page.on('pageerror',e=>errors.push(e.message))
const style='.site-header,.skip-link{visibility:hidden!important} *,*::before,*::after{caret-color:transparent!important}'
async function shot(key,title,selector){
 const el=page.locator(selector);await el.waitFor({state:'attached'});await el.scrollIntoViewIfNeeded()
 await page.evaluate(()=>document.fonts.ready);await page.mouse.move(0,0)
 const image=out+'/'+key+'.png';await el.screenshot({path:image,style,animations:'disabled'})
 const box=await el.boundingBox();figures.push({key,title,image,url:page.url(),css_width:box.width,css_height:box.height,device_scale_factor:4})
 console.log(JSON.stringify({captured:key,width:box.width,height:box.height}))
}
async function waitGenome(){
 await page.locator('.gene-genome-browser .igv-track-label').nth(3).waitFor({state:'attached'})
 await page.waitForFunction(()=>[...document.querySelectorAll('.gene-genome-browser .igv-spinner')].every(e=>getComputedStyle(e).display==='none'||!e.getBoundingClientRect().width),{},{timeout:30000})
 await page.waitForTimeout(1000)
}
try{
 await page.goto(origin+'/#/samples');await page.locator('.browse-profile-table tbody tr').nth(5).waitFor({state:'attached'})
 await page.locator('.browse-layout').evaluate(el=>{const rows=el.querySelectorAll('.browse-profile-table tbody tr'),height=rows[5].getBoundingClientRect().bottom-el.getBoundingClientRect().top;el.style.height=height+'px';el.style.overflow='hidden'})
 await shot('01_browse_main','Sample browse: filters, column labels and six representative rows','.browse-layout')
 await page.goto(origin+'/#/genes/SAMD11');await page.locator('.gene-pi-point').first().waitFor({state:'attached'});await waitGenome()
 await page.locator('.gene-pi-layout').evaluate(root=>{for(const region of root.querySelectorAll('.gene-pi-rows')){if(region.scrollHeight<=region.clientHeight)continue;const row=region.querySelector('[data-group]'),axis=region.querySelector('.gene-pi-axis-row');const rowHeight=row.getBoundingClientRect().height,axisHeight=axis.getBoundingClientRect().height,height=Math.floor((region.clientHeight-axisHeight)/rowHeight)*rowHeight+axisHeight;Object.assign(region.style,{height:height+'px',maxHeight:height+'px',minHeight:'0',flex:'none'})}})
 await shot('02_pi_charts_main','SAMD11: PI chart panels','.gene-pi-layout')
 await page.locator('.feature-charts').evaluate(e=>Object.assign(e.style,{maxHeight:'none',overflow:'visible'}))
 await shot('03_feature_charts_main','SAMD11: extracted feature charts, all 30 features','.feature-charts')
 await shot('05_sequence_statistics_main','SAMD11: ATGC, GC, CpG and regional GC profile','.region-stats')
 await page.goto(origin+'/#/genes/A1BG');await waitGenome()
 await page.locator('.gene-genome-browser').getByText('Crosshairs',{exact:true}).click()
 await shot('04_igv_tracks_main','A1BG: genome coordinate ruler, annotation and three sequencing tracks','.gene-genome-browser .genome-scroll')
 await page.goto(origin+'/#/samples/NP-7B10E7FF77');await page.locator('.profile-pair-heatmap td button').first().waitFor({state:'attached'})
 await page.getByLabel('Heatmap comparison sample',{exact:true}).selectOption('NP-410E006606');await page.waitForFunction(()=>document.querySelectorAll('.profile-pair-heatmap td button').length===4)
 await shot('06_correlation_main','Sample correlation matrix, scale and pair details','.pair-heatmap-layout')
 if(errors.length)throw Error(JSON.stringify(errors))
 figures.sort((a,b)=>a.key.localeCompare(b.key))
 await writeFile(out+'/captures.json',JSON.stringify({captured_on:'2026-10-25',origin,figures},null,2))
}finally{await browser.close()}
