"""Real saturation-mutagenesis plots, never synthetic effects or relabeled RNA."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator
import logomaker
import numpy as np
import pandas as pd

from deep_lift_core import REGIONS, OFFSETS, digest, write_json


def make_plots(root, config, plan):
    root=Path(root)
    out=root/'figures'
    out.mkdir()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'ps.fonttype':42})
    hypothetical=np.mean([np.load(root/'attributions'/f'seed_{s}'/'hypothetical.npy') for s in config['seeds']],0)
    metadata=pd.read_csv(Path(config['model_root'])/'data/genes_repartitioned.tsv',sep='\t')
    source=Path(config['source_root'])/'data/sequence_prepared'
    files=[]
    for item in plan:
        folder=root/'perturbations'/item['instance_id']
        heat=np.load(folder/'mean_delta.npy')
        sd=np.load(folder/'seed_sd_delta.npy')
        highlight=json.loads((folder/'highlight.json').read_text())
        left,right=item['window_start'],item['window_end']
        j=item['region_index']
        onehot=np.array(np.load(source/f'{item["region"]}.npy',mmap_mode='r')[item['gene_index']]).T
        observed=hypothetical[item['gene_row'],:,OFFSETS[j]+left:OFFSETS[j]+right]*onehot[:,left:right]
        assert np.isfinite(heat).all() and np.isfinite(observed).all()
        fig=plt.figure(figsize=(11,5.7))
        grid=fig.add_gridspec(3,2,width_ratios=(1,.025),height_ratios=(1.15,1.25,1),hspace=.15,wspace=.025)
        ax_logo=fig.add_subplot(grid[0,0])
        ax_line=fig.add_subplot(grid[1,0],sharex=ax_logo)
        ax_heat=fig.add_subplot(grid[2,0],sharex=ax_logo)
        cax=fig.add_subplot(grid[2,1])
        logomaker.Logo(pd.DataFrame(observed.T,columns=list('ACGT')),ax=ax_logo,
            color_scheme='classic',vpad=.05,fade_below=.5)
        ax_logo.axhline(0,color='.6',lw=.5)
        ax_logo.set_ylabel('DeepLIFT attribution\n(signed Y_PC1)',fontsize=9)
        alternatives=' / '.join(x['TF'] for x in item['all_TF_matches'][:4])
        if len(item['all_TF_matches'])>4:
            alternatives+=' / ...'
        fig.suptitle(f'{item["gene"]} | {item["region"].replace("_sequence", "")} | candidate motif: {alternatives}',fontsize=12,y=.985)

        increase=np.maximum(heat.max(0),0)
        decrease=np.maximum(-heat.min(0),0)
        x=np.arange(right-left)
        ax_line.plot(x,increase,color='#E43838',lw=1.2,label='Prediction increase (gain)')
        ax_line.plot(x,decrease,color='#286AC7',lw=1.2,label='Prediction decrease (loss magnitude)')
        ax_line.set_ylabel('Magnitude of prediction change\n(signed Y_PC1 units)',fontsize=9)
        ax_line.set_ylim(bottom=0)
        ax_line.legend(frameon=False,loc='upper right',fontsize=8)
        limit=max(float(np.abs(heat).max()),1e-6)
        im=ax_heat.imshow(heat,aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit,
            extent=(-.5,len(x)-.5,3.5,-.5),interpolation='nearest')
        ax_heat.set_yticks(range(4),labels=list('ACGT'))
        ax_heat.set_ylabel('Alternative base\n(input orientation)',fontsize=9)
        colorbar=fig.colorbar(im,cax=cax)
        colorbar.set_label('Mutant − WT\nΔ predicted signed Y_PC1',fontsize=9)
        for ax in (ax_logo,ax_line,ax_heat):
            ax.set_xlim(-.5,len(x)-.5)
        for ax in (ax_logo,ax_line):
            ax.tick_params(axis='x',bottom=False,labelbottom=False)
            ax.spines[['top','right']].set_visible(False)
        a,b=item['motif_start']-left,item['motif_end']-left
        for ax in (ax_logo,ax_heat):
            ymin,ymax=ax.get_ylim()
            ax.add_patch(Rectangle((a-.5,min(ymin,ymax)),b-a,abs(ymax-ymin),fill=False,edgecolor='#FF3991',lw=1))
        pos=highlight['position']-left
        annotation=f'{highlight["reference"]}>{highlight["alternative"]}: motif '
        annotation+='same-site loss' if highlight['FIMO_same_site_loss'] else 'score decrease'
        if highlight['FIMO_same_site_loss'] and highlight['overlapping_target_match_remaining']:
            annotation+='; overlapping match remains'
        ax_logo.text(.5,1.07,annotation,ha='center',va='bottom',transform=ax_logo.transAxes,fontsize=8,
            bbox={'facecolor':'white','edgecolor':'.4','pad':2})
        for ax in (ax_logo,ax_line,ax_heat):
            ax.axvline(pos,color='#FF3991',ls=':',lw=1)
        meta=metadata.iloc[item['gene_index']]
        requested_left=meta[f'{item["region"]}_requested_start0']
        requested_right=meta[f'{item["region"]}_requested_end0']
        def coordinate(value,_):
            oriented=left+value


            bp=requested_left+oriented+1 if meta.strand=='+' else requested_right-oriented
            return f'{int(round(bp)):,}'
        ax_heat.xaxis.set_major_locator(MaxNLocator(5,integer=True))
        ax_heat.xaxis.set_major_formatter(FuncFormatter(coordinate))
        ax_heat.set_xlabel(f'{item["chrom"]}, genomic position (1-based), strand {meta.strand}')
        fig.text(.08,.025,'Five frozen-seed prediction average; ChIP inputs unchanged. Motif annotation is sequence-based, not measured TF binding.',fontsize=8)
        fig.subplots_adjust(left=.09,right=.91,top=.84,bottom=.14)
        basename=f'{item["instance_id"]}_{item["gene"]}_{item["TF"]}'
        for ext in ('pdf','png'):
            path=out/f'{basename}.{ext}'
            fig.savefig(path,dpi=200)
            files.append(str(path.relative_to(root)))
        plt.close(fig)
    write_json(out/'plot_manifest.json',dict(files_sha256={name:digest(root/name) for name in files},
        heatmap='four alternative bases, mutant minus WT signed Y_PC1; WT cells exactly zero',
        curve_gain='max(0,max alternative delta) per position',
        curve_loss='max(0,-min alternative delta) per position; magnitude, not negative-valued curve',
        motif_labels='FIMO p<=1e-4 same frozen JASPAR matrix rescanned on WT and every mutant; alternate TF annotation retained',
        no_fold_numbers_visible=True,no_fabricated_values=True))
