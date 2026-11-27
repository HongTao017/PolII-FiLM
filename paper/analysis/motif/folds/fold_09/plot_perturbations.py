
from pathlib import Path
import matplotlib.pyplot as plt
from final_figure_small import apply_style
from draw_reference_style import draw, read_inputs, sha256


def make_plots(root, config, plan, output=None):
    root=Path(root)
    output=Path(output) if output else root/'figures'
    assert not output.exists(), 'Do not overwrite figures'
    registered={'frozen_perturbation_plan.json':sha256(root/'frozen_perturbation_plan.json')}
    for item in plan:
        for name in ('mean_delta.npy','highlight.json','variant_effects.tsv','WT_mutants.fa'):
            p=root/'perturbations'/item['instance_id']/name
            registered[p.relative_to(root).as_posix()]=sha256(p)
    data=[(item,read_inputs(root,item,registered)) for item in plan]
    apply_style()
    plt.rcParams.update({'axes.linewidth':.8,'savefig.facecolor':'white'})
    for i,(item,(seq,heat,coords,highlight,chosen,_)) in enumerate(data):
        draw(item,seq,heat,coords,highlight,chosen,None,chr(65+i) if i<26 else str(i+1))
    output.mkdir(parents=True)
    for i,(item,(seq,heat,coords,highlight,chosen,_)) in enumerate(data):
        draw(item,seq,heat,coords,highlight,chosen,output,chr(65+i) if i<26 else str(i+1))
