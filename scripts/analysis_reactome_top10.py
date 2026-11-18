"""Updated Reactome bubble plots, ranked by NEASE score after BH P <= 0.05."""
from pathlib import Path
import re
import textwrap
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import MultipleLocator, MaxNLocator

large_text='--large-text' in sys.argv
plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42,'ps.fonttype':42})
gene_col='Spliced genes (number of interactions affecting the pathway)'

def genes(value):
    if pd.isna(value):
        return []
    return sorted({re.sub(r'\s*\([^)]*\)','',v).strip() for v in value.split(',') if v.strip()})

for cohort in ['stroma', 'epithelium']:
    run = Path('results') / cohort / 'nease'
    (run / 'figures').mkdir(parents=True, exist_ok=True)
    epithelial = cohort == 'epithelium'
    frame=pd.read_csv(run/'results/Reactome_Python.tsv',sep='\t')
    p=frame.loc[frame['adj p_value']<=.05].sort_values(['Nease score','Pathway ID'],ascending=[False,True]).head(10).copy()
    assert len(p)==10 and p['Pathway ID'].is_unique and p['Nease score'].is_monotonic_decreasing
    p['Splicing-affected genes']=p[gene_col].map(lambda v:len(genes(v)))
    assert (p['Splicing-affected genes']>0).all() and (p['adj p_value']>0).all()
    p['-log10 adjusted P']=-np.log10(p['adj p_value'])
    p.insert(0,'Rank',range(1,11))
    p.to_csv(run/'results/Reactome_top10_for_figure.tsv',sep='\t',index=False)
    fig=plt.figure(figsize=(8.5,6.5 if large_text else 6.0),facecolor='white')
    ax=fig.add_axes([.515,.105,.25,.84])
    cmap=LinearSegmentedColormap.from_list('nease',['#FFFFFF','#CC79A7' if epithelial else '#FAA62B'])
    norm=Normalize(0,float(p['-log10 adjusted P'].max()))
    scale=40 if epithelial else 16
    dots=ax.scatter(p['Nease score'],range(10),s=p['Splicing-affected genes']*scale,
                    c=p['-log10 adjusted P'],cmap=cmap,norm=norm,edgecolors='white',linewidths=.7,zorder=3)
    labels=['\n'.join(textwrap.wrap(v,34 if large_text else 43,break_long_words=False,break_on_hyphens=False)) for v in p['Pathway name']]
    ax.set_yticks(range(10),labels)
    for label in ax.get_yticklabels():
        label.set_linespacing(1.03)
    ax.set_ylim(9.6,-.6)
    lo,hi=p['Nease score'].min(),p['Nease score'].max()
    span=hi-lo
    ax.set_xlim(lo-span*.15,hi+span*.20)
    ax.set_xlabel('NEASE score',fontsize=18 if large_text else 14,labelpad=8)
    ax.tick_params(axis='both',labelsize=16 if large_text else 12,length=3,pad=5)
    ax.xaxis.set_major_locator(MultipleLocator(2 if epithelial else 10))
    ax.xaxis.set_minor_locator(MultipleLocator(1 if epithelial else 5))
    ax.grid(True,axis='both',which='major',linestyle=(0,(1,3)),linewidth=.7,color='#999999')
    ax.grid(True,axis='x',which='minor',linestyle=(0,(1,3)),linewidth=.5,color='#BBBBBB')
    for spine in ax.spines.values():
        spine.set_color('#555555')
        spine.set_linewidth(.8)
    sizes=[2,7,14] if epithelial else [10,20,40]
    handles=[ax.scatter([],[],s=value*scale,color='black',edgecolors='none') for value in sizes]
    fig.legend(handles,[str(v) for v in sizes],title='Splicing-affected\ngenes',loc='upper left',
               bbox_to_anchor=(.79,.94),frameon=False,fontsize=14 if large_text else 11,title_fontsize=15 if large_text else 12,
               labelspacing=1.05,handletextpad=.9,borderpad=0)
    colorax=fig.add_axes([.80,.19,.025,.29])
    cb=fig.colorbar(dots,cax=colorax)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=14 if large_text else 11,length=0,pad=6)
    cb.locator=MaxNLocator(nbins=5)
    cb.update_ticks()
    fig.text(.80,.51,r'$-\mathrm{Log}_{10}$(P-adjust)',fontsize=15 if large_text else 12,ha='left')
    for ext in ['png','pdf']:
        suffix='_large_text' if large_text else ''
        target=run/'figures'/f'NEASE_Reactome_top10_updated{suffix}.{ext}'
        fig.savefig(target,dpi=300,facecolor='white',bbox_inches='tight',pad_inches=.08)
        print(target)
    plt.close(fig)
