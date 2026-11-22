"""Journal-size figures, generated only from exported numerical results."""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, PercentFormatter

ROOT=Path(__file__).resolve().parents[2]
C=json.loads((ROOT/"config/design.json").read_text())
metadata=json.loads((ROOT/"output/logs/run_metadata.json").read_text())
assert metadata["config_sha256"] == hashlib.sha256((ROOT/"config/design.json").read_bytes()).hexdigest(), "Simulation tables do not match config; run python run.py first."
TAB=ROOT/"output/tables"
FIG=ROOT/"output/figures"
FIG.mkdir(parents=True,exist_ok=True)
BLUE="#245A81"
ORANGE="#BB642B"
MAGENTA="#953D61"
GRAY="#66717A"
TAU=C["treatment_effect"]
SCALE=C["score_scale"]
ALPHA=C["quota_numerator"]/C["quota_denominator"]
plt.rcParams.update({
    "font.family":"DejaVu Serif", "font.size":9, "axes.titlesize":9.5,
    "axes.labelsize":9, "legend.fontsize":8, "xtick.labelsize":8,
    "ytick.labelsize":8, "axes.spines.top":False, "axes.spines.right":False,
    "axes.edgecolor":"#65717B", "axes.linewidth":.65, "lines.linewidth":1.5,
    "grid.color":"#DDE2E6", "grid.linewidth":.5, "pdf.fonttype":42,
    "ps.fonttype":42, "savefig.facecolor":"white", "svg.hashsalt":"scholarship-quota-v2",
    "svg.fonttype":"none"
})


def axes_pair():
    f,axs=plt.subplots(1,2,figsize=(7.5,3.15))
    f.subplots_adjust(left=.085,right=.985,bottom=.2,top=.87,wspace=.38)
    for a in axs:
        a.grid(axis="y",alpha=.7)
        a.set_axisbelow(True)
    return f,axs


def logx(ax, ticks):
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.minorticks_off()


def save(fig,name):
    fig.savefig(FIG/(name+".pdf"),metadata={"CreationDate":None,"ModDate":None,
        "Title":name,"Author":"Simulation replication package"})
    fig.savefig(FIG/(name+".png"),dpi=300)
    fig.savefig(FIG/(name+".svg"),metadata={"Date":None})
    plt.close(fig)


def cutoff_figure():
    d=pd.read_csv(TAB/"cutoff_and_spacing.csv")
    f,ax=axes_pair()
    f.subplots_adjust(top=.77)
    f.suptitle(f"Scores range from 0 to {SCALE:g}  |  Top {100*ALPHA:g}% receive scholarships",y=.985,fontsize=10.5,
               fontweight="bold",color=BLUE)
    ax[0].plot(d.m,SCALE*d.cutoff_rmse_exact,color=GRAY,label="Exact uniform-score benchmark")
    ax[0].scatter(d.m,SCALE*d.cutoff_rmse,s=23,color=BLUE,label="Monte Carlo RMSE",zorder=3)
    ax[0].set(title="(a) Cutoff uncertainty",ylabel="Cutoff RMSE (score points)",xlabel="Applicants per programme")
    ax[0].set_yscale("log")
    ax[0].legend(frameon=False,loc="upper right",fontsize=7.4)
    ax[1].plot(d.m,SCALE*d.gap_mean_exact,color=GRAY,label=f"Exact mean: {SCALE:g} / (size + 1)")
    ax[1].scatter(d.m,SCALE*d.gap_mean,s=23,color=ORANGE,label="Monte Carlo mean",zorder=3)
    ax[1].set(title="(b) Gap between adjacent ranks",ylabel="Score gap (points)",xlabel="Applicants per programme")
    ax[1].set_yscale("log")
    ax[1].legend(frameon=False,loc="upper right",fontsize=7.4)
    for a in ax: logx(a,[10,100,1000,10000])
    save(f,"fig01_cutoff_and_spacing")


def rd_figure():
    d=pd.read_csv(TAB/"large_cluster_rd.csv").query("dgp == 'linear'")
    f,ax=axes_pair()
    ax[0].axhline(TAU,color=MAGENTA,ls="--",label=f"True scholarship effect: {TAU:g} points")
    ax[0].fill_between(d.m,d.q025,d.q975,color=BLUE,alpha=.13,label="Central 95% Monte Carlo range")
    ax[0].plot(d.m,d.mc_mean,"o-",color=BLUE,ms=3.5,label="Monte Carlo mean")
    ax[0].set(title="(a) Local linear RD estimates",ylabel="Estimated effect (outcome points)",xlabel="Applicants per programme")
    ax[0].legend(frameon=False,loc="upper right",fontsize=7.2)
    ax[1].axhline(.95,color=GRAY,ls="--",label="Nominal coverage: 95%")
    ax[1].errorbar(d.m,d.coverage,yerr=1.96*d.coverage_mcse,fmt="o-",color=BLUE,
                   ms=3.5,capsize=2.5,label="Coverage with Monte Carlo error")
    ax[1].set(title="(b) Confidence interval coverage",ylabel="Coverage",xlabel="Applicants per programme")
    ax[1].set_ylim(min(.89,float((d.coverage-2.5*d.coverage_mcse).min())-.01),.99)
    ax[1].yaxis.set_major_formatter(PercentFormatter(1,decimals=0))
    ax[1].legend(frameon=False,loc="lower right",fontsize=7.2)
    for a in ax: logx(a,[100,1000,10000])
    save(f,"fig02_large_programme_rd")


def small_figure():
    d=pd.read_csv(TAB/"small_cluster_contrast.csv")
    f,ax=axes_pair()
    target=d.target.iloc[0]
    ax[0].fill_between(d.G,d.q025,d.q975,color=ORANGE,alpha=.16,label="Central 95% Monte Carlo range")
    ax[0].plot(d.G,d.mc_mean,"o-",color=ORANGE,ms=3.5,label="Monte Carlo mean")
    ax[0].axhline(target,color=GRAY,ls="-.",label=f"Outcome difference target: {target:.3f}")
    ax[0].axhline(TAU,color=MAGENTA,ls="--",label=f"True scholarship effect: {TAU:g}")
    ax[0].set(title="(a) More precise, same baseline gap",ylabel="Adjacent-rank outcome difference",xlabel="Independent 10-person programmes")
    ax[0].legend(frameon=False,loc="upper right",fontsize=7.2)
    ymin=min(float(d.q025.min()),TAU)-.35
    ymax=float(d.q975.max())+1.6
    ax[0].set_ylim(ymin,ymax)
    ax[1].axhline(.95,color=GRAY,ls="--",lw=1,label="Nominal coverage: 95%")
    ax[1].errorbar(d.G,d.coverage,yerr=1.96*d.coverage_mcse,fmt="o-",ms=3.3,capsize=2,
                   color=ORANGE,label="Coverage of outcome difference")
    ax[1].errorbar(d.G,d.causal_coverage,yerr=1.96*d.causal_coverage_mcse,fmt="s-",ms=3.3,capsize=2,
                   color=MAGENTA,label="Coverage of true scholarship effect")
    ax[1].set(title="(b) Coverage depends on the target",ylabel="Coverage",xlabel="Independent 10-person programmes")
    ax[1].set_ylim(-.035,1.06)
    ax[1].yaxis.set_major_formatter(PercentFormatter(1,decimals=0))
    ax[1].legend(frameon=False,loc="center right",fontsize=7.2)
    for a in ax: logx(a,[20,100,500,2000])
    save(f,"fig03_small_programmes")


def main():
    cutoff_figure();rd_figure();small_figure()
    print("Figures: 3 main figures; PDF/PNG/SVG",flush=True)


if __name__=="__main__":main()
