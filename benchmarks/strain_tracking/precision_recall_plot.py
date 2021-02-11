from pathlib import Path

import numpy
import pandas
import matplotlib.pyplot as plt
from matplotlib import ticker
from sklearn import metrics
import seaborn

seaborn.set(context="paper")

#%% Load eval data

BASE_PATH = Path(__file__).parent

strainge_df = pandas.read_csv(BASE_PATH / 'all_vs_all/results_strainge_c0.5.tsv', sep='\t')
# midas_df = pandas.read_csv(BASE_PATH / 'all_vs_all/results_midas_default.tsv', sep='\t')
midas_df = pandas.read_csv(BASE_PATH / 'all_vs_all/results_midas_all_snvs.tsv', sep='\t')
# strainphlan_df = pandas.read_csv(BASE_PATH / 'all_vs_all/results_strainphlan.tsv', sep='\t')
strainphlan_df = pandas.read_csv(BASE_PATH / 'all_vs_all/results_strainphlan_relaxed.tsv', sep='\t')

eval_df = pandas.concat([strainge_df, midas_df, strainphlan_df])

eval_df['timepoint1'] = eval_df['s1'].map(lambda s: int(s[s.rfind('.')+2:]))
eval_df['timepoint2'] = eval_df['s2'].map(lambda s: int(s[s.rfind('.')+2:]))

eval_df = eval_df.set_index(['tool', 'cov'])

COMPLETED = {}
for ix in eval_df.index.unique():
    num_completed = numpy.count_nonzero(eval_df.loc[ix, 'distance'].notnull())
    COMPLETED[ix] = num_completed / len(eval_df.loc[ix])


#%%

eval_finite = eval_df[eval_df['distance'].notnull()].copy()
color_palette = seaborn.color_palette()

TOOLS = ["strainge", "strainphlan", "midas"]

TEST_IX = [
    # Single strain samples only, no closely related strains
    (eval_finite['timepoint1'] <= 4) & (eval_finite['timepoint2'] <= 4),

    # With mixes
    (eval_finite['timepoint1'] <= 6) & (eval_finite['timepoint2'] <= 6),

    # all samples
    slice(None)
]

COVERAGES = ["0.1x", "0.5x", "1x", "10x"]
COV_HUE = dict(zip(COVERAGES, color_palette))
JITTER = dict(zip(COVERAGES, [-0.08, -0.04, 0, 0.04]))

TOOL_LABELS = ["StrainGE", "StrainPhlan", "MIDAS"]

NUM_TOOLS = len(TOOLS)
NUM_TESTS = len(TEST_IX)

# 3 tests, 3 plots each
num_rows = NUM_TESTS * 3
num_cols = NUM_TOOLS + 1  # one extra for color bar

gridspec_kw = {
    # 'height_ratios': [15, 3, 3] * NUM_TESTS,
    'width_ratios': [15] * NUM_TOOLS + [1],
    'top': 0.99,
    'bottom': 0.01,
    'hspace': 0.2,
    'wspace': 0.15
}

fig = plt.figure(figsize=(5, 7))
outer_grid = fig.add_gridspec(3, 4, **gridspec_kw)

inner_grid = outer_grid[-1, -1].subgridspec(3, 1, height_ratios=[9, 1.5, 1.5], hspace=0.5)
cbar_axis = fig.add_subplot(inner_grid[-2:])

for col, tool in enumerate(TOOLS):
    for row, test_ix in enumerate(TEST_IX):
        test_df = eval_finite.loc[test_ix, :].copy()

        inner_grid = outer_grid[row, col].subgridspec(3, 1, height_ratios=[9, 1.5, 1.5], hspace=0.5)
        pr_curve_ax, heatmap_ax, barplot_ax = inner_grid.subplots()
        heatmap_ax.get_shared_x_axes().join(heatmap_ax, barplot_ax)

        aucs = []
        completed = []

        for cov in COVERAGES:
            ix = (tool, cov)
            if ix not in test_df.index:
                aucs.append(numpy.nan)
                completed.append(0.0)
                continue

            truth = test_df.loc[ix, 'truth']
            score = test_df.loc[ix, 'distance']

            precision, recall, thresholds = metrics.precision_recall_curve(truth, -score)
            auc = metrics.auc(recall, precision)

            precision[precision == 1.0] += JITTER[cov]
            pr_curve_ax.plot(recall, precision, label=f'{cov}', marker='.', color=COV_HUE[cov], linewidth=1.5)

            aucs.append(auc)
            completed.append(COMPLETED[ix])

        if col == 0:
            pr_curve_ax.set_ylabel("Precision")

        pr_curve_ax.set_aspect("equal", adjustable="box")
        pr_curve_ax.set_xlabel("Recall")
        pr_curve_ax.set_xlim(-0.05, 1.05)
        pr_curve_ax.set_ylim(-0.05, 1.08)

        # AUC heatmap
        if row == 0 and col == 0:
            cbar_args = {
                'cbar': True,
                'cbar_ax': cbar_axis,
                'cbar_kws': {'ticks': [0, 0.5, 1]}
            }
        else:
            cbar_args = {
                'cbar': False
            }

        aucs = numpy.array(aucs).reshape((1, 4))
        seaborn.heatmap(aucs, vmin=0, vmax=1, center=0.0, cmap="coolwarm",
                        ax=heatmap_ax, **cbar_args)
        heatmap_ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
        heatmap_ax.grid(False)
        for i in range(1, 5):
            heatmap_ax.axvline(i, color='white')

        if col == 0:
            heatmap_ax.set_ylabel("AUC", rotation=0, horizontalalignment='right', va='center')

        # % successful runs bar plot
        barwidth = 0.98
        x = numpy.arange(0, len(COVERAGES)) + ((1 - barwidth) / 2)
        barplot_ax.bar(x, completed, width=barwidth, tick_label=COVERAGES, align='edge',
                       color=color_palette)
        barplot_ax.set_xticklabels(COVERAGES, ha='left')
        barplot_ax.set_yticks([0, 1])
        barplot_ax.set_xlim(0, 4.0)
        barplot_ax.set_ylim(0, 1.0)
        barplot_ax.grid(False)
        barplot_ax.yaxis.set_major_formatter(ticker.PercentFormatter(1.0))

        if col == 0:
            barplot_ax.set_ylabel("", rotation=0, horizontalalignment='right',
                                  verticalalignment='center')

        if col != 0:
            pr_curve_ax.tick_params(labelleft=False)
            barplot_ax.tick_params(labelleft=False)

fig.show()

fig.savefig('strain_comparisons_all.svg', bbox_inches='tight')
