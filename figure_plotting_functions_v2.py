import math

import seaborn as sns
import matplotlib.pyplot as plt
import cmcrameri.cm as cmc
import pandas as pd
import scipy
import numpy as np



def plot_variable_regions(ax, var_regs):
    for V in var_regs:
        ax.axvspan(V[0], V[1], facecolor='g', alpha=0.2)
    return


def set_params_for_plots(figsize, context):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': figsize}
    # cmap = np.append([cmc.batlow.colors[0]], [cmc.batlow.colors[-1]], axis=0)
    sns.set_theme(context=context, style="ticks", rc=custom_params, palette='viridis')
    return

# For graph 1a
def generate_summary_graphs(average_conserved: list, conserved_igs_true_perc_coverage: dict,
                            delta_composite_scores: dict, file_out: bool = False, file: str = ''):
    sns.set_theme(style='white', rc={"axes.spines.right": False, "axes.spines.top": False})
    fig, axs = plt.subplots(3, 1, layout='constrained')

    # sns.histplot(x=['No', 'Yes'], weights=is_u_conserved[1], ax=axs[0], binrange=[0, 1])
    sns.histplot(average_conserved, kde=True, ax=axs[0])
    axs[0].set_xlabel('Percent of background sequences with conserved Us to design')
    axs[0].set_title('Average conserved U sites between targeted designs and background sequences')

    sns.histplot(conserved_igs_true_perc_coverage, kde=True, ax=axs[1])
    axs[1].set_xlabel('True percent coverage in background sequences')
    axs[1].set_title('Distribution of conserved IGS sites in background sequences')
    axs[1].set_xlim(left=0, right=1)

    sns.histplot(delta_composite_scores, kde=True, ax=axs[2], binrange=[0, 1])
    axs[2].set_xlabel(u'Δ composite scores (design - background)')
    axs[2].set_title(u'Δ Composite scores of generated guides with conserved IGS on background sequences')
    axs[2].set_xlim(left=0, right=1)

    if file_out:
        plt.savefig(f'{file}/targeted designs vs background general stats.svg', transparent=False)
        plt.savefig(f'{file}/targeted designs vs background general stats.png')
    plt.show()

    # With more targets, make sure to graph location vs. scores (igs and guide) and the std dev per position

    return

# For graph 1b
def generate_igs_vs_guide_graph(conserved_igs_true_perc_coverage: dict, names_and_guide_scores: dict,
                                file_out: bool = False, file: str = '', title: str = 'targeted designs vs background'):
    igs_true_cov_vs_guide_comp_score_dict = {key: (int(key[5:]), conserved_igs_true_perc_coverage[key],
                                                   names_and_guide_scores[key]) for key in names_and_guide_scores}
    igs_vs_guide_df = pd.DataFrame(igs_true_cov_vs_guide_comp_score_dict,
                                   index=['Index', 'IGS true percent coverage',
                                          'Guide score at U site']).T

    slope, intercept, r, p, sterr = scipy.stats.linregress(x=igs_vs_guide_df['IGS true percent coverage'],
                                                           y=igs_vs_guide_df['Guide score at U site'])
    reg_plot = sns.jointplot(igs_vs_guide_df, x='IGS true percent coverage', y='Guide score at U site',
                             kind='reg', line_kws={'color': 'plum'}, xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))
    reg_plot.ax_joint.annotate(f'$r^2$={round(r, 3)}', xy=(0.1, 0.9), xycoords='axes fraction')

    if file_out:
        plt.savefig(f'{file}/{title} scoring correlations.svg', transparent=False)
        plt.savefig(f'{file}/{title} scoring correlations.png')
    plt.show()
    return

# For graph 2a
def generate_all_scores_vs_loc_graph(designs: np.ndarray, var_regs: list, group: str):
    # y-axis is the composite score, x-axis is the 16s rRNA gene, plot the universal, control, and each
    # selective for all designs in different panels (same data as above but order along gene)

    # Make figure (plt, ax):
    fig, axs = plt.subplots(3, 1, sharex=True, layout='constrained')

    # Set variable regions
    for ax in axs:
        plot_variable_regions(ax, var_regs)

    # Extract all relevant data
    u_conservation = [(ribodesign.ref_idx, ribodesign.u_conservation_background) for ribodesign in designs]
    igs_scoring = [(ribodesign.ref_idx, ribodesign.true_perc_cov) for ribodesign in designs]
    guide_scoring = [(ribodesign.ref_idx, ribodesign.score) for ribodesign in designs]

    # plot!
    sns.scatterplot(data=u_conservation, ax=axs[0], y='U conservation in background', alpha=0.8, legend=False)
    sns.scatterplot(data=igs_scoring, ax=axs[1], y='IGS true % coverage', alpha=0.8, legend=False)
    sns.scatterplot(data=guide_scoring, ax=axs[2], y=f'Guide score ({designs[0].score_type})', alpha=0.8, legend=False)
    sns.despine(offset=10, trim=False)
    ax[0].set_title(f'Scores for {group} designs')

    return fig, axs


def generate_axes(fig, graph: str, titles: list[str]):
    # Credit to https://tombohub.github.io/matplotlib-layout-generator/ for this amazing tool!!!
    axes = {}
    rows = 2
    cols = 3
    if graph == '1a' or graph == '2a':
        height = 2
        length = 5
        # 3 small graphs per condition, each graph takes 5 units of length + 1 of separation between graphs
        gridspec = fig.add_gridspec(nrows=rows * 3 + rows - 1, ncols=cols * 5 + cols - 1)
        for i, title in enumerate(titles):
            row_coord, col_coord = generate_coords(i, cols=cols, rows=rows, len=length)
            axes[f'U conservation in background {title}'] = fig.add_subplot(
                gridspec[row_coord:row_coord+height, col_coord:col_coord+length])
            axes[f'True percent coverage {title}'] = fig.add_subplot(
                gridspec[row_coord+2:row_coord+height*2, col_coord:col_coord+length])
            axes[f'Guide score {title}'] = fig.add_subplot(
                gridspec[row_coord+4:row_coord+height*3, col_coord:col_coord+5])
    elif graph == '1b':
        # plain old one ax per figure
        height = 4
        length = 4
        gridspec = fig.add_gridspec(nrows=rows * 3 + rows - 1, ncols=cols * 5 + cols - 1)
        for i, title in enumerate(titles):
            row_coord, col_coord = generate_coords(i, cols=cols, rows=rows, len=length, hei=height)
            axes[f'{title}'] = fig.add_subplot(gridspec[row_coord:row_coord+height, col_coord:col_coord+length])
    elif graph == '2b' or graph == '4':
        height = 2
        length = 5
        gridspec = fig.add_gridspec(nrows=len(titles)*height, ncols=length)
        for i, title in enumerate(titles):
            row_coord, col_coord = generate_coords(i, cols=1, rows=len(titles), len=length, hei=height)
            axes[f'{title}'] = fig.add_subplot(gridspec[row_coord:row_coord + height, col_coord:col_coord + length])
    elif graph == '5':
        height = 4
        length = 2 * len(titles)
        # 3 small graphs per condition, each graph takes 5 units of length + 1 of separation between graphs
        gridspec = fig.add_gridspec(nrows=len(titles)*height, ncols=length)
        row_coord_top, col_coord_top = generate_coords(0, cols=1, rows=len(titles), len=length, hei=height)
        row_coord_bot, col_coord_bot = generate_coords(1, cols=1, rows=len(titles), len=length, hei=height)
        axes[f'Fraction of order in each condition'] = fig.add_subplot(
                gridspec[row_coord_top:row_coord_top + height, col_coord_top:col_coord_top + length])
        axes[f'Count of order in each condition'] = fig.add_subplot(
            gridspec[row_coord_bot:row_coord_bot + height, col_coord_bot:col_coord_bot + length])

    return axes


def generate_coords(i, cols, rows, len, hei):
    row_coord = (hei + 1) * math.floor(i / cols)
    if math.floor(i / rows) == math.ceil(i / rows):
        col_coord = 0
    elif math.floor(i / rows) == round(i / rows):
        col_coord = 1 + len
    else:
        col_coord = 2 + len * 2

    return row_coord, col_coord



