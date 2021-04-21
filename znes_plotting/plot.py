"""Module with plot functions and tests based on sample data."""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from znes_plotting import shared


@shared.znes_style
def line(data=None, xlabel='', ylabel='', markersize=15, marker='None',
         **kwargs):
    """"Create line plot from given data structure."""
    colors = list(shared.znes_colors().values())[0:data.shape[0]]
    if hasattr(data, 'columns'):
        markers = shared.znes_markers(data.columns)
    else:
        markers = shared.znes_markers([0])

    # treat the case that colors are passed explicitly
    if 'color' in kwargs:
        colors = kwargs['color']
        kwargs.pop('color', None)

    ax = data.plot(kind='line', color=colors, grid=True,
                   markersize=markersize, **kwargs)

    if 'subplots' in kwargs.keys():
        if kwargs['subplots'] is True:
            one_axis = False
        else:
            one_axis = True
    else:
        one_axis = True

    if one_axis is True:
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        # draw vertical grid lines
        ax.grid('on', which='minor', axis='x')
        ax.grid('on', which='major', axis='x')

        # add markers
        if marker is not 'None':
            markers_ax1 = markers[0:len(ax.get_lines())]
            for i, line in enumerate(ax.get_lines()):
                line.set_marker(markers_ax1[i])

        # legend
        ax.legend(loc='upper right')
    else:
        # format labels
        [ax.set_xlabel(xlabel) for ax in plt.gcf().axes]
        [ax.set_ylabel(ylabel) for ax in plt.gcf().axes]

        # add markers
        if marker is not 'None':
            for ax in plt.gcf().axes:
                markers_ax1 = markers[0:len(ax.get_lines())]
                for i, line in enumerate(ax.get_lines()):
                    line.set_marker(markers_ax1[i])

        # add legends
        for ax in plt.gcf().axes:
            h1, l1 = ax.get_legend_handles_labels()
            ax.legend(h1, l1, loc='upper right')

    return ax


@shared.znes_style
def line_secondary_y(data=None, xlabel='', ylabel='', secondary_ylabel='',
                     split_legend=True, markersize=15, marker='o',
                     secondary_ylim=None, legend_fontsize=None,
                     color=list(shared.znes_colors().values()),
                     **kwargs):
    """Create line plot with secondary y axis from given data structure."""
    ax = data.plot(markersize=markersize, marker=marker, color=color, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.right_ax.set_ylabel(secondary_ylabel)

    # secondary ylim
    if secondary_ylim is not None:
        ax.right_ax.set_ylim(secondary_ylim[0], secondary_ylim[1])

    # draw vertical grid lines
    ax.grid(True)
    ax.set_axisbelow(True)
    ax.right_ax.grid(False)

    # add markers
    if marker is not 'None':
        markers = shared.znes_markers(data.columns)
        markers_ax1 = markers[0:len(ax.get_lines())]
        markers_ax2 = markers[len(markers_ax1)::]
        for i, line in enumerate(ax.get_lines()):
            line.set_marker(markers_ax1[i])
        for i, line in enumerate(ax.right_ax.get_lines()):
            line.set_marker(markers_ax2[i])

    # add legends
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax.right_ax.get_legend_handles_labels()
    if split_legend is True:
        if legend_fontsize is None:
            ax.legend(h1, l1, loc='upper left')
            ax.right_ax.legend(h2, l2, loc='upper right')
        else:
            ax.legend(h1, l1, loc='upper left', fontsize=legend_fontsize)
            ax.right_ax.legend(h2, l2, loc='upper right',
                               fontsize=legend_fontsize)
    else:
        if legend_fontsize is None:
            ax.legend(h1+h2, l1+l2, loc='upper right')
        else:
            ax.legend(h1+h2, l1+l2, loc='upper right',
                      fontsize=legend_fontsize)

    return ax


@shared.znes_style
def bar(data=None, xlabel='', ylabel='', legend_fontsize=None, **kwargs):
    """Create bar plot from given data structure."""
    if 'color' not in kwargs:
        kwargs['color'] = list(shared.znes_colors().values())[0:data.shape[1]]
    if 'rot' not in kwargs:
        kwargs['rot'] = 0
    ax = data.plot(kind='bar', grid=True, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if legend_fontsize is not None:
        ax.legend(fontsize=legend_fontsize)

    return ax


@shared.znes_style
def barh(data=None, xlabel='', ylabel='', legend_fontsize=None, **kwargs):
    """Create horizontal bar plot from given data structure."""
    if 'color' not in kwargs:
        kwargs['color'] = list(shared.znes_colors().values())[0:data.shape[1]]
    if 'rot' not in kwargs:
        kwargs['rot'] = 0
    ax = data.plot(kind='barh', grid=True, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if legend_fontsize is not None:
        ax.legend(fontsize=legend_fontsize)

    return ax


@shared.znes_style
def bar_secondary_y(data=None, xlabel='', ylabel='', secondary_ylabel='',
                    split_legend=True, secondary_ylim=None,
                    legend_fontsize=None,
                    color=list(shared.znes_colors().values()),
                    **kwargs):
    """Create bar plot with secondary y axis from given data structure."""
    ax = data.plot(kind='bar', color=color, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.right_ax.set_ylabel(secondary_ylabel)

    # secondary ylim
    if secondary_ylim is not None:
        ax.right_ax.set_ylim(secondary_ylim[0], secondary_ylim[1])

    # draw vertical grid lines
    ax.grid(True)
    ax.set_axisbelow(True)
    ax.right_ax.grid(False)

    # add legends
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax.right_ax.get_legend_handles_labels()
    if split_legend is True:
        if legend_fontsize is None:
            ax.legend(h1, l1, loc='upper left')
            ax.right_ax.legend(h2, l2, loc='upper right')
        else:
            ax.legend(h1, l1, loc='upper left', fontsize=legend_fontsize)
            ax.right_ax.legend(h2, l2, loc='upper right',
                               fontsize=legend_fontsize)
    else:
        if legend_fontsize is None:
            ax.legend(h1+h2, l1+l2, loc='upper right')
        else:
            ax.legend(h1+h2, l1+l2, loc='upper right',
                      fontsize=legend_fontsize)

    return ax


@shared.znes_style
def scatter(data=None, xlabel='', ylabel='', grid=True,
            markersize=100, marker='o', edgecolors='k', linewidth=1,
            **kwargs):
    """Create scatter plot from given data structure."""
    ax = data.plot(kind='scatter', grid=grid, rot=0, s=markersize,
                   marker=marker,
                   # edgecolors=edgecolors, linewidth=linewidth,
                   **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax


@shared.znes_style
def scatter_pareto_front(data=None, xlabel='', ylabel='', label='',
                         markersize=100, edgecolors='k', linewidth=1,
                         marker='o', color_df='darkgrey', color_ndf='darkblue',
                         alpha_df=1,
                         label_ndf='Non-dominated front',
                         label_df=None,
                         marker_df='o',
                         marker_ndf='o',
                         **kwargs):
    """Create grouped scatter plot from given data structure."""
    colors = shared.znes_colors()
    ax = data.plot(kind='scatter', color=colors[color_df], grid=True,
                   rot=0, s=markersize, label=label_df, alpha=alpha_df,
                   marker=marker_df,
                   # edgecolors=edgecolors, linewidth=linewidth,
                   **kwargs)
    data.loc[data['ndf'] == 1].plot.scatter(
        x=kwargs['x'], y=kwargs['y'], label=label_ndf,
        color=colors[color_ndf], s=markersize, marker=marker_ndf,
        # edgecolors=edgecolors, linewidth=linewidth,
        ax=ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax


@shared.znes_style
def box(data=None, xlabel='', ylabel='', **kwargs):
    """Create box plot from given data structure."""
    colors = shared.znes_colors()
    color = dict(boxes=colors['darkblue'], whiskers=colors['darkblue'],
                 medians=colors['darkblue'], caps=colors['darkblue'])
    props = shared.znes_boxprops()
    ax = data.plot(kind='box', color=color,
                   boxprops=props['boxprops'],
                   flierprops=props['flierprops'],
                   medianprops=props['medianprops'],
                   whiskerprops=props['whiskerprops'],
                   capprops=props['capprops'],
                   showmeans=False, showfliers=True,
                   grid=True, rot=0, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax


@shared.znes_style
def histogram(data=None, xlabel='', ylabel='', bins=30, legend=False,
              add_information=False, **kwargs):
    """Create histogram from given data structure."""
    colors = shared.znes_colors()
    ax = data.plot(kind='hist', color=colors['darkblue'], edgecolor='white',
                   linewidth=1.5, bins=bins, legend=legend,
                   grid=True, rot=0, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if add_information is True:
        if isinstance(bins, list):
            ax.set_title(
                'n='+str(data.shape[0])+', bin width='+str(bins[1]-bins[0]))
        else:
            ax.set_title('n='+str(data.shape[0])+', bins='+str(bins))

    return ax


@shared.znes_style
def surface_fit(x=None, y=None, z=None, scatter_data=None,
                xlabel='', ylabel='', zlabel='', markersize=50, alpha=0.5,
                **kwargs):
    """Create 3D plot which compares a fit to original data points."""
    colors = shared.znes_colors()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.plot_surface(x, y, z, color=colors['darkblue'], alpha=alpha)
    ax.scatter(scatter_data[:, 0], scatter_data[:, 1], scatter_data[:, 2],
               c=colors['red'], s=markersize)
    ax.set_xlabel('\n\n' + xlabel)
    ax.set_ylabel('\n\n' + ylabel)
    ax.set_zlabel('\n\n' + zlabel)

    return ax


if __name__ == '__main__':

    # Create sample data
    df = shared.znes_sample_dataframe(length=50)

    # Test bar plot
    bar(data=df[['A', 'B', 'C']][0:4].reset_index(drop=True), stacked=False,
        xlabel='Foo', ylabel='Bar')

    # Test line plot
    ax = line(data=df[['A', 'B', 'C']][0:4], xlabel='Foo', ylabel='Bar')

    # Test line plot with secondary y-axis
    ax = line_secondary_y(data=df[['A', 'B', 'C']][0:4],
                          secondary_y=['B', 'C'],
                          xlabel='Foo', ylabel='Bar',
                          secondary_ylabel='Foobar')

    # Test histogram
    histogram(data=df[['A']], xlabel='Foo', ylabel='Bar')

    # Test basic scatter plot
    scatter(data=df, xlabel='Foo',
            ylabel='Bar', x='Foo', y='Bar',
            color=shared.znes_colors()['darkblue'])

    # Test scatter plot with pareto front
    scatter_pareto_front(
        data=df, xlabel='Foo', ylabel='Bar',
        x='Foo', y='Bar')

    # Test grouping of scatter plots
    ax = scatter(data=df[df['ndf'] == 1], x='Foo', y='Bar', label='Scenario1',
                 color=shared.znes_colors()['darkblue'])
    scatter(data=df[df['ndf'] == 1]*1.5, xlabel='Foo',
            ylabel='Bar', x='Foo', y='Bar',
            label='Scenario2', color=shared.znes_colors()['orange'],
            ax=ax)

    # Test grouped scatter plot with color bar
    ax = scatter(data=df[df['ndf'] == 1], x='Foo', y='Bar',
                 colormap=shared.znes_linear_colormap(bins=25),
                 label='Non-dominated front',
                 # norm=mpl.colors.LogNorm(),
                 c='Foobar', zorder=2)
    ax2 = scatter(data=df[df['ndf'] == 0], x='Foo', y='Bar',
                  xlabel='Foo', ylabel='Bar',
                  label='Dominated fronts',
                  color='darkgrey', ax=ax, zorder=1)

    # Test ungrouped scatter plot with color bar
    ax = scatter(data=df, x='Foo', y='Bar',
                 xlabel='Foo', ylabel='Bar',
                 colormap=shared.znes_linear_colormap(bins=50),
                 norm=mpl.colors.LogNorm(),
                 c='Foobar')

    # Test line plot for two objectives depending on parameters
    df2 = df.loc[df['ndf'] == 1][['Foo', 'Bar', 'Foobar']].set_index('Foobar')
    ax = line_secondary_y(
        data=df2, xlabel='Foobar', ylabel='Foo',
        secondary_y=['Bar'], secondary_ylabel='Bar')

    # Test box plot with bins
    df2 = df[['A']].copy()
    bin_size = 0.25
    bins = [i for i in np.arange(0, 1, bin_size)] + [np.inf]
    labels = ['<= '+str(bin+bin_size) for bin in bins if bin is not np.inf]
    df2['bins'] = pd.cut(df['A'].values, bins, labels=labels)
    df2 = df2.pivot(columns='bins', values='A')
    box(data=df2, xlabel='Foo', ylabel='Bar')

    # Create multipage PDF from plots (set mute=False to write output)
    shared.create_multipage_pdf('plots.pdf', mute=True)
