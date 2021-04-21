"""Module with functionalities shared among energy systems."""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def create_multipage_pdf(file_name='plots.pdf', figs=None, dpi=300,
                         mute=False):
    """Save all open matplotlib figures into a multipage pdf-file.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>>
    >>> df1 = pd.DataFrame(np.random.randn(24, 2))
    >>> ax1 = df1.plot(kind='line')
    >>>
    >>> df2 = pd.DataFrame(np.random.randn(24, 2))
    >>> ax2 = df2.plot(kind='scatter', x=0, y=1)
    >>>
    >>> # mute is set to true to surpress writing a pdf file
    >>> create_multipage_pdf(file_name='plots.pdf', dpi=300, mute=True)
    False
    """
    if mute is True:
        # set return flag to false if no output is written
        flag = False
    else:
        pp = PdfPages(file_name)
        if figs is None:
            figs = [plt.figure(n) for n in plt.get_fignums()]
        for fig in figs:
            fig.savefig(pp, format='pdf')
        pp.close()

        # close all existing figures
        for fig in figs:
            plt.close(fig)

        # set return flag
        flag = True

    return flag


def znes_colors(n=None):
    """Return dict with ZNES colors.

    Examples
    --------
    >>> znes_colors().keys()  # doctest: +ELLIPSIS
    dict_keys(['darkblue', 'red', 'lightblue', 'orange', 'grey',...
    """
    colors = {
        'darkblue': '#00395B',
        'red': '#B54036',
        'lightblue': '#74ADC0',
        'orange': '#EC6707',
        'grey': '#BFBFBF',
        'dimgrey': 'dimgrey',
        'lightgrey': 'lightgrey',
        'slategrey': 'slategrey',
        'darkgrey': '#A9A9A9'
    }

    # allow for a dict of n colors
    if n is not None:
        return {k: colors[k] for k in list(colors)[:n]}
    else:
        return colors


def znes_boxprops():
    """Return dict with ZNES boxplot properties.

    Examples
    --------
    >>> znes_boxprops().keys()  # doctest: +ELLIPSIS
    dict_keys(['boxprops', 'flierprops', 'medianprops', 'whiskerprops',...
    """
    znes = znes_colors()

    props = {}
    # boxplot properties as dict since rcparams doesn't show effect with pd
    props['boxprops'] = dict(linewidth=3.5)
    props['flierprops'] = dict(
        linewidth=2, marker='D', markerfacecolor=znes['darkblue'])
    props['medianprops'] = dict(
        linewidth=3.5, markerfacecolor=znes['darkblue'])
    props['whiskerprops'] = dict(
        linewidth=3.5, markerfacecolor=znes['darkblue'])
    props['capprops'] = dict(
        linewidth=3.5, markerfacecolor=znes['darkblue'])
    # passing patch_artist=True to a pd.boxplot fills the box
    props['colorprops'] = dict(
        boxes=znes['darkblue'], whiskers=znes['darkblue'],
        medians=znes['darkblue'], caps=znes['darkblue'])

    return props


def znes_style(plotting_function):
    """Decorator to create basic matplotlib configuration with ZNES style.

    Set markers, lines and colors. Create a znes color palette.
    """
    mpl.rcParams.update(mpl.rcParamsDefault)  # reset to defaults
    mpl.style.use('default')
    znes = znes_colors()
    znes_palette = [znes['darkblue'], znes['orange'], znes['grey'],
                    znes['red'], znes['lightblue'], znes['dimgrey'],
                    znes['lightgrey'], znes['slategrey']]
    znes_palette_r = znes_palette[::-1]
    cmap_znes = mpl.colors.ListedColormap(znes_palette)
    cmap_znes_r = mpl.colors.ListedColormap(znes_palette_r)
    mpl.cm.register_cmap(name='znes', cmap=cmap_znes)
    mpl.cm.register_cmap(name='znes_r', cmap=cmap_znes_r)
    plt.rcParams['image.cmap'] = 'znes_r'
    # grid
    plt.rcParams['grid.color'] = 'k'  # grid color default b0b0b0
    plt.rcParams['grid.linestyle'] = 'dotted'  # solid
    plt.rcParams['grid.linewidth'] = 2.0  # in points default 0.8
    plt.rcParams['grid.alpha'] = 1.0  # transparency, between 0.0 and 1.
    # lines and markers
    plt.rcParams['lines.linewidth'] = 3.5
    plt.rcParams['scatter.marker'] = 'o'
    # axes, ticks
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['xtick.color'] = 'k'
    plt.rcParams['ytick.color'] = 'k'
    plt.rcParams['text.color'] = 'k'
    plt.rcParams['axes.labelcolor'] = 'k'
    plt.rcParams['axes.grid'] = True
    plt.rcParams['axes.axisbelow'] = True
    # legend
    plt.rcParams['legend.fontsize'] = 25
    plt.rcParams['legend.loc'] = 'upper right'
    plt.rcParams['legend.frameon'] = True
    # use normal text for math (non-italic)
    plt.rcParams.update({'mathtext.default':  'regular'})
    # figure
    plt.rcParams['figure.figsize'] = (20.0, 15.0)  # inches = px/dpi
    # font helvetica clone
    mpl.rcParams['font.size'] = 35
    mpl.rcParams['font.family'] = 'Carlito'  # Liberation Sans
    # mpl.rcParams['font.weight'] = 'light'  # only for ticks and legend
    mpl.rcParams['savefig.format'] = 'pdf'
    mpl.rcParams['savefig.bbox'] = 'tight'
    mpl.rcParams['savefig.pad_inches'] = 0.1  # padding when bbox is 'tight'

    # deactivate warning for high amount of figures
    mpl.rcParams.update({'figure.max_open_warning': 0})

    return plotting_function


def znes_style_plot():
    """Test function to test abovementioned decorator function.

    Checks if matplotlib figsize param is set locally.
    """
    param = plt.rcParams['figure.figsize']
    return param


def znes_style_test():
    """Test function to test abovementioned decorator function.

    Examples
    --------
    >>> plot = znes_style(znes_style_plot)
    >>> plot()
    [20.0, 15.0]
    """
    pass


def znes_linestyles(columns=None):
    """Return a list with line styles for a passed column list.

    Examples
    --------
    >>> znes_linestyles(['foo', 'bar'])
    ['-', '-']
    """
    linestyles = ['-', '-', '-', '-', '-', '-.', ':', '-.',
                  '-', '-', '-', '-', '-', '-.', ':', '-.']

    return linestyles[:len(columns)]


def znes_markers(columns=None):
    """Return a list with marker styles for a passed column list.

    Examples
    --------
    >>> znes_markers(['foo', 'bar'])
    ['o', 's']
    """
    markers = ['o', 's', 'v', 'x', 'H', '^', 'v', 's', '3', '.', '1', '_',
               'o', 's', 'v', 'x', 'H', '^', 'v', 's', '3', '.', '1', '_']

    return markers[:len(columns)]


def znes_linear_colormap(name='znes_linear', bins=256,
                         colors=['#00395B', '#FFFFFF', '#EC6707']):
    """Return a linear segmented colormap from three passed colors.

    Examples
    --------
    >>> cm = znes_linear_colormap(bins=127)
    >>> cm.N
    127
    >>> print(cm.name)
    znes_linear
    """
    cm = mpl.colors.LinearSegmentedColormap.from_list(name, colors, N=bins)

    return cm


def znes_linear_colormap2(name='znes_linear2', bins=20,
                          colors=['#00395B', '#EC6707', '#BFBFBF',
                                  '#B54036', '#74ADC0']):
    """Return a linear segmented colormap from three passed colors.

    Examples
    --------
    >>> cm = znes_linear_colormap2(bins=77)
    >>> cm.N
    77
    >>> print(cm.name)
    znes_linear2
    """
    cm = mpl.colors.LinearSegmentedColormap.from_list(name, colors, N=bins)

    return cm


def znes_sample_dataframe(length=25):
    """Return a sample DataFrame of defined length to test plots."""
    df = pd.DataFrame()
    df['Foo'] = np.arange(1, length)
    df['Bar'] = max(df['Foo']) - np.sqrt(df['Foo'])
    df['Foobar'] = df['Foo'] ** 2
    df = df.append(df * 2)
    df[['A', 'B', 'C']] = pd.DataFrame(np.random.rand(df.shape[0], 3))
    df['ndf'] = 0
    df['ndf'] = df['ndf'].where(df['Bar'] < (max(df['Bar'])/2)+1, 1)
    df = df.reset_index(drop=True)
    df.index = pd.DatetimeIndex(start='2018-01-01 00:00:00',
                                periods=df.shape[0], freq='h')

    return df


if __name__ == '__main__':
    import doctest
    doctest.testmod()
