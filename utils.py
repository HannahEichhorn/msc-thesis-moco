import numpy as np
import matplotlib.pyplot as plt
import os



def MakeBoxplot(metric, colors):
    '''
    Creating a box plot with preset colors

    Parameters
    ----------
    metric : array
        values to plot.
    colors : array or list
        colors for the different types.

    Returns
    -------
    None.

    '''
    
    small = dict(markersize=3)
    box1 = plt.boxplot(metric, flierprops=small, widths=0.5)
    for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors):
            patch.set(color=color)
            patch2.set(color=color)


def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.03, fs=None, maxasterix=None):
    """
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    lx, ly = center[num1], height[num1]
    rx, ry = center[num2], height[num2]

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh-0.035*barh)

    plt.plot(barx, bary, c='dimgrey')

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs, c='dimgrey')


def Show_Stars(p_cor, ind, bars, heights, arange_dh=False):
    '''
    Function to show brackets with asterisks indicating statistical 
    significance

    Parameters
    ----------
    p_cor : array or list
        corrected p-values.
    ind : list
        indices for the comparisons corresponding to the individual p-values.
    bars : list or array
        x position of boxplots.
    heights : list or array
        maximal value visualised in boxplots.
    arange_dh : float, optional
        offset above heights. The default is False.

    Returns
    -------
    int
        returns 0, when completed..

    '''
    
    star = p_cor < 0.05
    starstar = p_cor < 0.001
    
    all_nrs = []
    dh = .1
    if arange_dh==True:
        dhs = np.array([.1,.1,.2, .1, .1, .2, .1, .1, .2,.2, .2, .2, .2, .1, .1, .1, .1, .2, .3, .1, .2, .3, .2, .2])
    if arange_dh == 'REAC':
        dhs = np.array([.1, .1, .1, .2, .4, .4, .1, .1, .2, .4,.4])
    if arange_dh == 'OLD':
        dhs = np.array([.1, .1, .1, .2, .1, .1, .2 ])
    if arange_dh == '2D':
        dhs = np.array([.1, .1, .1, .2, .3, .3, .2 ])
    if arange_dh == 'diff':
        dhs = np.array([.1, .1, .3])
    for i in range(0,len(p_cor)):
        nr = ind[i]
        all_nrs.append(nr)
        if arange_dh:
            dh = dhs[i]
        if starstar[i]==True:
            barplot_annotate_brackets(nr[0], nr[1], '**', bars, heights, dh=dh)
        else:
            if star[i]==True:
                barplot_annotate_brackets(nr[0], nr[1], '*', bars, heights, dh=dh)

    return 0


def SortFiles(file):
    '''
    sort a number of files after date and pick the most recent file.

    Parameters
    ----------
    file : list
        list of filenames.

    Returns
    -------
    file : string
        filename corresponding to the most recent file.

    '''
    
    ind = []
    for i in range(len(file)):
        base = os.path.basename(file[i])
        if len(base)<30:
            ind.append(i)
    file = np.array(file)
    file = file[ind]
    file = np.sort(file)[::-1]
    return file


def Scale_Axis_Lim(data):
    '''
    Adjust the limits of the y axis.

    Parameters
    ----------
    data : array
        values to be plotted.

    Returns
    -------
    y_min : float
        minimal value for y-axis.
    y_max : float
        maixmal value for y-axis.

    '''
    
    maxim = np.amax(data)
    minim = np.amin(data)
    if np.isinf(maxim):
        temp = np.sort(data)
        maxim = temp[-2]

    y_min = minim - 0.1*(maxim-minim)
    y_max = maxim + 0.05*(maxim-minim)

    if np.isnan(y_min):
        y_min = 0

    if np.isnan(y_max):
        y_max = 1

    return y_min, y_max


