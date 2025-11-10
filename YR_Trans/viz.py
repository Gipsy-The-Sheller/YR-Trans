# Logistics of visualization datas
from plotnine import *
import numpy as np
import pandas as pd
import scipy.stats as st
import seaborn as sns
import logging

def VolcanoPlot(data, 
                GeneID=0, log2FC=1, pvalue=2, GeneName=None, 
                theme = 'linedraw', 
                color_scheme=['discrete', {'up': '#FF0000', 'down': '#0000FF', 'no-DEGs': '#ADAAAB'}], 
                color_distribution_convert=True,
                x_threshold=1, y_threshold=0.05, threshold_indicator=True, under_threshold_color=True, trimmode="none",
                alt = True):
    # This function is originally written in the SuperVolcano plugin in YRTools v0.0.1

    # trimmode: "99hpd" / "95hpd" / "none" NOTE: abandoned
    # color_scheme: ['discrete', {'up': up_color, 'down': down_color, 'no-DEGs': no_DEGs_color}]
    # color_scheme: ['continuous', cmap]
    # in this case the color will be set based on (log2FC-x_threshold)*(-log10(pvalue)-y_threshold), namely a inverse proportional relationship
    # extract log2FC and pvalue
    log2FC_index = data.columns[log2FC]
    pvalue_index = data.columns[pvalue]
    log2FC = data[log2FC_index]
    pvalue = data[pvalue_index]
    neg_log10_pvalue = -np.log10(pvalue)
    plot_data = pd.DataFrame({"log2FC": log2FC, "-log10(Pvalue)": neg_log10_pvalue})

    x_min = plot_data["log2FC"].min()
    x_max = plot_data["log2FC"].max()
    y_min = plot_data["-log10(Pvalue)"].min()
    y_max = plot_data["-log10(Pvalue)"].max()
    x_lim = max(abs(x_min), abs(x_max))
    y_limup = y_max

    # if trimmode == '99hpd':
    #     x_interval = st.norm.interval(0.99, loc=np.mean(plot_data["log2FC"]), scale=st.sem(plot_data["log2FC"]))
    #     x_lim = max(abs(x_interval[0]), abs(x_interval[1]))

    if color_scheme[0] == 'discrete':
        cmap = []
        stat = {"up":0, "down":0, "no-DEGs":0}
        for i in range(len(plot_data)):
            if data[pvalue_index][i] > y_threshold:
                cmap.append("no-DEGs")
                stat['no-DEGs'] += 1
            elif abs(data[log2FC_index][i]) < x_threshold:
                cmap.append("no-DEGs")
                stat['no-DEGs'] += 1
            elif data[log2FC_index][i] > 0:
                cmap.append("up")
                stat['up'] += 1
            else:
                cmap.append("down")
                stat['down'] += 1

        ups_txt = "up"
        downs_txt = "down"
        noDEGs_txt = "no-DEGs"

        if alt == True:
            ups_txt = f"up {str(stat['up'])}"
            downs_txt = f"down {str(stat['down'])}"
            noDEGs_txt = f"no-DEGs {str(stat['no-DEGs'])}"
        
        cmap = [ups_txt if i == 'up' else i for i in cmap]
        cmap = [downs_txt if i == 'down' else i for i in cmap]
        cmap = [noDEGs_txt if i == 'no-DEGs' else i for i in cmap]
        color_scheme[1][ups_txt] = color_scheme[1]['up']
        color_scheme[1][downs_txt] = color_scheme[1]['down']
        color_scheme[1][noDEGs_txt] = color_scheme[1]['no-DEGs']
        
        plot_data.insert(0, "cmap", cmap)
    if color_scheme[0] == 'cmap' or color_scheme[0] == 'gradient':
        cmap = [ ]
        stat = {"up":0, "down":0, "no-DEGs":0}
        for i in range(len(plot_data)):
            if under_threshold_color:
                x_div = plot_data["log2FC"][i]
                y_div = plot_data["-log10(Pvalue)"][i]
            else:
                if plot_data["log2FC"][i] > x_threshold or plot_data["log2FC"][i] < -x_threshold:
                    x_div = plot_data["log2FC"][i]
                else:
                    x_div = 0
                if plot_data["-log10(Pvalue)"][i] > -np.log10(y_threshold):
                    y_div = plot_data["-log10(Pvalue)"][i]
                else:
                    y_div = 0
            if y_div > -np.log10(y_threshold):
                if x_div < -x_threshold:
                    stat['down'] += 1
                if x_div > x_threshold:
                    stat['up'] += 1
                else:
                    stat['no-DEGs'] +=1
            else:
                stat['no-DEGs'] += 1
            cmap.append(x_div * y_div)
        plot_data.insert(0, "cmap", cmap)
        # find the max value of cmap
        color_scheme.append([-max(cmap), max(cmap)])

        ups_txt = "up"
        downs_txt = "down"
        noDEGs_txt = "no-DEGs"

        if alt == True:
            ups_txt = f"up {str(stat['up'])}"
            downs_txt = f"down {str(stat['down'])}"
            noDEGs_txt = f"no-DEGs {str(stat['no-DEGs'])}"

        if color_distribution_convert:
            # use normal distribution to determine the variance rate of cmap
            # 1. estimate distribution parameters
            cmap = np.array(cmap, dtype=float)
            mean = np.mean(cmap)
            std = np.std(cmap)
            # 2. use an approximate integral to estimate the area
            cmap = plot_data["cmap"]
            for i in range(len(cmap)):
                cmap[i] = st.norm.cdf(cmap[i], mean, std)
            plot_data["cmap"] = cmap
            maxx = max(cmap)
            minx = min(cmap)
            standard = st.norm.cdf(0, mean, std)
            color_scheme[2] =  [standard - max(maxx - standard, standard - minx), standard + max(maxx - standard, standard - minx)]
        
    def getcolorscheme(color_scheme):
        if color_scheme[0] == 'discrete':
            return scale_color_manual(values=color_scheme[1])
        elif color_scheme[0] == 'cmap':
            return scale_color_cmap(name='  ', cmap_name=color_scheme[1], limits=color_scheme[2], breaks = [0, 0.5, 1], labels = [downs_txt, noDEGs_txt, ups_txt])
        elif color_scheme[0] == 'gradient':
            return scale_color_gradientn(name='  ', colors=color_scheme[1], limits=color_scheme[2], breaks = [0, 0.5, 1], labels = [downs_txt, noDEGs_txt, ups_txt])

    themes = {'classic': theme_classic(), 'bw': theme_bw(), 'dark': theme_dark(), 'light': theme_light(), 'minimal': theme_minimal(), 'seaborn': theme_seaborn(), 'linedraw': theme_linedraw(), 'gray': theme_gray(), 'void': theme_void(), 'xkcd': theme_xkcd()}

    def decideindicator(threshold_indicator):
        if threshold_indicator:
            return [geom_vline(xintercept=x_threshold, color="black", linetype="dashed"), geom_vline(xintercept=-x_threshold, color="black", linetype="dashed"), geom_hline(yintercept=-np.log10(y_threshold), color="black", linetype="dashed")]
        else:
            return None

    return (
        ggplot(data=plot_data, mapping=aes(x="log2FC", y="-log10(Pvalue)", color="cmap")) 
        + geom_point() 
        + themes[theme]
        + xlab("log2FC") 
        + ylab("-log10(Pvalue)")
        + xlim(-x_lim, x_lim)
        + ylim(None, y_limup)
        + getcolorscheme(color_scheme)
        + decideindicator(threshold_indicator)
        # + theme(text=element_text(family='Times New Roman'))
    )

def EnrichmentScatter(data, color_scheme=['gradient', ['#f57f74', '#82cc5e']], sorting = {'column':2, 'ascending': False}, 
                vars = ['Rich Ratio', 'Pathway', "p-value", 'Gene Number'], theme_var=theme_bw()):
    """
    vars[4]: x, y, scatter_color, scatter_size
                ^  ^  ^              ^
    sorting  0  1  2              3
    """
    if sorting != None:
        data.sort_values(by=vars[sorting['column']], inplace=True, ascending=sorting['ascending'])
        ordered_categories = data.sort_values(by=vars[sorting['column']], ascending=sorting['ascending'])[vars[1]].unique()
        data[vars[1]] = pd.Categorical(
            data[vars[1]], 
            categories=ordered_categories,
            ordered=True
        )
    def getcolorscheme(color_scheme):
        # if color_scheme[0] == 'discrete':
        #     return scale_color_manual(values=color_scheme[1])
        if color_scheme[0] == 'cmap':
            return scale_color_cmap(name=vars[2], cmap_name=color_scheme[1])
        elif color_scheme[0] == 'gradient':
            return scale_color_gradientn(name=vars[2], colors=color_scheme[1])
    plot = (
        ggplot(data, aes(x=vars[0], y=vars[1], color=vars[2], size=vars[3])) + 
        geom_point() +
        # theme +
        getcolorscheme(color_scheme=color_scheme) +
        scale_size_continuous(range=[2,8]) +
        theme_var + 
        theme(figure_size=(9,5))
    )
    return plot

# def ClusteredHeatmap(data, theme_var=theme_bw(), color_scheme=['gradient', ['#4e6691','#eeeeee','#ba3e45']], cluster_columns=True, cluster_rows=True, cluster_type="cladogram", branch_width=1, standardization='Row', log_scale=False):


#     # log-scale transformation
#     if log_scale:
#         data = np.log10(data + 1)

#     # z-score standardization
#     if standardization == 'Row':
#         data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
#     elif standardization == 'Column':
#         data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
#     elif standardization == None or standardization == 'None':
#         pass
#     else:
#         # warning
#         logging.warning("Standardization method not found. Using None.")

#     # use seaborn to plot
#     if color_scheme[0] == 'cmap':
#         colors = color_scheme[1]
#     if color_scheme[0] == 'gradient':
#         # transit discrete colors to gradient colors
#         from matplotlib.colors import LinearSegmentedColormap
        
#         colors = LinearSegmentedColormap(color_scheme[1])
#         # transfer to list
#         colors = [colors(i) for i in range(colors.N)]
#         colors = sns.color_palette(color_scheme[1], as_cmap=True)

#     plot = sns.clustermap(data, figsize=(4,4), cmap=colors, row_cluster=cluster_rows, col_cluster=cluster_columns, tree_kws={'linewidth': branch_width})
#     return plot



# data = pd.DataFrame([
#     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
#     [2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
#     [3, 4, 5, 6, 7, 8, 9, 10, 11, 12.5],
# ])

# import matplotlib.pyplot as plt
# fig = ClusteredHeatmap(data)
# plt.show()