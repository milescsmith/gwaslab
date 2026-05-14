import gc as garbage_collect

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns
from adjustText import adjust_text
from allel import GenotypeArray, read_vcf, rogers_huff_r_between
from matplotlib import ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats

import gwaslab as gl
from gwaslab.bd.bd_common_data import get_chr_to_number, get_number_to_chr, get_recombination_rate
from gwaslab.g_Sumstats import Sumstats
from gwaslab.info.g_Log import Log
from gwaslab.info.g_version import _get_version
from gwaslab.util.util_in_get_sig import _anno_gene, _get_sig
from gwaslab.viz.viz_aux_annotate_plot import annotate_pair
from gwaslab.viz.viz_aux_quickfix import (
    _get_largenumber,
    _quick_add_tchrpos,
    _quick_assign_highlight_hue_pair,
    _quick_assign_i_with_rank,
    _quick_assign_marker_relative_size,
    _quick_extract_snp_in_region,
    _quick_fix,
    _quick_merge_sumstats,
)
from gwaslab.viz.viz_aux_reposition_text import adjust_text_position
from gwaslab.viz.viz_aux_save_figure import safefig, save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style
from gwaslab.viz.viz_plot_mqqplot import _mqqplot


@safefig
def plot_miami2(
          path1=None,
          path2=None,
          merged_sumstats=None,
          suffixes = None,
          cols=None,
          cols1=None,
          cols2=None,
          id0="TCHR+POS",
          id1=None,
          id2=None,
          mode="m",
          chr_dict=None,
          chr_dict1=None,
          chr_dict2=None,
          scaled=False,
          scaled1=False,
          scaled2=False,
          titles=None,
          titles_pad=None,
          use_rank=False,
          chrpad=0.03,
          region_hspace = 0.1,
          dpi=100,
          fontsize = 10,
          xtick_label_size = 10,
          font_family="Arial",
          xlabel_coords=None,
          xtick_label_pad=None,
          verbose=True,
          xtickpad=None,
          fig_kwargs=None,
          scatter_kwargs=None,
          same_ylim=False,
          save=None,
          save_kwargs=None,
          figax=None,
          build=None,
          log=Log(),
          **mqq_kwargs
          ):
    log.write(f"Start to create miami plot {_get_version()}:", verbose=verbose)

    # Auto-detect build from Sumstats objects if not provided
    if build is None:
        build1 = None
        build2 = None

        # Extract build from both Sumstats objects
        if isinstance(path1, Sumstats) and hasattr(path1, "build"):
            build1 = path1.build
        if isinstance(path2, Sumstats) and hasattr(path2, "build"):
            build2 = path2.build

        # Check consistency and set build
        if build1 is not None and build2 is not None:
            if build1 != build2:
                raise ValueError(f"Build version mismatch: sumstats1={build1}, sumstats2={build2}. Please ensure both sumstats use the same genome build.")
            build = build1
            log.write(f" -Auto-detected build from sumstats: {build}", verbose=verbose)
        elif build1 is not None:
            build = build1
            log.write(f" -Auto-detected build from sumstats1: {build}", verbose=verbose)
        elif build2 is not None:
            build = build2
            log.write(f" -Auto-detected build from sumstats2: {build}", verbose=verbose)

    ## figuring arguments ###########################################################################################################
    # figure columns to use (auto-detect scaled)
    if scaled is True:
        scaled1 = True
        scaled2 = True
    else:
        try:
            if isinstance(path1, Sumstats) and ("MLOG10P" in path1.data.columns):
                scaled1 = True
        except Exception:
            pass
        try:
            if isinstance(path2, Sumstats) and ("MLOG10P" in path2.data.columns):
                scaled2 = True
        except Exception:
            pass

    if cols is None:
        cols = ["CHR","POS","MLOG10P" if scaled or (scaled1 and scaled2) else "P"]

    if cols1 is None:
        cols1 = ["CHR","POS","MLOG10P" if scaled1 else "P"]

    if cols2 is None:
        cols2 = ["CHR","POS","MLOG10P" if scaled2 else "P"]

    # Auto-detect id1/id2 from Sumstats objects when not specified
    # This is needed for features like highlight, pinpoint, and annotation
    # that require SNP IDs to match variants by rsID/SNPID
    if id1 is None and isinstance(path1, Sumstats):
        for col in ["SNPID", "rsID"]:
            if col in path1.data.columns:
                id1 = col
                log.write(f" -Auto-detected id1='{col}' from sumstats1", verbose=verbose)
                break
    if id2 is None and isinstance(path2, Sumstats):
        for col in ["SNPID", "rsID"]:
            if col in path2.data.columns:
                id2 = col
                log.write(f" -Auto-detected id2='{col}' from sumstats2", verbose=verbose)
                break

    if id1 is not None:
        cols1.append(id1)
    if id2 is not None:
        cols2.append(id2)

    if chr_dict is None:
        chr_dict  = get_chr_to_number()
    if chr_dict1 is None:
        chr_dict1 = chr_dict
    if chr_dict2 is None:
        chr_dict2 = chr_dict



    if scatter_kwargs is None:
        scatter_kwargs = {}

    style = set_plot_style(
        plot="plot_miami",
        mode=mode,
        fig_kwargs=fig_kwargs,
        save_kwargs=save_kwargs,
        save=save,
        scatter_kwargs=scatter_kwargs,
        fontsize=fontsize,
        fontfamily=font_family,
        dpi=dpi,
        verbose=verbose,
        log=log,
    )
    fig_kwargs = style["fig_kwargs"]
    save_kwargs = style["save_kwargs"]
    scatter_kwargs = style["scatter_kwargs"]
    fontsize = style["fontsize"]
    font_family = style["font_family"]
    dpi = style["dpi"]

    if xlabel_coords is None:
        xlabel_coords = (-0.01,- region_hspace/2 )

    # figure out mqq args
    mqq_kwargs1,mqq_kwargs2 = _sort_kwargs_to_12(mqq_kwargs)

    # figure out args for pdf and svg
    from gwaslab.viz.viz_aux_style_options import figure_kwargs_for_vector_plot
    fig_kwargs, scatter_kwargs = figure_kwargs_for_vector_plot(save, fig_kwargs, scatter_kwargs)

    # add suffix if ids are the same
    id1_1, id2_2, mqq_kwargs1, mqq_kwargs2 = _solve_id_contradictory(id0, id1, id2, mqq_kwargs1, mqq_kwargs2)

    if titles is None:
        titles=["",""]

    titles_pad_adjusted=[1,0]
    if titles_pad is None:
        titles_pad=[0.2,0.2]
        if "anno" in mqq_kwargs.keys():
            if mqq_kwargs["anno"] is not None:
                titles_pad_adjusted= [1 + titles_pad[0], -titles_pad[1]]
        if "anno1" in mqq_kwargs.keys():
            if mqq_kwargs["anno1"] is not None:
                titles_pad_adjusted[0]= 1 + titles_pad[0]
        if "anno2" in mqq_kwargs.keys():
            if mqq_kwargs["anno2"] is not None:
                titles_pad_adjusted[1]=  - titles_pad[1]
    else:
        titles_pad_adjusted=[1 + titles_pad[0], -titles_pad[1]]

    if merged_sumstats is None:
        sumstats1 = _figure_type_load_sumstats(name="Sumstats1",
                                               path=path1,
                                               cols=cols1,
                                               log=log,
                                               verbose=verbose)
        sumstats2 = _figure_type_load_sumstats(name="Sumstats2",
                                               path=path2,
                                               cols=cols2,
                                               log=log,
                                               verbose=verbose)
    elif suffixes is not None:
        cols1[2] += suffixes[0]
        cols2[2] += suffixes[1]
        #sumstats1 = merged_sumstats[cols1].copy()
        #sumstats2 = merged_sumstats[cols2].copy()

    ## quick fix ###########################################################################################################
    if merged_sumstats is None:
        sumstats1 = _quick_fix(sumstats1,chr_dict=chr_dict1, scaled=scaled1, verbose=verbose, log=log)
        sumstats2 = _quick_fix(sumstats2,chr_dict=chr_dict2, scaled=scaled2, verbose=verbose, log=log)

        # get a large number ###########################################################################################################
        large_number = _get_largenumber(sumstats1["POS"].max(), sumstats2["POS"].max(),log=log)
    else:
        large_number = _get_largenumber(merged_sumstats["POS"].max(), merged_sumstats["POS"].max(),log=log)

    ## create merge index ###########################################################################################################
    if merged_sumstats is None:
        sumstats1 = _quick_add_tchrpos(sumstats1,large_number=large_number, dropchrpos=False, verbose=verbose, log=log)
        sumstats2 = _quick_add_tchrpos(sumstats2,large_number=large_number, dropchrpos=False, verbose=verbose, log=log)
        log.write(" -Merging sumstats using chr and pos...",verbose=verbose)

        ###### merge #####################################################################################################
        merged_sumstats = _quick_merge_sumstats(sumstats1=sumstats1,sumstats2=sumstats2)

    #assign index for x axis
    merged_sumstats, chrom_df = _quick_assign_i_with_rank(merged_sumstats,
                                                          chrpad=chrpad,
                                                          use_rank=use_rank,
                                                          chrom="CHR",
                                                          pos="POS",
                                                          drop_chr_start=False)

    # P_1  scaled_P_1  P_2  scaled_P_2  TCHR+POS CHR POS
    log.write(" -Columns in merged sumstats: {}".format(",".join(merged_sumstats.columns)), verbose=verbose)

    if merged_sumstats is None:
        del(sumstats1)
        del(sumstats2)
    garbage_collect.collect()
    #####################################################################################################################


    ##plotting
    if figax is None:
        #fig_kwargs["figsize"] = (15,10)
        fig, (ax1, ax5) = plt.subplots(2, 1, gridspec_kw={"height_ratios": [1, 1]},**fig_kwargs)
        plt.subplots_adjust(hspace=region_hspace)
    else:
        fig, ax1, ax5 = figax
    ##########################################################################################################################
    bbox1 = ax1.get_position()
    bbox5 = ax5.get_position()



    fig_height_inches = fig.get_figheight()
    ax_height_inches = (bbox1.height + bbox5.height) * fig_height_inches /2
    ax_height_points = ax_height_inches * 72
    tick = ax1.xaxis.get_major_ticks()[0].tick1line
    # Tick length is determined by its 'markersize' (in points)
    tick_length = tick.get_markersize()

    tick_points_to_pixels = tick_length * fig.dpi / 72.0
    ax_height_pixels = bbox1.height * fig.get_figheight() * fig.dpi
    tick_axes_length = tick_points_to_pixels / ax_height_pixels

    if xtick_label_pad is None:
        if "figsize" not in fig_kwargs.keys():
            fig_kwargs["figsize"] = (15,10)
        # (total hsapce - tick label font size) / 2
        xtick_label_pad = 0
        #xtick_label_pad =  ((ax_height_points * region_hspace) - 2*tick_length - xtick_label_size) / 2
    ########################################################################################################################
    #if same_ylim==True:
        #maxy = merged_sumstats[["scaled_P_1","scaled_P_2"]].max().max()

    log.write("Start to create Manhattan plot for sumstats1...", verbose=verbose)
    fig,log = _mqqplot(merged_sumstats,
                    chrom="CHR",
                    pos="POS",
                    p="P_1",
                    mlog10p="scaled_P_1",
                    snpid=id1_1,
                    scaled=scaled1,
                    log=log,
                    mode=mode,
                    figax=(fig,ax1),
                    scatter_kwargs=scatter_kwargs,
                    _chrom_df_for_i = chrom_df,
                    _invert=False,
                    _if_quick_qc=False,
                    font_family=font_family,
                    build=build,
                    **mqq_kwargs1
                    )
    log.write("Finished creating Manhattan plot for sumstats1", verbose=verbose)

    log.write("Start to create Manhattan plot for sumstats2...", verbose=verbose)
    fig,log = _mqqplot(merged_sumstats,
                      chrom="CHR",
                      pos="POS",
                      p="P_2",
                      mlog10p="scaled_P_2",
                      scaled=scaled2,
                      snpid=id2_2,
                      log=log,
                      mode=mode,
                      figax=(fig,ax5),
                      scatter_kwargs=scatter_kwargs,
                      _chrom_df_for_i = chrom_df,
                       _invert=True,
                      _if_quick_qc=False,
                      font_family=font_family,
                      build=build,
                     **mqq_kwargs2)
    log.write("Finished creating Manhattan plot for sumstats2", verbose=verbose)


    #####################################################################################################################
    ax1l, ax1r = ax5.get_xlim()
    ax5l, ax5r = ax1.get_xlim()
    ax1.set_xlim([min(ax1l,ax5l), max(ax1r,ax5r)])
    ax5.set_xlim([min(ax1l,ax5l), max(ax1r,ax5r)])

    #####################################################################################################################
    ax5.set_xlabel("")
    #ax5.set_xticks(chrom_df)
    ax5.set_xticklabels([])
    ax5.xaxis.set_ticks_position("top")
    ax5.tick_params(axis="x", which="major", pad=0)

    # Ad#just the visibility for spines #######################################################
    ax1, ax5 = _set_spine_visibility(ax1, ax5)
    ######################################################################################################################
#####################################################################################################################
    # set labels
    ax1.set_xlabel("Chromosome",fontsize=fontsize,family=font_family,labelpad=0, va="center",ha="center")
    ax1.xaxis.set_label_coords(xlabel_coords[0],xlabel_coords[1])

    #ax1.tick_params(axis='x', which='major', pad=xtick_label_pad, labelsize = xtick_label_size)

    for label in ax1.get_xticklabels():
        label.set_y( xlabel_coords[1] + tick_axes_length )
    ax1.tick_params(axis="x", which="major", pad=xtick_label_pad, labelsize = xtick_label_size)
    plt.setp(ax1.get_xticklabels(),  ha="center",va="center")

    ax1.set_ylabel(r"$\mathregular{-log_{10}(P)}$",fontsize=fontsize,family=font_family)
    ax5.set_ylabel(r"$\mathregular{-log_{10}(P)}$",fontsize=fontsize,family=font_family)

    ax1.set_title(titles[0],y=titles_pad_adjusted[0],family=font_family)
    ax5.set_title(titles[1],y=titles_pad_adjusted[1],family=font_family)

    ax5.invert_yaxis()

    save_figure(fig, save, keyword="miami",save_kwargs=save_kwargs, log=log, verbose=verbose)

    garbage_collect.collect()

    log.write("Finished creating miami plot successfully", verbose=verbose)
    #Return matplotlib figure object #######################################################################################
    return fig, log

def _sort_kwargs_to_12(mqq_kwargs):
    mqq_kwargs1={}
    mqq_kwargs2={}
    for key, value in mqq_kwargs.items():
        if key[-1]=="1":
            # for top mqq
            mqq_kwargs1[key[:-1]] = value
        elif key[-1]=="2":
            # for bottom mqq
            mqq_kwargs2[key[:-1]] = value
        else:
            # for both top and bottom
            mqq_kwargs1[key] = value
            mqq_kwargs2[key] = value
    return mqq_kwargs1, mqq_kwargs2

def _solve_id_contradictory(id0, id1, id2, mqq_kwargs1, mqq_kwargs2):
    if (id1 is not None) and (id2 is not None):
        if id1 == id2:
            id1_1 = id1 + "_1"
            id2_2 = id2 + "_2"
            if "anno" in mqq_kwargs1.keys():
                if mqq_kwargs1["anno"] == id1:
                    mqq_kwargs1["anno"] = id1_1
            if "anno" in mqq_kwargs2.keys():
                if mqq_kwargs2["anno"] == id2:
                    mqq_kwargs2["anno"] = id2_2
        else:
            id1_1 = id1
            id2_2 = id2

    if id1 is None:
        id1_1 = id0

    if id2 is None:
        id2_2 = id0

    return (id1_1, id2_2, mqq_kwargs1, mqq_kwargs2)

def _figure_kwargs_for_vector_plot(save, fig_kwargs, scatter_kwargs ):
    from gwaslab.viz.viz_aux_style_options import figure_kwargs_for_vector_plot
    return figure_kwargs_for_vector_plot(save, fig_kwargs, scatter_kwargs)

def _set_spine_visibility(ax1,ax5):
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_visible(True)
    ax1.spines["bottom"].set_visible(True)

    ax5.spines["top"].set_visible(True)
    ax5.spines["right"].set_visible(False)
    ax5.spines["left"].set_visible(True)
    ax5.spines["bottom"].set_visible(False)
    return ax1,ax5

def _figure_type_load_sumstats(name, path, cols, log, verbose):
    log.write(f" -Obtaining {name} CHR, POS, P and annotation from: {cols}", verbose=verbose)
    if isinstance(path, Sumstats):
        log.write(f" -Loading {name} from gwaslab.Sumstats Object", verbose=verbose)
        return path.data[cols].copy()
    raise TypeError(f"{name} must be a gwaslab.Sumstats")
