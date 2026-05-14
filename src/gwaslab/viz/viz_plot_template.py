from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import seaborn as sns

from gwaslab.info.g_Log import Log
from gwaslab.io.io_process_kwargs import _update_arg, _update_kwargs
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

def _plot(
    associations: Any,
    values: str = "Beta",
    sort: str = "P-value",
    mode: str = "gwaslab",
    fontsize: float = 12,
    font_family: str = "Arial",
    cmap: Any | None = None,
    ylabel: str = "Y",
    xlabel: str = "X",
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    yticks: list[float] | None = None,
    xticks: list[float] | None = None,
    ytick_labels: list[str] | None = None,
    xtick_labels: list[str] | None = None,
    title: str = "Title",
    title_pad: float = 1.08,
    title_fontsize: float = 13,
    title_kwargs: dict[str, Any] | None = None,
    linewidth: float | None = None,
    linestyle: str | None = None,
    linecolor: str | None = None,
    dpi: int = 200,
    anno_kwargs: dict[str, Any] | None = None,
    err_kwargs: dict[str, Any] | None = None,
    fig_kwargs: dict[str, Any] | None = None,
    scatter_kwargs: dict[str, Any] | None = None,
    save: Union[bool, str] | None = None,
    save_kwargs: dict[str, Any] | None = None,
    log: Log = Log(),
    verbose: bool = True,
    **args: Any
) -> tuple["Figure", "Axes"]:

    # update args

    cmap = _update_arg(cmap, "RdBu")

    style = set_plot_style(
        plot="plot_template",
        fig_kwargs=fig_kwargs or {"figsize":(10,10)},
        save_kwargs=save_kwargs,
        fontsize=fontsize,
        fontfamily=font_family,
        verbose=verbose,
        log=log,
    )
    fig,ax = plt.subplots(**(style.get("fig_kwargs", {})))

    # draw lines
    #horizontal_line = ax.axhline(y=1, linewidth = linewidth,
    #                        linestyle="--",
    #                        color=linecolor,zorder=1000)
    #vertical_line   = ax.axvline(x=1, linewidth = linewidth,
    #                        linestyle="--",
    #                        color=linecolor,zorder=1000)

    # ticks
    if xticks is not None:
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks,
                            fontsize=fontsize,
                            family=font_family)

    if yticks is not None:
        ax.set_xticks(yticks)
        ax.set_xticklabels(yticks,
                            fontsize=fontsize,
                            family=font_family)

    # labels
    ax.set_xlabel(xlabel,
                  fontsize=fontsize,
                  fontfamily=font_family)
    ax.set_ylabel(ylabel,
                  fontsize=fontsize,
                  fontfamily=font_family)

    ax.tick_params(axis="x",
                    labelsize=fontsize,
                    labelfontfamily=font_family)
    ax.tick_params(axis="y",
                    labelsize=fontsize,
                    labelfontfamily=font_family)

    # title
    title_pad = title_pad -0.05
    fig.suptitle(title ,
                 fontsize = title_fontsize,
                 x=0.5,
                 y=title_pad)

    # spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)

    save_kwargs = style.get("save_kwargs", style.get("save_kwargs", {}))
    save_figure(fig = fig, save = save, keyword=mode, save_kwargs=save_kwargs, log = log, verbose=verbose)

    log.write("Finished creating plots.", verbose=verbose)
    return fig, ax
