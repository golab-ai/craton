import base64
from itertools import zip_longest
from pathlib import Path
from typing import Dict, Union

import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import plotly.graph_objects as go
#from alchemlyb.estimators import BAR, MBAR, TI
from matplotlib.font_manager import FontProperties
from matplotlib.pyplot import MultipleLocator
from plotly.subplots import make_subplots

plt.style.use("ggplot")


def split_edge(edge_str):
    edge_list = edge_str.split("-")
    num = len(edge_list)
    if not num % 2:
        start = "-".join(edge_list[: num // 2])
        end = "-".join(edge_list[num // 2 :])
    else:
        raise RuntimeError("program cannot determine the node from edge string: {}".format(edge_str))
    return (start, end)


def plot_rmsd_rbfe(
    time,
    l0_ligand_rmsd,
    l0_protein_rmsd,
    l0_ligand_align_protein,
    l1_ligand_rmsd,
    l1_protein_rmsd,
    l1_ligand_align_protein,
    ax=None,
):
    """FEP result RMSD for ligand and protein"""
    if ax is None:
        fig, ax = plt.subplots(2, 1)

    ax[0].plot(time, l0_ligand_rmsd, label="Ligand")
    ax[0].plot(time, l0_protein_rmsd, label="Protein")
    ax[0].plot(time, l0_ligand_align_protein, label="Ligand after aligning")
    ax[0].set_xlim(left=0)
    ax[0].set_ylim(bottom=0)
    ax[0].set_xlabel("Time (ps)")
    ax[0].set_ylabel(r"RMSD ($\mathrm{\AA}$)")

    ax[1].plot(time, l1_ligand_rmsd, label="Ligand")
    ax[1].plot(time, l1_protein_rmsd, label="Protein")
    ax[1].plot(time, l1_ligand_align_protein, label="Ligand after aligning")
    ax[1].set_xlim(xmin=0)
    ax[1].set_ylim(ymin=0)
    ax[1].set_xlabel("Time (ps)")
    ax[1].set_ylabel(r"RMSD ($\mathrm{\AA}$)")

    plt.legend()
    return fig

def plot_exchange_rate(df, output):
    fig = go.Figure()
    row = len(df)
    for column in df.columns:
        fig.add_trace(go.Bar(x=np.arange(row) + 0.5, y=df[column], width=np.ones(row), name=column))
    fig.update_xaxes(tickvals=np.arange(row) + 0.5, ticktext=[f"lambda_{i}" for i in range(row)])
    # fig.update_yaxes(range=[0, self.count], visible=False)
    fig.update_layout(barmode="stack")
    fig.write_html(output)


def plot_rmsd_rhfe(l0_ligand_rmsd, l1_ligand_rmsd):
    """FEP result RMSD for ligand and protein

    Args:
        l0_ligand_rmsd (2d array): Liangd RMSD for lambda = 0 window
        l0_align_rmsd (2d array): Ligand RMSD and protein RMSD after aligning the protein for lambda = 0 window
        l1_ligand_rmsd (2d array): Liangd RMSD for lambda = 1 window
        l1_align_rmsd (2d array): Ligand RMSD and protein RMSD after aligning the protein for lambda = 1 window
    """

    fig, ax = plt.subplots(1, 1)

    ax.plot(l0_ligand_rmsd[:, 1], l0_ligand_rmsd[:, 2], label="lambda = 0")
    ax.plot(l1_ligand_rmsd[:, 1], l1_ligand_rmsd[:, 2], label="lambda = 1")
    ax.set_xlim(xmin=0)
    ax.set_ylim(ymin=0)

    plt.legend()
    return fig


def plot_correlation(calc, exp, calc_err=None, exp_err=None, system="", ax=None):

    """plot correlation between calcculation result and experiments result

    Args:
        calc_data : A dict object:
            key is the pair name
            value contains rbfe and std
        exp_data: A dict object has the same format as calc
    """

    import math

    from scipy import stats

    # calculate mue
    sumdev = 0.0
    for i in range(len(exp)):
        sumdev += abs(exp[i] - calc[i])
    mue = sumdev / len(exp)

    # calculate rmse
    sumdev = 0.0
    for i in range(len(exp)):
        sumdev += (exp[i] - calc[i]) * (exp[i] - calc[i])
    rmse = math.sqrt(sumdev / len(exp))

    # best fit
    b, a, r_value, p_value, std_err = stats.linregress(exp, calc)
    r2 = r_value * r_value

    # plot
    if ax is None:
        fig, ax = plt.subplots()

    ax.scatter(exp, calc)
    ax.errorbar(exp, calc, xerr=exp_err, yerr=calc_err, linestyle="None")

    # fit line
    t = np.linspace(min(exp), max(exp), 100)
    yfit = [a + b * i for i in t]
    ax.plot(t, yfit, color="C4", linestyle="-.")

    # set axis limit
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    min_lim = min(xlim[0], ylim[0])
    max_lim = max(xlim[1], ylim[1])

    ax.set_xlim(min_lim, max_lim)
    ax.set_ylim(min_lim, max_lim)

    # set axis interval
    x_major_locator = MultipleLocator(1)
    y_major_locator = MultipleLocator(1)
    ax.xaxis.set_major_locator(x_major_locator)
    ax.yaxis.set_major_locator(y_major_locator)

    # add y=x, x-2, x+2
    x = np.arange(min_lim, max_lim, 0.01)
    y_middle = x
    y_up = x + 2
    y_down = x - 2

    ax.plot(x, y_middle, "k--", alpha=0.7, zorder=0)
    ax.plot(x, y_up, "k--", alpha=0.7, zorder=0)
    ax.plot(x, y_down, "k--", alpha=0.7, zorder=0)
    ax.fill_between(x, y_middle, y_up, facecolor="mediumseagreen", alpha=0.1)
    ax.fill_between(x, y_down, y_middle, facecolor="mediumseagreen", alpha=0.1)
    ax.set_aspect("equal")

    # add text at the top left coner of the plot
    x_text = min_lim + 0.05 * (max_lim - min_lim)
    y_text = min_lim + 0.95 * (max_lim - min_lim)

    # t = "Mue:%3.2f\nRMSE:%3.2f"%(mue, rmse)
    t = "Mue:%3.2f\nRMSE:%3.2f\n$R^2:%3.2f$" % (mue, rmse, r2)
    ax.text(
        x_text,
        y_text,
        t,
        ha="left",
        va="top",
        bbox=dict(facecolor="lightgreen", alpha=0.5),
        fontsize=15,
    )
    ax.tick_params(axis="both", labelsize="large")

    # set axis label
    ax.set_xlabel(r"Exp. [kcal/mol]", size="large")
    ax.set_ylabel(r"Cal. [kcal/mol]", size="large")

    if system:
        ax.set_title(r"System: %s" % system, size="x-large")

    return ax


def show_cycle_closure(corr_data: Dict[str, Dict[str, Dict[str, float]]], ax: matplotlib.axes) -> matplotlib.axes:
    G = nx.DiGraph()
    edge_labels_dict = {}
    node_labels_dict = {}

    for edge in corr_data["edge"].items():
        key = edge[0]
        ddG = edge[1]["ddG"]
        ddG_std = edge[1]["std"]
        start, end = split_edge(key)
        G.add_edge(start, end)
        edge_labels_dict[(start, end)] = "{:.1f}\n({:.2f})".format(ddG, ddG_std)

    for node in G.nodes():
        node_labels_dict[node] = "{}\n{:.1f}\n({:.2f})".format(
            node, corr_data["node"][node]["dG"], corr_data["node"][node]["std"]
        )

    # pos = nx.spring_layout(G)
    pos = nx.nx_pydot.graphviz_layout(G, prog="dot")
    nx.draw(
        G,
        pos,
        edge_color="black",
        width=2,
        linewidths=1,
        node_size=1400,
        node_color="none",
        alpha=0.9,
        arrows=True,
        node_shape="s",
        ax=ax,
    )
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels_dict, font_color="red", font_size=18)
    nx.draw_networkx_labels(G, pos, labels=node_labels_dict, font_size=18)
    ax.axis("off")
    return ax


def show_ligand_dG(dG_dict: Dict[str, float], ax: matplotlib.axes = None) -> matplotlib.axes:
    dG_sorted_dict = dict(sorted(dG_dict.items(), key=lambda item: item[1]))
    ligand, dG = zip(*dG_sorted_dict.items())
    if ax is None:
        _, ax = plt.subplots(1, 1)
    num = len(ligand)
    ax.bar(range(num), dG)
    ax.set_ylabel("dG (kcal/mol)")
    ax.set_xticks(range(num))
    ax.set_xticklabels(ligand, rotation=90)
    return ax


def plot_convergence_block(forward, forward_error, ax=None, units="kcal/mol"):
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))
    n = len(forward)
    x = np.linspace(0, 1, n + 1)[1:]
    ax.errorbar(x, forward, yerr=forward_error, lw=3, marker="o")
    plt.xticks(x[::2], fontsize=10)
    plt.yticks(fontsize=10)

    ax.set_xlabel(r"Fraction of the molecule_dynamics time", fontsize=16)
    ax.set_ylabel(r"$\Delta G$ ({})".format(units), fontsize=16)
    plt.xticks(x, ["%.2f" % i for i in x])
    return ax


def plot_convergence(forward, forward_error, backward, backward_error, units="kcal/mol", ax=None):
    if ax is None:  # pragma: no cover
        _, ax = plt.subplots(figsize=(8, 6))

    plt.setp(ax.spines["bottom"], color="#D2B9D3", lw=3, zorder=-2)
    plt.setp(ax.spines["left"], color="#D2B9D3", lw=3, zorder=-2)

    for dire in ["top", "right"]:
        ax.spines[dire].set_color("none")

    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    f_ts = np.linspace(0, 1, len(forward) + 1)[1:]
    r_ts = np.linspace(0, 1, len(backward) + 1)[1:]

    ax.fill_between(
        [0, 1],
        backward.values[-1] - backward_error.values[-1],
        backward.values[-1] + backward_error.values[-1],
        color="#D2B9D3",
        zorder=1,
    )
    line1 = ax.errorbar(
        f_ts,
        forward,
        yerr=forward_error,
        color="#736AFF",
        lw=3,
        zorder=2,
        marker="o",
        mfc="w",
        mew=2.5,
        mec="#736AFF",
        ms=12,
    )
    line2 = ax.errorbar(
        r_ts,
        backward,
        yerr=backward_error,
        color="#C11B17",
        lw=3,
        zorder=3,
        marker="o",
        mfc="w",
        mew=2.5,
        mec="#C11B17",
        ms=12,
    )

    plt.xticks(r_ts[::2], fontsize=10)
    plt.yticks(fontsize=10)

    ax.legend(
        (line1[0], line2[0]),
        ("Forward", "Reverse"),
        loc=9,
        prop=FontProperties(size=18),
        frameon=False,
    )
    ax.set_xlabel(r"Fraction of the molecule_dynamics time", fontsize=16, color="#151B54")
    ax.set_ylabel(r"$\Delta G$ ({})".format(units), fontsize=16, color="#151B54")
    plt.xticks(f_ts, ["%.2f" % i for i in f_ts])
    plt.tick_params(axis="x", color="#D2B9D3")
    plt.tick_params(axis="y", color="#D2B9D3")
    return ax


def add_molecule_figure(fig, svg_str, row, col):
    src = "data:image/svg+xml;base64,{}".format(base64.b64encode(bytes(svg_str, "utf-8")).decode())
    fig.add_layout_image(
        dict(
            source=src,
            xref="x",
            yref="y",
            x=0,
            y=2,
            sizex=2,
            sizey=2,
            sizing="stretch",
            layer="above",
        ),
        row=row,
        col=col,
    )
    fig.update_xaxes(range=[0, 2], showticklabels=False, showgrid=False, row=row, col=col, visible=False)
    fig.update_yaxes(range=[0, 2], showticklabels=False, showgrid=False, row=row, col=col, visible=False)


def add_torsion_figure(fig, energy, rbfe, rhfe, row, col):
    x = list(energy.index)
    fig.add_trace(
        go.Scatter(
            x=x, y=energy - min(energy), name="energy", line_color="blue", opacity=0.6, line_width=1, marker_size=2
        ),
        secondary_y=True,
        row=row,
        col=col,
    )
    fig.add_trace(
        go.Histogram(
            x=rbfe,
            histnorm="probability",
            nbinsx=20,
            name="sampling",
            opacity=0.8,
            marker_color="red",
        ),
        secondary_y=False,
        row=row,
        col=col,
    )
    fig.add_trace(
        go.Histogram(
            x=rhfe,
            histnorm="probability",
            nbinsx=20,
            name="sampling",
            opacity=0.8,
            marker_color="pink",
        ),
        secondary_y=False,
        row=row,
        col=col,
    )
    fig.add_vline(rbfe[0], row=row, col=col)


def generate_torsion_figure_4col(
    energy,
    rbfe,
    rhfe,
    xvg,
    energy2,
    rbfe2,
    rhfe2,
    xvg2,
    show=False,
    title=None,
    filename: Union[str, Path] = "torsion_distriubiton.html",
):
    num_row = max(len(energy.columns), len(energy2.columns))
    spec_false, spec_true = {"secondary_y": False}, {"secondary_y": True}
    specs = [[spec_false, spec_true, spec_false, spec_true]] * num_row
    fig = make_subplots(num_row, 4, specs=specs, horizontal_spacing=0.05, vertical_spacing=0.05)
    for i, (key1, key2) in enumerate(zip_longest(energy.keys(), energy2.keys())):
        if key1 is not None:
            add_molecule_figure(fig, xvg[i], i + 1, 1)
            add_torsion_figure(fig, energy[key1], rbfe[key1], rhfe[key1], i + 1, 2)
        if key2 is not None:
            add_molecule_figure(fig, xvg2[i], i + 1, 3)
            add_torsion_figure(fig, energy2[key2], rbfe2[key2], rhfe2[key2], i + 1, 4)
    fig.update_layout(
        margin={"t": 10, "l": 0, "b": 0, "r": 0},
        showlegend=False,
        template="plotly_white",
        font=dict(size=5),
        title=title,
        height=num_row * 80,
        width=400,
    )
    if not show:
        fig.write_html(filename)
    else:
        fig.show()


def generate_torsion_figure(
    energy, rbfe, rhfe, xvg, show=False, title=None, filename: Union[str, Path] = "distribution.html"
):
    num_row = len(energy.columns)
    spec_false, spec_true = {"secondary_y": False}, {"secondary_y": True}
    specs = [[spec_false, spec_true]] * num_row
    fig = make_subplots(num_row, 2, specs=specs, horizontal_spacing=0.05, vertical_spacing=0.05)
    for i, key in enumerate(energy.keys()):
        add_molecule_figure(fig, xvg[i], i + 1, 1)
        add_torsion_figure(fig, energy[key], rbfe[key], rhfe[key], i + 1, 2)
    fig.update_layout(
        margin={"t": 10, "l": 0, "b": 0, "r": 0},
        showlegend=False,
        template="plotly_white",
        font=dict(size=5),
        title=title,
        height=num_row * 80,
        width=250,
    )
    if not show:
        fig.write_html(filename)
    else:
        fig.show()


def plot_df_state(estimators, ax=None, units="kcal/mol", colors=None, labels=None):

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))

    try:
        len(estimators)
    except TypeError:
        estimators = [
            estimators,
        ]

    formatted_data = []
    for dhdl in estimators:
        try:
            len(dhdl)
            formatted_data.append(dhdl)
        except TypeError:
            formatted_data.append(
                [
                    dhdl,
                ]
            )
    estimators = formatted_data

    # Get the dF
    dF_list = []
    error_list = []
    max_length = 0
    for dhdl_list in estimators:
        len_dF = sum([len(dhdl.delta_f_) - 1 for dhdl in dhdl_list])
        if len_dF > max_length:
            max_length = len_dF
        dF = []
        error = []
        for dhdl in dhdl_list:
            for i in range(len(dhdl.delta_f_) - 1):
                dF.append(dhdl.delta_f_.iloc[i, i + 1])
                error.append(dhdl.d_delta_f_.iloc[i, i + 1])

        dF_list.append(dF)
        error_list.append(error)

    # Sort out the colors
    if colors is None:
        colors_dict = {
            "TI": "#C45AEC",
            "TI-CUBIC": "#33CC33",
            "DEXP": "#F87431",
            "IEXP": "#FF3030",
            "GINS": "#EAC117",
            "GDEL": "#347235",
            "BAR": "#6698FF",
            "UBAR": "#817339",
            "RBAR": "#C11B17",
            "MBAR": "#F9B7FF",
        }
        colors = []
        for dhdl in estimators:
            dhdl = dhdl[0]
            if isinstance(dhdl, TI):
                colors.append(colors_dict["TI"])
            elif isinstance(dhdl, BAR):
                colors.append(colors_dict["BAR"])
            elif isinstance(dhdl, MBAR):
                colors.append(colors_dict["MBAR"])
    else:
        if len(colors) >= len(estimators):
            pass
        else:
            raise ValueError(
                "Number of colors ({}) should be larger than the number of data ({})".format(
                    len(colors), len(estimators)
                )
            )

    # Sort out the labels
    if labels is None:
        labels = []
        for dhdl in estimators:
            dhdl = dhdl[0]
            if isinstance(dhdl, TI):
                labels.append("TI")
            elif isinstance(dhdl, BAR):
                labels.append("BAR")
            elif isinstance(dhdl, MBAR):
                labels.append("MBAR")
    else:
        if len(labels) == len(estimators):
            pass
        else:
            raise ValueError(
                "Length of labels ({}) should be the same as the number of data ({})".format(
                    len(labels), len(estimators)
                )
            )

    width = 1.0 / (len(estimators) + 1)
    elw = 30 * width
    ndx = 1
    for x in range(max_length):
        lines = []
        for i, (dF, error) in enumerate(zip(dF_list, error_list)):
            y = [dF[j] for j in x]
            ye = [error[j] for j in x]
            lw = 0.1 * elw
            line = ax.bar(
                x + len(lines) * width,
                y,
                width,
                color=colors[i],
                yerr=ye,
                lw=lw,
                error_kw=dict(elinewidth=elw, ecolor="black", capsize=0.5 * elw),
            )
            lines += (line[0],)
        for dir in ["left", "right", "top", "bottom"]:
            if dir == "left":
                ax.yaxis.set_ticks_position(dir)
            else:
                ax.spines[dir].set_color("none")

        plt.yticks(fontsize=8)
        ax.set_xlim(x[0] - width, x[-1] + len(lines) * width)
        plt.xticks(x + 0.5 * width * len(estimators), tuple(["%d--%d" % (i, i + 1) for i in x]), fontsize=8)
        ndx += 1
    x = np.arange(max_length)

    ax = plt.gca()

    for tick in ax.get_xticklines():
        tick.set_visible(False)
    leg = plt.legend(lines, labels, loc=3, ncol=2, prop=FontProperties(size=10), fancybox=True)
    plt.title("The free energy change breakdown", fontsize=12)
    plt.xlabel("States", fontsize=12, color="#151B54")
    plt.ylabel(r"$\Delta G$ ({})".format(units), fontsize=12, color="#151B54")
    return ax


def draw_correlation_figure(graph):
    def draw_ddg_correlation(ax):
        calc_ddg_list = []
        calc_ddg_error_list = []
        exp_ddg_list = []
        exp_ddg_error_list = []
        for edge in graph.edges_iter():
            if edge.get_data("ddg") is not None and edge.get_data("exp_ddg") is not None:
                calc_ddg_list.append(edge.get_data("ddg"))
                calc_ddg_error_list.append(edge.get_data("ddg_error"))
                exp_ddg_list.append(edge.get_data("exp_ddg"))
                exp_ddg_error_list.append(0)
        plot_correlation(
            calc_ddg_list,
            exp_ddg_list,
            calc_ddg_error_list,
            exp_ddg_error_list,
            system="ddG",
            ax=ax,
        )

    def draw_dg_correlation(ax):
        calc_dg = []
        exp_dg = []
        calc_dg_error = []
        exp_dg_error = []
        for node in graph.nodes_iter():
            if node.get_data("dg") is not None and node.get_data("cc_dg") is not None:
                calc_dg.append(node.get_data("cc_dg"))
                calc_dg_error.append(node.get_data("cc_dg_error"))
                exp_dg.append(node.get_data("dg"))
                exp_dg_error.append(0)
        plot_correlation(calc_dg, exp_dg, calc_dg_error, exp_dg_error, system="dg", ax=ax)

    def draw_corrected_ddg_correlation(ax):
        calc_ddg_list = []
        calc_ddg_error_list = []
        exp_ddg_list = []
        exp_ddg_error_list = []
        for edge in graph.edges_iter():
            if edge.get_data("cc_ddg") is not None and edge.get_data("exp_ddg") is not None:
                calc_ddg_list.append(edge.get_data("cc_ddg"))
                calc_ddg_error_list.append(edge.get_data("cc_ddg_error"))
                exp_ddg_list.append(edge.get_data("exp_ddg"))
                exp_ddg_error_list.append(0)
        plot_correlation(
            calc_ddg_list,
            exp_ddg_list,
            calc_ddg_error_list,
            exp_ddg_error_list,
            system="corrected ddG",
            ax=ax,
        )

    def draw_ligand_dG_bar(ax):
        dg_dict = {}
        for node in graph.nodes_iter():
            if node.get_data("cc_dg") is not None:
                dg_dict[node.name] = node.get_data("cc_dg")
        show_ligand_dG(dg_dict, ax=ax)

    fig, ax = plt.subplots(2, 2, figsize=(20, 18))
    draw_ddg_correlation(ax=ax[0][0])
    draw_dg_correlation(ax=ax[0][1])
    draw_corrected_ddg_correlation(ax=ax[1][0])
    draw_ligand_dG_bar(ax=ax[1][1])
    return fig, ax


def generate_rdkit_torsion_figure():
    pass
