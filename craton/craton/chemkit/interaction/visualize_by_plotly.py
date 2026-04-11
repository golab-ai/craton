from pathlib import Path
from typing import Union

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

COLORS = {
    "hbond": "#636EFA",
    "saltbridge": "#EF553B",
    "pistacking": "#00CC96",
    "pication": "#AB63FA",
    "hydrophobic": "#FFA15A",
    "halogen": "#19D3F3",
    "waterbridge": "#FF6692",
    "metal": "#B6E880",
    "weakhbond": "#FF97FF",
    "chpi": "#FECB52",
    "energy": "#636EFA",
    "Coul-SR": "#EF553B",
    "LJ-SR": "#00CC96",
}


def fep_bar_plot(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    pair_name,
    filename: Union[str, Path] = "interaction_result.html",
    save_figure=True,
):
    interaction_all = df1.columns

    fig = make_subplots(
        rows=1,
        cols=2,
        specs=[[{}, {}]],
        shared_xaxes=True,
        shared_yaxes=True,
        vertical_spacing=0,
        subplot_titles=[f"{pair_name}(A)", f"{pair_name}(B)"],
    )

    legend_group_dict = {interaction: f"group{i}" for i, interaction in enumerate(interaction_all)}
    for i, interaction in enumerate(interaction_all):
        
        fig.append_trace(
            go.Bar(
                name=interaction,
                y=df1.index,
                x=df1[interaction],
                orientation="h",
                marker=dict(color=COLORS[interaction]),
                legendgroup=legend_group_dict[interaction],
            ),
            1,
            1,
        )
    fig.update_xaxes(autorange="reversed", row=1, col=1)
    fig.update_yaxes(side="right", row=1, col=1, tickfont_size=22)

    for i, interaction in enumerate(interaction_all):
        fig.append_trace(
            go.Bar(
                name=interaction,
                y=df2.index,
                x=df2[interaction],
                orientation="h",
                showlegend=False,
                marker=dict(color=COLORS[interaction]),
                legendgroup=legend_group_dict[interaction],
            ),
            1,
            2,
        )

    fig.update_layout(barmode="stack")
    if save_figure:
        fig.write_html(filename)
    else:
        return fig


def fep_bar_plot_with_energy(
    interaction_df1, interaction_df2, energy_df1, energy_df2, pair_name, filename="interaction.html"
):
    fig = fep_bar_plot(interaction_df1, interaction_df2, pair_name, save_figure=False)
    interaction_data_left = [interaction_df1[column].to_list() for column in interaction_df1]
    interaction_data_right = [interaction_df2[column].to_list() for column in interaction_df1]
    bar_number = len(interaction_df1.columns)
    interaction_name = list(interaction_df1.columns)

    energy_data_left = [energy_df1["total_mean"], energy_df1["Coul-SR_mean"], energy_df1["LJ-SR_mean"]] + [None] * (
        bar_number - 3
    )
    energy_data_right = [energy_df2["total_mean"], energy_df2["Coul-SR_mean"], energy_df2["LJ-SR_mean"]] + [None] * (
        bar_number - 3
    )
    energy_name = ["energy", "Coul-SR", "LJ-SR"] + ([None] * (bar_number - 3))
    energy_visible = ([True, "legendonly", "legendonly"] + [False] * (bar_number - 3)) * 2

    energy_total_left_std = energy_df1["total_std"]
    energy_coul_left_std = energy_df1["Coul-SR_std"]
    energy_lj_left_std = energy_df1["Coul-SR_std"]
    error_energy_left = [
        dict(type="data", array=energy_total_left_std),
        dict(type="data", array=energy_coul_left_std),
        dict(type="data", array=energy_lj_left_std),
    ] + [None] * (bar_number - 3)

    energy_total_right_std = energy_df2["total_std"]
    energy_coul_right_std = energy_df2["Coul-SR_std"]
    energy_lj_right_std = energy_df2["Coul-SR_std"]
    error_energy_right = [
        dict(type="data", array=energy_total_right_std),
        dict(type="data", array=energy_coul_right_std),
        dict(type="data", array=energy_lj_right_std),
    ] + [None] * (bar_number - 3)

    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                buttons=list(
                    [
                        dict(
                            args=[
                                {
                                    "x": interaction_data_left + interaction_data_right,
                                    "name": interaction_name + interaction_name,
                                    "visible": [True] * bar_number * 2,
                                    "error_x": [dict(visible=False)],
                                },
                                list(range(bar_number * 2)),
                            ],
                            label="Interaction",
                            method="restyle",
                        ),
                        dict(
                            args=[
                                {
                                    "x": energy_data_left + energy_data_right,
                                    "name": energy_name,
                                    "visible": energy_visible,
                                    "error_x": error_energy_left + error_energy_right,
                                },
                                list(range(bar_number * 2)),
                            ],
                            label="Energy",
                            method="update",
                        ),
                    ]
                ),
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.11,
                xanchor="center",
                y=1.1,
                yanchor="top",
            ),
        ],
    )
    fig.write_html(filename)


def normal_md_bar_plot(df, output_file="result.html"):
    fig = make_subplots(rows=1, cols=1)
    legend_group_dict = {interaction: f"group{i}" for i, interaction in enumerate(df.columns)}
    for i, interaction in enumerate(df.columns):
        fig.append_trace(
            go.Bar(
                name=interaction,
                y=df.index,
                x=df[interaction],
                orientation="h",
                marker=dict(color=COLORS[interaction]),
                legendgroup=legend_group_dict[interaction],
            ),
            1,
            1,
        )
    fig.update_layout(barmode="stack")
    fig.write_html(output_file)


def normal_md_bar_plot_with_energy(interaction_df, energy_df, output_file="interaction_with_energy.html"):
    fig = make_subplots(rows=1, cols=1)
    legend_group_dict = {interaction: f"group{i}" for i, interaction in enumerate(interaction_df.columns)}
    interaction_data, interaction_name = [], []
    for i, interaction in enumerate(interaction_df.columns):
        fig.append_trace(
            go.Bar(
                name=interaction,
                y=interaction_df.index,
                x=interaction_df[interaction],
                orientation="h",
                marker=dict(color=COLORS[interaction]),
                legendgroup=legend_group_dict[interaction],
            ),
            1,
            1,
        )
        interaction_data.append(interaction_df[interaction])
        interaction_name.append(interaction)

    bar_number = len(interaction_df.columns)

    energy_name = ["energy", "Coul-SR", "LJ-SR"]
    energy_data = [energy_df["total_mean"], energy_df["Coul-SR_mean"], energy_df["LJ-SR_mean"]]
    energy_total_std = energy_df["total_std"]
    energy_coul_std = energy_df["Coul-SR_std"]
    energy_lj_std = energy_df["Coul-SR_std"]
    energy_visible = [True, "legendonly", "legendonly"] + [False] * (bar_number - 3)
    energy_error = [
        dict(type="data", array=energy_total_std),
        dict(type="data", array=energy_coul_std),
        dict(type="data", array=energy_lj_std),
    ]

    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                buttons=list(
                    [
                        dict(
                            args=[
                                {
                                    "x": interaction_data,
                                    "name": interaction_name,
                                    "visible": [True] * bar_number,
                                    "error_x": [dict(visible=False)],
                                },
                                list(range(len(interaction_df.columns))),
                            ],
                            label="Interaction",
                            method="restyle",
                        ),
                        dict(
                            args=[
                                {
                                    "x": energy_data,
                                    "name": energy_name,
                                    "visible": energy_visible,
                                    "error_x": energy_error,
                                },
                                list(range(len(interaction_df.columns))),
                            ],
                            label="Energy",
                            method="restyle",
                        ),
                    ]
                ),
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.11,
                xanchor="left",
                y=1.1,
                yanchor="top",
            ),
        ],
        barmode="stack",
    )
    fig.write_html(output_file)


def ligand_water_bar_plot(df, output_file="result.html"):
    fig = px.bar(df, x=["donor", "acceptor"], orientation="h")
    fig.update_layout(xaxis_title="Fraction", yaxis_title="Atom")
    fig.write_html(output_file)
