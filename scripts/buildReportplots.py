#!/usr/bin/env python3

#<IntronSeeker searches introns by splice-realigning reads on contigs.>
#Copyright (C) <2019-2024> INRAE
#<Sarah Maman, Philippe Bardou, Emilien Lasguignes, Faustine Oudin, Floréal Cabanettes, Christophe Klopp>
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.


import pandas as pd
import plotly as py
import numpy as np 
import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.subplots as psp
import re
import plotly.express as px 


# Histogram   https://plotly.com/python/v3/histograms/
def plot_hist_contigs_len(fastaContigsLen, mFastaContigsLen):
    """
    plot_hist_contigs_len function plot an histogram of contigs length from fasta and modified fasta contigs length data.
    """
    contigs = go.Histogram(
        x=fastaContigsLen,
        name='Contigs',
        opacity=0.75
    )
    modifiedContigs = go.Histogram(
        x=mFastaContigsLen,
        name='Modified contigs',
        opacity=0.75
    )
    data = [contigs, modifiedContigs]
    layout = go.Layout(
        xaxis=dict(
            title='Contigs length'
        ),
        yaxis=dict(
            title='Number of contigs'
        ),
        bargap=0.2,
        bargroupgap=0.1
    )
    fig = go.Figure(data=data, layout=layout)
    fig.update_layout(
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=20,
            t=30,
            pad=0
        )
    )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')

def plot_hist(values:dict, title:str, xtitle:str, ytitle:str):
    """
    plot_hist function is a reusable function to plot an histogram.
    """
    hist, bin_edges = np.histogram(values, bins=500)
    myhist = go.Bar(
        x=bin_edges[:-1],
        y=hist,
        name=title,
        opacity=0.85
    )
    data = [myhist]
    layout = go.Layout(
        xaxis=dict(
            title=xtitle
        ),
        yaxis=dict(
            title=ytitle
        )
    )
    fig = go.Figure(data=data, layout=layout)
    fig.update_layout(
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=20,
            t=30,
            pad=0
        )
    )   
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')


def plot_dist_features_len(len_by_features, feature_names):
    """
    plot_dist_features_len function plot length by features as a distribution.
    """
    hist_data = len_by_features
    group_labels = feature_names
    colors = ['#333F44', '#37AA9C', '#94F3E4']
    # Create distplot with curve_type set to 'normal'
    fig = ff.create_distplot(hist_data, group_labels, show_hist=False, show_rug=True, colors=colors)
    fig.update_layout(
        xaxis=dict(
            title='Features length'
        ),
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=20,
            t=30,
            pad=0
        )
    )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')


def plot_insertion_in_contig(positions) :
    """
    plot_insertion_in_contig function plot introns position on a contig
    """
    hist = go.Histogram(
            x=positions,
            xbins=dict(
                start=0,
                end=100,
                size=2)
    )
    layout = go.Layout(xaxis=dict(title="% of contig length"),
                       yaxis=dict(title="Number of introns"))
    fig = go.Figure(data=[hist],layout=layout)
    fig.update_layout(
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=20,
            t=30,
            pad=0
        )
    )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')

def plot_abondance_model(df_fasta:dict) :
    """
    plot_abondance_model function plot abundance from df_fasta dataframe
    """
    fig = go.Figure()
    if 'waiting' in df_fasta.columns:
        fig.add_trace(
            go.Scatter(
                x = df_fasta.index,
                y = df_fasta['waiting'],
                mode = 'lines',
                name ='Waited abundance'
            )
        )
    fig.add_trace(
        go.Scatter(
            x = df_fasta.index,
            y = df_fasta['real'],
            mode = 'lines',
            name = 'Abundance'
        )
    )
    fig.add_trace(
        go.Scatter(
            x = df_fasta.index,                              
            y = df_fasta['norm'],
            mode = 'lines',
            name = 'Normalized abundance'
        )
    )
    fig.update_layout(
        xaxis=dict(title="Contig names"),
        yaxis=dict(title="Abundance percentage"#,
        # range=[-0.25,0.5]
        ),
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=20,
            t=30,
            pad=0
        )
    )
    return py.offline.plot(fig, include_plotlyjs=False, output_type='div')
