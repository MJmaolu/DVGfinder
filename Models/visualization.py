#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
## visualization.py contiene las funciones para generar los gráficos que 
## necesita el MÓDULO REPORTS 
##                  DVGfinder: Defective Viral Genome-finder
###############################################################################
## Genera los tipos de gráficos: scatterplot, histograma y diagrama de arcos
###############################################################################
## Author: Maria Jose Olmo-Uceda
## Version: 1.0
## Email: maolu@alumni.uv.es
## Date: 2021/07
###############################################################################

# third party imports
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go

# -----------------------------------------------------------------------------
# SCATTERPLOTS
# -----------------------------------------------------------------------------

# With matplotly.pyplot
# -----------------------------------------------------------------------------

def scatterplot_complete_color_program(df, sample_name):
    """
    Scatterplot de las coordenadas RI~BP de los eventos DVG del dataframe 
    introducido como argumento, coloreado por programa y 

    Args:
        df  (pd.DataFrame)  Dataframe con DVGs. Al menos debe tener las columnas:
                            'BP', 'RI', 'DVG_type' y 'read_counts_{program}'
        sample_name (str)                         
    """
    rgb_values = sns.color_palette("Set2", 3)

    plt.scatter(x='BP', y='RI', c=rgb_values[0], s='read_counts_virema', data=df, 
                alpha=0.6, label='virema')
    plt.scatter(x='BP', y='RI', c=rgb_values[2], s='read_counts_ditector', data=df, 
                alpha=0.6, label='ditector')
    plt.axis('square')
    plt.legend(loc='upper right', bbox_to_anchor=(0.65, -0.15), fancybox=True)
    plt.title("All events detected in sample '{}'".format(sample_name))
    plt.ylabel("RI position (nt)")
    plt.xlabel("BP position (nt)")
    


def scatterplot_complete_color_DVGtype(df, sample_name):
    """
    Scatterplot de todos los eventos coloreando por tipo de DVG. 
    Decidir si sized by abundance
    """
    # Get Unique continents
    color_labels = df['DVG_type'].unique()
    # List of colors in the color palettes
    rgb_values = sns.color_palette("Set2", len(color_labels))
    # separate in dfs by DVG_type categories
    dfr = df[df['DVG_type'] == 'Deletion_forward']
    drv = df[df['DVG_type'] == 'Deletion_reverse']
    inf = df[df['DVG_type'] == 'Insertion_forward']
    inr = df[df['DVG_type'] == 'Insertion_reverse']
    cb5 = df[df['DVG_type'] == '5cb/sb']
    cb3 = df[df['DVG_type'] == '3cb/sb']
    #plot
    plt.scatter(x='BP', y='RI', color=rgb_values[0], s=100*dfr['rpht_virema'], data=dfr, alpha=0.6, label='Deletion_forward')
    plt.scatter(x='BP', y='RI', color=rgb_values[1], s=100*drv['rpht_virema'], data=drv, alpha=0.6, label='Deletion_reverse')
    plt.scatter(x='BP', y='RI', color=rgb_values[2], s=100*inf['rpht_virema'], data=inf, alpha=0.6, label='Insertion_forward')
    plt.scatter(x='BP', y='RI', color=rgb_values[3], s=100*inr['rpht_virema'], data=inr, alpha=0.6, label='Insertion_reverse')
    plt.scatter(x='BP', y='RI', color=rgb_values[4], s=100*cb5['rpht_virema'], data=cb5, alpha=0.6, label='5cb/sb')
    plt.scatter(x='BP', y='RI', color=rgb_values[5], s=100*cb3['rpht_virema'], data=cb3, alpha=0.6, label='3cb/sb')

    plt.axis('square')
    plt.legend(bbox_to_anchor=(1, 1), fancybox=True)
    plt.title("Events detected in sample '{}'".format(sample_name))
    plt.ylabel("RI position (nt)")
    plt.xlabel("BP position (nt)")


# With plotly
# -----------------------------------------------------------------------------

def scatter_color_by_DVGtype(df, size_by, title=None, leader_coord1=69, 
                            leader_coord2=75):
    """
    Scatterplot de RI~BP coloreando por tipo de DVG y tamaño en función del 
    argumento 'sized_by' (lo habitual es que sea 'read_counts_{program}' o
    'mean_read_counts')
    Args:
        df  
        size_by (str)   Columna del df según la que queremos 
    """
    if title == None:
        title = "Relationship between BP and RI coordinates"
    
    scatter = px.scatter(df, 
            x="BP", y="RI",
	        size=size_by, 
            color="DVG_type",
            hover_name="cID_DI",
            title=title
            )
    """ TRS-L & TRS-B of SARS-CoV-2 plot indication
    scatter.add_vrect(x0=leader_coord1, x1=leader_coord2, 
        annotation_text="TRS_L",annotation_position="top",
        fillcolor="green", opacity=0.40
    )
    scatter.add_hrect(y0=leader_coord1, y1=leader_coord2, 
        annotation_text="TRS_L",
        fillcolor="green", opacity=0.40
    )
    """
    scatter.update_xaxes(showgrid=False)
    scatter.update_yaxes(showgrid=False)

    scatter.update_layout(width=700, height=700)

    return scatter

# -----------------------------------------------------------------------------
# HISTOGRAMS OF LENGTH DISTRIBUTION
# -----------------------------------------------------------------------------

def prepare_df_for_histogram_distribution(df, read_counts_column):
    """
    Prepara el df para que al plotear el histograma de frecuencias se tenga 
    en cuenta el número de conteos de los eventos
    """
    grouped = df.groupby(['DVG_type', 'length_dvg', read_counts_column]).\
        size().reset_index()
    grouped['counts'] = grouped[read_counts_column] * grouped[0]
    # Repetimos el número de filas del dataframe en función del valor de counts
    # así cuando la función histograma cuente filas de eventos que se repiten
    # el número de reads estará siendo tenido en cuenta
    df_frecs_to_histogram = grouped.loc[grouped.index.repeat(grouped['counts'])].\
        reset_index()[['DVG_type', 'length_dvg']]
    
    return df_frecs_to_histogram


def histogram_distribution_lengthDVG(df, read_counts_column, title=None):
    """
    Dibuja un histograma de frecuencias de las longitudes teóricas de los
    DVGs coloreados por tipos de DVGs. Para ello primero prepara la info
    del df para que el histograma tenga en cuenta el número de reads de cada
    tipo.

    Args:
        df   (pd.DataFrame)  Tabla con los eventos de cuya longitud se quiere 
                            representar la distribución 
        read_counts_column  (str)   Nombre de la columna que contiene los conteos
        title   (str)   Título que se le quiere dar al plot
    
    Returns:
        histogram
    """
    df_frecs_to_histogram = prepare_df_for_histogram_distribution(df, 
                                                        read_counts_column)

    if title == None:
        title = "Distribution of DVG lengths"
    histogram = px.histogram(df_frecs_to_histogram,
                    x="length_dvg",
                    color="DVG_type",
                    marginal="box",
                    nbins=100,
                    title=title
                        )

    return histogram

# -----------------------------------------------------------------------------
# BOXPLOTS
# -----------------------------------------------------------------------------

## NÚMBER OF EVENTS BY TYPE AND MODE 

# -----------------------------------------------------------------------------
# ARC DIAGRAMS ONLY FOR DELETIONS
# -----------------------------------------------------------------------------

def list_shapes(bp_list, ri_list, color, dvg_type, opacity=None):
    """
    Genera la lista de formas que se introduce al método update_layout

    Args:
        bp_list (list)  Lista con las coordenadas BP ordenadas 
        ri_list (list)  Lista con las coordenadas BP ordenadas 
        color   (str)   color de la linea
        opacity (float) opacidad de la linea. Por defecto será 1
    """
    ply_shapes = dict()

    if opacity == None:
        opacity = 1

    for i in range(len(bp_list)):
        ply_shapes['shape_' + str(i)]= go.layout.Shape(type="path",
                                    path="M {},0 Q {},{} {},0,".format(bp_list[i], 
                                    0.5*(bp_list[i]+ri_list[i]), 
                                    abs(bp_list[i]-ri_list[i])*0.8,
                                    ri_list[i]),
                                    line_color=color,
                                    #line_width=0.95,
                                    opacity=opacity,
                                    name=dvg_type
                )
    list_shapes = list(ply_shapes.values())

    return list_shapes


def arcs_diagram(df, len_wt, title=None):
    """
    Genera un diagrama de arcos de las 'Deletion_forward' y 'Deletion_reverse' 
    que existan en el df. Si no hay de alguna de las dos categorías también 
    funciona

    """
    if title == None:
        title = "Deletions detected"

    # añadimos 1 nt de margen
    xmax = len_wt + 100
    # separamos en df de deleciones forward y reverse
    df_for = df[df['DVG_type'] == 'Deletion_forward']
    df_rev = df[df['DVG_type'] == 'Deletion_reverse']

    # extraemos la lista de shapes de cada tipo de deleción
    list_shapes_for = list_shapes(list(df_for['BP']), list(df_for['RI']), \
        '#636EFA', 'Deletions_forward', 0.9)
    list_shapes_rev = list_shapes(list(df_rev['RI']), list(df_rev['BP']), \
        '#FFA15A', 'Deletions_reverse', 0.7)

    list_del_shapes = list_shapes_for + list_shapes_rev

    # Generamos la figura
    arcs = go.Figure()
    # axis
    arcs.update_xaxes(
    range=[0,xmax]
    )
    arcs.update_yaxes(
        range=[0, 15000],
        zeroline=False,
        showticklabels=False,
        showgrid=False
    )
    # arcs
    arcs.update_layout(shapes=list_del_shapes,
        hovermode="x unified",
        title=title,
        xaxis_title="Genome position (nts)",
        showlegend=True)

    return arcs