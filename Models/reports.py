#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
## Functions for REPORTS module 
##                  DVGfinder: Defective Viral Genome-finder
###############################################################################
## Generate an HTML report. Also, the tables shown are saved in the three modes:
## 'All', 'Consensus' y 'ML' in the 'FinalReports/{sample_name}/' directory.
###############################################################################
## Author: Maria Jose Olmo-Uceda
## Version: 3.0
## Email: maolu@alumni.uv.es
## Date: 2021/12
###############################################################################

# third party imports
import pandas as pd
import numpy as np
import plotly.express as px
import datapane as dp

# DVGfinder functions
from Models import visualization

# -----------------------------------------------------------------------------
# CREATION OF FINAL TABLES
# -----------------------------------------------------------------------------
def create_table_all(df, case):
    """
    Selecciona la información que se mostrará en la tabla completa de eventos 
    detectados
    Args:
        df  (pd.DataFrame)  Tabla accesible con toda la información de todos 
                            los eventos detectados por el módulo metabuscador.
    
    Return:
        df_all_show (pd.DataFrame)  Tabla con la selección de variables a 
                                    mostrar como informe
    """
    # Variables a mostrar según el caso de uso
    # caso 1: 
    selected_features = ['cID_DI', 'BP', 'RI', 'sense', 'DVG_type', 'length_dvg',
        'read_counts_virema', 'pBP_virema', 'pRI_virema' ,'rpht_virema',
        'read_counts_ditector','pBP_ditector', 'pRI_ditector', 'rpht_ditector']

    # caso 2:
    vir_selected = ['cID_DI', 'BP', 'RI', 'sense', 'DVG_type', 'length_dvg',
    'read_counts_virema', 'pBP_virema', 'pRI_virema' ,'rpht_virema']

    # caso 3:
    dit_selected = ['cID_DI', 'BP', 'RI', 'sense', 'DVG_type', 'length_dvg',
    'read_counts_ditector','pBP_ditector', 'pRI_ditector', 'rpht_ditector']
    
    if case == "case_1":
        df_all_show = df[selected_features]
    elif case == "case_2":
        df_all_show = df[vir_selected]
    elif case == "case_3":
        df_all_show = df[dit_selected]

    return df_all_show


def create_table_consensus(df_all):
    """
    Selecciona los eventos detectados por ambos programas. Y muestra las 
    variables seleccionadas para mostrar (selected_features) más la media
    y las desviaciones típicas de las variables 'read_counts_{program}' y
    'rpht_{program}'

    Args:
        df_all  (pd.DataFrame)  Tabla con todos los eventos detectados, 
                                Resultado del módulo metabuscador.
    
    Return:
        df_consensus    (pd.DataFrame)  Tabla con los DVGs consenso y las 
                                variables seleccionadas para mostrar en el 
                                informe HTML
    """
    # filtrar los eventos detectados por ambos programas
    df_consensus = df_all.loc[((df_all['read_counts_virema']>0) & 
                                (df_all['read_counts_ditector']>0))]
    # In case there are consensus DVGs, add features
    if not df_consensus.empty:
        df_consensus['mean_read_counts'] = np.mean(df_consensus[['read_counts_virema', 
                                            'read_counts_ditector']], axis=1)
        df_consensus['std_read_counts'] = np.std(df_consensus[['read_counts_virema', 
                                            'read_counts_ditector']], axis=1)
        df_consensus['mean_rpht'] = np.mean(df_consensus[['rpht_virema',
                                                        'rpht_ditector']], axis=1)
        df_consensus['std_rpht'] = np.std(df_consensus[['rpht_virema', 
                                            'rpht_ditector']], axis=1)                                               
                                                    
        df_consensus = df_consensus[['cID_DI', 'BP', 'RI', 
                'sense', 'DVG_type', 'length_dvg',
                'mean_read_counts', 'std_read_counts', 'mean_rpht', 'std_rpht',
                'read_counts_virema', 'pBP_virema', 'pRI_virema' ,'rpht_virema',
            'read_counts_ditector','pBP_ditector', 'pRI_ditector', 'rpht_ditector']]
    
    return df_consensus

def create_table_filteredML(df_filtered):
    """
    ################################ NOT USED ################################
    Selecciona las variables a mostrar en la tabla de filtrados ML 
    Args:
        df_filtered (pd.DataFrame)  Tabla con los eventos que quedan tras
                                filtrar con el modelo de ML
    
    Return:
        df_filtered_show    (pd.DataFrame)  Tabla que se mostrará
    """
    selected_features = ['cID_DI', 'BP', 'RI', 
        'sense', 'DVG_type', 'length_dvg',
        'read_counts_virema', 'pBP_virema', 'pRI_virema' ,'rpht_virema',
        'read_counts_ditector','pBP_ditector', 'pRI_ditector', 'rpht_ditector']
    df_filtered_show = df_filtered[selected_features]

    return df_filtered_show


# -----------------------------------------------------------------------------
# OTHER UTIL FUNCTIONS
# -----------------------------------------------------------------------------

def separate_df_by_program(df):
    """
    Separa el df introducido en dos, cada uno con los valores específicos 
    detectados por programa

    Args:
        df  (pd.DataFrame)
    Returns:
        df_virema   (pd.DataFrame)  Eventos detectados por virema
        df_ditector (pd.DataFrame)  Eventos detectados por ditector
    """
    # Establecemos las variables finales que queremos mostrar en cada df
    common_var = ['cID_DI', 'BP', 'RI', 'sense', 'DVG_type', 'length_dvg']
    virema_var = common_var + ['read_counts_virema', 'pBP_virema', 
                                'pRI_virema' ,'rpht_virema']
    ditector_var = common_var + ['read_counts_ditector','pBP_ditector', 
                                'pRI_ditector', 'rpht_ditector']
    
    df_virema = df[df['read_counts_virema'] > 0]
    df_ditector = df[df['read_counts_ditector'] > 0]

    df_virema = df_virema[df_virema.columns.intersection(virema_var)]
    df_ditector = df_ditector[df_ditector.columns.intersection(ditector_var)]

    return df_virema, df_ditector


def separate_df_all_by_program(df_all_show):
    """
    Separa el df con todos los eventos en dos nuevos dfs con la información
    específica de cada programa. 
    Args:
        df_all_show (pd.DataFrame)  Tabla accesible con todos los eventos. Ha
                                    sido previamente reducida a las columnas
                                    que se quieren mostrar. Output de la 
                                    función 'create_table_all'
    
    Returns:
        df_all_virema   (pd.DataFrame)  Todos los eventos detectados por ViReMA 
        df_all_ditector (pd.DataFrame)  Todos los eventos detectados por DItector 
    """

    df_all_virema, df_all_ditector = separate_df_by_program(df_all_show)
    
    return df_all_virema, df_all_ditector


def separate_filteredML_by_program(df_filteredML):
    """
    Separa el df con los eventos filtrados tras la aplicación del modelo ML en 
    dos nuevos dfs con la información específica de cada programa. 
    Args:
        df_filteredML (pd.DataFrame)  Tabla accesible con los eventos que quedan
                                    tras el filtro con el modelo de ML.
    
    Returns:
        df_filteredML_virema   (pd.DataFrame)  Eventos filtrados detectados por
                                            ViReMA 
        df_filteredML_ditector (pd.DataFrame)   Eventos filtrados detectados 
                                            por DItector 
    """
    df_filtered_virema, df_filteredML_ditector = separate_df_by_program(df_filteredML)
    
    return df_filtered_virema, df_filteredML_ditector


def separate_df_in_DFs_by_DVGtype(df):
    """
    Separa el df en otros DFs según el 'DVG_type'. Algunos pueden estar vacíos
    si esa categoría no tiene elementos.

    Args:
        df  (pd.DataFrame)  
    Returns:
        dict_df (dict)  Diccionario con los dfs separados:
                        keys:   'DF_3cbsb'          
                                'DF_5cbsb'
                                'DF_Deletion_forward'
                                'DF_Deletion_reverse'
                                'DF_Insertion_forward'
                                'DF_Insertion_reverse'
                        values: los pd.df correspondientes
    """
    # paso necesario hasta que cambiemos nomenclatura de cb quitando el '/'
    dvg_names = ['3cbsb', '5cbsb', 'Deletion_forward', 'Deletion_reverse', 
    'Insertion_forward', 'Insertion_reverse']
    dvg_types = ['3cb/sb', '5cb/sb', 'Deletion_forward', 'Deletion_reverse', 
    'Insertion_forward', 'Insertion_reverse']

    dict_dfs = dict()

    for dvgType in dvg_names:
        dict_dfs['DF_{}'.format(dvgType)] = df[df['DVG_type']== 
                                        dvg_types[dvg_names.index(dvgType)]]
        
    return dict_dfs


# -----------------------------------------------------------------------------
# CREATION OF FINAL REPORTS
# -----------------------------------------------------------------------------

def generate_report(df, df_predicted_reals_show, sample_name, len_wt):
    """
    Función completa para generar el informe HTML. LLama a las funciones 
    necesarias para crear las tablas de los 3 modos y todas las visualizaciones.
    Monta el resultado en un fichero HTML (utiliza la librería datapane)
    que se abre automáticamente en el navegador establecido por defecto. Guarda
    el HTML y las tablas de los tres modos (.csv) en el directorio 
    'FinalReports/{sample_name}/'

    Args:
        df          (pd.DataFrame)  Tabla completa con todos los eventos 
                                detectados con el módulo metabuscador
        df_predicted_reals_show (pd.DataFrame) Tabla con los eventos predichos
                                como reales por el modelo.
        sample_name     (str)   Nombre de la muestra.
        len_wt          (int)   Longitud del genoma de referencia
    
    """
    # file to save
    paths_csv = 'FinalReports' + "/" + sample_name + "/" + \
                sample_name
    path_html = 'FinalReports' + "/" + sample_name + "/" + \
                sample_name + "_report.html"
    # Generating tables to show up
    # -------------------------------------------------------------------------
    ## All the DVGs reducing the number of features to show
    df_all_show = create_table_all(df, "case_1")
    # Separating by search program
    df_all_virema, df_all_ditector = separate_df_all_by_program(df_all_show)
    ## Consensus DVGs (detected by both programs, sense sensitive)
    df_consensus = create_table_consensus(df_all_show)
    ## Predicted as reals by the model with a probability > t
    ## If there are events predicted as reals, separate by program    
    if not df_predicted_reals_show.empty:
        predicted = True
        df_filtered_virema, df_filtered_ditector = separate_df_all_by_program\
                                                (df_predicted_reals_show)
        # check which program has filtered
        if not df_filtered_virema.empty:
            virema_predicted = True
        else:
            virema_predicted = False
        if not df_filtered_ditector.empty:
            ditector_predicted = True
        else:
            ditector_predicted = False

    # Write tables
    # -------------------------------------------------------------------------
    df_all_show.to_csv(paths_csv + "_ALL.csv", index=False)
    df_consensus.to_csv(paths_csv + "_CONSENSUS.csv", index=False)
    df_predicted_reals_show.to_csv(paths_csv + "_FILTERED_ML.csv", index=False)

    # Generate visualizations
    # -------------------------------------------------------------------------
    ## COMPLETE
    ### scatterplots
    all_scatter_vir = visualization.scatter_color_by_DVGtype(df_all_virema, 
                "read_counts_virema", 
                title="Relationship between BP and RI coordinates in all " +  
                "the DVGs detected by ViReMa-a",    
                #leader_coord1=69, leader_coord2=75
                )
    all_scatter_dit = visualization.scatter_color_by_DVGtype(df_all_ditector, 
                "read_counts_ditector", 
                title="Relationship between BP and RI coordinates in all " +  
                "the DVGs detected by DI-tector",   
                #leader_coord1=69, leader_coord2=75
                )
    ### histograms
    all_hist_vir = visualization.histogram_distribution_lengthDVG(df_all_virema,
                 "read_counts_virema", 
                 title="Distribution of DVG lengths detected by ViReMa-a")
    all_hist_dit = visualization.histogram_distribution_lengthDVG(df_all_ditector,
                 "read_counts_ditector",
                 title="Distribution of DVG lengths detected by DI-tector")
    ### arc diagrams
    all_arcs_vir = visualization.arcs_diagram(df_all_virema, len_wt, 
                title="Deletions detected by ViReMa-a")
    all_arcs_dit = visualization.arcs_diagram(df_all_ditector, len_wt,
                title="Deletions detected by DI-tector")

    ## CONSENSUS. Si hay eventos consenso generamos las visualizaciones
    ### scatterplots
    if not df_consensus.empty:
        consensus = True
        consensus_scatter = visualization.scatter_color_by_DVGtype(df_consensus, 
                    "mean_read_counts", 
                    title="""
                    Relationship between BP and RI coordinates in consensus DVGs
                    """,    
                    leader_coord1=69, leader_coord2=75)
        ### histograms
        consensus_hist = visualization.histogram_distribution_lengthDVG(df_consensus,
                    "mean_read_counts", 
                    title="Distribution of DVG lengths in consensus DVGs")
        ### arc diagrams
        consensus_arcs = visualization.arcs_diagram(df_consensus, len_wt, 
                    title="Consensus deletions detected",
                )
    else:
        consensus = False

    ## ML MODEL FILTER. Only if there are filtered DVGs as potential reals
    ### scatterplots
    if predicted:
        # virema has filtered events
        if virema_predicted:
            # scatter plots
            ML_scatter_vir = visualization.scatter_color_by_DVGtype(df_filtered_virema, 
                        "read_counts_virema", 
                        title="Relationship between BP and RI coordinates in " +  
                        "DVGs detected by ViReMa-a filtered by ML",    
                        leader_coord1=69, leader_coord2=75)
            # histogram
            ML_hist_vir = visualization.histogram_distribution_lengthDVG(df_filtered_virema,
                    "read_counts_virema", 
                    title="Distribution of DVG lengths detected by ViReMa-a")
            # arc diagram
            ML_arcs_vir = visualization.arcs_diagram(df_filtered_virema, len_wt, 
                    title="Deletions detected by ViReMa-a")
        if ditector_predicted:
            # scatter plots
            ML_scatter_dit = visualization.scatter_color_by_DVGtype(df_filtered_ditector, 
                    "read_counts_ditector", 
                    title="Relationship between BP and RI coordinates in " +  
                    "the DVGs detected by DI-tector filtered by ML",   
                    leader_coord1=69, leader_coord2=75)
            # histograms
            ML_hist_dit = visualization.histogram_distribution_lengthDVG(df_filtered_ditector,
                    "read_counts_ditector",
                    title="Distribution of DVG lengths detected by DI-tector")
            # arc diagrams
            ML_arcs_dit = visualization.arcs_diagram(df_filtered_ditector, len_wt,
                    title="Deletions detected by DI-tector")

    # Montamos el report en función del tipo de caso 1 en el que estemos
    # -------------------------------------------------------------------------
    if consensus and predicted:
        if virema_predicted and ditector_predicted:
            # Informe completo: modos ALL, CONSENSUS & MLfiltered
            r = dp.Report(
            dp.Page(
                # COMPLETE
                #----------------------------------
                title="All DVGs detected",        
                blocks=[
                    dp.HTML(header_html(sample_name)),
                    dp.Text("## Complete table of DVG detected by ViReMa-a _OR_ DI-tector"),
                    dp.DataTable(df_all_show),
                    dp.Text("---"),
                    dp.Text("### Complete DVGs plots"),
                    ## scatterplots
                    dp.Select(blocks=[
                        dp.Plot(all_scatter_vir, label= "ViReMa-a detection",
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the  number of DVG readings " + 
                        "detected by ViReMa-a."),
                        dp.Plot(all_scatter_dit, label="DI-tector detection",
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the number of DVG readings " + 
                        "detected by DI-tector.")
                    ]),
                    ## histograms
                    dp.Select(blocks=[
                        dp.Plot(all_hist_vir, label= "ViReMa-a detection",
                        caption="Frequency histogram of the theoretical lengths of the DVGs."),
                        dp.Plot(all_hist_dit, label="DI-tector detection",
                        caption="Frequency histogram of the theoretical lengths of the DVGs.")
                    ]),
                    ## Arc diagrams
                    dp.Select(blocks=[
                        dp.Plot(all_arcs_vir, label= "ViReMa-a detection",
                        caption="Arc diagram of the deletions detected. The range " + 
                        " between the two 'x' coordinates represent the deletioned " +
                        " genome part."),
                        dp.Plot(all_arcs_dit, label="DI-tector detection",
                        caption="Arc diagram of the deletions detected. The range " + 
                        " between the two 'x' coordinates represent the deletioned " +
                        " genome part.")
                    ]),
                    dp.HTML(arcs_caption()),
                    dp.HTML(footer_html())
                    ]            
                ),

            dp.Page(
                # CONSENSUS
                #----------------------------------
                title="DVG consensus selected",
                blocks=[
                    dp.HTML(header_html(sample_name)),
                    dp.Text("## Consensus table of DVG detected by ViReMa-a _AND_ DI-tector"),
                    dp.DataTable(df_consensus),
                    dp.Text("---"),
                    dp.Text("### Consensus DVG plots"),
                    dp.Plot(consensus_scatter, 
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the **mean** number of DVG readings " + 
                        "detected by the programs."),
                    dp.Plot(consensus_hist,
                    caption="Frequency histogram of the theoretical lengths of the DVGs."),
                    dp.Plot(consensus_arcs,
                    caption="Arc diagram of the deletions detected. The range " + 
                        " between the two x coordinates represent the deletioned " +
                        " genome part."),
                    dp.HTML(arcs_caption()),
                    dp.HTML(footer_html())
                    ]
                ),
            dp.Page(
                # FILTERED ML
                #----------------------------------
                title="DVG ML filtered",
                blocks=[
                    dp.HTML(header_html(sample_name)),
                    dp.Text("## Table of DVG predicted as reals with the DVGfinder-MLmodel"),
                    dp.DataTable(df_predicted_reals_show),
                    dp.Text("---"),
                    dp.Text("### ML filtered DVG plots"),
                    ## scatterplots
                    dp.Select(blocks=[
                        dp.Plot(ML_scatter_vir, label= "ViReMa-a detection",
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the  number of DVG readings " + 
                        "detected by ViReMa-a."),
                        dp.Plot(ML_scatter_dit, label="DI-tector detection",
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the number of DVG readings " + 
                        "detected by DI-tector."),
                    ]),
                    ## histograms
                    dp.Select(blocks=[
                        dp.Plot(ML_hist_vir, label= "ViReMa-a detection",
                        caption="Frequency histogram of the theoretical lengths of the DVGs."),
                        dp.Plot(ML_hist_dit, label="DI-tector detection",
                        caption="Frequency histogram of the theoretical lengths of the DVGs.")
                    ]),
                    ## Arc diagrams
                    dp.Select(blocks=[
                        dp.Plot(ML_arcs_vir, label= "ViReMa-a detection",
                        caption="Arc diagram of the deletions detected. The range " + 
                        " between the two x coordinates represent the deletioned " +
                        " genome part."),
                        dp.Plot(ML_arcs_dit, label="DI-tector detection",
                        caption="Arc diagram of the deletions detected. The range " + 
                        " between the two x coordinates represent the deletioned " +
                        " genome part.")]
                    ),
                    dp.HTML(arcs_caption()),
                    dp.HTML(footer_html())
                ])
            )
        elif virema_predicted and not ditector_predicted:
            # Informe completo: modos ALL, CONSENSUS & MLfiltered only virema
            r = dp.Report(
            dp.Page(
                # COMPLETE
                #----------------------------------
                title="All DVGs detected",        
                blocks=[
                    dp.HTML(header_html(sample_name)),
                    dp.Text("## Complete table of DVG detected by ViReMa-a _OR_ DI-tector"),
                    dp.DataTable(df_all_show),
                    dp.Text("---"),
                    dp.Text("### Complete DVGs plots"),
                    ## scatterplots
                    dp.Select(blocks=[
                        dp.Plot(all_scatter_vir, label= "ViReMa-a detection",
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the  number of DVG readings " + 
                        "detected by ViReMa-a."),
                        dp.Plot(all_scatter_dit, label="DI-tector detection",
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the number of DVG readings " + 
                        "detected by DI-tector.")
                    ]),
                    ## histograms
                    dp.Select(blocks=[
                        dp.Plot(all_hist_vir, label= "ViReMa-a detection",
                        caption="Frequency histogram of the theoretical lengths of the DVGs."),
                        dp.Plot(all_hist_dit, label="DI-tector detection",
                        caption="Frequency histogram of the theoretical lengths of the DVGs.")
                    ]),
                    ## Arc diagrams
                    dp.Select(blocks=[
                        dp.Plot(all_arcs_vir, label= "ViReMa-a detection",
                        caption="Arc diagram of the deletions detected. The range " + 
                        " between the two 'x' coordinates represent the deletioned " +
                        " genome part."),
                        dp.Plot(all_arcs_dit, label="DI-tector detection",
                        caption="Arc diagram of the deletions detected. The range " + 
                        " between the two 'x' coordinates represent the deletioned " +
                        " genome part.")
                    ]),
                    dp.HTML(arcs_caption()),
                    dp.HTML(footer_html())
                    ]            
                ),

            dp.Page(
                # CONSENSUS
                #----------------------------------
                title="DVG consensus selected",
                blocks=[
                    dp.HTML(header_html(sample_name)),
                    dp.Text("## Consensus table of DVG detected by ViReMa-a _AND_ DI-tector"),
                    dp.DataTable(df_consensus),
                    dp.Text("---"),
                    dp.Text("### Consensus DVG plots"),
                    dp.Plot(consensus_scatter, 
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the **mean** number of DVG readings " + 
                        "detected by the programs."),
                    dp.Plot(consensus_hist,
                    caption="Frequency histogram of the theoretical lengths of the DVGs."),
                    dp.Plot(consensus_arcs,
                    caption="Arc diagram of the deletions detected. The range " + 
                        " between the two x coordinates represent the deletioned " +
                        " genome part."),
                    dp.HTML(arcs_caption()),
                    dp.HTML(footer_html())
                    ]
                ),
            dp.Page(
                # FILTERED ML
                #----------------------------------
                title="DVG ML filtered",
                blocks=[
                    dp.HTML(header_html(sample_name)),
                    dp.Text("## Table of DVG predicted as reals with the DVGfinder-MLmodel"),
                    dp.DataTable(df_predicted_reals_show),
                    dp.Text("---"),
                    dp.Text("### ML filtered DVG plots"),
                    ## scatterplot
                    dp.Plot(ML_scatter_vir, label= "ViReMa-a detection",
                    caption="Scatterplot of RI~BP coordinates. " + 
                    "Bubble size is relative to the  number of DVG readings " + 
                    "detected by ViReMa-a."),
                    ## histograms
                    dp.Plot(ML_hist_vir, label= "ViReMa-a detection",
                    caption="Frequency histogram of the theoretical lengths of the DVGs."),
                    ## Arc diagrams
                    dp.Plot(ML_arcs_vir, label= "ViReMa-a detection",
                    caption="Arc diagram of the deletions detected. The range " + 
                    " between the two x coordinates represent the deletioned " +
                    " genome part."),
                    dp.HTML(arcs_caption()),
                    dp.HTML(footer_html())
                ])
            )
        elif not virema_predicted and ditector_predicted:
            # Informe completo: modos ALL, CONSENSUS & MLfiltered only ditector
            r = dp.Report(
            dp.Page(
                # COMPLETE
                #----------------------------------
                title="All DVGs detected",        
                blocks=[
                    dp.HTML(header_html(sample_name)),
                    dp.Text("## Complete table of DVG detected by ViReMa-a _OR_ DI-tector"),
                    dp.DataTable(df_all_show),
                    dp.Text("---"),
                    dp.Text("### Complete DVGs plots"),
                    ## scatterplots
                    dp.Select(blocks=[
                        dp.Plot(all_scatter_vir, label= "ViReMa-a detection",
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the  number of DVG readings " + 
                        "detected by ViReMa-a."),
                        dp.Plot(all_scatter_dit, label="DI-tector detection",
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the number of DVG readings " + 
                        "detected by DI-tector.")
                    ]),
                    ## histograms
                    dp.Select(blocks=[
                        dp.Plot(all_hist_vir, label= "ViReMa-a detection",
                        caption="Frequency histogram of the theoretical lengths of the DVGs."),
                        dp.Plot(all_hist_dit, label="DI-tector detection",
                        caption="Frequency histogram of the theoretical lengths of the DVGs.")
                    ]),
                    ## Arc diagrams
                    dp.Select(blocks=[
                        dp.Plot(all_arcs_vir, label= "ViReMa-a detection",
                        caption="Arc diagram of the deletions detected. The range " + 
                        " between the two 'x' coordinates represent the deletioned " +
                        " genome part."),
                        dp.Plot(all_arcs_dit, label="DI-tector detection",
                        caption="Arc diagram of the deletions detected. The range " + 
                        " between the two 'x' coordinates represent the deletioned " +
                        " genome part.")
                    ]),
                    dp.HTML(arcs_caption()),
                    dp.HTML(footer_html())
                    ]            
                ),

            dp.Page(
                # CONSENSUS
                #----------------------------------
                title="DVG consensus selected",
                blocks=[
                    dp.HTML(header_html(sample_name)),
                    dp.Text("## Consensus table of DVG detected by ViReMa-a _AND_ DI-tector"),
                    dp.DataTable(df_consensus),
                    dp.Text("---"),
                    dp.Text("### Consensus DVG plots"),
                    dp.Plot(consensus_scatter, 
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the **mean** number of DVG readings " + 
                        "detected by the programs."),
                    dp.Plot(consensus_hist,
                    caption="Frequency histogram of the theoretical lengths of the DVGs."),
                    dp.Plot(consensus_arcs,
                    caption="Arc diagram of the deletions detected. The range " + 
                        " between the two x coordinates represent the deletioned " +
                        " genome part."),
                    dp.HTML(arcs_caption()),
                    dp.HTML(footer_html())
                    ]
                ),
            dp.Page(
                # FILTERED ML
                #----------------------------------
                title="DVG ML filtered",
                blocks=[
                    dp.HTML(header_html(sample_name)),
                    dp.Text("## Table of DVG predicted as reals with the DVGfinder-MLmodel"),
                    dp.DataTable(df_predicted_reals_show),
                    dp.Text("---"),
                    dp.Text("### ML filtered DVG plots"),
                    ## scatterplot
                    dp.Plot(ML_scatter_dit, label= "DI-tector detection",
                    caption="Scatterplot of RI~BP coordinates. " + 
                    "Bubble size is relative to the  number of DVG readings " + 
                    "detected by DI-tector."),
                    ## histograms
                    dp.Plot(ML_hist_dit, label= "DI-tector detection",
                    caption="Frequency histogram of the theoretical lengths of the DVGs."),
                    ## Arc diagrams
                    dp.Plot(ML_arcs_dit, label= "DI-tector detection",
                    caption="Arc diagram of the deletions detected. The range " + 
                    " between the two x coordinates represent the deletioned " +
                    " genome part."),
                    dp.HTML(arcs_caption()),
                    dp.HTML(footer_html())
                ])
            )
    
    elif consensus and not predicted:
        r = report_AllConsensus(sample_name, 
                        df_all_show, all_scatter_vir, all_scatter_dit,
                        all_hist_vir, all_hist_dit, all_arcs_vir, all_arcs_dit,
                        df_consensus, consensus_scatter, consensus_hist,
                        consensus_arcs)
    ## With no Consensus
    elif not consensus and predicted:
        if virema_predicted and ditector_predicted:
            r = report_AllML(sample_name, 
                df_all_show, all_scatter_vir, all_scatter_dit,
                all_hist_vir, all_hist_dit, all_arcs_vir, all_arcs_dit,
                df_predicted_reals_show, ML_scatter_vir, ML_scatter_dit, 
                ML_hist_vir, ML_hist_dit, ML_arcs_vir, ML_arcs_dit)
        elif virema_predicted and not ditector_predicted:
            r = dp.Report(
            dp.Page(
                # COMPLETE
                #----------------------------------
                title="All DVGs detected",        
                blocks=[
                    dp.HTML(header_html(sample_name)),
                    dp.Text("## Complete table of DVG detected by ViReMa-a _OR_ DI-tector"),
                    dp.DataTable(df_all_show),
                    dp.Text("---"),
                    dp.Text("### Complete DVGs plots"),
                    ## scatterplots
                    dp.Select(blocks=[
                        dp.Plot(all_scatter_vir, label= "ViReMa-a detection",
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the  number of DVG readings " + 
                        "detected by ViReMa-a."),
                        dp.Plot(all_scatter_dit, label="DI-tector detection",
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the number of DVG readings " + 
                        "detected by DI-tector.")
                    ]),
                    ## histograms
                    dp.Select(blocks=[
                        dp.Plot(all_hist_vir, label= "ViReMa-a detection",
                        caption="Frequency histogram of the theoretical lengths of the DVGs."),
                        dp.Plot(all_hist_dit, label="DI-tector detection",
                        caption="Frequency histogram of the theoretical lengths of the DVGs.")
                    ]),
                    ## Arc diagrams
                    dp.Select(blocks=[
                        dp.Plot(all_arcs_vir, label= "ViReMa-a detection",
                        caption="Arc diagram of the deletions detected. The range " + 
                        " between the two 'x' coordinates represent the deletioned " +
                        " genome part."),
                        dp.Plot(all_arcs_dit, label="DI-tector detection",
                        caption="Arc diagram of the deletions detected. The range " + 
                        " between the two 'x' coordinates represent the deletioned " +
                        " genome part.")
                    ]),
                    dp.HTML(arcs_caption()),
                    dp.HTML(footer_html())
                    ]            
                ),
            dp.Page(
                # FILTERED ML ViReMa-a
                #----------------------------------
                title="DVG ML filtered",
                blocks=[
                    dp.HTML(header_html(sample_name)),
                    dp.Text("## Table of DVG predicted as reals with the DVGfinder-MLmodel"),
                    dp.DataTable(df_predicted_reals_show),
                    dp.Text("---"),
                    dp.Text("### ML filtered DVG plots"),
                    ## scatterplot
                    dp.Plot(ML_scatter_vir, label= "ViReMa-a detection",
                    caption="Scatterplot of RI~BP coordinates. " + 
                    "Bubble size is relative to the  number of DVG readings " + 
                    "detected by ViReMa-a."),
                    ## histograms
                    dp.Plot(ML_hist_vir, label= "ViReMa-a detection",
                    caption="Frequency histogram of the theoretical lengths of the DVGs."),
                    ## Arc diagrams
                    dp.Plot(ML_arcs_vir, label= "ViReMa-a detection",
                    caption="Arc diagram of the deletions detected. The range " + 
                    " between the two x coordinates represent the deletioned " +
                    " genome part."),
                    dp.HTML(arcs_caption()),
                    dp.HTML(footer_html())
                ])
            )
        elif not virema_predicted and ditector_predicted:
            r = dp.Report(
            dp.Page(
                # COMPLETE
                #----------------------------------
                title="All DVGs detected",        
                blocks=[
                    dp.HTML(header_html(sample_name)),
                    dp.Text("## Complete table of DVG detected by ViReMa-a _OR_ DI-tector"),
                    dp.DataTable(df_all_show),
                    dp.Text("---"),
                    dp.Text("### Complete DVGs plots"),
                    ## scatterplots
                    dp.Select(blocks=[
                        dp.Plot(all_scatter_vir, label= "ViReMa-a detection",
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the  number of DVG readings " + 
                        "detected by ViReMa-a."),
                        dp.Plot(all_scatter_dit, label="DI-tector detection",
                        caption="Scatterplot of RI~BP coordinates. " + 
                        "Bubble size is relative to the number of DVG readings " + 
                        "detected by DI-tector.")
                    ]),
                    ## histograms
                    dp.Select(blocks=[
                        dp.Plot(all_hist_vir, label= "ViReMa-a detection",
                        caption="Frequency histogram of the theoretical lengths of the DVGs."),
                        dp.Plot(all_hist_dit, label="DI-tector detection",
                        caption="Frequency histogram of the theoretical lengths of the DVGs.")
                    ]),
                    ## Arc diagrams
                    dp.Select(blocks=[
                        dp.Plot(all_arcs_vir, label= "ViReMa-a detection",
                        caption="Arc diagram of the deletions detected. The range " + 
                        " between the two 'x' coordinates represent the deletioned " +
                        " genome part."),
                        dp.Plot(all_arcs_dit, label="DI-tector detection",
                        caption="Arc diagram of the deletions detected. The range " + 
                        " between the two 'x' coordinates represent the deletioned " +
                        " genome part.")
                    ]),
                    dp.HTML(arcs_caption()),
                    dp.HTML(footer_html())
                    ]            
                ),
            dp.Page(
                # FILTERED ML
                #----------------------------------
                title="DVG ML filtered",
                blocks=[
                    dp.HTML(header_html(sample_name)),
                    dp.Text("## Table of DVG predicted as reals with the DVGfinder-MLmodel"),
                    dp.DataTable(df_predicted_reals_show),
                    dp.Text("---"),
                    dp.Text("### ML filtered DVG plots"),
                    ## scatterplot
                    dp.Plot(ML_scatter_dit, label= "DI-tector detection",
                    caption="Scatterplot of RI~BP coordinates. " + 
                    "Bubble size is relative to the  number of DVG readings " + 
                    "detected by DI-tector."),
                    ## histograms
                    dp.Plot(ML_hist_dit, label= "DI-tector detection",
                    caption="Frequency histogram of the theoretical lengths of the DVGs."),
                    ## Arc diagrams
                    dp.Plot(ML_arcs_dit, label= "DI-tector detection",
                    caption="Arc diagram of the deletions detected. The range " + 
                    " between the two x coordinates represent the deletioned " +
                    " genome part."),
                    dp.HTML(arcs_caption()),
                    dp.HTML(footer_html())
                ])
            )
    
    else:
        r = report_onlyAll(sample_name, 
                df_all_show, all_scatter_vir, all_scatter_dit,
                all_hist_vir, all_hist_dit, all_arcs_vir, all_arcs_dit)

    r.save(path=path_html, open=True)


def report_AllConsensus(sample_name, 
                        df_all_show, all_scatter_vir, all_scatter_dit,
                        all_hist_vir, all_hist_dit, all_arcs_vir, all_arcs_dit,
                        df_consensus, consensus_scatter, consensus_hist,
                        consensus_arcs):
    r = dp.Report(
    dp.Page(
        # COMPLETE
        #----------------------------------
        title="All DVGs detected",        
        blocks=[
            dp.HTML(header_html(sample_name)),
            dp.Text("## Complete table of DVG detected by ViReMa-a _OR_ DI-tector"),
            dp.DataTable(df_all_show),
            dp.Text("---"),
            dp.Text("### Complete DVGs plots"),
            ## scatterplots
            dp.Select(blocks=[
                dp.Plot(all_scatter_vir, label= "ViReMa-a detection",
                caption="Scatterplot of RI~BP coordinates. " + 
                "Bubble size is relative to the  number of DVG readings " + 
                "detected by ViReMa-a."),
                dp.Plot(all_scatter_dit, label="DI-tector detection",
                caption="Scatterplot of RI~BP coordinates. " + 
                "Bubble size is relative to the number of DVG readings " + 
                "detected by DI-tector.")
            ]),
            ## histograms
            dp.Select(blocks=[
                dp.Plot(all_hist_vir, label= "ViReMa-a detection",
                caption="Frequency histogram of the theoretical lengths of the DVGs."),
                dp.Plot(all_hist_dit, label="DI-tector detection",
                caption="Frequency histogram of the theoretical lengths of the DVGs.")
            ]),
            ## Arc diagrams
            dp.Select(blocks=[
                dp.Plot(all_arcs_vir, label= "ViReMa-a detection",
                caption="Arc diagram of the deletions detected. The range " + 
                " between the two 'x' coordinates represent the deletioned " +
                " genome part."),
                dp.Plot(all_arcs_dit, label="DI-tector detection",
                caption="Arc diagram of the deletions detected. The range " + 
                " between the two 'x' coordinates represent the deletioned " +
                " genome part.")
            ]),
            dp.HTML(arcs_caption()),
            dp.HTML(footer_html())
            ]            
        ),
    dp.Page(
        # CONSENSUS
        #----------------------------------
        title="DVG consensus selected",
        blocks=[
            dp.HTML(header_html(sample_name)),
            dp.Text("## Consensus table of DVG detected by ViReMa-a _AND_ DI-tector"),
            dp.DataTable(df_consensus),
            dp.Text("---"),
            dp.Text("### Consensus DVG plots"),
            dp.Plot(consensus_scatter, 
                caption="Scatterplot of RI~BP coordinates. " + 
                "Bubble size is relative to the **mean** number of DVG readings " + 
                "detected by the programs."),
            dp.Plot(consensus_hist,
            caption="Frequency histogram of the theoretical lengths of the DVGs."),
            dp.Plot(consensus_arcs,
            caption="Arc diagram of the deletions detected. The range " + 
                " between the two x coordinates represent the deletioned " +
                " genome part."),
            dp.HTML(arcs_caption()),
            dp.HTML(footer_html())
            ]
        )
    )
    return r

def report_AllML(sample_name, 
                df_all_show, all_scatter_vir, all_scatter_dit,
                all_hist_vir, all_hist_dit, all_arcs_vir, all_arcs_dit,
                df_predicted_reals_show, ML_scatter_vir, ML_scatter_dit, 
                ML_hist_vir, ML_hist_dit, ML_arcs_vir, ML_arcs_dit):
    r = dp.Report(
    dp.Page(
        # COMPLETE
        #----------------------------------
        title="All DVGs detected",        
        blocks=[
            dp.HTML(header_html(sample_name)),
            dp.Text("## Complete table of DVG detected by ViReMa-a _OR_ DI-tector"),
            dp.DataTable(df_all_show),
            dp.Text("---"),
            dp.Text("### Complete DVGs plots"),
            ## scatterplots
            dp.Select(blocks=[
                dp.Plot(all_scatter_vir, label= "ViReMa-a detection",
                caption="Scatterplot of RI~BP coordinates. " + 
                "Bubble size is relative to the  number of DVG readings " + 
                "detected by ViReMa-a."),
                dp.Plot(all_scatter_dit, label="DI-tector detection",
                caption="Scatterplot of RI~BP coordinates. " + 
                "Bubble size is relative to the number of DVG readings " + 
                "detected by DI-tector.")
            ]),
            ## histograms
            dp.Select(blocks=[
                dp.Plot(all_hist_vir, label= "ViReMa-a detection",
                caption="Frequency histogram of the theoretical lengths of the DVGs."),
                dp.Plot(all_hist_dit, label="DI-tector detection",
                caption="Frequency histogram of the theoretical lengths of the DVGs.")
            ]),
            ## Arc diagrams
            dp.Select(blocks=[
                dp.Plot(all_arcs_vir, label= "ViReMa-a detection",
                caption="Arc diagram of the deletions detected. The range " + 
                " between the two 'x' coordinates represent the deletioned " +
                " genome part."),
                dp.Plot(all_arcs_dit, label="DI-tector detection",
                caption="Arc diagram of the deletions detected. The range " + 
                " between the two 'x' coordinates represent the deletioned " +
                " genome part.")
            ]),
            dp.HTML(arcs_caption()),
            dp.HTML(footer_html())
            ]            
        ),
        dp.Page(
        # FILTERED ML
        #----------------------------------
        title="DVG ML filtered",
        blocks=[
            dp.HTML(header_html(sample_name)),
            dp.Text("## Table of DVG predicted as reals with the DVGfinder-MLmodel"),
            dp.DataTable(df_predicted_reals_show),
            dp.Text("---"),
            dp.Text("### ML filtered DVG plots"),
            ## scatterplots
            dp.Select(blocks=[
                dp.Plot(ML_scatter_vir, label= "ViReMa-a detection",
                caption="Scatterplot of RI~BP coordinates. " + 
                "Bubble size is relative to the  number of DVG readings " + 
                "detected by ViReMa-a."),
                dp.Plot(ML_scatter_dit, label="DI-tector detection",
                caption="Scatterplot of RI~BP coordinates. " + 
                "Bubble size is relative to the number of DVG readings " + 
                "detected by DI-tector."),
            ]),
            ## histograms
            dp.Select(blocks=[
                dp.Plot(ML_hist_vir, label= "ViReMa-a detection",
                caption="Frequency histogram of the theoretical lengths of the DVGs."),
                dp.Plot(ML_hist_dit, label="DI-tector detection",
                caption="Frequency histogram of the theoretical lengths of the DVGs.")
            ]),
            ## Arc diagrams
            dp.Select(blocks=[
                dp.Plot(ML_arcs_vir, label= "ViReMa-a detection",
                caption="Arc diagram of the deletions detected. The range " + 
                " between the two x coordinates represent the deletioned " +
                " genome part."),
                dp.Plot(ML_arcs_dit, label="DI-tector detection",
                caption="Arc diagram of the deletions detected. The range " + 
                " between the two x coordinates represent the deletioned " +
                " genome part.")]
            ),
            dp.HTML(arcs_caption()),
            dp.HTML(footer_html())
        ])
    )
    return r

def report_onlyAll(sample_name, 
                df_all_show, all_scatter_vir, all_scatter_dit,
                all_hist_vir, all_hist_dit, all_arcs_vir, all_arcs_dit):
    r = dp.Report(
    dp.Page(
        # COMPLETE
        #----------------------------------
        title="All DVGs detected",        
        blocks=[
            dp.HTML(header_html(sample_name)),
            dp.Text("## Complete table of DVG detected by ViReMa-a _OR_ DI-tector"),
            dp.DataTable(df_all_show),
            dp.Text("---"),
            dp.Text("### Complete DVGs plots"),
            ## scatterplots
            dp.Select(blocks=[
                dp.Plot(all_scatter_vir, label= "ViReMa-a detection",
                caption="Scatterplot of RI~BP coordinates. " + 
                "Bubble size is relative to the  number of DVG readings " + 
                "detected by ViReMa-a."),
                dp.Plot(all_scatter_dit, label="DI-tector detection",
                caption="Scatterplot of RI~BP coordinates. " + 
                "Bubble size is relative to the number of DVG readings " + 
                "detected by DI-tector.")
            ]),
            ## histograms
            dp.Select(blocks=[
                dp.Plot(all_hist_vir, label= "ViReMa-a detection",
                caption="Frequency histogram of the theoretical lengths of the DVGs."),
                dp.Plot(all_hist_dit, label="DI-tector detection",
                caption="Frequency histogram of the theoretical lengths of the DVGs.")
            ]),
            ## Arc diagrams
            dp.Select(blocks=[
                dp.Plot(all_arcs_vir, label= "ViReMa-a detection",
                caption="Arc diagram of the deletions detected. The range " + 
                " between the two 'x' coordinates represent the deletioned " +
                " genome part."),
                dp.Plot(all_arcs_dit, label="DI-tector detection",
                caption="Arc diagram of the deletions detected. The range " + 
                " between the two 'x' coordinates represent the deletioned " +
                " genome part.")
            ]),
            dp.HTML(arcs_caption()),
            dp.HTML(footer_html())
            ]            
        )
    )
    return r

def generate_incomplete_report(df, sample_name, len_wt, case):
    """
    Función completa para generar el informe HTML para los casos incompletos,
    ie, casos en los que solo un programa haya encontrado DVGs en los datos.
    - LLama a las funciones necesarias para crear la tabla de DVGs a mostrar
    y las visualizaciones, y 
    - Monta el resultado en un fichero HTML (utiliza la librería datapane)
    que se abre automáticamente en el navegador establecido por defecto. Guarda
    el HTML y las tablas de los tres modos (.csv) en el directorio 
    'FinalReports/{sample_name}/'

    Args:
        df          (pd.DataFrame)  Tabla completa con todos los eventos 
                                detectados con el módulo metabuscador
        sample_name     (str)   Nombre de la muestra.
        len_wt          (int)   Longitud del genoma de referencia
        case    (str)   Caso de uso. 
                        case_2: solo output de ViReMa-a
                        case_3: solo output de DI-tector
    
    """
    # Para extraer el nombre del programa según lo necesitemos
    case_dict = {"case_2" : ["ViReMa-a", "virema"], 
                "case_3" : ["DI-tector", "ditector"]}
    
    # file to save
    paths_csv = 'FinalReports' + "/" + sample_name + "/" + \
                sample_name
    path_html = 'FinalReports' + "/" + sample_name + "/" + \
                sample_name + "_report.html"
    # Generamos las tablas a mostrar
    # -------------------------------------------------------------------------
    ## Todos los eventos reduciendo a las variables de interés
    df_all_show = create_table_all(df, case)
    

    # Escribimos la tabla
    # -------------------------------------------------------------------------
    df_all_show.to_csv(paths_csv + "_{}.csv".format(case_dict[case][1]), index=False)
    

    # Generamos las visualizaciones
    # -------------------------------------------------------------------------
    ## COMPLETE
    ### scatterplots
    scatter_plot = visualization.scatter_color_by_DVGtype(df_all_show, 
                "read_counts_{}".format(case_dict[case][1]), 
                title="Relationship between BP and RI coordinates in all " +  
                "the DVGs detected by {}".format(case_dict[case][0]))   
    
    ### histograms
    hist_plot = visualization.histogram_distribution_lengthDVG(df_all_show,
                 "read_counts_{}".format(case_dict[case][1]), 
                 title="Distribution of DVG lengths detected by {}"\
                     .format(case_dict[case][0]))
    
    ### arc diagrams
    arcs_plot = visualization.arcs_diagram(df_all_show, len_wt, 
                title="Deletions detected by {}".format(case_dict[case][0]))
    

    # Montamos el report
    # -------------------------------------------------------------------------
    r = dp.Report(
    dp.Page(
        # CONSENSUS
        #----------------------------------
        title="DVG detected by {}".format(case_dict[case][0]),
        blocks=[
            dp.HTML(header_html(sample_name)),
            dp.Text("## Table of DVGs detected by {}".format(case_dict[case][0])),
            dp.DataTable(df_all_show),
            dp.Text("---"),
            dp.Text("### DVG plots"),
            dp.Plot(scatter_plot, 
                caption="Scatterplot of RI~BP coordinates. " + 
                "Bubble size is relative to the number of DVG readings " + 
                "detected by {}.".format(case_dict[case][0])),
            dp.Plot(hist_plot,
            caption="Frequency histogram of the theoretical lengths of the DVGs."),
            dp.Plot(arcs_plot,
            caption="Arc diagram of the deletions detected. The range " + 
                " between the two x coordinates represent the deletioned " +
                " genome part."),
            dp.HTML(arcs_caption()),
            dp.HTML(footer_html())
            ]
        )
    )
    r.save(path=path_html, open=True)


def header_html(sample_name):
    header = """
    <html>
        <style type='text/css'>
        @keyframes colorines {
            0%   {color: rgb(193, 58, 211);}
            25%  {color: #EC4899;}
            50%  {color: #8B5CF6;}
            100% {color: #EF4444;}
        }
        #container {
            background: #303844;
            padding: 1.2em;
        }
        img {
            /*text-align: center*/
        }
        h2 {
            /*color:#deef44;
            color:rgb(193, 58, 211);
            text-align: center;*/
            color:#eee;
            animation-name: colorines;
            animation-duration: 6s;
            animation-iteration-count: infinite;
            font-size: 25px;
            font-weight: 200;
            margin-left: 9px;
        }
    </style>
    <div id="container">
      <img src="http://147.156.206.144/appweb/LOGODVGfinder_simple_white.png" width=300>
      <h2> Sample """ + sample_name + """ </h2>
    </div>
    </html>"""

    return header

def footer_html():
    footer = """
    <footer>
        <style type='text/css'>
            #container {
                background: #303844;
                padding: 0.3em;
                text-align: right;
            }
            h1 {
                color:#eee;
                font-size: 14px;
                text-align: right;
                font-style: italic;
            }
            a {
                color:#eee;
                font-size: 14px;
                text-align: right;
            }
            h2 {
                color:#deef44;
                font-size: 22px;
                text-align: right;
            }
        </style>
        <div id="container">
        <h1> María José Olmo-Uceda, PhD student</h1>
        <a href = "mailto: mariajose.olmo@csic.es"> mariajose.olmo@csic.es</a>
        <h2> EvolSysVir Group, I2SysBio (CSIC-UV) <\h2>
        </div>
    </footer>"""

    return footer

def arcs_caption():
    arc_caption = """
    <!DOCTYPE html>
    <html>
    <body>
        <style type='text/css'>
            #container {
                background: #ffffff;
                padding: 0.3em;
                margin-left: auto;
                width: 1000px;
                }
            h1 {
                color:rgb(47, 63, 92);
                font-size: 12px;
                text-align: left;
                font-family: "Open Sans", verdana, arial, sans-serif;      
                font-weight: lighter; 
                margin-bottom: 0;
            }
            h2 {
                color:rgb(47, 63, 92);
                font-size: 12px;
                text-align: left;
                width:100px;
                margin-left: 23px;
                font-family: "Open Sans", verdana, arial, sans-serif;
                font-weight: lighter;            
            }
            #cuadrado_azul {
                height:12px ;
                width: 12px;
                background-color: #636EFA;
                color:#636EFA;
            }
            
            #cuadrado_marron {
                height:12px ;
                width: 12px;
                background-color: #b88f6f;
                color:#b88f6f;
            }
            
            #cuadrado_naranja {
                height:12px ;
                width: 12px;
                background-color: #FFA15A;
                color:#FFA15A;
            }
            #titulo {
                width: 100px;
            }
            #leyendita div{
                float:left;
                margin-top: 5px
        
            }
            #leyendita {
                width: 15%;	
                float: left;
                top: 0;
            }
            #leyendita h2{
                margin-top: 4px;
            }
        </style>
        <div id="container">
            <div id="titulo"> 
                <h1> DVG_type</h1>
            </div>
            <div id="leyendita">
                <div id="cuadrado_azul"> &nbsp; </div>
                <h2> Deletion_forward </h2>
            </div>
            <div id="leyendita">
            <div id="cuadrado_naranja"> &nbsp;</div>
                <h2> Deletion_reverse </h2>
            </div>
            <div id="leyendita">
            <div id="cuadrado_marron"> &nbsp;</div>
                <h2> both </h2>
            </div>
        </div>

    </body>
    </html> """
    return arc_caption


# Versión anterior sin logo
"""
    <html>
        <style type='text/css'>
            @keyframes example {
                0%   {color: rgb(193, 58, 211);}
                25%  {color: #EC4899;}
                50%  {color: #8B5CF6;}
                100% {color: #EF4444;}
            }
            #container {
                background: #1F2937;
                padding: 0.6em;
                /*position: sticky;
                top: 0;*/
            }
            h1 {
                color:rgb(193, 58, 211);
                animation-name: example;
                animation-duration: 6s;
                animation-iteration-count: infinite;
            }
            h2 {
                color: #deef44
            }
        </style>
        <div id="container">
        <h1> DVGfinder report </h1>
        <h2> Sample """ #+ sample_name + """ <\h2>
        #</div>"""
