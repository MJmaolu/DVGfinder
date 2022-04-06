#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
## Functions for the PREDICTION module
##                  DVGfinder: Defective Viral Genome-finder
###############################################################################
## Genera la predicción y selecciona los DVGs predichos como reales de entre
## todos los detectados por el módulo METABUSCADOR. Devuelve el resultado en 
## un pd.DataFrame con las variables que se quieren mostrar en el informe
## HTML
###############################################################################
## Author: Maria Jose Olmo-Uceda
## Version: 3.0
## Email: mariajose.olmo@csic.es
## Date: 2021/12
###############################################################################

# third party imports
import pandas as pd
import numpy as np
import pickle

# DVGfinder functions
from Models import prediction
from Models import reports


selected_features = ['rpht_virema', 'rpht_ditector',
    'pBP_virema', 'pRI_virema', 'pBP_ditector', 'pRI_ditector',
    'sdrm_junction_bt2_mapped', 'sdrm_junction_bm_mapped', 
    'sdrm_junction_bt2_H', 'sdrm_junction_bm_H',
    'sdrm_neighborsBP_bt2_mapped', 'sdrm_neighborsRI_bt2_mapped',
    'sdrm_neighborsBP_bm_mapped', 'sdrm_neighborsRI_bm_mapped',
    'sdrm_neighborsBP_bt2_H', 'sdrm_neighborsRI_bt2_H', 
    'sdrm_neighborsBP_bm_H', 'sdrm_neighborsRI_bm_H']

# Preparación del df para dejarlo en el formato que acepta el modelo
# -----------------------------------------------------------------------------

# Obtenemos el df_metrics
def calculate_sdrm(a,b):
    """
    Calcula la semi-diferencia relativa a la media de a y b
    """
    if a + b == 0:
        sdrm = 1
    else:
        sdrm = abs(a - b)/(a + b)
    return sdrm

def df_to_metrics_df(df):
    """
    Genera la tabla de métricas a partir del df_resume. 
    Las métricas generadas son las sdrm (semi-diferencia relativa a la media) 
    de cada par de valores. Se tiene en cuenta que pueden haber 0s.
    Además se mantienen las métricas generadas anteriormente de 'rpht' y pCoordinates.

    Args:
        df  (pd.DataFrame)  Tabla con la información de todos los eventos.
    Return:
        df_metrics  (pd.DataFrame)  Tabla con todas las métricas generadas de 
                            los eventos del df
    """

    metrics = ['rpht_virema', 'rpht_ditector',
    'pBP_virema', 'pRI_virema', 'pBP_ditector', 'pRI_ditector',
    'sdrm_junction_bt2_mapped', 'sdrm_junction_bm_mapped', 
    'sdrm_junction_bt2_H', 'sdrm_junction_bm_H',
    'sdrm_neighborsBP_bt2_mapped', 'sdrm_neighborsRI_bt2_mapped',
    'sdrm_neighborsBP_bm_mapped', 'sdrm_neighborsRI_bm_mapped',
    'sdrm_neighborsBP_bt2_H', 'sdrm_neighborsRI_bt2_H', 
    'sdrm_neighborsBP_bm_H', 'sdrm_neighborsRI_bm_H', 'labels']

    df_w_metrics = df[:] # unlinked copy

    # semidiferencia de las profundidades en la 'junction' relativa a su media
    df_w_metrics['sdrm_junction_bt2_mapped'] = df_w_metrics.apply(lambda row:
                        calculate_sdrm(row['depthBP_bt2_mapped'], 
                                    row['depthRI_bt2_mapped']), 
                                    axis=1)
    df_w_metrics['sdrm_junction_bt2_H'] = df_w_metrics.apply(lambda row:
                        calculate_sdrm(row['depthBP_bt2_H'], 
                                    row['depthRI_bt2_H']), 
                                    axis=1)
    df_w_metrics['sdrm_junction_bm_mapped'] = df_w_metrics.apply(lambda row:
                        calculate_sdrm(row['depthBP_bm_mapped'], 
                                    row['depthRI_bm_mapped']), 
                                    axis=1)
    df_w_metrics['sdrm_junction_bm_H'] = df_w_metrics.apply(lambda row:
                        calculate_sdrm(row['depthBP_bm_H'], 
                                    row['depthRI_bm_H']), 
                                    axis=1)
    # semidiferencia de las profundidades en los alrededores de la 'junction'
    # relativa a su media
    ## En los alineamientos completos
    df_w_metrics['sdrm_neighborsBP_bt2_mapped'] = df_w_metrics.apply(lambda row:
                        calculate_sdrm(row['mean_depth_pre5BP_bt2_mapped'], 
                                    row['mean_depth_post5BP_bt2_mapped']), 
                                    axis=1)
    df_w_metrics['sdrm_neighborsRI_bt2_mapped'] = df_w_metrics.apply(lambda row:
                        calculate_sdrm(row['mean_depth_pre5RI_bt2_mapped'], 
                                    row['mean_depth_post5RI_bt2_mapped']), 
                                    axis=1)
    df_w_metrics['sdrm_neighborsBP_bm_mapped'] = df_w_metrics.apply(lambda row:
                        calculate_sdrm(row['mean_depth_pre5BP_bm_mapped'], 
                                    row['mean_depth_post5BP_bm_mapped']), 
                                    axis=1)
    df_w_metrics['sdrm_neighborsRI_bm_mapped'] = df_w_metrics.apply(lambda row:
                        calculate_sdrm(row['mean_depth_pre5RI_bm_mapped'], 
                                    row['mean_depth_post5RI_bm_mapped']), 
                                    axis=1)
    ## En los alineamientos con Hard-clipping   
    df_w_metrics['sdrm_neighborsBP_bt2_H'] = df_w_metrics.apply(lambda row:
                        calculate_sdrm(row['mean_depth_pre5BP_bt2_H'], 
                                    row['mean_depth_post5BP_bt2_H']), 
                                    axis=1)
    df_w_metrics['sdrm_neighborsRI_bt2_H'] = df_w_metrics.apply(lambda row:
                        calculate_sdrm(row['mean_depth_pre5RI_bt2_H'], 
                                    row['mean_depth_post5RI_bt2_H']), 
                                    axis=1)
    df_w_metrics['sdrm_neighborsBP_bm_H'] = df_w_metrics.apply(lambda row:
                        calculate_sdrm(row['mean_depth_pre5BP_bm_H'], 
                                    row['mean_depth_post5BP_bm_H']), 
                                    axis=1)
    df_w_metrics['sdrm_neighborsRI_bm_H'] = df_w_metrics.apply(lambda row:
                        calculate_sdrm(row['mean_depth_pre5RI_bm_H'], 
                                    row['mean_depth_post5RI_bm_H']), 
                                    axis=1)

    # Generamos un nuevo df solo con las métricas
    df_metrics = df_w_metrics[df_w_metrics.columns.intersection(metrics)]

    return df_metrics

# Reducimos el df_metrics a las variables que según el método RFE son más 
# informativas para nuestro modelo

def df_to_df_metrics_reduced(df):

    # Los 17 predictores seleccionados por RFE (durante fase de entrenamiento 
    # del modelo) 
    selected_metrics = ['rpht_virema', 'rpht_ditector',
    'pBP_virema', 'pRI_virema', 'pBP_ditector', 'pRI_ditector',
    'sdrm_junction_bt2_mapped', 'sdrm_junction_bm_mapped', 
    'sdrm_junction_bt2_H', 'sdrm_junction_bm_H',
    'sdrm_neighborsBP_bt2_mapped', 'sdrm_neighborsRI_bt2_mapped',
    'sdrm_neighborsBP_bm_mapped', 'sdrm_neighborsRI_bm_mapped',
    'sdrm_neighborsRI_bt2_H', 
    'sdrm_neighborsBP_bm_H', 'sdrm_neighborsRI_bm_H']

    # tabla de métricas
    df_metrics = df_to_metrics_df(df)

    # reducida a las variables que acepta el modelo
    df_metrics_reduced = df_metrics[df_metrics.columns.\
                                    intersection(selected_metrics)]
    df_metrics_reduced.columns = ['rpht_virema', 'rpht_ditector',
    'pBP_virema', 'pRI_virema', 'pBP_ditector', 'pRI_ditector',
    'sdrm_junction_bt_mapped', 'sdrm_junction_bm_mapped', 
    'sdrm_junction_bt_H', 'sdrm_junction_bm_H',
    'sdrm_neighborsBP_bt_mapped', 'sdrm_neighborsRI_bt_mapped',
    'sdrm_neighborsBP_bm_mapped', 'sdrm_neighborsRI_bm_mapped',
    'sdrm_neighborsRI_bt_H', 
    'sdrm_neighborsBP_bm_H', 'sdrm_neighborsRI_bm_H']

    return df_metrics_reduced

def df_to_df_metrics_reduced_v3(df):

    # Los 17 predictores seleccionados por RFE (durante fase de entrenamiento 
    # del modelo) 
    selected_metrics = ['rpht_virema', 'rpht_ditector',
    'pBP_virema', 'pRI_virema', 'pBP_ditector', 'pRI_ditector',
    'sdrm_junction_bt2_mapped', 'sdrm_junction_bm_mapped', 
    'sdrm_junction_bt2_H', 'sdrm_junction_bm_H',
    'sdrm_neighborsBP_bt2_mapped', 'sdrm_neighborsRI_bt2_mapped',
    'sdrm_neighborsBP_bm_mapped', 'sdrm_neighborsRI_bm_mapped',
    'sdrm_neighborsBP_bt2_H', 'sdrm_neighborsRI_bt2_H', 
    'sdrm_neighborsBP_bm_H', 'sdrm_neighborsRI_bm_H']

    # tabla de métricas
    df_metrics = df_to_metrics_df(df)

    # reducida a las variables que acepta el modelo
    df_metrics_reduced = df_metrics[df_metrics.columns.\
                                    intersection(selected_metrics)]
    """
    df_metrics_reduced.columns = ['rpht_virema', 'rpht_ditector',
    'pBP_virema', 'pRI_virema', 'pBP_ditector', 'pRI_ditector',
    'sdrm_junction_bt_mapped', 'sdrm_junction_bm_mapped', 
    'sdrm_junction_bt_H', 'sdrm_junction_bm_H',
    'sdrm_neighborsBP_bt_mapped', 'sdrm_neighborsRI_bt_mapped',
    'sdrm_neighborsBP_bm_mapped', 'sdrm_neighborsRI_bm_mapped',
    'sdrm_neighborsBP_bt_H', 'sdrm_neighborsRI_bt_H', 
    'sdrm_neighborsBP_bm_H', 'sdrm_neighborsRI_bm_H']
    """
    return df_metrics_reduced


# Aplicación del modelo
# -----------------------------------------------------------------------------

def generate_prediction_and_add_to_df_resume(df, model_file):
    """
    Carga el modelo y genera la predicción

    Args:
        df  (pd.DataFrame)  Toda la información extraída con el módulo 
                            metabuscador
        model_file  (str)   Fichero .sav con el modelo de predicción 

    Return:
        df_w_pred   (pd.DataFrame)  df (completo) con la predicción (1/0) 
                                    añadida como nueva columna 'prediction'
    """
    # load model
    model = pickle.load(open(model_file, 'rb'))
    # generación del df con las métricas que acepta el modelo
    #df_metrics = df_to_metrics_df(df)
    df_metrics = df_to_df_metrics_reduced(df)

    # predicción
    y_pred = model.fit(df_metrics, model_file) # vector predicción

    df_w_pred = df[:] # unlinked copy
    df_w_pred['prediction'] = y_pred # añadimos la predicción al df completo

    return df_w_pred

def filter_predicted_as_reals_v3(df, model_file, threshold):
    """
    Función completa a la que se llama desde el controlador (DVGfinder.py).
    Realiza todo el proceso de: 
        - generar las métricas a partir del df completo,
        - reducir los predictores a los que acepta el modelo
        - Carga el modelo y genera la predicción
        - Extrae los DVGs predichos como reales y devuelve la tabla a mostrar
            (features_to_show) en el informe HTML
    
    Args:
        df  (pd.DataFrame)  Tabla completa resultado del módulo metabuscador
        model_file  (str)   Fichero .sav con el modelo de predicción
    
    Return:
        df_ML_show  (pd.DataFrame)  Tabla con los DVGs predichos como reales
                                por el modelo y las columnas que se quieren
                                mostrar (features_to_show)
    """

    features_to_show = ['cID_DI', 'p(real)', 'BP', 'RI', 'sense', 'DVG_type', 
        'length_dvg', 'read_counts_virema', 'pBP_virema', 'pRI_virema',
        'rpht_virema', 'read_counts_ditector','pBP_ditector', 'pRI_ditector', 
        'rpht_ditector']

    # Cargamos el modelo
    model = pickle.load(open(model_file, 'rb'))

    # generación del df con las métricas que acepta el modelo
    #df_metrics = df_to_metrics_df(df)  # all metrics
    df_metrics = df_to_df_metrics_reduced_v3(df)   # 18 features

    # predicción
    y_pred = model.predict_proba(df_metrics)[:,1] # vector with p(label1)

    df_w_pred = df[:] # unlinked copy
    df_w_pred['p(real)'] = y_pred # añadimos la predicción al df

    # filtramos DVGs predichos como reales    
    df_ML = df_w_pred[df_w_pred['p(real)'] >= threshold]
    # Dejamos solo las variables a mostrar
    df_ML_show = df_ML[df_ML.columns.intersection(features_to_show)]

    return df_ML_show

def filter_predicted_as_reals(df, model_file):
    """
    Función completa a la que se llama desde el controlador (DVGfinder.py).
    Realiza todo el proceso de: 
        - generar las métricas a partir del df completo,
        - reducir los predictores a los que acepta el modelo
        - Carga el modelo y genera la predicción
        - Extrae los DVGs predichos como reales y devuelve la tabla a mostrar
            (features_to_show) en el informe HTML
    
    Args:
        df  (pd.DataFrame)  Tabla completa resultado del módulo metabuscador
        model_file  (str)   Fichero .sav con el modelo de predicción
    
    Return:
        df_ML_show  (pd.DataFrame)  Tabla con los DVGs predichos como reales
                                por el modelo y las columnas que se quieren
                                mostrar (features_to_show)
    """

    features_to_show = ['cID_DI', 'BP', 'RI', 'sense', 'DVG_type', 'length_dvg',
        'read_counts_virema', 'pBP_virema', 'pRI_virema' ,'rpht_virema',
        'read_counts_ditector','pBP_ditector', 'pRI_ditector', 'rpht_ditector']

    # Cargamos el modelo
    model = pickle.load(open(model_file, 'rb'))

    # generación del df con las métricas que acepta el modelo
    #df_metrics = df_to_metrics_df(df)  # all metrics
    df_metrics = df_to_df_metrics_reduced(df)   # 17 features

    # predicción
    y_pred = model.predict(df_metrics) # vector predicción

    df_w_pred = df[:] # unlinked copy
    df_w_pred['prediction'] = y_pred # añadimos la predicción al df

    # filtramos DVGs predichos como reales    
    df_ML = df_w_pred[df_w_pred['prediction'] == 1]
    # Dejamos solo las variables a mostrar
    df_ML_show = df_ML[df_ML.columns.intersection(features_to_show)]

    return df_ML_show
