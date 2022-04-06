#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Funciones útiles para la fase de generación del modelo de ML para el programa 
DVGfinder
Basado en el homólogo que tengo en choose+trainingML

MJ
2021/05/07
"""
import os
from numpy.core.fromnumeric import mean
import pandas as pd
import numpy as np
from pandas.core.tools.numeric import to_numeric
from pandas.io.parsers import TextParser


def concatenate_event(DVG_type, BP, RI):
    """
    Create a string with the form DVG_type + _ + BP + _ + RI.
    """
    real_event =  DVG_type + '_' + str(BP)+ '_' + str(RI)
    return real_event
    

def extract_real_events(input_reals_csv):
    """
    Extrae la lista con los eventos reales
    Los introduce como cadena concatenada de DVGtype_BP_RI

    Args:
        input_reals_csv (str)   Path to the input csv file with the real events
                                It has to content at least the columns:
                                BP, RI and DVG_type
    
    Returns:
        real_dvgs   (list)  List with the real dvg events where each event is
                            a str formed by 'DVG_type' + '_' + 'BP' + '_' + 'RI'
    """
    df = pd.read_csv(input_reals_csv, header=0)
    real_dvgs = list(df.apply(lambda row : concatenate_event(row['DVG_type'], 
                                        row['BP'], row['RI']), axis=1))
    # add also the complementary dvgs, ie, the product of the replication of 
    # the 'original' events.
    real_dvgs = real_dvgs + list(df.apply(lambda row : complement_dvg(row['DVG_type'], 
                                        row['BP'], row['RI']), axis=1))

    return real_dvgs



def asign_status_event(DVG_type, BP, RI, real_dvgs):
    """
    1º concatena el evento en la forma 'DVG_type' + '_' + 'BP' + '_' + 'RI' y 
    consulta si está en la lista de dvgs reales. En caso afirmativo devuelve
    True y si no False. 
    Si la lista 'real_dvgs' está hecha con la función 
    extrac_real_events también tiene en cuenta los eventos complementarios

    Args:
        DVG_type
        BP
        RI
        real_dvgs   (list)  Lista con los eventos considerados reales. Pueden
                            (y deben) estar incluídos también los complementarios
    Returns:
        status  (int)  [out]   True for the real events and False for the 
                                artefacts
    """

    event = DVG_type + '_' + str(BP)+ '_' + str(RI)
    #print('* ', event)
    if event in real_dvgs:
        status = 1
    else:
        status = 0
    
    return status


def create_Y_vector(input_reals_csv, unificated_table):
    """
    Generate the vector with the labels: real (1), artifact(0). The order of
    the events is exactly the same that in unificate_table.

    Args:
        input_reals_csv (str)   Path to the csv file with the real events. At
                                least has the columns 'DVG_type', 'BP' and 'RI'
        unificated_table    (str)   Path to the csv file with all the events
                                    identified by a program (optimezed for
                                    DVGfinder) at least has the columns 
                                    'DVG_type', 'BP' and 'RI'
    Returns:
        Y   (list)  Inform about the reals and artefactual events
    """

    real_dvgs = extract_real_events(input_reals_csv)

    df = pd.read_csv(unificated_table, header=0)
    
    Y = list(df.apply(lambda row : asign_status_event(row['DVG_type'], row['BP'],
                                row['RI'], real_dvgs), axis=1))
    
        
    return Y


def complement_dvg(DVG_type, BP, RI):
    """
    Returns the name of the theoretical transcription of the real event

    Args:
        DVG_type
        BP
        RI
    """
    complements = {'Deletion_forward' : 'Deletion_reverse',
                'Deletion_reverse' : 'Deletion_forward',
                'Insertion_forward' : 'Insertion_reverse',
                'Insertion_reverse' : 'Insertion_forward', 
                '5cb/sb' : '5cb/sb',
                '3cb/sb' : '3cb/sb'}
            
    complement = complements[DVG_type] + '_' + str(RI) + '_' + str(BP)

    return complement

def complement_dvg_senseID(sense, BP, RI):
    """
    Returns the name of the theoretical transcription of the real event

    Args:
        DVG_type
        BP
        RI
    """
    complements = {'++' : '--',
                '--' : '++',
                '+-' : '+-',
                '-+' : '-+'}
            
    complement = complements[sense] + '_' + str(RI) + '_' + str(BP)

    return complement

def generate_df_resume_labeled(csv_reals_composition, csv_resume, csv_labels):
    """
    Takes the csv files with the composition of the dataset (or input of 
    SDgenerator) and of the resume generated by metabuscador step of DVGfinder
    and generate a df with the selected features (this part can be modified)
    and the labels

    Args:
        csv_reals_composition   (str)   [in]    csv file with the composition 
                                        of the dataset. At least has the 
                                        columns BP, RI, DVG_type, proportion 
                                        and basename_files
        csv_resume  (str)   [in]    csv with the info generated in DVGfinder
        csv_labels  (str)   [in]    csv with the labels (0/1). For the moment
                                    this csv is produced by write_reals.py

    Returns:
        df_resume_labeled   (pd.DataFrame)  [out]   Table with the selected
                                    features and the status (artifact/real)
                                    labeled
    """
    # composition of the dataset (real events introduced)
    df_composition = pd.read_csv(csv_reals_composition, header=0)
    # table with all the events detected by DVGfinder
    df_resume = pd.read_csv(csv_resume, header=0)
    # corresponding labels to the events identified by DVGfinder
    y = pd.read_csv(csv_labels, header=None)
    y.columns = ['labels']

    # list of features considering informatives at the moment
    #---------># the idea is to reduce in metabuscador step 
    features_interest = ['cID_DI', 'sense', 'BP', 'RI', 'DVG_type', 
        'read_counts_virema', 'read_counts_ditector', 'depthBP_bt2_mapped', 
        'depthRI_bt2_mapped', 'depthBP_bt2_H', 'depthRI_bt2_H', 
        'depthBP_bm_mapped', 'depthRI_bm_mapped', 'depthBP_bm_H', 
        'depthRI_bm_H', 'pBP_virema', 'pRI_virema', 'pBP_ditector', 
        'pRI_ditector', 'rpht_virema', 'rpht_ditector', 
        'mean_depth_pre5BP_bt2_mapped', 'mean_depth_post5BP_bt2_mapped',
        'mean_depth_pre5RI_bt2_mapped', 'mean_depth_post5RI_bt2_mapped',
        'mean_depth_pre5BP_bt2_H', 'mean_depth_post5BP_bt2_H',
        'mean_depth_pre5RI_bt2_H', 'mean_depth_post5RI_bt2_H', 
        'mean_depth_pre5BP_bm_mapped', 'mean_depth_post5BP_bm_mapped', 
        'mean_depth_pre5RI_bm_mapped', 'mean_depth_post5RI_bm_mapped', 
        'mean_depth_pre5BP_bm_H', 'mean_depth_post5BP_bm_H', 
        'mean_depth_pre5RI_bm_H', 'mean_depth_post5RI_bm_H']
    # df with only the selected features
    df_interest = df_resume[features_interest]
    # df with the selected features labeled
    df_resume_labeled = pd.concat([df_interest, y], axis=1)

    return df_resume_labeled


def reduce_resume_df_to_reals(df_resume_labeled):
    """
    Returns:
        df_resume_reals (pd.DataFrame)  df_resume_labeled where labels=1, ie,
                                with the real events
    """
    # simplify the df_resume table to the theoretical real events only
    df_resume_reals = df_resume_labeled[df_resume_labeled.labels == 1]

    return df_resume_reals


def extract_program_counts(df_resume_reals, program, event):
    """
    Extract the number of reads of the event detected by the program
    Args:
        df_resume_reals (pd.DataFrame)  [in]    tabla con los eventos anotados
                                como reales a la que se le ha añadido la 
                                columna event (DVGtype_BP_RI)
        program (str)   [in]    virema or ditector
        event   (str)   [in]    complete name of the event (DVGtype_BP_RI)
    
    Returns:
        number_reads_program  (int) [out]   Number of reads detected

    """
    name_column = "read_counts_" + program

    if event in list(df_resume_reals.event):
        number_reads_program = df_resume_reals[df_resume_reals.event == event][name_column].values[0]
    else:
        number_reads_program = 0
    
    return number_reads_program


def generate_info_th_real_events(csv_reals_composition, csv_resume, csv_labels,
                                length_reads):
    """
    Generate a df with info about the teorethical real events and the detection
    of each one of them.
    
    Args: 
        df_reals_composition (str)   Path to the input csv file with the real 
                                events
                                It has to content at least the columns:
                                BP, RI, DVG_type, proportion and basename_files
        df_resume_labeled (pd.DataFrame)    Accesible table with the info generated
                                through metabuscador step, with labels and 
                                preferibly with number of variables simplified.
    Return:
        detected_th_real_events (pd.DataFrame)  df columns: th_real_event, 
                                direct/complementary, read_counts_virema, 
                                read_counts_ditector
    """
    df_resume_labeled = generate_df_resume_labeled(csv_reals_composition, \
                        csv_resume, csv_labels)

    # Only theoretical real events
    df_resume_reals = reduce_resume_df_to_reals(df_resume_labeled)

    # if there are real events detected:
    if not df_resume_reals.empty:
        # add event column (complete name of the event: DVGtype_BP_RI)
        df_resume_reals['event'] = df_resume_reals.apply(lambda row : 
                            concatenate_event(row['DVG_type'], row['BP'], 
                            row['RI']), axis=1)
        # df with the composition of real events
        df_composition = pd.read_csv(csv_reals_composition, header=0)
        # generate new df with events
        ## first add the direct events with the proportion. We also add the number
        ## of reads that has been generated the dvg in the SD and the lenght of its
        ## genome to do posterior calcs
        detected_th_real_events = df_composition[['basename_files', 'proportion',\
            'length_dvg', 'N_dvg']]
        detected_th_real_events.columns = ['direct_event', 'proportion', \
            'length_dvg', 'N_dvg']

        ## add complementary events
        detected_th_real_events['complementary_event'] = \
                            list(df_composition.apply(lambda row :\
                                complement_dvg(row['DVG_type'], row['BP'], 
                                row['RI']), axis=1))
        
        # write if the events has been detected by each program    
        ## add virema reads
        program = 'virema'
        detected_th_real_events['direct_virema(reads)'] = \
            detected_th_real_events.apply(lambda row : extract_program_counts\
                (df_resume_reals, program, row['direct_event']), axis=1)
        detected_th_real_events['complementary_virema(reads)'] = \
            detected_th_real_events.apply(lambda row : extract_program_counts\
                (df_resume_reals, program, row['complementary_event']), axis=1)

        ## add ditector reads
        program = 'ditector'
        name_column = "read_counts_" + program
        detected_th_real_events['direct_ditector(reads)'] = \
            detected_th_real_events.apply(lambda row : extract_program_counts\
                (df_resume_reals, program, row['direct_event']), axis=1)
        detected_th_real_events['complementary_ditector(reads)'] = \
            detected_th_real_events.apply(lambda row : extract_program_counts\
                (df_resume_reals, program, row['complementary_event']), axis=1)

        ## add total reads per program
        detected_th_real_events['total_virema(reads)'] = \
            detected_th_real_events.apply(lambda row : \
            row['direct_virema(reads)'] + row['complementary_virema(reads)'],
            axis=1)
        detected_th_real_events['total_ditector(reads)'] = \
            detected_th_real_events.apply(lambda row : \
            row['direct_ditector(reads)'] + row['complementary_ditector(reads)'],
            axis=1)
        
        ## add theoretical number of reads spanning the junction
        detected_th_real_events['N_reads(theoretical)'] = detected_th_real_events.\
            apply(lambda row : th_mean_depth_dvg(row['N_dvg'], length_reads, \
                row['length_dvg']), axis=1)
        
        ## add percentage of detection of each software
        detected_th_real_events['detection_virema(%)'] = detected_th_real_events.\
            apply(lambda row : percent_detected(row['total_virema(reads)'], 
                row['N_reads(theoretical)']), axis=1)
        detected_th_real_events['detection_ditector(%)'] = detected_th_real_events.\
            apply(lambda row : percent_detected(row['total_ditector(reads)'], 
                row['N_reads(theoretical)']), axis=1)
    else:
        print("There are not real DVG events detected")
        detected_th_real_events = "No real DVG events have been detected"
    
    return detected_th_real_events

def th_mean_depth_dvg(N_dvg, length_reads, length_dvg):
    """
    Returns the theoretical mean depth which the dvg has been generated in the
    synthetic dataset (SD). We take this value as the aproximate number of 
    reads at the junction.

    Args:
        N_dvg   (int)   [in]    Number of total reads of the dvg in the SD
        length_reads    (int)   [in]    Number of bases of each read
        length_dvg  (int)   [int]   Length of the dvg genome
    
    Returns:
        mean_depth  (float -> int) [out]   Theoretical mean depth, i.e., number of 
                                reads mapping the junction (BP-RI)
    """
    mean_depth = (N_dvg * length_reads)/length_dvg

    return round(mean_depth)

def percent_detected(total_program, th_reads):
    """
    Calculate the percentage of theoretical real reads detected 

    Args:
        total_program   (int)   [in]    Number of total reads with the junction
                                of the event detected by the program
        th_reads    (int)   [in]    Theoretical number of reads spanning the
                                junction in the SD
    
    Returns:
        pc  (int)   [out]   Percentage of reads detected
    """
    pc = (total_program / th_reads) * 100
    return pc


def write_detected_th_real_events(csv_reals_composition, csv_resume, 
                                    csv_labels, length_reads, out_name):
    """
    Write in a csv file the info about the real events detected by the 
    programs 
    """
    #dir_tables = "Output/tables"
    out_name = out_name + ".csv"
    detected_th_real_events = generate_info_th_real_events\
                            (csv_reals_composition, csv_resume, csv_labels, 
                            length_reads)
    if not type(detected_th_real_events) == str:
        detected_th_real_events.to_csv(out_name, index=None)

###############################################################################
## FUNCIONES PARA GENERAR LAS TABLAS DE EVALUACIÓN DE DETECCIÓN
##
## Algunas funciones son muy parecidas a las anteriores
###############################################################################

def trad_dvgType_to_sense(dvg_type):
    """
    Traduce dvg_type a sense
    """
    dvg_types = {"++" : ('Deletion_forward', 'Insertion_forward'), 
    '--' : ('Deletion_reverse', 'Insertion_reverse'), 
    '+-' : '5cb/sb', '-+' : '3cb/sb'}

    for key, value in dvg_types.items():
        if dvg_type in value:
            sense = key
    
    return sense
    

def list_events_without_filter_complementaries(df):
    """
    Genera una lista con los eventos del dataframe en formato sense_BP_RI
    Columnas: DVG_type, BP, RI
    La función ya se encarga de traducir a sense

    Args:
        df  (pd.DataFrame)  Tabla accesible con al menos las columnas
                            sense, BP, RI
    
    Return:
        list_events  (list)  Lista con todos los eventos en formato
                            cID_DI que hay en el dataframe
                            Columnas que necesita: BP, RI, DVG_type
                            La fu
    """
    
    df_ids = df[['DVG_type', 'BP', 'RI']]
    df_ids['sense'] =  df_ids.apply(lambda row : 
                        trad_dvgType_to_sense(row['DVG_type']), 
                        axis=1)

    df_ids['cID_DI'] = df_ids.apply(lambda row : complement_dvg_senseID(
                row['sense'], row['BP'], row['RI']), axis = 1)

    list_events = list(df_ids['cID_DI'])
    return list_events     

def list_events_filtering_complementaries(df):
    """
    Genera la lista de eventos dvg teniendo en cuenta que no haya duplicaciones
    por complementariedad. Va introduciendo los eventos tras comprobar que no
    son complementarios a ninguno de los existentes en la lista

    Necesita las columnas sense, BP, RI
    """
    # lista con todos los eventos en formato sense_BP_RI
    list_dvgs = list_events_without_filter_complementaries(df)


    list_events_filtered_complementarities = list()

    for dvg in list_dvgs:
        sense, BP, RI = dvg.split("_")
        complement = complement_dvg_senseID(sense, BP, RI)
        if complement not in list_dvgs:
    
            list_events_filtered_complementarities.append(dvg)
    
    return list_events_filtered_complementarities


def generate_list_events_from_composition(df_composition):
    """
    Extrae la lista de eventos reales del fichero de composición de la muestra.
    No debería, pero por si acaso comprueba que no haya eventos complementarios
    entre ellos.
    Args:
        df_composition  (pd.Dataframe)  
    
    Returns:
        list_reals  (list)  Lista con los cIDDI en formato sense_BP_RI reales
    """

    list_reals = list_events_filtering_complementaries(df_composition)

    return list_reals

def generate_list_artifacts_from_resume(df_resume_labeled):

    # solo artefactos
    df_artifacts = df_resume_labeled.loc[df_resume_labeled['labels'] == 0]
    list_artifacts = list_events_filtering_complementaries(df_artifacts)

    return list_artifacts

def create_df_with_cIDDI_and_class(df_composition, df_resume_labeled):
    """
    Genera a partir de los df de composicion y resumen etiquetado el nuevo 
    df con todos los eventos reales posibles (hayan sido o no detectados) y 
    los artefactos generados con los algoritmos. Tanto unos como en otros se 
    comprueba que no hayan duplicaciones por complementarios (es decir, que 
    aparezca el mismo evento en un sentido y en el opuesto: ++_3_5 y --_5_3).

    Args:
        df_composition
        df_resume_labeled
    Returns:
        new_df  (pd.DataFrame)  Tabla con todos los eventos positivos y negativos
                                Columnas:
                                    direct_event    (cID_DI)
                                    complementary_event (cID_DI)
                                    class   (1: real, 0: artefacto)
    """
    
    list_reals = generate_list_events_from_composition(df_composition)
    
    list_artifacts = generate_list_artifacts_from_resume(df_resume_labeled)
    # concateno las listas de cID_DI
    list_events = list_reals + list_artifacts
    # generamos la lista con los 1 y 0 ordenados por evento
    reals = list('1' * len(list_reals) + '0' * len(list_artifacts))

    # generamos lista de complementarios
    complements = list()
    for event in list_events:
        sense, BP, RI = event.split("_")
        complement = complement_dvg_senseID(sense, BP, RI)
        complements.append(complement)

    d = {'direct_event' : list_events, 'complementary_event' : complements,
        'real_class' : reals}

    # creamos el dataframe
    new_df = pd.DataFrame(d)

    return new_df


def label_resume_table(composition_csv, resume_csv, write=None):
    """
    Añade a la tabla resumen la clase predicha ('label') :
        1: El evento está en el fichero de composición (directamente o su complementario reverso) 
        0: artefacto

    Args:
        composition_csv
        resume_csv
        write   (bool)  Si es True, escribe el df en un fichero con el mismo 
                        nombre que el basename de resume_csv + _labeled.csv
    Returns:
        df_labeled  (pd.DataFrame) 
    """
    if write == None:
        write = False
    
    labeled_name = os.path.basename(resume_csv).split(".")[0] + "_labeled.csv"

    composition = pd.read_csv(composition_csv, header=0)
    resume = pd.read_csv(resume_csv, header=0)

    y_list = create_Y_vector(composition_csv, resume_csv)
    y = pd.DataFrame(y_list, columns=['labels'])
    # concatenamos labels al df resumen
    resume_labeled = pd.concat([resume, y], axis=1)
    # si se ha indicado que se escriba 
    
    if write:
        resume_labeled.to_csv(labeled_name, index=False)
        print("The resume labeled has been write in {}".format(labeled_name))
    

    return resume_labeled

def split_resume_labeled_in_reals_and_artifacts(resume_labeled):
    """
    Separa la tabla resumen en reales y artefactos en función del valor de la
    columna 'labels'.

    Args:
        resume_labeled
    Returns:
        reals_detected  (pd.DataFrame)  
        artifacts_detected  (pd.DataFrame)
    """

    reals_detected = resume_labeled.loc[resume_labeled['labels']==1]
    artifacts_detected = resume_labeled.loc[resume_labeled['labels']==0]

    return reals_detected, artifacts_detected

def generar_red_resume_labeled(composition_csv, resume_csv):
    """
    Genera el fichero resumen con la primera selección de variables de interes
    etiquetado.
    Args:
        composition_csv
        resume_csv
    
    Returns:
        df_labeled  (pd.DataFrame)  Resumen con las variables originales de 
                                interés + etiqueta de la clase a la que 
                                pertenece el evento (1: real, 0: artefacto)
    """
    
    composition = pd.read_csv(composition_csv, header=0)
    resume = pd.read_csv(resume_csv, header=0)

    features_interest = ['cID_DI', 'sense', 'BP', 'RI', 'DVG_type', 
        'read_counts_virema', 'read_counts_ditector', 'depthBP_bt2_mapped', 
        'depthRI_bt2_mapped', 'depthBP_bt2_H', 'depthRI_bt2_H', 
        'depthBP_bm_mapped', 'depthRI_bm_mapped', 'depthBP_bm_H', 
        'depthRI_bm_H', 'pBP_virema', 'pRI_virema', 'pBP_ditector', 
        'pRI_ditector', 'rpht_virema', 'rpht_ditector', 
        'mean_depth_pre5BP_bt2_mapped', 'mean_depth_post5BP_bt2_mapped',
        'mean_depth_pre5RI_bt2_mapped', 'mean_depth_post5RI_bt2_mapped',
        'mean_depth_pre5BP_bt2_H', 'mean_depth_post5BP_bt2_H',
        'mean_depth_pre5RI_bt2_H', 'mean_depth_post5RI_bt2_H', 
        'mean_depth_pre5BP_bm_mapped', 'mean_depth_post5BP_bm_mapped', 
        'mean_depth_pre5RI_bm_mapped', 'mean_depth_post5RI_bm_mapped', 
        'mean_depth_pre5BP_bm_H', 'mean_depth_post5BP_bm_H', 
        'mean_depth_pre5RI_bm_H', 'mean_depth_post5RI_bm_H']
    
    # reducción al primer conjunto de variables de interés
    resume_reduced = resume[resume.columns.intersection(features_interest)]

    y_list = create_Y_vector(composition_csv, resume_csv)
    y = pd.DataFrame(y_list)
    y.columns = ['labels']

    # Generación del df compuesto por todas las variables + la etiqueta y 
    df_labeled = pd.concat([resume, y], axis=1)

    return df_labeled

def create_2vectors_reads_virema(new_df, df_labeled):
    """
    Creamos el vector con los conteos del programa virema en el df_labeled
    
    """

    vector_direct_reads = list()
    vector_complementary_reads = list()
    
    # Extraemos los read_counts de virema que aparecen de forma directa
    for event in new_df['direct_event']:
        v = df_labeled.loc[df_labeled['cID_DI'] == event]['read_counts_virema'].values
        if len(v) > 0:
            vector_direct_reads.append(v[0])
        else:
            vector_direct_reads.append(0)

    # Extraemos los read_counts de virema que aparecen de forma complementaria
    for event in new_df['complementary_event']:
        v = df_labeled.loc[df_labeled['cID_DI'] == event]['read_counts_virema'].values
        if len(v) > 0:
            vector_complementary_reads.append(v[0])
        else:
            vector_complementary_reads.append(0)
    
    return vector_direct_reads, vector_complementary_reads
    
def add_virema_counts(df_mode, df_labeled):

    direct, complementary = create_2vectors_reads_virema(df_mode, df_labeled)

    direct = pd.DataFrame(direct, columns=['read_counts_direct_virema'])
    complementary = pd.DataFrame(complementary, columns=['read_counts_complement_virema'])

    df_new = pd.concat([df_mode, direct, complementary], axis=1)
    
    df_new['predicted_class_virema'] = df_new.apply(lambda row: \
                is_predicted(row['read_counts_direct_virema'], 
                            row['read_counts_complement_virema']), axis=1)

    return df_new

def is_predicted(direct, complement):
    """
    Retorna 1 o 0 en función de si el evento ha sido detecta en alguno de sus
    sentidos

    Args:
        df_mode (pd.df) Al menos: read_counts_direct_{program} y
                                  read_counts_complement_{program}
        program (str)   virema o ditector
    Return:
        is_predicted    (int)   0: no, 1: sí
    """
    if direct > 0 or complement > 0:
        predicted = 1
    else:
        predicted = 0

    return predicted


def create_2vectors_reads_ditector(new_df, df_labeled):
    """
    Creamos el vector con los conteos del programa virema en el df_labeled
    
    """

    vector_direct_reads = list()
    vector_complementary_reads = list()
    
    # Extraemos los read_counts de virema que aparecen de forma directa
    for event in new_df['direct_event']:
        v = df_labeled.loc[df_labeled['cID_DI'] == event]['read_counts_ditector'].values
        if len(v) > 0:
            vector_direct_reads.append(v[0])
        else:
            vector_direct_reads.append(0)

    # Extraemos los read_counts de virema que aparecen de forma complementaria
    for event in new_df['complementary_event']:
        v = df_labeled.loc[df_labeled['cID_DI'] == event]['read_counts_ditector'].values
        if len(v) > 0:
            vector_complementary_reads.append(v[0])
        else:
            vector_complementary_reads.append(0)
    
    return vector_direct_reads, vector_complementary_reads


def add_ditector_counts(df_mode, df_labeled):

    direct, complementary = create_2vectors_reads_ditector(df_mode, df_labeled)

    direct = pd.DataFrame(direct, columns=['read_counts_direct_ditector'])
    complementary = pd.DataFrame(complementary, columns=['read_counts_complement_ditector'])

    df_new = pd.concat([df_mode, direct, complementary], axis=1)
    # añadimos columna predicted_ditector (1/0)
    df_new['predicted_class_ditector'] = df_new.apply(lambda row: \
                is_predicted(row['read_counts_direct_ditector'], 
                            row['read_counts_complement_ditector']), axis=1)

    return df_new


def metabuscador_mode_df(df_composition, df_labeled):
    """
    """
    df = create_df_with_cIDDI_and_class(df_composition, df_labeled)
    meta_df = add_virema_counts(df, df_labeled)
    meta_df = add_ditector_counts(meta_df, df_labeled)

    return meta_df

def complete_metabuscador_mode(df_composition, df_labeled):
    """
    Genera la tabla de evaluación del modo metabuscador, de esta se podrán 
    hacer las subselecciones de los distintos modos (excepto del ML)
    """
    df_metabuscador = metabuscador_mode_df(df_composition, df_labeled)
    df_metabuscador['predicted_class_metabuscador'] = df_metabuscador.apply(
        lambda row : is_predicted(row['predicted_class_virema'], 
                                row['predicted_class_ditector']), axis=1)
    
    df_metabuscador['is_metabuscador'] = df_metabuscador.apply(lambda row:
           eval_category_tpfptn(row['real_class'], row['predicted_class_metabuscador']),
           axis=1)

    return df_metabuscador

def is_consensus(in_virema, in_ditector):
    """
    Devuelve 1 en el caso de que el evento haya sido detectado por los dos 
    programas. Si no, devuelve 0
    """
    in_virema = int(in_virema)
    in_ditector = int(in_ditector)
    
    if in_virema == 1 and in_ditector == 1:
        is_consensus = 1
    else:
        is_consensus = 0
    
    return is_consensus

def is_ml(direct_event, complementary_event, ml_list):
    """
    Retorna 1 o 0 en función de si el evento está dentro de la lista ml_list
    Función necesaria para generar la columna 'predicted_ML' en el DF general
    """
    if direct_event in ml_list or complementary_event in ml_list:
        is_ml = 1
    else: 
        is_ml = 0

    return is_ml

def complete_mode(df_composition, df_labeled, df_ml):
    """
    Genera la tabla de evaluación del modo metabuscador, de esta se podrán 
    hacer las subselecciones de los distintos modos (excepto del ML)
    """
    # lista con los dvgs filtrados con el modelo
    ml_list = list(df_ml['cID_DI'])

    df_metabuscador = metabuscador_mode_df(df_composition, df_labeled)
    df_metabuscador['predicted_class_metabuscador'] = df_metabuscador.apply(
        lambda row : is_predicted(row['predicted_class_virema'], 
                                row['predicted_class_ditector']), axis=1)
    
    df_metabuscador['predicted_class_consensus'] = df_metabuscador.apply(
        lambda row : is_consensus(row['predicted_class_virema'], 
                                row['predicted_class_ditector']), axis=1)
    df_metabuscador['predicted_class_ML'] = df_metabuscador.apply(
        lambda row : is_ml(row['direct_event'], row['complementary_event'],
                            ml_list), axis=1)                          

    # Añadimos el tipo de predicción que es teniendo en cuenta la clase real
    df_metabuscador['is_metabuscador'] = df_metabuscador.apply(lambda row:
           eval_category_tpfptn(row['real_class'], row['predicted_class_metabuscador']),
           axis=1)
    df_metabuscador['is_virema'] = df_metabuscador.apply(lambda row:
           eval_category_tpfptn(row['real_class'], row['predicted_class_virema']),
           axis=1)
    df_metabuscador['is_ditector'] = df_metabuscador.apply(lambda row:
           eval_category_tpfptn(row['real_class'], row['predicted_class_ditector']),
           axis=1)
    df_metabuscador['is_consensus'] = df_metabuscador.apply(lambda row:
           eval_category_tpfptn(row['real_class'], row['predicted_class_consensus']),
           axis=1)
    df_metabuscador['is_ML'] = df_metabuscador.apply(lambda row :
            eval_category_tpfptn(row['real_class'], row['predicted_class_ML']), 
            axis=1)

    return df_metabuscador

def eval_category_tpfptn(real_class, predicted_class):
    """
    Clasifica en TP, TN, FP o FN
    """
    real_class = int(real_class)
    predicted_class = int(predicted_class)

    if real_class == 1 and predicted_class == 1:
        contingency_class = 'TP'
    elif real_class == 1 and predicted_class == 0:
        contingency_class = 'FN'
    elif real_class == 0 and predicted_class == 1:
        contingency_class = 'FP'
    elif real_class == 0 and predicted_class == 0:
        contingency_class = 'TN'
    
    return contingency_class

def calc_sensitivity(TP, FN):
    if (TP + FN) > 0:
        return TP/(TP + FN)
    else:
        return 0

def calc_precision(TP, FP):
    if (TP + FP) > 0:
        return TP/(TP + FP)
    else:
        return 0

def f1(TP, FP, FN):
    ppv = calc_precision(TP, FP)
    tpr = calc_sensitivity(TP, FN)

    if (ppv+tpr) > 0:
        return (2*ppv*tpr)/(ppv+tpr)
    else:
        return 0


def count_TP_FP_FN(is_class_vector):
    """
    Retorna los valores
    """
    TP = FP = FN = TN = 0

    for c in is_class_vector:
        if c == 'TP':
            TP += 1
        elif c == 'FP':
            FP += 1
        elif c == 'FN':
            FN += 1
        elif c == 'TN':
            TN += 1
    contingency_dict = {'TP' : TP, 'FP' : FP, 'FN' : FN}
        
    return contingency_dict

def count_TP_TN_FP_FN(is_class_vector):
    """
    Retorna los valores TP, TN, FP, FN en un diccionario

    Returns:
        contingency_dict    (dict)  Valores TP, TN, FP y FN
    """
    TP = FP = FN = TN = 0

    for c in is_class_vector:
        if c == 'TP':
            TP += 1
        elif c == 'FP':
            FP += 1
        elif c == 'FN':
            FN += 1
        elif c == 'TN':
            TN += 1
    contingency_dict = {'TP' : TP, 'TN': TN, 'FP' : FP, 'FN' : FN}
        
    return contingency_dict
