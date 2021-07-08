#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
## 
## 
##                  DVGfinder: Defective Viral Genome-finder
###############################################################################
## PROGRAMA PRINCIPAL (controlador):
##       1.1. Metabúsqueda de DVGs con ViReMa-a (0.23) y DI-tector (v0.6)
##       1.2. Procesamiento de los alineamientos
##       1.3. Generación de una tabla conjunta con todos los DVGs detectados y 
##          extracción de variables
##       2.   Predicción con el modelo de Random Forest entrenado de los 
##            eventos reales
##       3.   Generación del informe HTML con las tablas a mostrar y las 
##            visualizaciones. 3 modos:
##          - COMPLETO: DVGs detectados por ViReMa-a + DI-tector
##          - CONSENSO: solos DVGs detectados por ambos 
##          - PREDICHO (ML): DVGs predichos como reales con el modelo de ML       
##                               
###############################################################################
## Author: Maria Jose Olmo-Uceda
## Version: 1.0
## Email: maolu@alumni.uv.es
## Date: 2021/07
###############################################################################

# system imports
import os
import time
import subprocess
from multiprocessing import Pool
import argparse

# third party imports
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import datapane as dp
import plotly.express as px

# DVGfinder functions
from Models import metabuscador, prediction, reports, visualization, headerName


def main():
    #print(os.listdir())
    t0 = time.time()

    # Lee los parámetros introducidos por el usuario
    fq, virus_ref, virus_ref_path, margin, n_processes = metabuscador.read_arguments()

    # Extraemos el nombre de la muestra
    sample_name = os.path.basename(fq).split(".")[-2]

    #ref = "NC_045512.2"
    ### nombre de la referencia sin la extensión
    virus_ref_no_extension = os.path.basename(virus_ref).split(".")[-2]

    # Lee el genoma de referencia del fichero fasta para extraer su longitud
    wt_sequence = SeqIO.read(virus_ref_path, "fasta")
    len_wt = len(wt_sequence)

    ## Fichero donde tenemos guardado el modelo
    # model_file = "rf_noTunning_rfe17.sav"
    model_file = "rf_tunning_x17.sav"

    ###########################################################################
    # PARTE 1: METABUSCADOR
    ###########################################################################
    
    # 1.1: IDENTIFICACIÓN DE LOS DVGS
    # -------------------------------------------------------------------------
    t1 = time.time() 
    
    # Con ViReMa
    metabuscador.run_virema_v023(fq, sample_name, virus_ref_no_extension, 
                            n_processes)
    t2 = time.time()
    # Con DItector
    metabuscador.run_ditector(fq, sample_name, virus_ref, n_processes)
    t3 = time.time()

    
    # 1.2: PROCESAMIENTO DE LOS ALINEAMIENTOS 
    # -------------------------------------------------------------------------
    # ALINEAMIENTOS
    ## Si no existe el directorio Outputs/alignments lo creamos
    alignment_directory = "Outputs/alignments"
    dir_outputs = os.listdir("Outputs/")
    if alignment_directory not in dir_outputs:
        os.system("mkdir {}".format(alignment_directory))

    ## Movemos los dos ficheros de alinemiento completos (bowtie y bwa mem)
    metabuscador.move_alignments(sample_name)

    ## Ordenamos los ficheros de alineamiento originales con `samtools sort`
    metabuscador.sort_complete_alignments(sample_name)

    ## Extraemos las reads que mapean bien (las consideraremos 'reads nativas')
    metabuscador.extract_mapped_reads(sample_name)

    ## generamos la selección de lecturas con CIGAR H (hard clipping)
    metabuscador.extract_H_reads_and_sort(sample_name)

    ## Escribimos los ficheros con la profundidad de cada posición
    metabuscador.write_depth_file(sample_name)

    ## Escribimos los ficheros con la información sobre los alineamientos
    metabuscador.generate_stats_alignment(sample_name)

    t4 = time.time()

    # 1.3: GENERACIÓN DE LA TABLA RESUMEN CON LAS MÉTRICAS
    # -------------------------------------------------------------------------

    # EXTRACCIÓN INFORMACIÓN RAW_OUTPUTS
    ## info base: BP, RI, read_counts_{program}, sense/DVG_type
    metabuscador.extract_recombination_events_virema(sample_name)
    metabuscador.extract_recombination_events_ditector(sample_name)
    ## Añadimos DVG_type o sense + identificadores ID_DI y cID_DI
    metabuscador.add_sense_and_IDs(sample_name) 
    metabuscador.add_dvgtype_and_IDs(sample_name)
        
    # Creamos una primera versión de la tabla unificada con la información 
    # extraída directamente de los raw_{program}
    metabuscador.create_unificate_dvg_table(sample_name)
    
    ## Borrar los ficheros con el output separado por programa
    os.system("rm Outputs/*_from_raw_*")
    
    ##Cargamos el dataframe en el que añadiremos las variables
    df, out_file = metabuscador.generate_df(sample_name, fq)
    
    ## Añadimos la información sobre profundidad de las coordenadas BP y RI
    ## de cada evento para los diferentes ficheros de alineamiento y escribimos
    df = metabuscador.parallelize_add_features(df, sample_name, fq, margin, 
                                        len_wt, 
                                        metabuscador.add_features, 
                                        n_processes)
    df.to_csv(out_file, index=False)

    t5 = time.time()

    ###########################################################################
    # PARTE 2: APLICACIÓN DEL MODELO DE ML A LA TABLA DE EVENTOS DETECTADOS
    ###########################################################################
    
    df_ML_show = prediction.filter_predicted_as_reals(df, model_file)

    t6 = time.time()
    ###########################################################################
    # PARTE 3: GENERACIÓN DEL INFORME FINAL
    ###########################################################################
    
    # 3.1: GENERACIÓN DEL DIRECTORIO FINALREPORTS
    # -------------------------------------------------------------------------
    # Si no existe el directorio FinalReports lo creamos
    dir_reports = 'FinalReports'
    if dir_reports not in os.listdir():
        os.mkdir(dir_reports)
    
    # Si no existe el directorio FinalReports/{sample_name} lo creamos
    report_directory = dir_reports + "/" + sample_name
    try:
        os.mkdir(report_directory)
    except OSError as error:
        print(error) 

    # 3.2: GENERACIÓN DE LAS TABLAS, LAS VISUALIZACIONES Y EL REPORT EN HTML
    # -------------------------------------------------------------------------
    # La siguente función:
    #       - crea las tablas ALL, CONSENSUS y MLfiltered, 
    #       - las escribe en el directorio /FinalReports/{sample} y 
    #       - genera las visualizaciones
    #       - crea el report en HTML que abre automáticamente

    reports.generate_report(df, df_ML_show, sample_name, len_wt)

    ###########################################################################
    # PARTE 4: LIMPIEZA DE DIRECTORIOS
    ###########################################################################
    # Borramos el directorio /Outputs/alignments
    os.system('rm -r Outputs/alignments') 

    # Movemos los resultados individuales de los programas a un nuevo 
    # directorio
    dir_olds = "OldOutputs"   
    dir_olds_virema = "OldOutputs/virema/"
    dir_olds_ditector = "OldOutputs/ditector"

    if dir_olds not in os.listdir():
        os.mkdir(dir_olds)
    if dir_olds_virema in os.listdir():
        os.mkdir(dir_olds_virema)
    if dir_olds_ditector in os.listdir():
        os.mkdir(dir_olds_ditector)

    os.system('mv Outputs/virema/* {}'.format(dir_olds_virema))
    os.system('mv Outputs/ditector/* {}'.format(dir_olds_ditector))
    os.system('mv Outputs/* {}'.format(dir_olds))
    # Borramos el directorio Outputs
    os.system('rm -r Outputs')

    tf = time.time()
    
    # Informamos sobre los tiempos de ejecución
    print(headerName.s5)
    print("-" * 79)
    print("- Virema DVG search: {} s".format(t2-t1))
    print("- Ditector DVG search: {} s".format(t3-t2))
    print("- Alignments inspection: {} s".format(t4-t3))
    print("- Metrics and table generation: {} s".format(t5-t4))
    print("- DVG filter with ML algorithm: {} s".format(t6-t5))
    print("- Visualizations and HTML report generation: {} s".format(tf-t6))
    print()
    print("- Total time: {} s".format(tf-t0))
    print("-" * 79)

    
if __name__=='__main__':
    main()  
