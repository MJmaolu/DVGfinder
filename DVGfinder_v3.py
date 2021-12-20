#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
## 
## 
##                  DVGfinder: Defective Viral Genome-finder
###############################################################################
##  STRUCTURE:
##       1.1. DVGS metasearch with ViReMa-a (0.23) and DI-tector (v0.6)
##       1.2. Alignment processing
##       1.3. Generation of a unify table with all the detected DVGs and 
##          
##       2.   Prediction of the real events with a Gradient Boosting classifier
##            algorithm trained model
##       3.   Generation of an HTML report with interactive tables and plots 
##            3 modes:
##          - ALL: all the DVGs detected by ViReMa-a and DI-tector
##          - CONSENSUS: intersection 
##          - FILTERED: DVGs predicted as reals (p(TP) >= defined threshold)      
##                               
###############################################################################
## Author: Maria Jose Olmo-Uceda
## Version: 3.0     
## - Change in the prediction model --> Gradient Boosting Classifier: 
##                                      "gbc_randomOpt.sav"
## Email: mariajose.olmo@csic.es
## Date: 2021/12/01
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

    # User parameters
    fq, virus_ref, virus_ref_path, margin, threshold, n_processes = \
        metabuscador.read_arguments_v3()

    sample_name = os.path.basename(fq).split(".")[-2]

    #ref = "NC_045512.2"
    ### Reference name without extension
    virus_ref_no_extension = os.path.basename(virus_ref).split(".")[-2]

    # Genome length
    wt_sequence = SeqIO.read(virus_ref_path, "fasta")
    len_wt = len(wt_sequence)
    ref = wt_sequence.id # extraemos el identificador utilizado como referencia

    ## Predictive model file
    model_file = "gbc_randomOpt_train.sav"

    ###########################################################################
    # PARTE 1: METASEARCH
    ###########################################################################
    
    # 1.1: DVG SEARCH
    # -------------------------------------------------------------------------
    t1 = time.time() 
    
    # With ViReMa-a
    metabuscador.run_virema_v023(fq, sample_name, virus_ref_no_extension, 
                            n_processes)
    t2 = time.time()
    # With DItector
    metabuscador.run_ditector(fq, sample_name, virus_ref_path, n_processes)
    t3 = time.time()

    
    # 1.2: ALIGNMENT PROCESSING
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
    metabuscador.extract_recombination_events_virema(sample_name, ref)
    metabuscador.extract_recombination_events_ditector(sample_name)
    
    ## Comprobamos qué programas han encontrado eventos --> booleano
    bool_vir, bool_dit = metabuscador.which_program_found(sample_name)
    
    #bool_vir = bool_dit = True 
    ## Definimos caso de uso en el que nos encontramos (4 posibles)
    case = metabuscador.define_case(bool_vir, bool_dit)
    print()
    print("~~~"*34)
    print(case)
    # Definimos variable que determinará si seguimos o no con el programa
    go = True

    # Si no se encuentran DVGs (caso 4) se finaliza el flujo de trabajo
    if case == "case_4":
        go = False

    ## Añadimos DVG_type o sense + identificadores ID_DI y cID_DI
    while go:
        # Si Virema encuentra DVGs
        if case in ["case_1", "case_2"]:
            metabuscador.add_sense_and_IDs(sample_name)
            print("ViReMa-a found DVGs in the sample")
        # Si DI-tector encuentra DVGs 
        if case in ["case_1", "case_3"]:
            metabuscador.add_dvgtype_and_IDs(sample_name)
            print("DI-tector found DVGs in the sample")
            
        # Creamos una primera versión de la tabla unificada con la información 
        # extraída directamente de los raw_{program}
        # Excluye cualquier DVG que esté caracterizado fuera del rango del 
        # genoma de referencia
        metabuscador.create_unificate_dvg_table(sample_name, bool_vir, bool_dit, len_wt)
        
        ## Borrar los ficheros con el output separado por programa
        #os.system("rm Outputs/*_from_raw_*")
        
        ##Cargamos el dataframe en el que añadiremos las variables
        df, out_file = metabuscador.generate_df(sample_name)
        
        ## Añadimos la información sobre profundidad de las coordenadas BP y RI
        ## de cada evento para los diferentes ficheros de alineamiento y escribimos
        #df = metabuscador.add_features(sample_name, margin, len_wt, bool_vir, bool_dit, df)
        df = metabuscador.parallelize_add_features(df, sample_name, margin, len_wt, 
                                            bool_vir, bool_dit,
                                            metabuscador.add_features, 
                                            n_processes)
        df.to_csv(out_file, index=False)

        t5 = time.time()

        ###########################################################################
        # PARTE 2: APLICACIÓN DEL MODELO DE ML A LA TABLA DE EVENTOS DETECTADOS
        ###########################################################################
        # Solo en el CASO 1 --> DVGs obtenidos con ambos programas
        if case == "case_1":
            df_ML_show = prediction.filter_predicted_as_reals_v3(df, model_file, threshold)
            print("~"*36)
            print("df_ML generated")
            df_ML_show.to_csv("{}_df_ML.csv".format(sample_name), index=False)

        t6 = time.time()
        
        go = False
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
    # El contenido del Informe dependerá del caso de uso:
    # Caso 1: report completo (modos ALL, CONSENSUS y MLfiltered)
    #       - crea las tablas ALL, CONSENSUS y MLfiltered, 
    #       - las escribe en el directorio /FinalReports/{sample} y 
    #       - genera las visualizaciones
    #       - crea el report en HTML que abre automáticamente
    # Caso 2 y 3: modo {programa}
    #       - crea la tabla PROGRAMA
    #       - la escribe en el directorio /FinalReports/{sample} y 
    #       - genera las visualizaciones
    #       - crea el report en HTML que abre automáticamente
    # Caso 4: Muestra el mensaje: "Ningún DVG detectado en estos datos"
    if case == "case_1":
        reports.generate_report(df, df_ML_show, sample_name, len_wt)
    elif case in ["case_2", "case_3"]:
        reports.generate_incomplete_report(df, sample_name, len_wt, case)
    elif case == "case_4":
        print("")
        print("~~~ No DVGs detected in the data! ~~~")
        print("")
    
    ###########################################################################
    # PARTE 4: LIMPIEZA DE DIRECTORIOS
    ###########################################################################
    # Borramos el directorio /Outputs/alignments
    #os.system('rm -r Outputs/alignments') 

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
    #os.system('rm -r Outputs')

    tf = time.time()
    
    # Informamos sobre los tiempos de ejecución
    print(headerName.s5)
    print("-" * 79)
    """
    print("- Virema DVG search: {} s".format(t2-t1))
    print("- Ditector DVG search: {} s".format(t3-t2))
    print("- Alignments inspection: {} s".format(t4-t3))
    print("- DVGfinder's processing: {} s".format(tf-t4))
    print("- Metrics and table generation: {} s".format(t5-t4))
    print("- DVG filter with ML algorithm: {} s".format(t6-t5))
    print("- Visualizations and HTML report generation: {} s".format(tf-t6))
    """
    print()
    print("- Total time: {} s".format(tf-t0))
    print("-" * 79)

    
if __name__=='__main__':
    main()
  
