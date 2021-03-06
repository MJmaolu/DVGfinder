#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
## 
## 
##                  DVGfinder: Defective Viral Genome-finder
###############################################################################
## MAIN PROGRAM (controller):
##       1.1. DVG Metasearch with ViReMa-a (0.23) and DI-tector (v0.6)
##       1.2. Alignment process
##       1.3. Resume in a unique table of all the DVGs detected and generation
##          of informative features
##       2.   Prediction of the TP DVGs with the Gradient Boosting Classifier
##       3.   Generation of the final outputs, including the rendering of the
##           HTML report. Levels:
##          - ALL: All DVGs detected by ViReMa-a + DI-tector
##          - CONSENSUS: only DVGs detected for both programs. Sense sensitive 
##          - FILTERED (ML): DVGs predicted as reals by our ML model     
###############################################################################
## Author: Maria Jose Olmo-Uceda
## Version: 3.1
## - Take into account the polarity of the reference viral genome
## Email: mariajose.olmo@csic.es; maolu@alumni.uv.es
## Last update: 2022/04/12
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

    # Read user parameters
    fq, virus_ref, virus_ref_path, margin, threshold, n_processes, polarity = \
        metabuscador.read_arguments_v3()

    # Extract sample name
    sample_name = os.path.basename(fq).split(".")[-2]

    #ref = "NC_045512.2"
    ### Rerence name without extension
    virus_ref_no_extension = os.path.basename(virus_ref).split(".")[-2]

    # Lee el genoma de referencia del fichero fasta para extraer su longitud
    wt_sequence = SeqIO.read(virus_ref_path, "fasta")
    len_wt = len(wt_sequence)
    ref = wt_sequence.id # extraemos el identificador utilizado como referencia

    ## Fichero donde tenemos guardado el modelo
    # model_file = "rf_noTunning_rfe17.sav"
    model_file = "gbc_randomOpt.sav"

    ###########################################################################
    # PARTE 1: METABUSCADOR
    ###########################################################################
    
    # 1.1: IDENTIFICACI??N DE LOS DVGS
    # -------------------------------------------------------------------------
    t1 = time.time() 
    
    # Con ViReMa
    metabuscador.run_virema_v023(fq, sample_name, virus_ref_no_extension, 
                            n_processes)
    t2 = time.time()
    # Con DItector
    metabuscador.run_ditector(fq, sample_name, virus_ref_path, n_processes)
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

    ##??Ordenamos los ficheros de alineamiento originales con `samtools sort`
    metabuscador.sort_complete_alignments(sample_name)

    ##??Extraemos las reads que mapean bien (las consideraremos 'reads nativas')
    metabuscador.extract_mapped_reads(sample_name)

    ## generamos la selecci??n de lecturas con CIGAR H (hard clipping)
    metabuscador.extract_H_reads_and_sort(sample_name)

    ##??Escribimos los ficheros con la profundidad de cada posici??n
    metabuscador.write_depth_file(sample_name)

    ## Escribimos los ficheros con la informaci??n sobre los alineamientos
    metabuscador.generate_stats_alignment(sample_name)

    t4 = time.time()

    # 1.3: GENERACI??N DE LA TABLA RESUMEN CON LAS M??TRICAS
    # -------------------------------------------------------------------------

    #??EXTRACCI??N INFORMACI??N RAW_OUTPUTS
    ##??info base: BP, RI, read_counts_{program}, sense/DVG_type
    metabuscador.extract_recombination_events_virema(sample_name, ref)
    metabuscador.extract_recombination_events_ditector(sample_name)
    
    ##??Comprobamos qu?? programas han encontrado eventos --> booleano
    bool_vir, bool_dit = metabuscador.which_program_found(sample_name)
    
    #bool_vir = bool_dit = True 
    ## Definimos caso de uso en el que nos encontramos (4 posibles)
    case = metabuscador.define_case(bool_vir, bool_dit)
    print()
    print("~~~"*34)
    print(case)
    # Definimos variable que determinar?? si seguimos o no con el programa
    go = True

    #??Si no se encuentran DVGs (caso 4) se finaliza el flujo de trabajo
    if case == "case_4":
        go = False

    ## Add DVG_type or sense + identifiers ID_DI (BP_RI) y cID_DI (sense_BP_RI)
    while go:
        #??Si Virema encuentra DVGs
        if case in ["case_1", "case_2"]:
            metabuscador.add_sense_and_IDs(sample_name)
            print("ViReMa-a found DVGs in the sample")
        # Si DI-tector encuentra DVGs 
        if case in ["case_1", "case_3"]:
            metabuscador.add_dvgtype_and_IDs(sample_name)
            print("DI-tector found DVGs in the sample")
            
        # Creamos una primera versi??n de la tabla unificada con la informaci??n 
        # extra??da directamente de los raw_{program}
        #??Excluye cualquier DVG que est?? caracterizado fuera del rango del 
        #??genoma de referencia
        metabuscador.create_unificate_dvg_table(sample_name, bool_vir, bool_dit, len_wt)
        
        ##??Borrar los ficheros con el output separado por programa
        #os.system("rm Outputs/*_from_raw_*")
        
        ##Cargamos el dataframe en el que a??adiremos las variables
        df, out_file = metabuscador.generate_df(sample_name)
        
        ## A??adimos la informaci??n sobre profundidad de las coordenadas BP y RI
        ## de cada evento para los diferentes ficheros de alineamiento y escribimos
        #df = metabuscador.add_features(sample_name, margin, len_wt, bool_vir, bool_dit, df)
        df = metabuscador.parallelize_add_features(df, sample_name, margin, len_wt, 
                                            bool_vir, bool_dit,
                                            metabuscador.add_features, 
                                            n_processes)
        # Write resume file, in case virus has (-) polarity, invert the terms
        if polarity == 0:
            df_negative = metabuscador.invert_polarity(df)
            df_negative.to_csv(out_file, index=False)
        # if (+) write directly
        else:
            df.to_csv(out_file, index=False)

        t5 = time.time()

        ###########################################################################
        # PARTE 2: APLICACI??N DEL MODELO DE ML A LA TABLA DE EVENTOS DETECTADOS
        ###########################################################################
        #??Solo en el CASO 1 --> DVGs obtenidos con ambos programas
        if case == "case_1":
            df_ML_show = prediction.filter_predicted_as_reals_v3(df, model_file, threshold)
            print("~"*36)
            print("df_ML generated")
            # in case polarity (-)
            if polarity == 0:
                df_ML_show_negative = metabuscador.invert_polarity(df_ML_show)
                df_ML_show_negative.to_csv("{}_df_ML.csv".format(sample_name), 
                                            index=False)
            else:
                df_ML_show.to_csv("{}_df_ML.csv".format(sample_name), index=False)

        t6 = time.time()
        
        go = False
    ###########################################################################
    # PARTE 3: GENERACI??N DEL INFORME FINAL
    ###########################################################################
    
    # 3.1: GENERACI??N DEL DIRECTORIO FINALREPORTS
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
    
    
    # 3.2: TABLES GENERATION, VISUALIZATIONS AND HTML REPORT
    # -------------------------------------------------------------------------
    #??The report content depends on the usage case:
    # Case 1: report (levels ALL, CONSENSUS y MLfiltered)
    #       - crea las tablas ALL, CONSENSUS y MLfiltered, 
    #       - las escribe en el directorio /FinalReports/{sample} y 
    #       - genera las visualizaciones
    #       - crea el report en HTML que abre autom??ticamente
    # Caso 2 y 3: modo {programa}
    #       - crea la tabla PROGRAMA
    #       - la escribe en el directorio /FinalReports/{sample} y 
    #       - genera las visualizaciones
    #       - crea el report en HTML que abre autom??ticamente
    # Caso 4: Muestra el mensaje: "Ning??n DVG detectado en estos datos"
    if case == "case_1":
        #??if (-)sRNA virus
        if polarity == 0:
            reports.generate_report(df_negative, df_ML_show_negative, sample_name, len_wt)
        else:
            reports.generate_report(df, df_ML_show, sample_name, len_wt)
    elif case in ["case_2", "case_3"]:
        if polarity == 0:
            reports.generate_incomplete_report(df_negative, sample_name, len_wt, case)
        else:
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
    
    # Informamos sobre los tiempos de ejecuci??n
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
  
