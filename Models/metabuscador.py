#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
## Functions for the METASEARCH module (v3) 
##                  DVGfinder: Defective Viral Genome-finder
###############################################################################
## Run the programs ViReMa-a (v0.23) and DI-tector (v0.6). Process the 
## alignments, generate the depth files and extra information about the 
## recombinations detected by the DVGs search algorithms.
## From this point, generate an unified table with all the DVGs detected by
## the programs and add informative features.
##
## The result of each search algorithm is stored on 'oldReports' when the 
## whole process is finished.
###############################################################################
## Author: Maria Jose Olmo-Uceda
## Version: 3.0
## Email: mariajose.olmo@csic.es
## Date: 2021/10/10
###############################################################################

# system imports
import os
import time
import subprocess
from multiprocessing import  Pool
import argparse

# third party imports
import pandas as pd
import numpy as np
import pyfastx

# READ USER PARAMETERS
# -----------------------------------------------------------------------------

def read_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument('-fq', '--fastq_file', 
            help=
            """Complete path to the fastq file we want to process.
            [Default is SD2_0h]""")
    parser.add_argument('-r', '--virus_reference', 
            help="""Genome reference in fasta format. Needs the indexed genome
            for bwa and for bowtie in the same directory.
            [Default is SARSCoV2_Wuhan1.fasta]""")
    parser.add_argument('-m', '--margin', 
            help=""" Number of positions to take in account for the average
            post and pre coordinate metrics. [Default is 5].""")
    parser.add_argument('-n', '--n_processes', 
            help="Number of processes. [Default is 1].")
    #parser.print_help()
    args = parser.parse_args()

    # Definimos las variables de entrada y salida
    ## fq file with single ends to process
    if args.fastq_file:
        fq = args.fastq_file
    else:
        fq = "SD2_0h.fq"

    ## Reference genome
    if args.virus_reference:
        virus_ref = args.virus_reference
        virus_ref_path = args.virus_reference
    else:
        virus_ref = "SARSCoV2WuhanHu1.fasta"
        virus_ref_path = "ExternalNeeds/references/" + virus_ref
    
    ## margin positions
    if args.margin:
        margin = int(args.margin)
    else:
        margin = 5
    
    ## number of processes
    if args.n_processes:
        n_processes = int(args.n_processes)
    else:
        n_processes = 1
    
    return fq, virus_ref, virus_ref_path, margin, n_processes

def read_arguments_v3():

    parser = argparse.ArgumentParser()
    parser.add_argument('-fq', '--fastq_file', 
            help=
            """Complete path to the fastq file we want to process.
            [Default is SD2_0h]""")
    parser.add_argument('-r', '--virus_reference', 
            help="""Genome reference in fasta format. Needs the indexed genome
            for bwa and for bowtie in the same directory.
            [Default is SARSCoV2_Wuhan1.fasta]""")
    parser.add_argument('-m', '--margin', 
            help=""" Number of positions to take in account for the average
            post and pre coordinate metrics. [Default is 5].""")
    parser.add_argument('-t', '--threshold',
            help="""Probability threshold (float). The prediction module will show only
            events with p(real) >= threshold. [Default is 0.5].""")
    parser.add_argument('-n', '--n_processes', 
            help="Number of processes. [Default is 1].")
    #parser.print_help()
    args = parser.parse_args()

    # Definimos las variables de entrada y salida
    ## fq file with single ends to process
    if args.fastq_file:
        fq = args.fastq_file
    else:
        fq = "SD2_0h.fq"

    ## Reference genome
    if args.virus_reference:
        virus_ref = args.virus_reference
        virus_ref_path = args.virus_reference
    else:
        virus_ref = "SARSCoV2WuhanHu1.fasta"
        virus_ref_path = "ExternalNeeds/references/" + virus_ref
    
    ## margin positions
    if args.margin:
        margin = int(args.margin)
    else:
        margin = 5
    
    ## threshold
    if args.threshold:
        threshold = float(args.threshold)
    else:
        threshold = 0.5
    
    ## number of processes
    if args.n_processes:
        n_processes = int(args.n_processes)
    else:
        n_processes = 1
    
    return fq, virus_ref, virus_ref_path, margin, threshold, n_processes

# 1.1: DVG SEARCH
# -----------------------------------------------------------------------------
def run_virema_v023(fq, sample_name, virus_ref_no_extension, n_processes):
    """
    Lanza ViReMa_v0.23 (funciona con python3).
    Genera el output propio de este programa en la carpeta: 
        /Output/virema/{sample}

    Args:
        fq  (str)   [in]    fastq sobre el que se realiza la identificación
                            de los DVGs
        sample_name (str)   [in]    Nombre de la muestra. Será igual al nombre
                                    del fq sin la extensión.
        virus_ref_no_extension (str)   [in]  Genoma de referencia. ViReMa usa
                             bowtie por defecto por lo que el nombre de la 
                            referencia será solo el basename. 
                            Ej: SARSCoV2WuhanHu1
                            Si no encuentra el genoma indexado en el path
                            usará Bowtie-Build para generarlo.
        n_processes (int)   [in]   Number of process  
    """
    # Creamos las variables necesarias
    ## Guardamos el nombre de la muestra
    sample_name = os.path.basename(fq).split(".")[-2]
    ## Rutas para el programa y la referencia
    virema = "ExternalNeeds/thirdPrograms/ViReMa_0.23/ViReMa.py"
    ref = "ExternalNeeds/references/" + virus_ref_no_extension
    out_dir_virema = "Outputs/virema/" + sample_name
    sam_bt = sample_name + "_bt2.sam"   

    # Lanzamos el proceso de búsqueda
    subprocess.run(['python3', virema, '--Output_Tag', sample_name, 
                    '--Output_Dir', out_dir_virema, '--P', str(n_processes), 
                    ref, fq, sam_bt])

def run_ditector(fq, sample_name, virus_ref_path, n_processes):
    """
    Lanza ViReMa_v0.23 (funciona con python3).
    Genera el output propio de este programa en la carpeta: 
        /Output/virema/{sample}

    Args:
        fq  (str)   [in]    fastq sobre el que se realiza la identificación
                            de los DVGs
        sample_name (str)   [in]    Nombre de la muestra. Será igual al nombre
                                    del fq sin la extensión.
        virus_ref (str)   [in]  Genoma de referencia. DItector usa bwa mem por 
                                lo que le damos el nombre del fichero en formato
                                fasta y deberá estár indexado en el mismo path
                                Ej SARSCoV2WuhanHu1.fasta
        n_processes (int)   [in]   Number of process  
    """

    # Creamos las variables necesarias
    ## Guardamos el nombre de la muestra
    #sample_name = os.path.basename(fq).split(".")[-2]
    ## Rutas para el programa y la referencia
    ditector = "ExternalNeeds/thirdPrograms/DI-tector_06.py"
    ref = virus_ref_path
    out_dir_ditector = "Outputs/ditector/" + sample_name 
    #sam_bm = sample_name + "_bm.sam"

    # Lanzamos proceso de búsqueda guardando alineamientos intermedios
    subprocess.run(['python3', ditector, '-n', '1', '-o', out_dir_ditector,
                    '-t', sample_name, '-k', '-x', str(n_processes), ref, fq])
    
    # Borramos todos los ficheros temporales 
    ## Lecturas que no han alineado o lo han hecho con H|S
    wo_virus = out_dir_ditector + "/" + sample_name + "_temp_file_woVirus.sam"
    ## Fichero con las reads no mapeadas segmentadas
    seqment = out_dir_ditector + "/" + sample_name + "_temp_seqment.fq"
    ## alineamiento con bwa aln del fichero seqment (.fq)
    aln = out_dir_ditector + "/" + sample_name + "_temp_aln.sam"
    ## segmentos alineados (.txt)
    ali = out_dir_ditector + "/" + sample_name + "_Ali.txt"
    ## lecturas que mapean (samtools view -F4)
    well_aln = out_dir_ditector + "/" + sample_name + "_temp_file_Virus.bam"
    ## profundidad de todas las posiciones del well_aln
    depth_well_aln = out_dir_ditector + "/" + sample_name + "_temp_file_Virus_Coverage.txt"
    
    subprocess.run(['rm', wo_virus, seqment, aln, ali, well_aln, depth_well_aln])


# 1.2: EXTRACTION OF DATA AND METRICS GENERATION
# -----------------------------------------------------------------------------

def move_alignments(sample_name):
    """
    Mueve los alineamientos {sample}_bt2.sam y {sample}_temp_file_onVirus.sam
    a una nueva carpeta /Outputs/alignments/sample. Primero comprueba si la 
    carpeta existe y si no la crea
    Args:
        sample_name (str)   [in]    Nombre de la muestra
    """
    list_dir_alignments = os.listdir("Outputs/alignments")
    dir_sample_al = "Outputs/alignments" + "/" + sample_name + "/"
   
    if sample_name not in list_dir_alignments:
        os.system("mkdir {}".format(dir_sample_al))
    

    bt2_alignment = "Outputs/virema/" + sample_name + "/" + sample_name \
        + "_bt2.sam"
    bm_alignment = "Outputs/ditector/" + sample_name + "/" + sample_name \
        + "_temp_file_onVirus.sam"
    bm_oldname = "Outputs/alignments/" + sample_name + "/" + sample_name \
        + "_temp_file_onVirus.sam"
    bm_newname = "Outputs/alignments/" + sample_name + "/" + sample_name \
        + "_bm.sam"
    
    # movemos al directorio especifico
    subprocess.run(["mv", bt2_alignment, dir_sample_al])
    subprocess.run(["mv", bm_alignment, dir_sample_al])
    # cambio nombre del bm
    subprocess.run(["mv", bm_oldname, bm_newname])

def sort_complete_alignments(sample_name):
    """
    Ordena (`samtools sort`) los alineamientos con el nombre de la muestra 
    sample_name que hay en el directorio Outputs/alignments.
    Args:
        sample_name (str)   [in]    Nombre de la muestra
    
    Returns:
        sorted_alignment    (bam file) [out] Alineamiento binario ordenado 
    """
    dir_alignments = "Outputs/alignments/" + sample_name + "/"
    alignments = ['_bm.sam', '_bt2.sam']
    sorted_tails = ['_bm_sorted.bam', '_bt2_sorted.bam']

    current_files = os.listdir(dir_alignments)

    for i in range(len(alignments)):
        alignment = alignments[i] 
        tail = sorted_tails[i]      
        unsorted_sam = dir_alignments + sample_name + alignment
        sorted_bam = dir_alignments + sample_name + tail
        # si el alineamiento ya está ordenado avisamos por pantalla
        if sorted_bam in current_files:
            print("{} is already sorted".format(sample_name + tail))
        else:
            os.system("samtools sort {} > {}".format(unsorted_sam, sorted_bam))
        

def extract_mapped_reads(sample_name):
    """
    Usa samtools view para extraer las reads que mapean contra la referencia
    Excribe la selección en un nuevo fichero de alineamiento

    Args:
        sample_name (str)   [in]    Nombre de la muestra
    Returns:
        out_file    (file.sam)  [out]
    """
    dir_alignments = "Outputs/alignments/" + sample_name + "/"
    sorted_alignments = ['_bm_sorted.bam', '_bt2_sorted.bam']
    mapped_tails = ['_bm_sorted_mapped.bam', '_bt2_sorted_mapped.bam']
    
    current_files = os.listdir(dir_alignments)

    for i in range(len(sorted_alignments)):
        sorted_alignment = sorted_alignments[i] 
        tail = mapped_tails[i]     

        sorted_bam = dir_alignments + sample_name + sorted_alignment
        mapped_bam = dir_alignments + sample_name + tail
        # si el alineamiento ya está ordenado avisamos por pantalla
        if mapped_bam in current_files:
            print("{} already exist".format(sample_name + tail))
        else:
            os.system("samtools view -Sbh -F4 {} > {}".format(sorted_bam, \
                mapped_bam))


def extract_H_reads_and_sort(sample_name):
    """
    Extrae las reads que contienen H en su CIGAR. Para ello llama al script 
    'extract_H_reads.sh' que toma como argumentos cada alineamiento completo
    de los existentes en el directorio /Outputs/alignments y el nombre que
    se le da al fichero resultante del filtro. 
    Args:
        sample_name (str)   [in]    Nombre de la muestra

    Returns:
        sample_name_{alignment}_H.sam (file)  Fichero con las reads que 
                                    contienen H en su CIGAR
        sample_name_{alignment}_sorted_H.bam  (file)  Fichero anterior en 
                                    formato binario y ordenado
    """
    dir_alignments = "Outputs/alignments/" + sample_name + "/"
    alignments = ['_bm.sam', '_bt2.sam']
    H_tails = ['_bm_H.sam', '_bt2_H.sam']
    sorted_H_tails = ['_bm_sorted_H.bam', '_bt2_sorted_H.bam']

    # Para cada alineamiento completo
    for i in range(len(alignments)):
        alignment = alignments[i] 
        tail = H_tails[i]      
        sorted_tail = sorted_H_tails[i]

        complete_alignment = dir_alignments + sample_name + alignment
        reads_H = dir_alignments + sample_name + tail
        sorted_H = dir_alignments + sample_name + sorted_tail
        
        # extraemos las reads con CIGAR H
        os.system("./Models/extract_H_reads.sh {} {}".format\
            (complete_alignment, reads_H))
        os.system("samtools sort {} > {}".format(reads_H, sorted_H))

def write_depth_file(sample_name):
    """
    Para los tres tipos de alineamientos (completo, mapped y H) de cada
    uno de los algoritmos usados (bm y bt2), escribe los ficheros que contienen
    la profundidad de cada una de las posiciones del genoma.

    Args:
        sample_name (str)   [in]    Nombre de la muestra

    Returns:
        sample_name_{alignment}_depth.txt   (file)  Profundidad de todas las 
                                posiciones del alineamiento completo
        sample_name_{alignment}_sorted_mapped_depth.txt (file)  Profundidad
                                del alineamiento con las lecturas que podemos
                                considerara nativas
        sample_name_{alignment}_sorted_H_depth.txt  (file)  Profundidad del
                                alineamiento con las lecturas que que contienen
                                algún evento de recombinación 
    """
    dir_alignments = "Outputs/alignments/" + sample_name + "/"
    #'_bm_sorted.bam', '_bt2_sorted.bam', 
    alignments = ['_bm_sorted_mapped.bam', '_bt2_sorted_mapped.bam', 
                '_bm_sorted_H.bam', '_bt2_sorted_H.bam']
    tail = "_depth.txt"

    for alignment in alignments:
        alignment_file = dir_alignments + sample_name + alignment 
        base_name = os.path.basename(alignment_file).split(".")[-2]
        depth_file = dir_alignments + base_name + tail

        os.system("samtools depth -a {} > {}".format(alignment_file, depth_file)) 

def extract_recombination_events_virema(sample_name, ref):
    """
    Genera un fichero tsv con la información: BP, RI, read_counts_virema, sense
    Para ello parsea el fichero generado por ViReMa 
                '{sample_name}_Virus_Recombination_Results.txt'
    Args: 
        sample_name (str)   [in]    Nombre de la muestra
        ref (str)   [in]    Identificador del genoma de referencia. Extraemos
                            del id del fasta de referencia

    Returns:
        {sample_name}_from_raw_virema   (file)  tsv con la información de cada 
                            evento de recombinación detectado por ViReMa.
                            BP, RI, read_counts_virema, sense
                            
    """
    dir_output_virema = "Outputs/virema/"
    recombinations_file = dir_output_virema + sample_name + "/" + sample_name \
        + "_Virus_Recombination_Results.txt"
    
    # extracción de la información
    os.system("./Models/extract_recombination_events_virema.sh {} {} {}".format\
        (sample_name, ref, recombinations_file))

def extract_recombination_events_ditector(sample_name):
    """
    Genera un fichero tsv con la información: BP, RI, read_counts_virema, DVG_type
    Para ello parsea el fichero generado por DItector 
                '{sample_name}_counts.txt'
    Args: 
        sample_name (str)   [in]    Nombre de la muestra
        
    Returns:
        {sample_name}_from_raw_ditector   (file)  tsv con la información de cada 
                            evento de recombinación detectado por DItector.
                            BP, RI, read_counts_ditector, DVG_type
    """
    dir_output_ditector = "Outputs/ditector/"
    recombinations_file = dir_output_ditector + sample_name + "/" + sample_name \
        + "_counts.txt"
        
    # extracción de la información
    os.system("Models/extract_recombination_events_ditector.sh {} {}".format\
        (sample_name, recombinations_file))

def which_program_found(sample_name):
    """
    Interroga los ficheros donde se extraen los eventos detectados por cada
    programa para establecer con los resultados de cuál continuar con el flujo
    de trabajo

    Args: 
        sample_name (str)
    Return:
        bool_vir    (bool)  ¿Virema ha encontrado DVGs? True/False
        bool_dit    (bool)  ¿Ditector ha encontrado DVGs? True/False
    """

    vir = "Outputs/{}_from_raw_virema.tsv".format(sample_name)
    dit = "Outputs/{}_from_raw_ditector.tsv".format(sample_name)

    if os.path.getsize(vir) == 0:
        bool_vir = False
    else:
        bool_vir = True
    
    if os.path.getsize(dit) == 0:
        bool_dit = False
    else:
        bool_dit = True

    return bool_vir, bool_dit

def define_case(bool_vir=bool, bool_dit=bool):
    """
    Define el caso de uso en el que está la búsqueda. Existen 4 posibilidades
    en función de qué programas hayan encontrado DVGs en la muestra:
        case_1 = Virema-a + DI-tector
        case_2 = Virema-a
        case_3 = DI-tetor
        case_4 = None
    En función del caso de uso, el programa determinará el consecuente flujo de 
    trabajo.

    Args:
        bool_vir    (bool)  ¿Virema ha encontrado DVGs? True/False
        bool_dit    (bool)  ¿Ditector ha encontrado DVGs? True/False

    Return:
        case    (str)   Caso de uso en el que se encuentra el programa
    """
    if bool_vir:
        if bool_dit:
            case = "case_1"
        else: 
            case = "case_2"
    else:
        if bool_dit:
            case = "case_3"
        else:
            case = "case_4"

    return case

def asign_dvg_type(sense, BP, RI):
    """
    Clasifica el dvg en Deleción, Inserción o copy/snap back con su dirección
    en función de los parámetros sense (++, --, +-, -+) y las coordenadas 
    BP (breakpoint) y RI (reinitiation site) 
    El criterio es:
        sense       BP ~ RI     Clasificación       DVG_type   
        ---------------------------------------------------- 
        ++          <           Deletion Forward    Deletion_forward
        ++          >=          Insertion Forward   Insertion_forward
        --          >           Deletion Reverse    Deletion_reverse
        --          <=          Insertion Reverse   Insertion_reverse
        +-          ><=         5' cb/sb            5cb/sb
        -+          ><=         3' cb/sb            3cb/sb          

    Args:

        sense (str) [in]    Strand sense of each segment (preBP & postRI)
        BP  (int)   [in]   Breakpoint: Última posición que mapea con la referencia
                            en el 1er fragmento del dvg
        RI  (int)   [in]   Reinitiation site: primera posición que mapea con
                            la referencia en el 2º fragmento del dvg 

    Returns:

        dvg_type    (str)   [out]   tipo de dvg asignado en función de los 
                                    criterios definidos
    """
    # Nos aseguramos de que los argumentos son del tipo que necesitamos
    BP = int(BP)
    RI = int(RI)
    sense = str(sense)

    if sense not in ('++', '--', '+-', '-+'):
        print("Sense must be: '++', '--', '+-' or '-+'")
        dvg_type = "None"
        ## esto es feo, añadir TypeError
    else:
        # Asignamos los dvg_type en función del criterio detallado
        if sense == '++':
            if BP < RI:
                dvg_type = 'Deletion_forward'
            else:
                dvg_type = 'Insertion_forward'
        elif sense == '--':
            if BP > RI:
                dvg_type = 'Deletion_reverse'
            else:
                dvg_type = 'Insertion_reverse'
        elif sense == '+-':
            dvg_type = '5cb/sb'
        else:
            dvg_type = '3cb/sb'
    
        return dvg_type

def asign_sense(dvg_type):
    """
    Indica el sentido de los dos fragmentos de la lectura en función del tipo
    de dvg en el que haya sido clasificado. 
    Importante: DI-tector podría no haber asignado bien el tipo.
    El criterio es:
        DVG_type            sense       
        ----------------------------
        Deletion_forward     ++        
        Insertion_forward    ++        
        Deletion_reverse     --        
        Insertion_reverse    --            
        5cb/sb               +-        
        3cb/sb               -+          

    Args:
        dvg_type    (str)   [in]   Tipo de dvg con el que DItector ha 
                                clasificado el evento de recombinación
    
    Returns:
        sense   (str)   [out]   Strand sense de los dos fragmentos. El orden  
                                es siempre 5' -> 3'
    """

    dvg_type = str(dvg_type)
    # Diccionario con las conversiones sense --> dvg_type
    dvg_types = {"++" : ('Deletion_forward', 'Insertion_forward'), 
    '--' : ('Deletion_reverse', 'Insertion_reverse'), 
    '+-' : '5cb/sb', '-+' : '3cb/sb'}

    sense = "dvg_type unrecognized"
    for key, value in dvg_types.items():
        if dvg_type in value:
            sense = key
    
    return sense

def create_IDDI(BP, RI):
    """
    Genera el ID_DI corto como BP_RI

    Args:
        BP  (int)   [in]   Coordenada BP
        RI  (int)   [in]   Coordenada RI
    
    Returns:
        ID_DI   (str)   [out]   Identificador corto = BP_RI
    """
    ID_DI = str(BP) + "_" + str(RI)

    return ID_DI

def create_cIDDI(sense, ID_DI):
    """
    Genera el cID_DI

    Args:
        sense   (str)   [in]    Strand sense of the first and second part of 
                                the DVG
        ID_DI   (str)   [out]   Identificador corto = BP_RI
    Return:
        cID_DI  (str)   [out]   Identificador completo = sense_BP_RI
    """
    cDI_ID = str(sense) + "_" + ID_DI

    return cDI_ID

def add_sense_and_IDs(sample_name):
    """
    Para el output de ViReMa-a
    Lee el fichero {sample_name}_from_raw_virema.tsv, lo convierte en un 
    pd.DataFrame, añade el nombre de las columnas existente y le añade:
     'DVG_type', 'ID_DI' y 'cID_DI'. Por último escribe el df en un nuevo 
    fichero: {sample}_from_raw_virema_wID.csv

    Args:
        sample_name (str)   [in]    Nombre de la muestra
    """

    input_file = "Outputs/" + sample_name + "_from_raw_virema.tsv"
    output_file = "Outputs/" + sample_name + "_from_raw_virema_wIDs.tsv"

    df = pd.read_csv(input_file, sep="\t", header=None, low_memory=False)
    # Nombre de las columnas existentes
    df.columns = ['BP', 'RI', 'read_counts_virema', 'sense']
    df[['DVG_type']] = df.apply(lambda row :
                    asign_dvg_type(row['sense'], row['BP'], row['RI']), axis=1)
    df[['ID_DI']] = df.apply(lambda row : 
                    create_IDDI(row['BP'], row['RI']), axis=1)
    df[['cID_DI']] = df.apply(lambda row : 
                    create_cIDDI(row['sense'], row['ID_DI']), axis=1)

    df.to_csv(output_file, sep="\t", index=False)

def add_dvgtype_and_IDs(sample_name):
    """
    Para el output de DI-tector
    Lee el fichero {sample_name}_from_raw_ditector.tsv, lo convierte en un 
    pd.DataFrame, añade el nombre de las columnas existente y le añade:
     'sense', 'ID_DI' y 'cID_DI'. Por último escribe el df en un nuevo 
    fichero: {sample}_from_raw_ditector_wID.csv

    Args:
        sample_name (str)   [in]    Nombre de la muestra
    """

    input_file = "Outputs/" + sample_name + "_from_raw_ditector.tsv"
    output_file = "Outputs/" + sample_name + "_from_raw_ditector_wIDs.tsv"

    df = pd.read_csv(input_file, sep="\t", header=None, low_memory=False)
    # Nombre de las columnas existentes
    df.columns = ['DVG_type', 'BP', 'RI', 'read_counts_ditector']
    df[['sense']] = df.apply(lambda row :
                    asign_sense(row['DVG_type']), axis=1)
    df[['ID_DI']] = df.apply(lambda row : 
                    create_IDDI(row['BP'], row['RI']), axis=1)
    df[['cID_DI']] = df.apply(lambda row : 
                    create_cIDDI(row['sense'], row['ID_DI']), axis=1)

    df.to_csv(output_file, sep="\t", index=False)


def lists_dvgs(sample_name, bool_vir, bool_dit):
    """
    Extrae todos los dvgs detectados en la muestra (sample_name). 
    Parsea los

    Args:
        sample_name (str)   [in]    Nombre de la muestra
        bool_vir    (bool)  ¿Virema ha encontrado DVGs? True/False
        bool_dit    (bool)  ¿Ditector ha encontrado DVGs? True/False    
    
    Returns:
        dvg_list    (list)  [out]   Lista con el total de cID_DI no repetidos
                                 detectados con todos los programas
        dvg_list_vir    (list)  [out]   Lista con los cID_DI de ViReMa
        dvg_list_dit    (list)  [out]   Lista con los cID_DI de DItector
    """

    infile_vir ="Outputs/" + sample_name + "_from_raw_virema_wIDs.tsv"
    infile_dit = "Outputs/" + sample_name + "_from_raw_ditector_wIDs.tsv"

    dvg_list_vir = list()
    dvg_list_dit = list()

    if bool_vir:
        df1 = pd.read_csv(infile_vir, header=0, sep="\t")
        dvg_list_vir = list(df1['cID_DI'])
    if bool_dit:
        df2 = pd.read_csv(infile_dit, header=0, sep="\t")
        dvg_list_dit = list(df2['cID_DI'])
        
    # Todos los DVGs, sin duplicaciones 
    dvg_list = sorted(set(dvg_list_vir) | set(dvg_list_dit))

    return dvg_list, dvg_list_vir, dvg_list_dit

# 1.3: GENERATION OF RESUME TABLE 
# -----------------------------------------------------------------------------

def create_unificate_dvg_table(sample_name, bool_vir, bool_dit, len_wt):
    """
    Generamos el df unificado con la información de salida de los programas
    """

    output_file = "Outputs/" + sample_name + "_unificated_table.csv"

    # generamos la lista completa de dvgs
    dvg_list, dvg_list_vir, dvg_list_dit = lists_dvgs(sample_name, 
                                                    bool_vir, bool_dit)
    #print(dvg_list, dvg_list_vir, dvg_list_dit)
    
    df = pd.DataFrame()
    df[['cID_DI']] = dvg_list
    df[['sense', 'BP', 'RI']] = df.cID_DI.str.split("_", expand=True,)
    # convertimos BP y RI a tipo integer
    df = df.astype({'BP' : int, 'RI' : int}) 
    df[['DVG_type']] = df.apply(lambda row : 
                             asign_dvg_type(row['sense'], row['BP'], 
                             row['RI']), axis=1)
    if bool_vir:
        df[['read_counts_virema']] = df.apply(lambda row : 
                            extract_counts_virema(row['cID_DI'], sample_name, 
                            dvg_list_vir), axis=1)
    if bool_dit:
        df[['read_counts_ditector']] = df.apply(lambda row : 
                            extract_counts_ditector(row['cID_DI'], sample_name,
                             dvg_list_dit), axis=1)
    # Antes de escribir revisamos que todos los DVGs tengan un rango de coordenadas
    # adecuado. Es decir, que BP o RI no puedan ser > len_wt ni < 1 (este tipo
    # de errores puede ser cometido por DI-tector. Los eliminamos del df.
    df_possible = df.loc[~((df.BP > len_wt) | (df.RI > len_wt) | (df.BP < 1) | (df.RI < 1))]
    df_possible.to_csv(output_file, index=False)
    
    #return df_possible

#–---------------------    
# COUNTS NUMBER
#–---------------------    
def extract_counts_virema(cID_DI, sample_name, dvg_list_vir):
    """
    Comprueba si el cID_DI ha sido detectado por ViReMa, si es así extrae del 
    fichero {sample_name}_from_raw_virema_wIDs.tsv el número, si no devuelve 0
    """

    infile ="Outputs/" + sample_name + "_from_raw_virema_wIDs.tsv"
    df = pd.read_csv(infile, header=0, sep="\t")


    if cID_DI in dvg_list_vir:
        read_count = df.loc[df['cID_DI'] == cID_DI,'read_counts_virema'].values[0]
    else:
        read_count = 0

    return read_count

def extract_counts_ditector(cID_DI, sample_name, dvg_list_dit):
    """
    Comprueba si el cID_DI ha sido detectado por DItector, si es así extrae del 
    fichero {sample_name}_from_raw_ditector_wIDs.tsv el valor, si no devuelve 0
    """

    infile = "Outputs/" + sample_name + "_from_raw_ditector_wIDs.tsv"
    df = pd.read_csv(infile, header=0, sep="\t")

    if cID_DI in dvg_list_dit:
        read_count = df.loc[df['cID_DI'] == cID_DI,'read_counts_ditector'].values[0]
    else:
        read_count = 0

    return read_count

#–---------------------    
# DEPTHS
#–--------------------- 
def extract_depth_from_coordinate(coordinate, sample_name, depth_txt, len_wt):
    """
    Obtenemos la profundidad de la posición (base 1) establecida como coordenada.
    Se tiene en cuenta que el fichero pueda estar vacío, en este caso devolverá
    0 en cualquier coordenada. Esto puedde pasar principalmente con el filtro
    de solo reads con CIGAR H.

    Args:
        coordinate  (int)   [in]    Posición (base 1) de la que queremos 
                                    extraer la profundidad 
        depth_txt   (str)   [in]    Fichero con la profundidad de todas las
                                    posiciones
        
    """
    depth_file = "Outputs/alignments/" + sample_name + "/" + depth_txt

    if (os.path.getsize(depth_file) != 0 and 1 < coordinate < len_wt):
        header_names =  ("chr", "position", "depth")
        depth = pd.read_csv(depth_file, sep="\t", header=None, \
            names=header_names)
        depth_value = depth.loc[depth['position'] == \
            coordinate, 'depth'].values[0]
    else:
        depth_value = 0

    return depth_value

def write_depth_coordinate(sample_name, df, len_wt):
    """
    Añade las columnas depth_BP y depth_RI de:
    los alineamientos mapped y soloH generados con cada algoritmo (bt2 y bwa mem)

    Args:
        sample_name (str)   [in]    Nombre de la muestra  
        df  (pd.DataFrame)  [in]    Tabla unificada cargada en memoria
    Returns:
        df  (pd.Dataframe)  [out]   Tabla con todos los eventos DVGs detectados
                                    más la información sobre la profundidad
                                    de cada coordenada de cada uno de los 
                                    ficheros de alineamiento    
    """

    #bt2 = sample_name + "_bt2_sorted_depth.txt"
    #bm = sample_name + "_bm_sorted_depth.txt"
    bt2_mapped = sample_name + "_bt2_sorted_mapped_depth.txt"
    bm_mapped = sample_name + "_bm_sorted_mapped_depth.txt"
    bt2_H = sample_name + "_bt2_sorted_H_depth.txt"
    bm_H = sample_name + "_bm_sorted_H_depth.txt"

    #df = pd.read_csv("Outputs/unificated_table.csv", header=0)
    
    # Añadimos columnas con información sobre la profundidad
    """ TENDRÁ SENTIDO CUANDO LAS MUESTRAS ESTÉN CONTAMINADAS CON GENOMAS EXÓGENOS
    df[['depthBP_bt2']] = df.apply(lambda row : 
                        extract_depth_from_coordinate(row['BP'], sample_name,
                                                     bt2), axis=1)
    df[['depthRI_bt2']] = df.apply(lambda row : 
                        extract_depth_from_coordinate(row['RI'], sample_name,
                                                     bt2), axis=1)
    """
    df[['depthBP_bt2_mapped']] = df.apply(lambda row : 
                        extract_depth_from_coordinate(row['BP'],  sample_name,
                                                    bt2_mapped, len_wt), axis=1)
    df[['depthRI_bt2_mapped']] = df.apply(lambda row : 
                        extract_depth_from_coordinate(row['RI'],  sample_name,
                                                    bt2_mapped, len_wt), axis=1)
    df[['depthBP_bt2_H']] = df.apply(lambda row : 
                        extract_depth_from_coordinate(row['BP'],  sample_name,
                                                    bt2_H, len_wt), axis=1)
    df[['depthRI_bt2_H']] = df.apply(lambda row : 
                        extract_depth_from_coordinate(row['RI'],  sample_name,
                                                    bt2_H, len_wt), axis=1)
    """
    df[['depthBP_bm']] = df.apply(lambda row : 
                        extract_depth_from_coordinate(row['BP'],  sample_name,
                                                    bm), axis=1)                   
    df[['depthRI_bm']] = df.apply(lambda row : 
                        extract_depth_from_coordinate(row['RI'],  sample_name,
                                                    bm), axis=1)
    """
    df[['depthBP_bm_mapped']] = df.apply(lambda row : 
                        extract_depth_from_coordinate(row['BP'],  sample_name,
                                                    bm_mapped, len_wt), axis=1)
    df[['depthRI_bm_mapped']] = df.apply(lambda row : 
                        extract_depth_from_coordinate(row['RI'],  sample_name,
                                                    bm_mapped, len_wt), axis=1)
    df[['depthBP_bm_H']] = df.apply(lambda row : 
                        extract_depth_from_coordinate(row['BP'],  sample_name,
                                                    bm_H, len_wt), axis=1)
    df[['depthRI_bm_H']] = df.apply(lambda row :  
                        extract_depth_from_coordinate(row['RI'],  sample_name,
                                                    bm_H, len_wt), axis=1)
    # write
    #df.to_csv("Outputs/unificate__temp_table_wDepths.csv", index=False)

    return df


def mean_depth_neighborhood_pre_coordinate(coordinate, sample_name, depth_txt,
                                            margin):
    """
    Calcula la profundidad media de las 'margin' posiciones antes y después 
    de la coordenada
    Args:
        coordinate  (int)   [in]    Posición (base 1) de la que vamos a evaluar
                                    la profundidad media de 
        depth_txt   (str)   [in]    Fichero con la profundidad de todas las
                                    posiciones
        margin      (int)   [in]    Número de posiciones de distancia con 
                                    respecto a la coordenada con las que se 
                                    quiere calcular la media de profundidad
    
    Returns:
        meanPre (float) [out]   Profundidad media de las posiciones 'margin'
                                anteriores a la coordenada
        meanPost (float) [out]   Profundidad media de las posiciones 'margin'
                                posteriores a la coordenada
    """
    
    # cargamos el fichero con la profundidad como pd.df

    depth_file = "Outputs/alignments/" + sample_name + "/" + depth_txt
    header_names =  ("chr", "position", "depth")
    df = pd.read_csv(depth_file, sep="\t", header=None, names=header_names)
    
    # comprobamos que el punto de inicio no está fuera del rango del genoma
    start = coordinate - margin

    if start < 1:
        start = 1
    # nos aseguramos de que no va a dividir por 0.
    if coordinate == 1 or df.empty:
        meanPre = 0
    else:
        positions_checked = coordinate - start
        if positions_checked == 0:
            positions_checked = 1 # aseguramos que no pueda haber /0
        sum_depth_values = 0
        for pos in range(start, coordinate):
            value = df[df['position'] == pos]['depth'].values[0]
            sum_depth_values += value
    
        # calculamos la media Pre-coordenada
        meanPre = sum_depth_values / positions_checked

    return meanPre


def mean_depth_neighborhood_post_coordinate(coordinate, sample_name, depth_txt,
                                            margin, len_wt):
    """
    Calcula la profundidad media de las 'margin' posiciones después 
    de la coordenada
    Args:
        coordinate  (int)   [in]    Posición (base 1) de la que vamos a evaluar
                                    la profundidad media de 
        depth_txt   (str)   [in]    Fichero con la profundidad de todas las
                                    posiciones
        margin      (int)   [in]    Número de posiciones de distancia con 
                                    respecto a la coordenada con las que se 
                                    quiere calcular la media de profundidad
    
    Returns:
        meanPre (float) [out]   Profundidad media de las posiciones 'margin'
                                anteriores a la coordenada
        meanPost (float) [out]   Profundidad media de las posiciones 'margin'
                                posteriores a la coordenada
    """
    
    # cargamos el fichero con la profundidad como pd.df

    depth_file = "Outputs/alignments/" + sample_name + "/" + depth_txt
    header_names =  ("chr", "position", "depth")
    df = pd.read_csv(depth_file, sep="\t", header=None, names=header_names)
    
    # comprobamos que el punto de inicio no está fuera del rango del genoma
    end = coordinate + margin + 1
    
    # si el final proyectado sobrepasa la longitud del genoma llegamos hasta la última posición	
    if end >= len_wt:
        end = len_wt
    # si el final proyectado es <= al margen (significa que la coordenada es < 0), establecemos la coordenada en 1
    elif end <= margin:
        coordinate = 1
    # nos aseguramos de que no va a dividir por 0.
    if coordinate == len_wt or df.empty:
        meanPost = 0
    else:
        positions_checked = end - (coordinate + 1)
        if positions_checked <= 0:
            positions_checked = 1 # aseguramos que no pueda haber /0
        sum_depth_values = 0
        for pos in range(coordinate + 1, end):
            value = df[df['position'] == pos]['depth'].values[0]
            sum_depth_values += value
                    
        # calculamos la media Pre-coordenada
        meanPost = sum_depth_values / positions_checked

    return meanPost
    
def write_depth_neighborhood_coordinate(sample_name, df, margin, len_wt):
    """
    Escribe el fichero con

    Args:
        sample_name (str)   [in]    Nombre de la muestra
        df  (pd.DataFrame)  [in]    Tabla con los datos resumen ya cargada en 
                                    memoria y preparada para modificarse
        margin  (int)   [in]    Número de posiciones antes y después de las 
                                coordenadas BP y RI que vamos a tener en cuenta
                                para el cálculo de la profundidad media.

    """

    #bt2 = sample_name + "_bt2_sorted_depth.txt"
    #bm = sample_name + "_bm_sorted_depth.txt"
    bt2_mapped = sample_name + "_bt2_sorted_mapped_depth.txt"
    bm_mapped = sample_name + "_bm_sorted_mapped_depth.txt"
    bt2_H = sample_name + "_bt2_sorted_H_depth.txt"
    bm_H = sample_name + "_bm_sorted_H_depth.txt"

    #df = pd.read_csv("Outputs/unificate_table_wDepths.csv", header=0)
    
    # Añadimos columnas con información sobre la profundidad
    """
    df[['mean_depth_pre'+str(margin)+'BP_bt2']] = df.apply(lambda row : 
        mean_depth_neighborhood_pre_coordinate(row['BP'], sample_name, bt2, 
                                            margin), axis=1)
    df[['mean_depth_post'+str(margin)+'BP_bt2']] = df.apply(lambda row : 
        mean_depth_neighborhood_post_coordinate(row['BP'], sample_name, bt2, 
                                            margin), axis=1)    
    df[['mean_depth_pre'+str(margin)+'RI_bt2']] = df.apply(lambda row : 
        mean_depth_neighborhood_pre_coordinate(row['RI'], sample_name, bt2, 
                                            margin), axis=1)
    df[['mean_depth_post'+str(margin)+'RI_bt2']] = df.apply(lambda row : 
        mean_depth_neighborhood_post_coordinate(row['RI'], sample_name, bt2, 
                                            margin), axis=1) 
    """
    df[['mean_depth_pre'+str(margin)+'BP_bt2_mapped']] = df.apply(lambda row : 
        mean_depth_neighborhood_pre_coordinate(row['BP'], sample_name, bt2_mapped,
                                            margin), axis=1)
    df[['mean_depth_post'+str(margin)+'BP_bt2_mapped']] = df.apply(lambda row : 
        mean_depth_neighborhood_post_coordinate(row['BP'], sample_name, 
                                            bt2_mapped, margin, len_wt), axis=1)    
    df[['mean_depth_pre'+str(margin)+'RI_bt2_mapped']] = df.apply(lambda row : 
        mean_depth_neighborhood_pre_coordinate(row['RI'], sample_name, bt2_mapped, 
                                            margin), axis=1)
    df[['mean_depth_post'+str(margin)+'RI_bt2_mapped']] = df.apply(lambda row : 
        mean_depth_neighborhood_post_coordinate(row['RI'], sample_name, bt2_mapped,
                                            margin, len_wt), axis=1)    
    
    df[['mean_depth_pre'+str(margin)+'BP_bt2_H']] = df.apply(lambda row : 
        mean_depth_neighborhood_pre_coordinate(row['BP'], sample_name, bt2_H, 
                                            margin), axis=1)
    df[['mean_depth_post'+str(margin)+'BP_bt2_H']] = df.apply(lambda row : 
        mean_depth_neighborhood_post_coordinate(row['BP'], sample_name, bt2_H, 
                                            margin, len_wt), axis=1)    
    df[['mean_depth_pre'+str(margin)+'RI_bt2_H']] = df.apply(lambda row : 
        mean_depth_neighborhood_pre_coordinate(row['RI'], sample_name, bt2_H, 
                                            margin), axis=1)
    df[['mean_depth_post'+str(margin)+'RI_bt2_H']] = df.apply(lambda row : 
        mean_depth_neighborhood_post_coordinate(row['RI'], sample_name, bt2_H, 
                                            margin, len_wt), axis=1)    

    """
    df[['mean_depth_pre'+str(margin)+'BP_bm']] = df.apply(lambda row : 
        mean_depth_neighborhood_pre_coordinate(row['BP'], sample_name, bm, margin), 
        axis=1)
    df[['mean_depth_post'+str(margin)+'BP_bm']] = df.apply(lambda row : 
        mean_depth_neighborhood_post_coordinate(row['BP'], sample_name, bm, margin), 
        axis=1)    
    df[['mean_depth_pre'+str(margin)+'RI_bm']] = df.apply(lambda row : 
        mean_depth_neighborhood_pre_coordinate(row['RI'], sample_name, bm, margin), 
        axis=1)
    df[['mean_depth_post'+str(margin)+'RI_bm']] = df.apply(lambda row : 
        mean_depth_neighborhood_post_coordinate(row['RI'], sample_name, bm, margin), 
        axis=1) 
    """
    df[['mean_depth_pre'+str(margin)+'BP_bm_mapped']] = df.apply(lambda row : 
        mean_depth_neighborhood_pre_coordinate(row['BP'], sample_name, bm_mapped, margin), 
        axis=1)
    df[['mean_depth_post'+str(margin)+'BP_bm_mapped']] = df.apply(lambda row : 
        mean_depth_neighborhood_post_coordinate(row['BP'], sample_name, bm_mapped, 
                                                margin, len_wt), axis=1)    
    df[['mean_depth_pre'+str(margin)+'RI_bm_mapped']] = df.apply(lambda row : 
        mean_depth_neighborhood_pre_coordinate(row['RI'], sample_name, bm_mapped, margin), 
        axis=1)
    df[['mean_depth_post'+str(margin)+'RI_bm_mapped']] = df.apply(lambda row : 
        mean_depth_neighborhood_post_coordinate(row['RI'], sample_name, bm_mapped, 
                                            margin, len_wt), axis=1)    
    
    df[['mean_depth_pre'+str(margin)+'BP_bm_H']] = df.apply(lambda row : 
        mean_depth_neighborhood_pre_coordinate(row['BP'], sample_name, bm_H, margin), 
        axis=1)
    df[['mean_depth_post'+str(margin)+'BP_bm_H']] = df.apply(lambda row : 
        mean_depth_neighborhood_post_coordinate(row['BP'], sample_name, bm_H, 
                                            margin, len_wt), axis=1)    
    df[['mean_depth_pre'+str(margin)+'RI_bm_H']] = df.apply(lambda row : 
        mean_depth_neighborhood_pre_coordinate(row['RI'], sample_name, bm_H, margin), 
        axis=1)
    df[['mean_depth_post'+str(margin)+'RI_bm_H']] = df.apply(lambda row : 
        mean_depth_neighborhood_post_coordinate(row['RI'], sample_name, bm_H, 
                                            margin, len_wt), axis=1)    
    
    # write
    #df.to_csv("Outputs/unificate_temp_table_wMeanDepths.csv", index=False)
    return df

def proportion_coordinate(read_counts, depth_coordinate_mapped):
    """
    Calcula la ratio entre el dvg detectado y la población total. Considerando
    como total la suma de población "nativa" (=mapped) y defectiva.

     pCoordinate = read_counts / depthCoordinate_mapped
    
    Args:
        read_counts (int)   [in]    Número de lecturas en las que se ha 
                                    encontrado el dvg
        depth_coordinate_mapped  (int)   [in] Reads que han mapeado en la 
                                    coordenada (BP o RI)

    Returns:
        ratio   (float) [out]   Proporción de reads defectivas vs totales 
    """
    if depth_coordinate_mapped != 0:
        ratio = read_counts / depth_coordinate_mapped
    else:
        ratio = 0

    return ratio


def write_proportion_reads(df, bool_vir, bool_dit):
    """
    Añade los ratios pBP y pRI de cada dvg en el pd.DataFrame df

    Args:
        df  (pd.DataFrame)  [in]    Tabla con los datos resumen ya cargada en 
                                    memoria y preparada para modificarse
        bool_vir    (bool)  Indica si virema ha encontrado DVGs
        bool_dit    (bool)  Indica si ditector ha encontrado DVGs
    Returns:
        df  (pd.DataFrame)  [in]    df al que se le han añadido los datos sobre
                                    ratio
    """

    #df = pd.read_csv("Outputs/unificate_table_wDepths.csv", header=0)

    # Proporciones según lo detectado por ViReMa
    if bool_vir:
        df[['pBP_virema']] = df.apply(lambda row : 
            proportion_coordinate(row['read_counts_virema'],
            row['depthBP_bt2_mapped']), axis=1)
        df[['pRI_virema']] = df.apply(lambda row : 
            proportion_coordinate(row['read_counts_virema'],
            row['depthRI_bt2_mapped']), axis=1)

    # Proporciones según lo detectado por DItector
    if bool_dit:
        df[['pBP_ditector']] = df.apply(lambda row : 
            proportion_coordinate(row['read_counts_ditector'],
            row['depthBP_bm_mapped']), axis=1)
        df[['pRI_ditector']] = df.apply(lambda row : 
            proportion_coordinate(row['read_counts_virema'],
            row['depthRI_bm_mapped']), axis=1)
    #df.to_csv("Outputs/unificate_temp_table_wProportions.csv", index=False)

    return df

#–--------------------------   
# LENGTH DVG RECONSTRUCTION
#–--------------------------
def first_seqment(len_wt, DVG_type, BP):
    """
    Longitud del primer segmento del DVG (hasta BP incluido) en función del tipo
    de DVG.
    Dos procedimientos diferentes para averiguar la longitud:
        a) Deletion_forward, Insertion_forward, 5cb/sb --> 1:BP
        b) Deletion_reverse, Insertion_reverse, 3cb/sb --> len_wt - (RI - 1)


    Args:
        len_wt   (int)   [in]    Longitud del genoma tomado como referencia
        DVG_type        (str)   [in]    Tipo de DVG
        BP              (int)   [in]    Break point: último nucleótido del 
                                        primer segmento
    Returns:
        len_first_seqment   (int)   [out]   Longitud del 1er segmento del DVG
    """

    a = ['Deletion_forward', 'Insertion_forward', '5cb/sb']
    b = ['Deletion_reverse', 'Insertion_reverse', '3cb/sb']

    if DVG_type in a:
        len_first_seqment = BP
    
    elif DVG_type in b:
        len_first_seqment = len_wt - BP + 1

    else:
        len_first_seqment = 'Error'
    
    return len_first_seqment


def second_seqment(len_wt, DVG_type, RI):
    """
    Longitud del segundo segmento del DVG (desde RI incluido) en función del tipo
    de DVG.
    Dos procedimientos diferentes para averiguar la longitud:
        a) Deletion_forward, Insertion_forward, 5cb/sb --> 1:BP
        b) Deletion_reverse, Insertion_reverse, 3cb/sb --> len_wt - (RI - 1)


    Args:
        len_wt   (int)   [in]    Longitud del genoma tomado como referencia
        DVG_type        (str)   [in]    Tipo de DVG
        BP              (int)   [in]    Break point: último nucleótido del 
                                        primer segmento
    Returns:
        len_second_seqment   (int)   [out]   Longitud del 1er segmento del DVG
    """

    a = ['Deletion_reverse', 'Insertion_reverse', '5cb/sb']
    b = ['Deletion_forward', 'Insertion_forward', '3cb/sb']

    if DVG_type in a:
        len_second_seqment = RI
    
    elif DVG_type in b:
        len_second_seqment = len_wt - RI + 1
    else:
        len_second_seqment = 'Error'
    
    return len_second_seqment

def len_dvg(len_wt, DVG_type, BP, RI):
    """
    Longitud completa del dvg: len_first_seqment + len_second_seqment
    """
    len_first_seqment = first_seqment(len_wt, DVG_type, BP)
    len_second_seqment = second_seqment(len_wt, DVG_type, RI)

    len_dvg = len_first_seqment + len_second_seqment

    return len_dvg

def add_lengths_to_df(df, len_wt):
    """
    Añade la longitud teórica del DVG al df

    Args:
        df  (pd.DataFrame)  [in]   Tabla accesible con la información de todos
                                    los DVGs encontrados
        len_wt  (int)   [in]    Longitud del genoma de referencia. Default is
                                29903 (SARS-CoV-2 Wuhan1)    
    Returns:
        df  (pd.DataFrame)  [out]   Table with the theoric length of the DVG
                                    added
    """
    df[['length_dvg']] = df.apply(lambda row : len_dvg(len_wt,
                        row['DVG_type'], row['BP'], row['RI']), axis=1)
    
    return df


def generate_stats_alignment(sample_name):
    """
    Generamos los estadísticos de todos los alineamientos
    Args:
        sample_name (str)   [in]    Nombre de la muestra
    """

    dir_alignments = "Outputs/alignments/" + sample_name + "/"

    #sorted_alignments = ['_bm_sorted.bam', '_bt2_sorted.bam']
    mapped_tails = ['_bm_sorted_mapped.bam', '_bt2_sorted_mapped.bam']
    H_tails = ['_bm_sorted_H.bam', '_bt2_sorted_H.bam']
    #sorted_alignments + 
    all_tails = mapped_tails + H_tails

    for tail in all_tails:
        alignment = dir_alignments + sample_name + tail
        base_file = os.path.basename(alignment).split(".")[-2]
        out_file = dir_alignments + base_file + "_stats.txt"
        # generamos el fichero con los datos

        os.system('samtools coverage {} > {}'.format(alignment, out_file))
        os.system('')


def read_stats_alignments(df_stats):
    """
    Asigna las variables de `samtools coverage`que nos interesan, es decir:
    numreads,	covbases,	coverage,	meandepth,	meanbaseq y	meanmapq

    Args:
        df_stats  (pd.DataFrame)  [in]    Fichero con la información sobre los 
                                    estadísticos del alineamiento en un formato
                                    preparado para asignar los valores de las
                                    variables de interés
    Returns:
        numreads    (int)   [out]   Número de lecturas del alineamiento
        covbases    (int)   [out]   Número de bases cubiertas con una profundidad
                                    >= 1
        coverage    (int)   [out]   Porcentaje que representa: covbases/totalBases
        meandepth   (float) [out]   Profundidad media del alineamiento
        meanbaseq   (int)   [out]   Mean baseQ in covered region
        meanmapq    (int)   [out]   Mean mapQ of selected reads
    """
    
    numreads = df_stats['numreads'][0]
    #covbases = df_stats['covbases'][0]
    #coverage = df_stats['coverage'][0]
    #meandepth = df_stats['meandepth'][0]
    #meanbaseq = df_stats['meanbaseq'][0]
    #meanmapq = df_stats['meanmapq'][0]

    return numreads 
    #,covbases, coverage, meandepth, meanbaseq, meanmapq

# Hay que contar con que puede estar vacío
def add_stats_alignment(df, alignment_stats): 
    """
    Añade al df los valores de estádisticos extraídos con `samtools coverage` 
    del alineamiento

    Args:
        df  (pd.DataFrame)  [in]    Tabla con los datos resumen ya cargada en 
                                    memoria y preparada para modificarse
        alignment_stats   (str)  [in] Path relativo del fichero del que se van a 
                                    extraer los valores
    Return:
        df  (pd.DataFrame)  [out]   Tabla con los datos resumen y los datos 
                                    sobre estadísticos del alineamiento 
                                    añadidos
    """

    #df = pd.read_csv("Outputs/unificated_table.csv", header=0)
    #file_stats = "Outputs/alignments/" + alignment_stats

    # Comprobamos que no está vacío
    if os.path.getsize(alignment_stats) != 0:
        # cargamos el df con los datos sobre los estadísticos
        df_stats = pd.read_csv(alignment_stats, header=0, sep="\t")
        numreads = read_stats_alignments(df_stats)
        #, covbases, coverage, meandepth, meanbaseq, meanmapq = \
        #    read_stats_alignments(df_stats)
    else:
        numreads = 0
    
    # asignar el tail_name que añadir en las columnas
    ## algoritm
    if "bt2" in alignment_stats:
        al = "_bt2"
    elif "bm" in alignment_stats:
        al = "_bm"
    ## filter
    if "H" in alignment_stats:
        tail = "_H" 
    elif "mapped" in alignment_stats:
        tail = "_mapped"
    else:
        tail = "_all" 
    
    df[['numReads' + al + tail]] = numreads
    #df[['covBases' + al  + tail]] = covbases
    #df[['coverage' + al  + tail]] = coverage
    #df[['meanDepth' + al  + tail]] = meandepth
    #df[['meanBaseQ' + al  + tail]] = meanbaseq
    #df[['meanMapQ' + al  + tail]] = meanmapq

    return df


def add_stats_all_alignments(sample_name, df):
    """
    Añadimos los datos sobre numreads, covbases, coverage, meandepth, meanbaseq
    y meanmapq de cada alineamiento a la tabla resumen
    Args:
        sample_name (str)   [in]    Nombre de la muestra
        df  (pd.DataFrame)  [in]    Tabla con los datos resumen ya cargada en 
                                    memoria y preparada para modificarse
    Returns:
        df  (pd.DataFrame)  [out]   Tabla con los datos resumen y los datos 
                                    sobre estadísticos de todos los alineamientos
                                    ya añadidos        

    """
    dir_alignments =  "Outputs/alignments/" + sample_name + "/"

    #bt2 = sample_name + "_bt2_sorted_stats.txt"
    #bm = sample_name + "_bm_sorted_stats.txt"
    bt2_mapped = sample_name + "_bt2_sorted_mapped_stats.txt"
    bm_mapped = sample_name + "_bm_sorted_mapped_stats.txt"
    bt2_H = sample_name + "_bt2_sorted_H_stats.txt"
    bm_H = sample_name + "_bm_sorted_H_stats.txt"
    #bt2, bm, 
    stat_files = [bt2_mapped, bm_mapped, bt2_H, bm_H]

    for stat_file in stat_files:
        alignment_stats = dir_alignments + stat_file
        df = add_stats_alignment(df, alignment_stats)

    return df


#–------------------------------------------------     
# RPHT Reads Per Hundred Thousand
#–------------------------------------------------  

def calculate_rpht(read_counts, mapped_reads):
    """
    El rpht o 'reads per hundred thousand' es una medida rudimentaria de 
    normalización que tiene en cuenta el número de reads con DVG con respecto
    al total de reads que han mapeado por 10^5.
        rpht = (read_counts / mapped_reads) * 100000

    Args:
        read_counts_{program}   (int)   [in]    Número de reads en las que el 
                                        programa {program} ha detectado el
                                        evento de fusión
                                        Columna 'read_counts_{virema/ditector}'
        mapped_reads_{al}   (int)   [in]    Número de reads totales del 
                                    alineamiento {al}_mapped 
                                    (tras samtools -F 4 al alineamiento original)
                                    Columna 'numReads_{bt2/bm}_mapped'
    Returns
        rpht    (float) [out]   Lecturas con dvg por cada 10000 lecturas
                                mapeadas. Métrica normalizada de abundancia del
                                evento.
    """

    return (read_counts / mapped_reads) * 100000 

def add_rpht(df, bool_vir, bool_dit):
    """
    Añade los valores de rpht al df

    """
    #import pandas as pd
    # valores de rpht según la detección de ViReMa
    if bool_vir:
        df[['rpht_virema']] = df.apply(lambda row : 
            calculate_rpht(row['read_counts_virema'], row['numReads_bt2_mapped']), 
            axis=1)

    # valores de rpht según la detección de DItector
    if bool_dit:
        df[['rpht_ditector']] = df.apply(lambda row : 
            calculate_rpht(row['read_counts_ditector'], row['numReads_bm_mapped']),
            axis=1)
    
    #df.to_csv("Outputs/unificate_temp_table_wRpht.csv", index=False)

    return df


def extract_info_fq(fq):
    """
    Usa el módulo pyfastx para extraer la info del fichero fastq
    Args:
        fq  (str)   [in]    Ruta completa al fichero fastq que se va a 
                            procesar
    Returns:
        fq_reads
        gc_content
        average_read_length
        max_read_length
        min_read_length
        max_quality
        min_quality
    """

    fq = pyfastx.Fastq(fq)

    fq_reads = len(fq)
    gc_content = fq.gc_content
    average_read_length = fq.avglen
    max_read_length = fq.maxlen
    min_read_length = fq.minlen
    max_quality = fq.maxqual
    min_quality = fq.minqual

    return fq_reads, gc_content, average_read_length, max_read_length, \
    min_read_length, max_quality, min_quality


def add_info_fq(sample_name, df, fq):
    """
    # REPENSAR, ES UNA INFO QUE NO VAMOS A USAR PARA EL ENTRENAMIENTO DEL MODELO
    # Y PARA QUÉ ESCRIBIRLA EN UN FICHERO. A NO SER QUE SEA INFO QUE QUERAMOS
    # AÑADIR A ALGUNA DE LAS VISUALIZACIONES DE DESPUÉS. DE MOMENTO YA NO LA
    # VAMOS A ESCRIBIR EN LA TABLA RESUMEN.
    De momento añade la información extraída del fichero fq a analizar en el
    df tabla_resumen. 
    Args:
        sample_name (str)   [in]    Nombre de la muestra
        df  (pd.DataFrame)  [in]    Tabla con los datos resumen ya cargada en 
                                    memoria y preparada para modificarse
        fq  (str)   [in]    Ruta completa al fichero fq que queremos procesar

    Returns:
        df  (pd.DataFrame)  [out]   Tabla con los datos resumen y los datos 
                                    con la infromación del fichero fastq 
                                    añadidos
    """

    info_fq = extract_info_fq(fq)
    
    df[['fq_reads_' + sample_name]] = info_fq[0]
    df[['gc_content_' + sample_name]] = info_fq[1]
    df[['average_read_length_' + sample_name]] = info_fq[2]
    df[['max_read_length_' + sample_name]] = info_fq[3]
    df[['min_read_length_' + sample_name]] = info_fq[4]
    df[['max_quality_' + sample_name]] = info_fq[5]
    df[['min_quality_' + sample_name]] = info_fq[6]

    return df


#–------------------------------------------------     
# GENERATING AND ADDING THE INFO TO THE RESUME DF
#–------------------------------------------------  

def generate_df(sample_name):
    """
    Genera la estructura pd.DataFrame donde carga los datos del fichero 
    a procesar

    Args:    
        sample_name (str)   [in]    Nombre de la muestra
        fq  (str)   [in]    Ruta completa al fichero fastq que se va a procesar
    
    Returns:
        df  (pd.DataFrame)  [out]   Tabla cargada en memoria y accesible
    """
    dir_name =  "Outputs/"
    table_name = sample_name + "_unificated_table.csv"
    table_txt = dir_name + table_name

    out_file = dir_name + sample_name + "_resume_table.csv"

    df = pd.read_csv(table_txt, header=0)

    return df, out_file

def add_features(sample_name, margin, len_wt, bool_vir, bool_dit, df):
    """
    Añade las variables a la tabla unificada de DVGs y las escribe en el 
    fichero 
    Args:
        df  (pd.DataFrame)  [in]    Tabla unoficada con la información básica
                                de todos los eventos identificados por los
                                programas.
        margin  (int)   [in]    Número de posiciones antes y después de las 
                                coordenadas que van a ser tenidos en cuenta 
                                para el cálculo de la media. Default is 5.

    Return:
        df  (pd.DataFrame)  [out]   Tabla unificada con todos los eventos 
                                identificados y las variables añadidas
    """

    # añadimos la profundidad de cada coordenada BP y RI de los diferentes 
    # alineamientos 
    df = write_depth_coordinate(sample_name, df, len_wt)
    # añadimos las ratios
    df = write_proportion_reads(df, bool_vir, bool_dit)
    # añadimos los stats propios de cada alineamiento
    df = add_stats_all_alignments(sample_name, df)
    # añadimos los valores de rpht ((read_counts_dvg / mapped_reads)*100000)
    df = add_rpht(df, bool_vir, bool_dit)
    ## añadimos la información sobre el fastq
    #df = add_info_fq(sample_name, df, fq)
    # longitud de los dvgs
    df = add_lengths_to_df(df, len_wt)

    # añadimos la profundidad media en las fronteras de BP y RI de todos los 
    # tipos de alineamiento de que disponemos.
    df = write_depth_neighborhood_coordinate(sample_name, df, margin, len_wt)

    #df.to_csv(out_file, index=False)
    return df


def parallelize_add_features(df, sample_name, margin, len_wt, 
                            bool_vir, bool_dit, add_features, 
                            n_cores):
    """
    This function parallelize the process of adding features to the final df.
    First, check the number of events is at least equal to the n_cores, if not
    reduce the n_cores to the n_events

    Args:
        df  (pd.DataFrame)  [in]    
        sample_name (str)   [in]
        margin  (int)   [in]
        add_features    (function)  [in]
        n_cores (int)   [in]
    
    Returns:
        df  (pd.DataFrame)  [out]  
    """
    # check the number of events is >= n_cores. If not, reduce the n_cores to 
    # the n_events
    if len(df) < n_cores:
        n_cores = len(df)
    # limit the number of cores to 10, because in some serves could cause 
    # problems
    elif n_cores > 10:
        n_cores = 10
    # iterable argument
    df_split = np.array_split(df, n_cores)
    # constant arguments
    constant_args = [sample_name, margin, len_wt, bool_vir, bool_dit]

    pool = Pool(n_cores)

    results = [pool.apply_async(add_features, constant_args + [i]) for i in df_split]
    pool.close()
    pool.join()

    df = pd.concat(r.get() for r in results)
    
    return df
