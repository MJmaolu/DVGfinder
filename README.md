<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/MJmaolu/DVGfinder">
    <img src="LOGO DVGfinder_marron.png" alt="Logo" width="700">
  </a>  
  <h3 align="center">DVGfinder_v3</h3>

<p align="center">
    <b>DVGfinder</b> is a tool that integrates the most used DVGs search algorithms, <b>ViReMa-a</b> and <b>DI-tector</b>, in a unique workflow, making the analysis of a sample easier and more intuitive. 
  
  Also, DVGfinder_v3 implements a Gradient Boosting classifier to try to reduce the number of false positives introduced by the search algorithms and generates an HTML report with interactive tables and plots that facilitates a first exploration of the results.
    <br />
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a hrf="#usage">Usage</a></li>
    <li><a href="#tutorial">Tutorial</a></li>   
    <li><a href="#labeled-dataset">Labeled Dataset</a></li>
    <li><a href="#reference">Reference</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<!-- BUILT WITH -->
## Built with
### DVGs search algorithms

* **ViReMa-a (0.23)**:

Routh A, Johnson JE. Discovery of functional genomic motifs in viruses with ViReMa-a Virus Recombination Mapper-for analysis of next-generation sequencing data. Nucleic Acids Res. 2014 Jan;42(2):e11. [doi: 10.1093/nar/gkt916](https://academic.oup.com/nar/article/42/2/e11/1024459). Epub 2013 Oct 16. PMID: 24137010; PMCID: PMC3902915.
  

* **DI-tector_06.py**: 

Beauclair G, Mura M, Combredet C, Tangy F, Jouvenet N, Komarova AV. DI-tector: defective interfering viral genomes' detector for next-generation sequencing data. RNA. 2018 Oct;24(10):1285-1296. [doi: 10.1261/rna.066910.118](https://pubmed.ncbi.nlm.nih.gov/30012569/). Epub 2018 Jul 16. PMID: 30012569; PMCID: PMC6140465.
 
  

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Prerequisites

DVGfinder uses the ViReMa-a (v0.23) and DI-tector_06.py programs. 

This third party scripts are in the ExternalNeeds directory so you only have to follow the nexts steps to run DVGfinder.


### Installation

1. Clone the repo in the directory of your choice
   ```sh
   git clone https://github.com/MJmaolu/DVGfinder.git
   ```

2. Go to the DVGfinder directory
   ```sh
   cd DVGfinder
   ```
   
3. Give execution permission to all the scripts in the DVGfinder directory
   ```sh
   chmod -R +x .
   ```
4. Create a new environment with conda with all the dependencies needed to run DVGfinder
   ```sh
   conda env create -f dvgfinder_env.yaml
   ```
   
5. Activate DVGfinder environment 
   ```sh
   conda activate dvgfinder_env
   ```

<!-- USAGE -->
## Usage

```python3 DVGfinder_v3.py -fq path_to_fastq_file [-r path_to_fasta_virus_reference] [-m margin] [-t threshold] [-n number_processes]```

`-r` The fasta of the viral reference and its indexed files by `bwa` and `bowtie` should be all in the path ExternalNeeds/references.
 
<!-- TUTORIAL -->
## Tutorial

You can explore an example of results in the directory 'tumvas72_N100K_l100'. 

To test the program follow the next steps:

1. Activate the environment 

2. Run DVGfinder on the example sample

```
python3 DVGfinder_v3.py -fq tumvas72_N100K_l100/tumvas72_N100K_l100.fq -r ExternalNeeds/references/TuMV-AS.fasta -t probability_threshold_to_filter_as_realsDVGs -n number_of_process
```

3. Wait and your results will appear in the 'FinalReports' directory. In addition, an html report will open in your default browser.

[Link to an example HTML report](http://147.156.206.144/appweb/tumvas72_N100K_l100_report.html)

#### About the HTML report: 

- The results are displayed at three levels (tabs at the top): ALL, CONSENSUS and FILTERED ML. 

- CONSENSUS and FILTERED ML only appear if both search algorithms have identified DVGs in the sample.

- The same information displayed in the interactive tables is written in csv files.

<!-- LABELED DATASET -->
## Labeled Dataset

The dataset used to generate the classificator is avalaible as '630N5Ml100_v2_metrics_labeledDataset.csv'. 

For this version of DVGfinder we have used a Gradient Boosting Classifier algorithm to generate the model but I encorage you to play directly with the data and try to improve it. 

<!-- REFERENCE -->
## Reference

*Preprint version*:

- Olmo-Uceda, M.J.; Muñoz-Sánchez, J.C.; Lasso-Giraldo, W.; Arnau, V.; Díaz-Villanueva, W.; Elena, S.F. DVGfinder: A Metasearch Engine for Identifying Defective Viral Genomes in RNA-Seq Data. *Preprints* 2022, 2022030110 [doi: 10.20944/preprints202203.0110.v1](https://www.preprints.org/manuscript/202203.0110/v1).
 
<!-- CONTACT -->
## Contact

María José Olmo-Uceda - mariajose.olmo@csic.es

PhD student

*EvolSysVir Group*, I<sup>2</sup>SysBio (CSIC-UV) 

---

Project Link: [https://github.com/MJmaolu/DVGfinder](https://github.com/MJmaolu/DVGfinder)

Page Link: [https://mjmaolu.github.io/DVGfinder/](https://mjmaolu.github.io/DVGfinder/)

<p align='right'> 
  <b>Under Construction</b> 
</p> 
<p align='right'> 
  Any suggestions will be welcome :hugs:
</p>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/MJmaolu/DVGfinder.svg?style=for-the-badge
[contributors-url]: https://github.com/MJmaolu/DVGfinder/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/MJmaolu/DVGfinder.svg?style=for-the-badge
[forks-url]: https://github.com/MJmaolu/DVGfinder/network/members
[stars-shield]: https://img.shields.io/github/stars/MJmaolu/DVGfinder.svg?style=for-the-badge
[stars-url]: https://github.com/MJmaolu/DVGfinder/stargazers
[issues-shield]: https://img.shields.io/github/issues/MJmaolu/DVGfinder.svg?style=for-the-badge
[issues-url]: https://github.com/gMJmaolu/DVGfinder/issues
[license-shield]: https://img.shields.io/github/license/MJmaolu/DVGfinder.svg?style=for-the-badge
[license-url]: https://github.com/MJmaolu/DVGfinder/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: www.linkedin.com/in/maria-jose-olmo-uceda
