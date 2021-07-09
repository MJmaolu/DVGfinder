<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/MJmaolu/DVGfinder">
    <img src="LOGO%20DVGfinder_marron.png" alt="Logo" width="800" /*height="80"*/>
  </a>

  <h3 align="center">DVGfinder</h3>

  <p align="center">
    DVGfinder is a program that integrates the DVGs ViReMa-a and DI-tector search algorithms, making the analysis of a sample easier and more intuitive. 
    <br />
    <a href="https://github.com/MJmaolu/DVGfinder"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/MJmaolu/DVGfinder">View Demo</a>
    ·
    <a href="https://github.com/MJmaolu/DVGfinder/issues">Report Bug</a>
    ·
    <a href="https://github.com/MJmaolu/DVGfinder/issues">Request Feature</a>
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
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>

### DVGs identifying algorithms

* **ViReMa-a (0.23)**:

Routh A, Johnson JE. Discovery of functional genomic motifs in viruses with ViReMa-a Virus Recombination Mapper-for analysis of next-generation sequencing data. Nucleic Acids Res. 2014 Jan;42(2):e11. <doi: 10.1093/nar/gkt916>. Epub 2013 Oct 16. PMID: 24137010; PMCID: PMC3902915.
  

* **DI-tector_06.py**: 

Beauclair G, Mura M, Combredet C, Tangy F, Jouvenet N, Komarova AV. DI-tector: defective interfering viral genomes' detector for next-generation sequencing data. RNA. 2018 Oct;24(10):1285-1296. <doi: 10.1261/rna.066910.118>. Epub 2018 Jul 16. PMID: 30012569; PMCID: PMC6140465.
 
  

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Prerequisites

DVGfinder uses the ViReMa-a (v0.23) and DI-tector_06.py programs. 

This third party scripts are in the ExternalNeeds directory so you only have to follow the nexts steps to run DVGfinder.


### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/MJmaolu/DVGfinder.git
   ```
2. Unzip model file

  ```sh
  unzip -d rf_tunning_x17.sav.zip
  ```

   
2. Create a new environment with conda with all the dependencies needed to run DVGfinder
   ```sh
   conda env create -f env env/DVGfinder_env.yaml
   ```
   
3. Activate DVGfinder environment 
   ```sh
   conda activate dvgfinder_env
   ```

<!-- USAGE EXAMPLES -->
## Usage

```python3 DVGfinder.py -fq fastq_file [-r virus_reference] [-m margin] [-n number_processes]```


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

María José Olmo-Uceda - maolu@alumni.uv.es

Project Link: [https://github.com/MJmaolu/DVGfinder](https://github.com/MJmaolu/DVGfinder)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

* []()
* []()
* []()





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
