ADD ZENODO LINK XXX

Table of contents
=================

<!--ts-->
   * [Purpose of the codes](#purpose-of-the-codes)
   * [How to cite](#how-to-cite)
   * [How to Install the Project](#how-to-install-the-project)
      * [Running in a local installation](#running-in-a-local-installation)
      * [Running on Google drive](#running-on-google-drive)
   * [Output data](#output-data)
   * [Examples used as an input](#examples-used-as-an-input)
   * [Acknowledgements](#acknowledgements)
   * [License](#license)
<!--te-->

# Purpose of the codes

The AINTBAD (Analysis of Ice nucleation Temperature for $B$ and $A$ Determination) includes the numerical integration in Equation 7 and the analytical model of Equation 8 in an python code to estimate $A$, $B$ and $J$ from experimental drop-freezing data. The code outputs the parameters $A$ and $B$. These are used to compute the nucleation barriers $\Delta G$, the temperature that corresponds to 50\% of frozen droplets $T_{50}$, and the homogeneous nucleation rate evaluated at $T_{50}$.

The code IPA (Inhomogeneous Poisson Analysis) is capable of taking various parametrizations for the homogeneous nucleation rate $J_{hom}(T)$, the droplet size distributions (Gaussian, Gamma, uniform, exponential, etc.), and cooling rates to compute the survival probability or fraction of frozen droplets. We use the nucleation rate data vs temperature as the input to compute the survival probability.

# How to cite

This code is the active version of the code archived here and documented in Addula, de Almeida Ribeiro, Molinero & Peters, Ice nucleation from drop-freezing experiments: Impact of droplet volume dispersion and cooling rates, DOI: XXX

You are welcome to use and distribute this code as you see fit, but it remains the intellectual property of the authors and must be cited appropriately (please cite the paper). Direct any questions about this code to: Ravi Addula (raddula@illinois.edu) and Ingrid de A. Ribeiro (ingrid.ribeiro@utah.edu), or create an issue in this Github repository.

-----
# How to Install the Project
## Running in a local installation

Launch with:
```
python AINTBAD_CODE.py
```
or
```
python IPA_CODE.py
```

## Running on Google drive

Colaboratory, or [Colab](https://colab.research.google.com/?utm_source=scs-index) for short, is a product from Google Research. It allows anybody to write and execute arbitrary Python code through the browser (Firefox, Safari, Chrome etc).

To use it, follow the steps:

- download the directory AINTBAD (or IPA), and then upload it in a **Google drive** directory;

- inside the directory, click on "add more apps" if you do not see Colaboratory,

<img width="1528" alt="image" src="https://github.com/Molinero-Group/volume-dispersion/assets/60315074/30be9f33-4fe8-4254-84da-53a7c0ae1d5a">

- type *colaboratory* and install it

<img width="969" alt="image" src="https://user-images.githubusercontent.com/60315074/184701819-8e4baaf3-f2b7-47d2-b067-bb3947151fba.png">

- open the Python notebook AINTBAD_code.ipynb (or IPA_code.ipynb). You need to check the path to your files.
  Run it by clicking on the play icon and allow the notebook to access your Google Drive files.
  While the code is executing, a series of questions will appear.

![image](https://github.com/Molinero-Group/volume-dispersion/assets/60315074/6bb6fc50-c672-4299-9abe-c311abdf79bd)

# Output data

After answering all of these questions, the AINTBAD and IPA code will generate in the same directory a few files in .pdf and .txt extensions.

# Examples used as an input

AINTBAD code requires an input file located in the directory called input that has to be located in the same directory as the code. The file has to be a text format .txt.

Examples are located inside the input directory. \
The input data used as examples for the codes are not our own, and were obtained from the following sources:\
Nm_bacteria.txt (Schwidetzky et al., 2021) DOI: [10.1021/acs.jpclett.1c03118](https://doi.org/10.1021/acs.jpclett.1c03118)\

# Acknowledgements 

The authors gratefully acknowledge support by AFOSR through MURI Award No. FA9550-20-1-0351. 

# License

The details on the license are shown in [LICENSE](https://github.com/Molinero-Group/underlying-distribution/blob/main/LICENSE).

