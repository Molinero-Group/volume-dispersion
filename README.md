ADD ZENODO LINK
<a href="https://zenodo.org/badge/latestdoi/559044289"><img src="https://zenodo.org/badge/559044289.svg" alt="DOI"></a>

Table of contents
=================

<!--ts-->
   * [Purpose of the codes](#purpose-of-the-codes)
   * [How to cite](#how-to-cite)
   * [How to Install the Project](#how-to-install-the-project)
      * [Running in a local installation](#running-in-a-local-installation)
      * [Running on Google drive](#running-on-google-drive)
   * [How to Run the Project](#how-to-run-the-project)
   * [Output data](#output-data)
   * [Examples used as an input](#examples-used-as-an-input)
   * [Acknowledgements](#acknowledgements)
   * [License](#license)
<!--te-->

# Purpose of the codes

The AINTBAD (Analysis of Ice nucleation Temperature for $B$ and $A$ Determination) includes the numerical integration in Equation 7 and the analytical model of Equation 8 in an python code to estimate $A$, $B$ and $J$ from experimental drop-freezing data. The code outputs the parameters $A$ and $B$. These are used to compute the nucleation barriers $\Delta G$, the temperature that corresponds to 50\% of frozen droplets $T_{50}$, and the homogeneous nucleation rate evaluated at $T_{50}$.

The code IPA (Inhomogeneous Poisson Analysis) is capable of taking various parametrizations for the homogeneous nucleation rate $J_{hom}(T)$, the droplet size distributions (Gaussian, Gamma, uniform, exponential, etc.), and cooling rates to compute the survival probability or fraction of frozen droplets. We use the nucleation rate data vs temperature as the input to compute the survival probability.

# How to cite

This code is the active version of the code archived here and documented in Addula, de Almeida Ribeiro, Molinero & Peters, ["Ice nucleation from drop-freezing experiments: Impact of droplet volume dispersion and cooling rates"](https://chemrxiv.org/engage/chemrxiv/article-details/63691c92b58850396f407923), DOI: 10.26434/chemrxiv-2022-ddzv8 

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

- download the directory HUB-forward (or HUB-backward), and then upload it in a **Google drive** directory;

- inside the directory, click on "add more apps",

![image](https://user-images.githubusercontent.com/60315074/199114124-7ffb328d-dd1d-44f4-9d04-246bba6f7538.png)

- type *colaboratory* and install it

<img width="969" alt="image" src="https://user-images.githubusercontent.com/60315074/184701819-8e4baaf3-f2b7-47d2-b067-bb3947151fba.png">

- open the Python notebook HUB-forward.ipynb (or HUB-backward.ipynb)

<img width="848" alt="image" src="https://user-images.githubusercontent.com/60315074/199114807-10b687f2-4221-444b-a48f-0b299eb307fb.png">

# How to Run the Project

- run it by clicking on the play icon and allow the notebook to access your Google Drive files. While the code is executing, a series of questions will appear,

<img width="613" alt="image" src="https://user-images.githubusercontent.com/60315074/199116848-1bd67052-2f33-456f-8005-598beb18a1c0.png">

# Output data

After answering all of these questions, the HUB-forward code will generate in the same directory a few files in .pdf and .txt extensions: the underlying and the concentration-dependent distributions, the fraction of frozen droplets $f_{ice} (T)$ and the cumulative freezing spectrum $N_m (T)$.

After answering a series of questions, HUB-backward will generate a few files in .pdf and .txt extensions: a spline fit of the input data, the optimized cumulative freezing spectrum $N_m (T)$ data, and the optimized differential freezing spectrum $n_m (T)$.

# Examples used as an input

HUB-backward requires an input file located in the directory called input that has to be located in the same directory as the code. The file has to be a text format .txt, the first column is temperature, and the second column is $N_m (T)$ or $f_{ice} (T)$. 

Examples are located inside the input directory. \
The input data used as examples for the HUB-backward code are not our own, and were obtained from the following sources:\
Nm_bacteria.txt (Schwidetzky et al., 2021) DOI: [10.1021/acs.jpclett.1c03118](https://doi.org/10.1021/acs.jpclett.1c03118)\
Nm_fusarium_kunert2019_strain_3-68.txt (Kunert et al., 2019) DOI: [10.5194/bg-16-4647-2019](https://doi.org/10.5194/bg-16-4647-2019)\
Nm_pollen_thesis.txt (Dreischmeier, 2019) DOI: [10.4119/unibi/2907691](https://doi.org/10.4119/unibi/2907691)\
Nm_pH_6p2.txt, Nm_pH_5p6.txt and Nm_pH_4p4.txt (Lukas et al., 2020) DOI: [10.1021/jacs.9b13069](https://doi.org/10.1021/jacs.9b13069)\
fractionofice_cholesterol_fig10.txt and fractionofice_cholesterol_fig13.txt (Zhang and Maeda, 2022) DOI: [10.1016/j.ces.2022.118017](https://doi.org/10.1016/j.ces.2022.118017)

# Acknowledgements 

I. de A. R and V. M. gratefully acknowledge support by AFOSR through MURI Award No. FA9550-20-1-0351. K. M. acknowledges support by the National Science Foundation under Grant No. (NSF 2116528) and from the Institutional Development Awards (IDeA) from the National Institute of General Medical Sciences of the National Institutes of Health under Grants #P20GM103408, P20GM109095.

# License

The details on the license are shown in [LICENSE](https://github.com/Molinero-Group/underlying-distribution/blob/main/LICENSE).

