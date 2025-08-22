---
title: 'Lumabi: a Python package to streamline the computation of phonon-resolved luminescence spectra of defects and dopants in solids'
tags:
  - Python
  - materials science
  - defect physics
  - luminescence
  - density functional theory
  - vibronic peaks
authors:
  - name: Julien Bouquiaux
    orcid: 0000-0003-1982-052X
    corresponding: true
    affiliation: "1, 5"
  - name: Matteo Giantomassi
    orcid: 0000-0002-7007-9813
    affiliation: "1"
  - name: Samuel Poncé
    orcid: 0000-0003-1159-8389
    affiliation: "1, 2"
  - name: Yongchao Jia
    orcid: 0000-0002-7319-9546
    affiliation: "3"
  - name: Masayoshi Mikami
    orcid: 0000-0002-5598-4882
    affiliation: "4"
  - name: Xavier Gonze
    orcid: 0000-0002-8377-6829
    affiliation: "1"

affiliations:
 - name:  Institute of Condensed Matter and Nanosciences, Université catholique de Louvain, B-1348 Louvain-la-Neuve, Belgium
   index: 1
 - name:  WEL Research Institute, avenue Pasteur 6, 1300 Wavre, Belgium.
   index: 2
 - name: Yanshan University, Hebei Key Laboratory of Applied Chemistry, Yanshan University, 066004 Qinhuangdao, P. R. China
   index: 3
 - name: Science & Innovation Center, Mitsubishi Chemical Corporation, Yokohama 227-8502, Japan
   index: 4
 - name: Matgenix, A6K Advanced Engineering Centre, Charleroi, Belgium.
   index: 5

date: 07 August 2025
bibliography: paper.bib
---

# Summary

![Lumabi logo](Lumabi_logo.pdf){width=50%}

Lumabi is a Python package within the AbiPy framework [@gonze2020abinit] that automates the computation of phonon-resolved luminescence spectra of point defects and dopants in inorganic solids using the ABINIT first-principles software [@gonze2020abinit]. The package provides an end-to-end workflow: from $\Delta$SCF density-functional theory calculations with constrained occupations, to the generation of defect phonon modes in large supercells, right through the generation of luminescence spectra based on Huang–Rhys theory [@huang1950theory;@jin2021photoluminescence].  

Lumabi addresses the growing need for reproducible, automated workflows in defect physics [@lejaeghere2016reproducibility;@bosoni2024verify], with applications ranging from quantum technologies [@wolfowicz2021quantum;@dreyer2018first] to phosphors for solid-state lighting [@pust2015revolution;@lin2017inorganic;@fang2022evolutionary]. Tutorials and examples are available in the [AbiPy Book](https://abinit.github.io/abipy_book/lumabi/intro/intro.html).

# Statement of need

Defect-induced luminescence plays a key role in materials design for optoelectronics, quantum information, and phosphor technologies. Accurate predictions require ground- and excited-state calculations, phonon computations, and multiple post-processing steps, which are typically laborious to set up.  

Existing tools focus either on defect energetics [@naik2018coffee;@pean2017presentation;@goyal2017computational;@broberg2018pycdt;@kumagai2021insights;@neilson2022defap;@arrigoni2021spinney;@shen2024pymatgen] or luminescence post-processing [@Kavanagh2024;@turiansky2021nonrad;@cavignac2024], and most are tied to the commercial VASP software [@kresse1996efficiency]. To our knowledge, no open-source package has provided a fully automated pipeline for computing defect phonon modes in large supercells together with luminescence spectra.  

Lumabi aims at filling this gap. Built on ABINIT and AbiPy, with interfaces to Phonopy [@togo2015first;@togo2023first] and Pymatgen [@ong2013python], it streamlines the entire workflow. It enables reproducible simulations with limited intervention and produces structured data suitable for data-driven searches [@hariyani2023guide] and machine learning [@lee2025machine].

# Software Description, Features, and Computational Workflow

The code is organized into four modules that can be combined into a complete workflow or used independently. We describe here the overall working principles of each module. For a more practical approach, we provide [online tutorials](https://abinit.github.io/abipy_book/lumabi/lumiwork/lesson_lumiwork.html).

## LumiWork Module

![The LumiWork module, an AbiPy Workflow that automates ABINIT DFT tasks with $\Delta$SCF constrained occupations.](LumiWork.pdf)

A computational workflow for calculating phonon-resolved photoluminescence (PL) spectra of defect systems starts with the LumiWork module, which automates ABINIT DFT tasks with $\Delta$SCF constrained occupations [@jones1989density;@hellman2004potential]. Users provide the defect supercell structure, the DFT input parameters, and constrained occupations of the Kohn-Sham states designed to mimic the excited state of the system under study. This module manages two structural relaxations for
the ground- and the excited-state, and offers optional static SCF computations followed by non-SCF band structure calculations. As the relaxed excited state is not known in advance, input files are generated dynamically.

## $\Delta$SCF Post-Processing Module

![The $\Delta$SCF module, designed to post-process $\Delta$SCF constrained-occupation calculations using a one-dimensional configuration-coordinate model.](dSCF_post_process.pdf)

The next step in the workflow is handled by the $\Delta$SCF post-processing module. This tool takes the NetCDF output files generated by the previous LumiWork module, and processes them following a one-dimensional configuration-coordinate model [@jia2017first;@bouquiaux2021importance]. This analysis provides insights into the luminescence characteristics of the defect under study by computing properties such as transition energies, Huang-Rhys factors, effective phonon frequencies, and lineshapes following this 1D model or within a semi-classical approximation.
It also facilitates the analysis of atomic relaxations by, for example, automatically generating VESTA [@momma2011vesta] files that include relaxation vectors.

## IFCs Embedding Module

![The IFCs embedding module, allowing to calculate defect phonons in large supercells.](IFCs_embedding.pdf)

The Interatomic Force Constants (IFCs) Embedding module enables defect phonon calculations in large supercells, which are otherwise computationally prohibitive with standard density-functional perturbation theory or finite differences approach. The method combines short-range defect force constants, obtained in a small supercell, with pristine host force constants computed from the bulk and folded into a large supercell. The resulting embedded IFC matrix captures both localized defect modes and host phonons, allowing accurate spectral simulations at dilute defect concentrations. The implementation interfaces with Phonopy and produces phonon objects compatible with later analysis. First employed in the context of the luminescence of the NV center by Alkauskas et al. [@alkauskas2014], this embedding approach has been then used in various materials [@jin2021photoluminescence;@bouquiaux2023first;@razinkovas2021vibrational;@maciaszek2023application;@jin2022vibrationally]. For the mathematical details and the technical implementation, we refer to the accompanying [Jupyter Book](https://abinit.github.io/abipy_book/lumabi/theory/lesson_theory.html).

## Lineshape Calculation Module

![The lineshape module, allowing to compute the temperature-dependent spectra.](lineshape.pdf)

As a final step, the Lineshape module is used. The main task of this module is to compute the Huang-Rhys spectral function $S(\hbar\omega)$ [@alkauskas2014; @bouquiaux2023first] and to generate temperature-dependent PL spectra using the efficient generating function approach [@jin2021photoluminescence].

The code takes as input the zero-phonon line energy, the atomic displacements $\Delta R_{\kappa\alpha}$ or forces $\Delta F_{\kappa\alpha}$ induced by the electronic transition (obtained from the $\Delta$SCF post-processing step), and the phonon modes provided as a Phonopy object (potentially obtained from the IFCs embedding module). Notice that the use of the displacements is only compatible if the phonon supercell is of the same size as the $\Delta$SCF supercell. The use of the forces (equivalent under the harmonic approximation) allows one to use efficiently the previous block and enlarge the supercell size, ensuring a good convergence of the Huang-Rhys spectral function [@alkauskas2014;@jin2021photoluminescence;@bouquiaux2023first]. An analysis of the different phonon mode localization can also be performed.

# Examples and Applications

This computational workflow has been used for inorganic phosphors activated with Eu$^{2+}$ dopants and has also been tested on a variety of other systems including  F-centers (oxygen vacancy) in CaO and the NV center in diamond. Its versatility allows for any kind of point defect.
These developments have been particularly useful in understanding the luminescence properties of technologically significant red-emitting Eu-doped phosphor materials. Notably, the workflow has been applied to SrAl$_2$Li$_2$O$_2$N$_2$:Eu$^{2+}$ and SrLiAl$_3$N$_4$:Eu$^{2+}$, shedding new light on their phonon sideband [@bouquiaux2021importance;@bouquiaux2023first]. We refer the reader to the accompanying notebook tutorials for practical examples demonstrating the application of this workflow.


# Acknowledgements

J.B. acknowledge support from the F.R.S.-FNRS.
S.P. is a Research Associate of the Fonds de la Recherche Scientifique - FNRS.
X.G. also acknowledges support from the F.R.S.-FNRS, under grant n°T.0103.19 (PDR - ALPS).
S.P. also acknowledges support from the Walloon Region in the strategic axe FRFS-WEL-T.
This work was supported by the Communauté française de Belgique through the SURFASCOPE project (ARC 19/24-102).
Y.J acknowledges the funding support from National Key R\&D Program (No.2022YFB3503800), Natural Science Foundation of Hebei Province (No.E2021203126) and Cultivation Project for Basic Research and Innovation of Yanshan University (No.2021LGQN033).
Computational resources have been provided by the Consortium des Équipements de Calcul Intensif (CÉCI), funded by the FRS-FNRS under Grant No. 2.5020.11 and by the Walloon Region.

# References