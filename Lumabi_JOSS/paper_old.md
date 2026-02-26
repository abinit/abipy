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
 - name:  Institute of Condensed Matter and Nanosciences, Université catholique de Louvain, B-1348 Louvain-la-Neuve, Belgium
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

Lumabi is a Python package integrated within the AbiPy framework [@gonze2020abinit] designed to automate and streamline the computation of phonon-resolved luminescence spectra of point defects (including the case of dopants) in inorganic solids using the ABINIT first-principles software application [@gonze2002first;@gonze2009abinit;@gonze2016recent;@gonze2020abinit]. The package addresses the growing need for efficient, reproducible workflows in materials science [@lejaeghere2016reproducibility;@bosoni2024verify], particularly in the study of defect-related luminescent properties, which are critical for applications ranging from quantum technologies [@wolfowicz2021quantum;@dreyer2018first] to down-conversion phosphors materials used in white light LEDs [@pust2015revolution;@lin2017inorganic;@fang2022evolutionary]. Lumabi automates key steps in the computational workflow, from initial $\Delta$SCF density-functional theory calculations with constrained occupations, to the generation of defect phonon modes in large supercells, right through the final generation of luminescence spectra based on the Huang-Rhys theory [@huang1950theory;@jin2021photoluminescence].
A collection of tutorials and examples, organized as Jupyter Books,
is available on the [Abipy Book project webpage](https://abinit.github.io/abipy_book/lumabi/intro/intro.html).

# Statement of need

The study of defect-induced luminescence in materials is crucial for understanding and designing materials with specific optical properties. However, the computational workflow required to accurately predict these properties is complex and involves ground- and excited-state calculations,
phonon computations, as well as pre- and post-processing steps.

Therefore, a number of software packages have been developed to manage the pre- or post-processing of defect calculations [@naik2018coffee;@pean2017presentation;@goyal2017computational;@broberg2018pycdt;@kumagai2021insights;@neilson2022defap;@arrigoni2021spinney;@shen2024pymatgen]. The few that focus on luminescent properties [@Kavanagh2024;@turiansky2021nonrad;@cavignac2024] are all interfaced with the commercial software VASP [@kresse1996efficiency], and only provide the post-processing part of DFT computations following the formalism proposed by Alkauskas et al. [@alkauskas2014]. To our knowledge, the generation of defect phonon modes in large supercells following an embedding of the interatomic force constants is currently not available.

Divided in four main Python modules, Lumabi provides a solution to automate all the necessary tasks required to compute phonon-resolved luminescence spectra of defects. The density functional theory (DFT) computations are performed with the open-source ABINIT software [@gonze2002first;@gonze2009abinit;@gonze2016recent;@gonze2020abinit]. The generation of ABINIT input files, automatic workflow management, and post-processing  are done with the AbiPy package, which is interfaced with Pymatgen [@ong2013python]. Lumabi is also interfaced with the Phonopy package [@togo2015first;@togo2023first] for the calculations of phonons.
This tool is intended for researchers aiming to perform accurate defect luminescence simulations with minimal human intervention.
Moreover, the systematic generation of well-structured data through this automated workflow could enable data-driven searches for new phosphors [@hariyani2023guide] and could provide a robust foundation for training machine-learning models [@lee2025machine] to accelerate the discovery of materials with targeted luminescent properties.

# Software Description, Features, and Computational Workflow

The code is structured around four Python modules, each handling a different stage of the computational workflow.
Although these modules are designed to work seamlessly together --where the output of one serves as the input to the next--
each individual building block can be also used independently with ease. In this article, we describe here the overall working principles of each block. For a more practical approach, we provide [online tutorials](https://abinit.github.io/abipy_book/intro.html) for the different modules.

## LumiWork Module

![The LumiWork module, an AbiPy Workflow that automates ABINIT DFT tasks with $\Delta$SCF constrained occupations.](LumiWork.pdf)

A computational workflow for calculating phonon-resolved photoluminescence (PL) spectra of defect systems starts with the LumiWork module, which automates ABINIT DFT tasks with $\Delta$SCF constrained occupations [@jones1989density;@hellman2004potential]. Users provide the defect supercell structure (e.g. in CIF format), the DFT input parameters (e.g., the plane-wave kinetic energy cutoff, the stopping SCF criterion, etc.), and constrained occupations of the Kohn-Sham states designed to mimic the excited state of the system under study. This module manages two structural relaxations for
the ground- and the excited-state, and offers optional static SCF computations followed by non-SCF band structure calculations. Note that the workflow is dynamic: since the relaxed excited state is not known beforehand, input files must be generated at run-time.

The `LumiWork` class inherits from the Abipy's `Work` class.
Consequently, it leverages error detection and automatic restart handling capabilities available for all ABINIT-related Works,
ensuring that most computations proceed smoothly without interruption.
The majority of outputs are in NetCDF format, serving as the foundation for subsequent modules

## $\Delta$SCF Post-Processing Module

![The $\Delta$SCF module, designed to post-process $\Delta$SCF constrained-occupation calculations using a one-dimensional configuration-coordinate model.](dSCF_post_process.pdf)

The next step in the workflow is handled by the $\Delta$SCF post-processing module. This tool takes the NetCDF output files generated by the previous LumiWork module, and processes them following a one-dimensional configuration-coordinate model [@jia2017first;@bouquiaux2021importance]. This analysis provides detailed insights into the luminescence characteristics of the defect under study by computing properties such as transition energies, Huang-Rhys factors, effective phonon frequencies, and lineshapes following this 1D model or within a semi-classical approximation.
It also facilitates the analysis of atomic relaxations by, for example, automatically generating VESTA [@momma2011vesta] files that include relaxation vectors.

In cases where experimental lineshapes exhibit resolved phonon peaks,
it is important to note that the lineshapes produced by this tool *are not* supposed to closely match experimental data at the level of the phonon peaks.
Indeed, the model accounts for only a single effective phonon frequency $\omega_{\rm{eff}}$,
and improved agreement with experiment requires the use of subsequent modules.
Nevertheless, the 1D model provides the overall spectrum shape, or the total Huang-Rhys factor, which can be compared to experiment.

## IFCs Embedding Module

![The IFCs embedding module, allowing to calculate defect phonons in large supercells.](IFCs_embedding.pdf)

The Interatomic Force Constants (IFCs) Embedding module computes the phonon of defected materials in large supercells at the $\mathbf{q}=(0,0,0)$ ($\Gamma$ point) wave vector. This module addresses the computational challenges associated with capturing both the coupling with long-wavelength phonons and the localized phonon modes introduced by defects. Direct methods, such as Density Functional Perturbation Theory (DFPT) [@gonze1997dynamical] and finite difference approaches [@togo2015first], are impractical for large supercells due to their high computational costs. First employed in the context of the luminescence of the NV center by Alkauskas et al. [@alkauskas2014], the embedding approach has been then used in various materials [@jin2021photoluminescence;@bouquiaux2023first;@razinkovas2021vibrational;@maciaszek2023application;@jin2022vibrationally]. The procedure relies on the short-range nature of interatomic forces to build IFCs matrices for defect supercells containing thousands of atoms, targeting the dilute limit.

First, the real-space IFCs for the defect system $C_{\kappa\alpha,\kappa'\alpha'}^{\mathrm{defect}}$ are computed within a relatively small supercell. Typically, this is achieved using a finite-difference approach at the $\Gamma$ point using Phonopy.

Second, the IFCs of the pristine system $C_{\kappa\alpha,\kappa'\alpha'}^{\mathrm{pristine}}$ are computed in a much larger supercell.
A practical approach consists in using DFPT on the unit cell to compute the dynamical matrix
on a regular $\mathbf{q}$-point mesh. This is often followed by Fourier interpolation to obtain a denser $\mathbf{q}$-mesh, equivalent to a large supercell.
This folding procedure, which yields $C_{\kappa\alpha,\kappa'\alpha}^{\mathrm{pristine}}$ has been implemented.
It takes as input an ABINIT DDB file (or a phonopy object) and produces as output a Phonopy object containing the folded pristine phonons.

Third, an embedded IFC matrix $C_{\kappa\alpha,\kappa'\alpha'}^{\mathrm{emb}}$ is constructed using the following rules: if both atoms $\kappa$ and $\kappa'$ are within a sphere centered around the defect, with a cut-off radius $R_c$, then $C_{\kappa\alpha,\kappa'\alpha'}^{\mathrm{emb}}$ is set to $C_{\kappa\alpha,\kappa'\alpha'}^{\mathrm{defect}}$. For all other atomic pairs, the embedded IFC is set to the pristine value: $C_{\kappa\alpha,\kappa'\alpha'}^{\mathrm{emb}} = C_{\kappa\alpha,\kappa'\alpha'}^{\mathrm{pristine}}$. The cut-off radius $R_c$ is typically determined by the size of the initial supercell used to compute the defect IFCs. Additionally, one can use a second cut-off radius $R_b$. If two atoms are not within a sphere of radius $R_b$, the IFC is set to zero. This allows one to obtain a sparse matrix, for which faster diagonalization techniques may be used, although the current implementation does not natively support such technique. Finally, the acoustic sum rule is reinforced following the conventions of Ref. [@jin2021photoluminescence].

The code performs a mapping between the structures of the pristine and defect computations. Such mapping is not trivial to automate
because the atomic positions of the defect supercell might not be the same as the pristine positions. Also, vacancies or interstitials should be handled with care during the mapping. For this reason, the mapping is done in a semi-automatic way, meaning that the user is responsible for the input of the different structure manipulations needed to go from the pristine cell to the defect cell. For details about this mapping, several examples are provided with the code. At the end, a Phonopy object containing the embedded IFCs and the new large defect supercell structure is created.
For the treatment of Born effective charges (BECs), if the BEC tensor of the substituted atom was computed in the defect phonon calculation, it is replaced by the BEC of the dopant atom. If not, one can choose to use the most common oxidation state of the dopant atom as its new BEC tensor that is now assumed to be diagonal. The same approach applies to interstitials, but in this case, the dopant atom and its corresponding BEC are simply added. For vacancies, the atom and its BEC are removed.

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
