# Awesome Python Chemistry [![Awesome](https://cdn.rawgit.com/sindresorhus/awesome/d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)](https://github.com/sindresorhus/awesome)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

A curated list of awesome Python frameworks, libraries, software and resources related to Chemistry.

Inspired by [awesome-python](https://awesome-python.com).

## Table of contents

- [Awesome Python Chemistry ![Awesome](https://github.com/sindresorhus/awesome)](#awesome-python-chemistry-)
  - [General Chemistry](#general-chemistry)
  - [Machine Learning](#machine-learning)
  - [Generative Molecular Design](#generative-molecular-design)
  - [Simulations](#simulations)
  - [Molecular Visualization](#molecular-visualization)
  - [Database Wrappers](#database-wrappers)
  - [Learning Resources](#learning-resources)
  - [See Also](#see-also)

---

## General Chemistry

*Packages and tools for general chemistry.*

- [AQME](https://github.com/jvalegre/aqme) - Ensemble of automated QM workflows that can be run through jupyter notebooks, command lines and yaml files.
- [aizynthfinder](https://github.com/MolecularAI/aizynthfinder) - A tool for retrosynthetic planning.
- [batchcalculator](http://lukaszmentel.com/batchcalculator/) - A GUI app based on wxPython for calculating the correct amount of reactants (batch) for a particular composition given by the molar ratio of its components.
- [cctbx](https://cctbx.github.io/) - The Computational Crystallography Toolbox.
- [ChemFormula](https://github.com/molshape/ChemFormula) - ChemFormula provides a class for working with chemical formulas. It allows parsing chemical formulas, calculating formula weights, and generating formatted output strings (e.g. in HTML, LaTeX, or Unicode).
- [chemlib](https://chemlib.readthedocs.io/en/latest/) - A robust and easy-to-use package that solves a variety of chemistry problems.
- [chempy](http://pythonhosted.org/chempy/) - ChemPy is a package useful for chemistry (mainly physical/inorganic/analytical chemistry).
- [datamol](https://github.com/datamol-org/datamol): - Molecular Manipulation Made Easy. A light wrapper build on top of RDKit.
- [GoodVibes](https://github.com/bobbypaton/GoodVibes) - A Python program to compute quasi-harmonic thermochemical data from Gaussian frequency calculations.
- [hgraph2graph](https://github.com/wengong-jin/hgraph2graph) - Hierarchical Generation of Molecular Graphs using Structural Motifs.
- [ionize](http://lewisamarshall.github.io/ionize/) - Calculates the properties of individual ionic species in aqueous solution, as well as aqueous solutions containing arbitrary sets of ions.
- [LModeA-nano](https://lmodea-nano.readthedocs.io/en/latest/) - Calculates the intrinsic chemical bond strength based on local vibrational mode theory in solids and molecules.
- [mendeleev](http://mendeleev.readthedocs.io/en/stable/) - A package that provides a python API for accessing various properties of elements from the periodic table of elements.
- [nmrglue](https://github.com/jjhelmus/nmrglue) - A package for working with nuclear magnetic resonance (NMR) data including functions for reading common binary file formats and processing NMR data.
- [NistChemPy](https://github.com/IvanChernyshov/NistChemPy) - A package for accessing data from the NIST webbook. API includes access to thermodynamic properties, molecular structures and more.
- [Open Babel](http://openbabel.org/) - A chemical toolbox designed to speak the many languages of chemical data.
- [periodictable](https://github.com/pkienzle/periodictable) - This package provides a periodic table of the elements with support for mass, density and xray/neutron scattering information.
- [propka](https://github.com/jensengroup/propka) - Predicts the pKa values of ionizable groups in proteins and protein-ligand complexes based in the 3D structure.
- [pybaselines](https://pybaselines.readthedocs.io/en/latest/) - A package for fitting baselines of spectra for baseline correction.
- [pybel](https://openbabel.org/docs/UseTheLibrary/Python_Pybel.html) - Pybel provides convenience functions and classes that make it simpler to use the Open Babel libraries from Python.
- [pycroscopy](https://pycroscopy.github.io/pycroscopy/index.html) - Scientific analysis of nanoscale materials imaging data.
- [pyEQL](https://pyeql.readthedocs.io/en/latest/index.html) - A set of tools for conventional calculations involving solutions (mixtures) and electrolytes.
- [pyiron](http://pyiron.org/) - pyiron - an integrated development environment (IDE) for computational materials science.
- [pymatgen](http://pymatgen.org) - Python Materials Genomics is a robust, open-source library for materials analysis.
- [pymatviz](https://github.com/janosh/pymatviz) - A toolkit for visualizations in materials informatics.
- [symfit](https://symfit.readthedocs.io/en/stable/) - a curve-fitting library ideally suited to chemistry problems, including fitting experimental kinetics data.
- [symmetry](http://pythonhosted.org/symmetry/) - Symmetry is a library for materials symmetry analysis.
- [stk](https://github.com/lukasturcani/stk) - A library for building, manipulating, analyzing and automatic design of molecules, including a genetic algorithm.
- [spectrochempy](https://github.com/spectrochempy/spectrochempy) - A library for processing, analyzing and modeling spectroscopic data.

## Machine Learning

*Packages and tools for employing machine learning and data science in chemistry.*

- [amp](http://amp.readthedocs.io/en/latest/) - Is an open-source package designed to easily bring machine-learning to atomistic calculations.
- [atom3d](https://github.com/drorlab/atom3d) - Enables machine learning on three-dimensional molecular structure.
- [chainer-chemistry](https://github.com/chainer/chainer-chemistry) - A deep learning framework (based on Chainer) with applications in Biology and Chemistry.
- [chemml](https://hachmannlab.github.io/chemml/) - A machine learning and informatics program suite for the analysis, mining, and modeling of chemical and materials data.
- [chemprop](https://github.com/chemprop/chemprop) - Message Passing Neural Networks for Molecule Property Prediction .
- [cgcnn](https://github.com/txie-93/cgcnn) - Crystal graph convolutional neural networks for predicting material properties.
- [deepchem](http://deepchem.io/) - Deep-learning models for Drug Discovery and Quantum Chemistry.
- [DeepPurpose](https://github.com/kexinhuang12345/DeepPurpose) - A Deep Learning Library for Compound and Protein Modeling DTI, Drug Property, PPI, DDI, Protein Function Prediction.
- [DescriptaStorus](https://github.com/bp-kelley/descriptastorus) - Descriptor computation (chemistry) and (optional) storage for machine learning.
- [DScribe](https://github.com/SINGROUP/dscribe) - Descriptor library containing a variety of fingerprinting techniques, including the Smooth Overlap of Atomic Positions (SOAP).
- [graphein](https://github.com/a-r-j/graphein) - Provides functionality for producing geometric representations of protein and RNA structures, and biological interaction networks.
- [Matminer](https://github.com/hackingmaterials/matminer) - Library of descriptors to aid in the data-mining of materials properties, created by the Lawrence Berkeley National Laboratory.
- [MoleOOD](https://github.com/yangnianzu0515/MoleOOD) - a robust molecular representation learning framework against distribution shifts.
- [megnet](https://github.com/materialsvirtuallab/megnet) - Graph Networks as a Universal Machine Learning Framework for Molecules and Crystals.
- [MAML](https://github.com/materialsvirtuallab/maml) - Aims to provide useful high-level interfaces that make ML for materials science as easy as possible.
- [MORFEUS](https://github.com/kjelljorner/morfeus) - Library for fast calculations of **mo**lecula**r** **fe**at**u**re**s** from 3D structures for machine learning with a focus on steric descriptors.
- [olorenchemengine](https://github.com/Oloren-AI/olorenchemengine) - Molecular property prediction with unified API for diverse models and respresentations,
  with integrated uncertainty quantification, interpretability, and hyperparameter/architecture tuning.
- [ROBERT](https://github.com/jvalegre/robert) - Ensemble of automated machine learning protocols that can be run sequentially through a single command line. The program works for regression and classification problems.
- [schnetpack](https://github.com/atomistic-machine-learning/schnetpack) - Deep Neural Networks for Atomistic Systems.
- [selfies](https://github.com/aspuru-guzik-group/selfies) - Self-Referencing Embedded Strings (SELFIES): A 100% robust molecular string representation.
- [Summit](https://github.com/sustainable-processes/summit) - Package for optimizing chemical reactions using machine learning (contains 10 algorithms + several benchmarks).
- [TDC](https://github.com/mims-harvard/TDC) - Therapeutics Data Commons (TDC) is the first unifying framework to systematically access and evaluate machine learning across the entire range of therapeutics.
- [XenonPy](https://github.com/yoshida-lab/XenonPy) - Library with several compositional and structural material descriptors, along with a few pre-trained neural network models of material properties.

## Generative Molecular Design

*Packages and tools for generating molecular species*

- [GraphINVENT](https://github.com/MolecularAI/GraphINVENT) - A platform for graph-based molecular generation using graph neural networks.
- [GuacaMol](https://github.com/BenevolentAI/guacamol) - A package for benchmarking of models for _de novo_ molecular design.
- [moses](https://github.com/molecularsets/moses) - A benchmarking platform for molecular generation models.
- [perses](https://github.com/choderalab/perses) - Experiments with expanded ensembles to explore chemical space.

## Simulations

*Packages for atomistic simulations and computational chemistry.*

- [alchemlyb](https://github.com/alchemistry/alchemlyb) - Makes alchemical free energy calculations easier by leveraging the full power and flexibility of the PyData stack.
- [atomate2](https://github.com/materialsproject/atomate2) - atomate2 is a library of computational materials science workflows.
- [Atomic Silumation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/index.html) - Is a set of tools and modules for setting up, manipulating, running, visualizing and analyzing atomistic simulations.
- [basis_set_exchange](https://github.com/MolSSI-BSE/basis_set_exchange) - A library containing basis sets for use in quantum chemistry calculations. In addition, this library has functionality for manipulation of basis set data.
- [CACTVS](https://www.xemistry.com/academic/) - Cactvs is a universal, scriptable cheminformatics toolkit, with a large collection of modules for property computation, chemistry data file I/O and other tasks.
- [CalcUS](https://github.com/cyllab/CalcUS) - Quantum chemisttry web platform that brings all the necessary tools to perform quantum chemistry in a user-friendly web interface.
- [cantera](https://github.com/Cantera/cantera) - A collection of object-oriented software tools for problems involving chemical kinetics, thermodynamics, and transport processes.
- [CatKit](https://github.com/SUNCAT-Center/CatKit) - General purpose tools for high-throughput catalysis.
- [ccinput](https://github.com/cyllab/ccinput/) - A tool and library for creating quantum chemistry input files.
- [cclib](https://cclib.github.io/) - A library for parsing output files various quantum chemical programs.
- [cinfony](http://cinfony.github.io/) - A common API to several cheminformatics toolkits (Open Babel, RDKit, the CDK, Indigo, JChem, OPSIN and cheminformatics webservices).
- [chemlab](http://chemlab.readthedocs.io/en/latest/index.html) - Is a library that can help the user with chemistry-relevant calculations.
- [emmet](https://github.com/materialsproject/emmet) - A package to 'build' collections of materials properties from the output of computational materials calculations.
- [fromage](https://github.com/Crespo-Otero-group/fromage/) - The "FRamewOrk for Molecular AGgregate Excitations" enables localised QM/QM' excited state calculations in a solid state environment.
- [GPAW](https://wiki.fysik.dtu.dk/gpaw/) - Is a density-functional theory (DFT) Python code based on the projector-augmented wave (PAW) method and the atomic simulation environment (ASE).
- [horton](http://theochem.github.io/horton/2.0.1/index.html) - Helpful Open-source Research TOol for N-fermion system, a quantum-chemistry program that can perform computations involving model Hamiltonians.
- [HTMD](https://github.com/Acellera/htmd) - High-Throughput Molecular Dynamics: Programming Environment for Molecular Discovery.
- [Indigo](https://github.com/epam/Indigo) - Universal cheminformatics libraries, utilities and database search tools.
- [IoData](https://github.com/theochem/iodata) - File parser/converter for QM, MD and plane-wave DFT programs.
- [Jarvis-tools](https://github.com/usnistgov/jarvis) - An open-access software package for atomistic data-driven materials design
- [mathchem](http://mathchem.iam.upr.si/) - Is a free open source package for calculating topological indices and other invariants of molecular graphs.
- [MDAnalysis](http://www.mdanalysis.org/) - Is an object-oriented library to analyze trajectories from molecular dynamics (MD) simulations in many popular formats.
- [MDTraj](http://mdtraj.org) - Package for manipulating molecular dynamics trajectories with support for multiple formats.
- [MMTK](http://dirac.cnrs-orleans.fr/MMTK/) - The Molecular Modeling Toolkit is an Open Source program library for molecular simulation applications.
- [MolMod](http://molmod.github.io/molmod/index.html) - A library with many components that are useful to write molecular modeling programs.
- [nmrsim](https://nmrsim.readthedocs.io/en/latest/index.html) - A library for simulating first- or second-order NMR spectra and dynamic NMR resonances.
- [oddt](https://github.com/oddt/oddt) - Open Drug Discovery Toolkit, a modular and comprehensive toolkit for use in cheminformatics, molecular modeling etc.
- [OPEM](https://github.com/ECSIM/opem) - Open source PEM (Proton Exchange Membrane) fuel cell simulation tool.
- [openmmtools](https://github.com/choderalab/openmmtools) - A batteries-included toolkit for the GPU-accelerated OpenMM molecular simulation engine.
- [overreact](https://github.com/geem-lab/overreact) - A library and command-line tool for building and analyzing complex homogeneous microkinetic models from quantum chemistry calculations, with support for quasi-harmonic thermochemistry, quantum tunnelling corrections, molecular symmetries and more.
- [ParmEd](https://github.com/ParmEd/ParmEd) -  Parameter/topology editor and molecular simulator with visualization capability.
- [pGrAdd](https://github.com/VlachosGroup/PythonGroupAdditivity) - A library for estimating thermochemical properties of molecules and adsorbates using group additivity.
- [phonopy](http://atztogo.github.io/phonopy/) - An open source package for phonon calculations at harmonic and quasi-harmonic levels.
- [PLAMS](https://github.com/SCM-NV/PLAMS) - Python Library for Automating Molecular Simulation: input preparation, job execution, file management, output processing and building data workflows.
- [pMuTT](https://vlachosgroup.github.io/pMuTT/) - A library for ab-initio thermodynamic and kinetic parameter estimation.
- [PorePy](https://github.com/pmgbergen/porepy) - A Simulation Tool for Fractured and Deformable Porous Media.
- [ProDy](http://prody.csb.pitt.edu/) - An open source package for protein structural dynamics analysis with a flexible and responsive API.
- [ProLIF](https://github.com/chemosim-lab/ProLIF) - Interaction Fingerprints for protein-ligand complexes and more.
- [Psi4](http://psicode.org) - A hybrid Python/C++ open-source package for quantum chemistry.
- [Psi4NumPy](https://github.com/psi4/psi4numpy/) - Psi4-based reference implementations and Jupyter notebook-based tutorials for foundational quantum chemistry methods.
- [pyEMMA](http://www.emma-project.org/latest/) - Library for the estimation, validation and analysis Markov models of molecular kinetics and other kinetic and thermodynamic models from molecular dynamics data.
- [pygauss](https://pygauss.readthedocs.io/en/stable/index.html) - An interactive tool for supporting the life cycle of a computational molecular chemistry investigations.
- [PyQuante](http://pyquante.sourceforge.net/) - Is an open-source suite of programs for developing quantum chemistry methods.
- [pysic](https://github.com/thynnine/pysic) - A calculator incorporating various empirical pair and many-body potentials.
- [Pyscf](https://github.com/sunqm/pyscf) - A quantum chemistry package written in Python.
- [pyvib2](http://pyvib2.sourceforge.net/) - A program for analyzing vibrational motion and vibrational spectra.
- [RDKit](http://www.rdkit.org/) - Open-Source Cheminformatics Software.
- [ReNView](https://github.com/VlachosGroup/ReNView/wiki/Reaction-Network-Viewer-(ReNView)-Usage-Instructions) - A program to visualize reaction networks.
- [stk](https://github.com/lukasturcani/stk) - A library for building, manipulating, analyzing and automatic design of molecules.
- [QMsolve](https://github.com/quantum-visualizations/qmsolve) - A module for solving and visualizing the Schrödinger equation.
- [QUIP](http://libatoms.github.io/QUIP/) - A collection of software tools to carry out molecular dynamics simulations.
- [torchmd](https://github.com/torchmd/torchmd) - End-To-End Molecular Dynamics (MD) Engine using PyTorch.
- [tsase](http://theory.cm.utexas.edu/tsase/) - The library which depends on ASE to tackle transition state calculations.
- [yank](https://github.com/choderalab/yank) - An open, extensible Python framework for GPU-accelerated alchemical free energy calculations.

### Force Fields

*Packages related to force fields*

- [acpype](https://github.com/alanwilter/acpype) - Convert AMBER forcefields from ANTECHAMBER to GROMACS format.
- [CHGNet](https://github.com/CederGroupHub/chgnet) - Pretrained universal neural network potential for charge-informed atomistic modeling.
- [FitSNAP](https://github.com/FitSNAP/FitSNAP) - A Package For Training SNAP Interatomic Potentials for use in the LAMMPS molecular dynamics package.
- [fftool](https://github.com/paduagroup/fftool) - Tool to build force field input files for molecular simulation.
- [FLARE](https://github.com/mir-group/flare) - A package for creating fast and accurate interatomic potentials.
- [global-chem](https://github.com/Sulstice/global-chem) - A Chemical Knowledge Graph and Toolkit, writting in IUPAC/SMILES/SMARTS, for common small molecules from diverse communities to aid users in selecting compounds for forcefield parametirization.
- [matbench-discovery](https://github.com/janosh/matbench-discovery) - A benchmark for ML-guided high-throughput materials discovery.
- [NeuralForceField](https://github.com/learningmatter-mit/NeuralForceField) - Neural Network Force Field based on PyTorch.
- [openff-toolkit](https://github.com/openforcefield/openff-toolkit) - The Open Forcefield Toolkit provides implementations of the SMIRNOFF format, parameterization engine, and other tools.

## Molecular Visualization

*Packages for viewing molecular structures.*

- [ase-gui](https://wiki.fysik.dtu.dk/ase/ase/gui/gui.html#module-ase.gui) - The graphical user-interface allows users to visualize, manipulate, and render molecular systems and atoms objects.
- [chemiscope](https://github.com/lab-cosmo/chemiscope) - An interactive structure/property explorer for materials and molecules.
- [chemview](http://chemview.readthedocs.io/en/latest/) - An interactive molecular viewer designed for the IPython notebook.
- [imolecule](http://patrickfuller.github.io/imolecule/) - An embeddable webGL molecule viewer and file format converter.
- [moleculekit](https://github.com/Acellera/moleculekit) - A molecule manipulation library.
- [nglview](https://github.com/arose/nglview) - A [Jupyter](https://jupyter.org/) widget to interactively view molecular structures and trajectories.
- [PyMOL](https://pymol.org/) - A user-sponsored molecular visualization system on an open-source foundation, maintained and distributed by Schrödinger.
- [pymoldyn](https://pgi-jcns.fz-juelich.de/portal/pages/pymoldyn-main.html) - A viewer for atomic clusters, crystalline and amorphous materials in a unit cell corresponding to one of the seven 3D Bravais lattices.
- [rdeditor](https://github.com/EBjerrum/rdeditor) - Simple RDKit molecule editor GUI using PySide.
- [sumo](http://sumo.readthedocs.io/en/latest/) - A toolkit for plotting and analysis of ab initio solid-state calculation data.
- [surfinpy](https://surfinpy.readthedocs.io/en/latest/) - A library for the analysis, plotting and visualisation of ab initio surface calculation data.
- [trident-chemwidgets](https://github.com/tridentbio/trident-chemwidgets) - Jupyter Widgets to interact with molecular datasets.

## Database Wrappers

*Providing a python layer for accessing chemical databases*

- [ccdc](https://downloads.ccdc.cam.ac.uk/documentation/API/index.html) - An API for the Cambridge Structural Database System.
- [ChemSpiPy](http://chemspipy.readthedocs.io/en/latest/) - [ChemSpider](http://www.chemspider.com/) wrapper, that allows chemical searches, chemical file downloads, depiction and retrieval of chemical properties.
- [CIRpy](http://cirpy.readthedocs.io/en/latest/) - An interface for the Chemical Identifier Resolver (CIR) by the CADD Group at the NCI/NIH.
- [pubchempy](http://pubchempy.readthedocs.io/en/latest/) - PubChemPy provides a way to interact with PubChem in Python.
- [chembl-downloader](https://github.com/cthoyt/chembl-downloader) - Automate downloading and querying the latest (or a given) version of [ChEMBL](https://www.ebi.ac.uk/chembl/)
- [drugbank-downloader](https://github.com/cthoyt/drugbank_downloader) - Automate downloading, opening, and parsing [DrugBank](https://www.drugbank.com/)

## Learning Resources

*Resources for learning to apply python to chemistry.*

- [An Introduction to Applied Bioinformatics](https://github.com/applied-bioinformatics/iab2) - A Jupyter book demonstrating working with biochemical data using the scikit-bio library for tasks such as sequence alignment and calculating Hamming distances.
- [Computational Thermodynamics](https://kyleniemeyer.github.io/computational-thermo/content/intro.html) - This collection of Jupyter notebooks demonstrates solutions to a range of thermodynamic problems including solving chemical equilibria, comparing real versus ideal gas behavior, and calculating the temperature and composition of a combustion reaction.
- [SciCompforChemists](https://github.com/weisscharlesj/SciCompforChemists) - Scientific Computing for Chemists with Python is a Jupyter book teaching basic python in chemistry skills, including relevant libraries, and applies them to solving chemical problems.

## Miscellaneous Awesome

- [Colorful Nuclide Chart](https://people.physics.anu.edu.au/~ecs103/chart/) - A beatuful, interactive visualization of nuclides with access to a varirty of nuclear properties and allows saving high quality images for publications, presentations and outreach.

## See Also

- [awesome-cheminformatics](https://github.com/hsiaoyi0504/awesome-cheminformatics) Another list focuses on Cheminformatics, including tools not only in Python.
- [awesome-small-molecule-ml](https://github.com/benb111/awesome-small-molecule-ml) A collection of papers, datasets, and packages for small-molecule drug discovery. Most links to code are in Python.
- [awesome-molecular-docking](https://github.com/yangnianzu0515/awesome-molecular-docking) A curated list of molecular docking software, datasets, and papers.
- [jarvis](https://jarvis.nist.gov/) Joint Automated Repository for Various Integrated Simulations is a repository designed to automate materials discovery and optimization using classical force-field, density functional theory, machine learning calculations and experiments.
- [polypharmacy-ddi-synergy-survey](https://github.com/AstraZeneca/polypharmacy-ddi-synergy-survey) A collection of research papers (with Python implementations) focusing on drug-drug interactions, synergy and polypharmacy.
