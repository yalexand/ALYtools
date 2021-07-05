# ALYtools

Set of various Bioimaging analysis applications written in Matlab

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation and Running](#installation-and-running)
- [Demo](#demo)
- [CIDR Problem](#cidr-problem)
- [ic_OPTtools](#ic_opttools)
- [OPT_ZFish_Embryo Problem](#opt_zfish_embryo_problem)
- [MicroscopyImageFormatter](#microscopyimageformatter)
- [License](#license)

# Overview

ALYtools massively uses Open Source code designed by other researchers. 
References are not kept systematically, but License section should cover all.  
If you find any of ALYtools solutions/ideas worth using, - please take without reference either.
Software uses Bioformats for image I/O, Icy for visualization and Excel format for quantification outputs; it can be run headless on HPC, and in batch processing mode.
"ALYtools.m" is the "Dispatcher"-style GUI capable of switching between different (unrelated) image analysis "Problems". Majority of "Problems" are short attempts or exercises aimed at analyzing Biomicroscopy data.
Several applications have proved to be useful; some of them were reported in the literature, and some are mentioned below.

# System Requirements

## Software
ALYtools applications require at least Matlab 2018a and Toolboxes, and Icy (http://icy.bioimageanalysis.org/).
They can be run under Windows 10, Mac OSX, and RedHat and CentOS linux.
## Hardware
Hardware requirements depend on "Problem", but in most cases average 2020-ish desktop specs suffice. NVIDIA GPU may be optionally used.

# Installation and Running

Nothing specific is needed beyond "git clone" or equivalent option from GUI. To use Icy visualization, one needs installing Icy with 2 plugins: "matlabcommuncator" and "matlabxserver". Icy should be installed somewhere on C:\ in common folder (not "Program Files").
When first running under Win10, ALYtools tries to find Icy directory; that can take couple of minutes. In general, it is better to run Icy first for all software mentioned below.

# Demo

In the "TestData" directory, one can find collection of Demo images, - a folder for "Problem". Not all "Problems" are covered (some demand too many images).

# CIDR Problem

CIDR ("collagen induced DDR1 redistribution") is the texture characterization method reported in the article
https://www.nature.com/articles/s41598-019-53176-4
Briefly, the algorithm counts average number of zero-crossings along random line profile on a segmented "patch" on the LF-subtracted texture image. "TestData" folder contains 2 images with different textures (rough and fine) acquired with SIM microscopy. Demo is run with default ALYtools options.

# ic_OPTtools

ic_OPTtools is the separate (not part of ALYtools) GUI-ed software for reconstructing OPT datasets. It uses basic FBP ("iradon") reconstruction, but also can apply TwIST algorithm. One can apply downsampling, choose FOV in axial direction, use GPU and multicore acceleration, use batch mode, visualize in Icy, and save resulting 3D images. Software also can process OPT-FLIM images.
Also, it can deduce and optionally apply simple (translation) axial mismatch correction.


# OPT_ZFish_Embryo Problem

This application was used in the article "Smad4 controls signaling robustness and embryo morphogenesis by playing distinct roles in the Nodal and BMP pathways" by Guglielmi et al (in press). In order to run Demo version, - 
```
- run Matlab file ALYtools.m
- choose "Settings->Problem->OPT_ZFish_Embryo"
- use "File->Load Single" to load the multi-channel OPT image "..your_path\ALYtools\TestData\ZFish_Embryo\WT_2.ome.tif" 
- click on Settings->Microns per pixel and enter 6.0
- click "Analysis->Run Current"
```
The input image contains 2 embryos. Results of analysis (quantification xls file and 2 images, one per embryo) will be saved in newly created directory with name consisting of "ALYtools Analysis" prefix, and time stamp postfix.  


# MicroscopyImageFormatter

MicroscopyImageFormatter is the application providing the model-based intensity artefacts correction in multi-channel fluorescence images acquired with time-lapse in Widefield mode. 
Detailed user guide will be added soon.

# License
This project is covered under GNU General Public License as published by the Free Software Foundation.







  









