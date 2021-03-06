# Thoron Detection Project 
This project was done by Stephen Styles, Edward Timko, and Liam Wrubleski for the Pacific Institute of Mathematical Sciences (PIMS) Math^Industry workshop. The project was given by Environmental Instruments Canada (EIC), and the industry contact and mentor during this project was Kai Kaletsch.
## Description
The goal of this project was to find a method of determining the amount of Radon-220 and Radon-222 in a sample of air that may contain any amount of either. This works by using the EIC Radon Sniffer to detect alpha particles from the decay of these isotopes and their progeny. As the source of these alpha particles cannot by directly determined by the sniffer, an algorithm must be developed that allows for the amount of each isotope to be recovered from the decay signal. This repository contains the code used during the course of this project.

## Repository
The contents of this repository are divided into folders, by language and purpose.
### Data
This folder contains all data files used or generated by any of the code. The file `August 20, 2020 19_31.xlsx` is sample data provided by EIC. The files labelled with `Trial` are specific trial runs, isolated from the sample data. The other files record the result of the most recent grid search.
### Images
This folder contains some of the images generated for use in the final paper and presentation. Other images are stored with the specific document.
### LaTeX
This folder contains the LaTeX files, images, and compiled documents for the final paper, final presentation, and check-in presentation from the workshop.
### Matlab
This folder contains the Matlab code used to generate some of the images used in various documents.
### Python
This folder contains the Python code used for the primary analysis and algorithm development. To get the parameters for the resulting model, run `model_parameters.py`. To run a gridsearch, run `grid_search.py`. To compare the simulated data against real-world data, and determine the activity and amount from that real world data, run `realworld_compare.py`. These files were created after development was completed, for the clarity of this repository. The code used during development is in `old code`. Each of the above files should be edited to change parameters or files.
