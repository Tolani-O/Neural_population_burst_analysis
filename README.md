# Population firing rate and peak time estimation method

## Description

This is a method to obtain statistical estimates of trial peak times, as well as trial to trial correlations between pairs of brain areas. It is based on the work of [Olarinre et al. (2024)]. It was applied to data from the Allen Brain Observatory Neuropixels Visual Coding dataset, which can be accessed via the [AllenSDK](https://allensdk.readthedocs.io/en/latest/visual_coding_neuropixels.html), the [DANDI Archive](https://dandiarchive.org/dandiset/000021), and through [AWS Registry of Open Data](https://registry.opendata.aws/allen-brain-observatory/).

## Getting Started

### Dependencies

* Python 3.7.9
* R 4.0.5

For Python scripts, the following packages are required:
```
numpy==1.19.5
pandas==1.1.5
allensdk==2.11.2
```

For R scripts, the following packages are required:
```
rstan_2.21.5         
gridExtra_2.3        
gratia_0.7.3         
mgcv_1.8-39          
dplyr_1.0.9           
ggplot2_3.3.6        
tidyverse_1.3.1  
```     

### Executing program

* To run the Python scripts, you will need to have python and the AllenSDK installed. You can find instructions for installing the AllenSDK [here](https://allensdk.readthedocs.io/en/latest/install.html).
* To run the R scripts, you will need to have R installed. You can find instructions for installing R [here](https://www.r-project.org/).
* To run the R scripts, you will need to have RStan installed. You can find instructions for installing RStan [here](www.mc-stan.org/users/interfaces/rstan).
* First run the Python script `download_and_format_allen_data.py` to download the data from the Allen Brain Observatory Neuropixels Visual Coding dataset. This will create a folder called `RStudioProjects/rNeuroPixel/rDataset` in the home directory. For each mouse id specified in the script, it will create a folder called `units_data_[mouse_id]`, where for each stimulus configuration and visual region, it will save the mouses spike train data for all trials of each as a datatable, which is a format readable to the R script. All mouse ids are specified in the script.
* Next, run the R script `generate_single_mouse_data.R`. This can be run from commandline with a single command line argument being a mouse id. The script will read the mouse data from the corresponding `units_data_[mouse_id]` folder, perform subset preselection for each visual area, perform naive estimation of peak times, with standard errors, compute posterior estimates of peak times and their trial to trial correlations, and finally it will write posterior summaries of the peaktimes, correlations and partial correlaitons to their corresponding folders, together with the naive correlation estimates.     

## Author

Motolani Olarinre

## Version History

* 0.1
    * Initial Release
