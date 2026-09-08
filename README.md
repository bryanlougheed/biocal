# biocal: Using sedimentological priors for more accurate calibration of 14C determinations from bioturbated sediment archives.

## Short description
A function for creating an accurate credible calibrated age distribution for radiocarbon-dated multi-specimen samples sourced from bioturbated sediment samples such as deep-sea sediment cores. Priors relating to sediment accumulation rate (SAR), bioturbation depth (BD), fraction broken microfossils and temporal changes in species abundance can be included. Traditional priors such as reservoir effect can also be included.

See published paper (open-access): 

B.C. Lougheed, 2022. "Using sedimentological priors to improve 14C calibration of bioturbated sediment archives." _Radiocarbon_, vol 64(1), pp 135-151. https://doi.org/10.1017/RDC.2021.116 

## Matlab and Python versions
The Python version is the version that is currently receiving active updates and has been heavily optimised using machine code. The Matlab version is no longer updated but still included as `biocal_matlab.zip` for legacy purposes. Note that because the Matlab version has not been similarly optimised, it will be significantly slower and use many more GB of RAM. I no longer have a Matlab license so you are on your own if you run into problems.

## Install as a Python package
Install (or upgrade) the package in your python environment using the `pip` terminal command (you may need to [install](https://github.com/git-guides/install-git) `git` first):

`pip install git+https://github.com/bryanlougheed/biocal.git`

If the install has been successful then you should be able to import the package in python in the usual way, e.g.:

`import biocal`

## Licence
Please consult the licence file in the repo.

## Tutorial

There is a functional example in the jupyter notebook `biocal_tutorial.ipynb`
