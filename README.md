# NetREX
A python tool to reconstruct a gene regulatory network given context-specific expression data and a prior network.

## Pre-request
1. Install Anaconda (python3 version) based on https://docs.anaconda.com/anaconda/install/
2. Install package progressbar2
```bash
   $ pip install progressbar2
```

## How to use
1. Default
```bash
$ python ./NetREX.py -e express_file -p prior_file 
```   
-e [expression file name] (Required)[Expression file format](#exp-file-format) is explained below. 

-p [prior file name] (Required)[Prior file format](#prior-file-format) is explained below. 

2. Advanced
```bash
$ python ./NetREX.py -e express_file -p prior_file -k 0.6 0.7 0.8 -t 1.2
``` 
-k [a list of percentage] (Optional)"0.6 0.7 0.8" means that NetREX would keep 60%, 70%, and 80% respectively in its parameters space.

-t [a ratio] (Optional)"1.2" means # edges in the output network is 1.2 times to the # edges in the prior network 


## Expression file format

## Prior file format

## Format of the prior network file for NetREX

# PriorBoost

