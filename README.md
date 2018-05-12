# [NetREX](#netrex) and [PriorBoost](#priorboost)

# NetREX
A python tool to reconstruct a gene regulatory network given context-specific expression data and a prior network.

## Pre-request
1. Install Anaconda (python3 version) based on https://docs.anaconda.com/anaconda/install/
2. Install package progressbar2
```bash
   $ pip install progressbar2
```

## How to use?
### Default
```bash
$ python ./NetREX.py -e express_file -p prior_file 

-e [expression file name] <Required> expression file format is explained below. 
-p [prior file name] <Required> prior network file format is explained below. 
```   


### Advanced
```bash
$ python ./NetREX.py -e express_file -p prior_file -k 0.6 0.7 0.8 -t 1.2

-k [a list of percentages] (Optional)"0.6 0.7 0.8" means that NetREX would keep 60%, 70%, and 80% edges in the prior respectively. The final predicted network is the consensus based on networks predicted from those percentages.
-t [a ratio] (Optional)"1.2" means the number of  edges in the output network is 1.2 times to the edges in the prior network 
``` 

## Output of NetREX
Two files will be outputed: 1) The predicted network with edge pair format; 2) The predicted network with adjacency matrix format.

## Example of use NetREX
Put NetREX.py and prior.txt and expression.txt in the NexREX_Example_Data folder in the same folder. Then run:
```bash
$ python ./NetREX.py -e expression.txt -p prior.txt
```
The rank of each edges will be outputed in "NetREX_PredictedEdgeList.txt" and "NetREX_PredictedNetwork.tsv" as shown in this repo.


## Expression file format
A Tab sepeated file with the first colomn stroing the gene names. For example:

| E1 | 0.2508 | 0.2684 | 0.2786 | 0.2878  | ... |
|----|--------|--------|--------|---------|-----|
| E2 | 0.3149 | 0.3323 | 0.3427 | 0.3538  | ... |
| E3 | 0.2361 | 0.2526 | 0.2616 | 0.2698  | ... |
| E4 | 0.2530 | 0.2699 | 0.2791 | 0.2910  | ... |
| E5 | 0.2521 | 0.2693 | 0.2797 | 0.2885  | ... |
| E6 | 0.3162 | 0.3291 | 0.3390 | 0.3482  | ... |
| E7 | 0.2866 | 0.3021 | 0.3105 | 00.3175 | ... |

## Prior file format
A Tab sepeated file with the first colomn stroing the gene names and the first row stroing the TF names. For example:

|    | M1  | M2  | M3  |
|----|-----|-----|-----|
| E1 | 0.0 | 1.0 | 1.0 |
| E2 | 1.0 | 1.0 | 1.0 |
| E3 | 1.0 | 1.0 | 0.0 |
| E4 | 1.0 | 1.0 | 0.0 |
| E5 | 1.0 | 1.0 | 1.0 |
| E6 | 1.0 | 1.0 | 1.0 |
| E7 | 1.0 | 0.0 | 1.0 |



# PriorBoost
A python tool to compare the explaination power of two different networks. Assume Net1 is the network obtained from a prior-based method and Net2 is a network obtained from a expression-based method. 
## How to use?
```bash
$ python ./PriorBoost.py -e express_file -p Net1_file -b Net2_file 

-e [expression file name] <Required> expression file format is explained above. 
-p [Net1 file name] <Required> Net1 network file format is explained above. Elements in the Network are their ranks. rank 1 is the best.
-b [Net2 file name] <Required> Net2 network file format is explained above. Elements is the weights of each edges. The larger the better.
```   
