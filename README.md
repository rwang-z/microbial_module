
# Deciphering associations between gut microbiota and clinical factors using microbial modules

## Introduction

The method is used to identify underlying microbial modules, i.e., groups of taxa with similar abundance patterns, from longitudinal human gut microbiota. The identified modules reflected both the compositional and dynamic characteristics of the gut microbiota.

The samples of the longitudinal data at species level was modeled as a third-order tensor and the compositional relative abundance was factorized into the contributions of microbial modules. On account of that the gut microbiota composition of a subject should change smoothly over time, we assumed correlated activities of the microbial modules at neighbor time points in the model.

For more details, please read the paper: Deciphering associations between gut microbiota and clinical factors using microbial modules

Workflow of the mehtod:

![figures](https://user-images.githubusercontent.com/57746198/152830171-2332eb42-3bf6-4fd6-a1fa-a5b2bcf351cf.png)


## Platform

R version 4.0.3

## Usage

Save the complete tensor of dimension (subject * taxa * visits) constructed from your samples as a variable 'data_tensor' in a .RData file.
The tensor will be factorized into three factor matrices in terms of subjects (A), taxa (X) and visits (B), respectively, which share a common dimension of C microbial modules. 

The rows of X indicate the members of the modules, whereas the entries in A and B represent the activities of the modules in the subjects and visits, respectively.

Execute the program by running:

```
> source(‘main.r’)
> res = main(num_components, tensor_file, seed_flag, seed)
```

- num_components: the maximum number of microbial modules. The algorithm will automatically determine an appropriate number of modules less than the maximum limit. Default 300.
 
- tensor_file: the path of the .RData file containing the tensor. The tensor file constructed from the “10-visit-set” in the paper is provided as default ('data/iHMP_10_visit_tensor.RData').

- seed_flag: how to initialize the variables.

  - ‘random’: do not set seed and initialize the variables randomly. Default value.
  
  - ‘set’: using the seed provided to the parameter ‘seed’.
	
- seed: the seed used for initialization when ‘seed_flag’ is set to ‘set’. Default 12345.

### Output

The algorithm returns a variable containing the factorization results of A, B and X. 

The activities of the taxa in the modules can be extracted by

```
> X = res$X$mom1 * round(res$X$gamma)
```

The activities of the modules in the subjects:

```
> A = res$A$mu
```

The activities of the modules in the subjects:

```
> B = res$B$nu
```

Save results to .RData file

```
> save(res, file = ‘output_path’)
```


