# Identifying microbial modules from logitudinal human gut microbiota 

source('noweight_tensor_factorization_time_dependency.r')

run_factorization <- function(params, data_list, max_iters=2000, track=10){
    # params: parameters of factorization
    # data, a list containing the data_tensor, the indicator (matrix of ones for complete tensor), and the size of the tensor
    # max_iters: maximum iteration for tensor factorization
    # track: the number of iterations to consider in checking convergence

    if(any(is.na(data_list$data_tensor))){    # check NA in tensor data
      stop('ERROR: the data tensor contains missing data')
    }
    if(any(is.infinite(data_list$data_tensor))){   # check infinite values 
      stop('ERROR: the data tensor contains infinite data')
    }
    if(!is.numeric(data_list$data_tensor)){   # check numeric values
      stop('ERROR: the data tensor contains non-numeric data')
    }

    ########### run the method and return variables ############
    res<-tensor_factorization(params, data_list, maxiter=max_iters,track=track,debugging=FALSE)
    print("Estimation finished!")
    
    # output parameters
    print('Parameters used')
    print(paste0('e, f, g, h, u, v, r, z: ', params$e, ',', params$f, ',', params$g, ',', params$h, ',', 
                    params$u, ',', params$v, ',', params$r, ',', params$z))
    print(paste0('Data dimension: ', params$N, ',', params$L, ',', params$T))
    
    print(paste0('Number of componets: ', params$C))
    if(params$seed_flag == 'set'){
        print(paste0('seed: ', params$seed))
    }
    
    if(res$Error==1){
        print('Finished but Negative Free Energy shows signs of a decrease, this should be checked.')
    }else{
        print('Finished')
    }
    print(paste(res$maximumiteration,' Iterations were carried out.'))
    print('returning output')
    return(res)
}

main = function(num_components = 300, tensor_file = 'data/iHMP_10_visit_tensor.RData', seed_flag = 'random', seed = 12345){
    ##### num_components: number of components
    ##### tensor_file: a RData file containing a variable 'data_tensor' of dimension (subject * taxa * visits), a complete tensor
    ##### seed_flag: the way to set seed
    #####   -- 'set': to use the provided seed 
    #####   -- 'random': to initialize variables randomly
    ##### seed: seed used for initialization when seed_flag is 'set'

    # loading data tensor
    load(tensor_file)
    num_patient = dim(data_tensor)[1]
    num_taxa = dim(data_tensor)[2]
    num_visit = dim(data_tensor)[3]
    print(paste0('Data dimension: ', num_patient, ',', num_taxa, ',', num_visit))
    print(paste0('Number of componets: ', num_components))

    data <- list(num_patient = num_patient, num_taxa = num_taxa, week_num = num_visit, 
                 data_tensor = data_tensor,
                 indicator = matrix(1, num_patient, num_visit))

    params=list(N=num_patient, L=num_taxa, T=num_visit, C=num_components, seed = seed, seed_flag = seed_flag,
                a=1e-6, b=1e6, e=1e-6, f=1e6, g=0, h=0, u=1e-6, v=1e6, r=1, z=1, c0=1e-6, d0=1e6)

    # run tensor factorization
    res <- run_factorization(params, data)
    # save(res, file = 'result/factorized_results.RData')

    return(res)
}

# source('main.r')
# res = main()

