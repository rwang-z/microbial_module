# tensor factorization with dependencies on the loadings of factor matrix X
#
# The code is modified from SDA4D (https://github.com/marchinilab/SDA4D)

outer_product <- function(A,B,C, num_components){
    dim_A = dim(A)[1]
    dim_B = dim(B)[1]
    dim_C = dim(C)[1]
    product <- array(0,c(dim_A, dim_B, dim_C))
    for(compo in 1:num_components){
        product = product + A[,compo] %o% B[,compo] %o% C[,compo]
    }

    return(product)
}

# calculate reconstruction RMSE
cal_rmse = function(x, y){
    se = (x-y)^2
    rmse = sqrt(mean(se))

    return(rmse)
}


tensor_factorization<-function(params, data, maxiter, track=10, debugging=FALSE){

    initialization <- function(params){
        # initialise the variables
        list_of_vars <- list()
        list_of_vars$Error=0
        list_of_vars$Neg_FE=c()

        # A: N by C, first dimension, normal distribution
        list_of_vars$A <- list(mu = matrix(rnorm(params$N * params$C),params$N,params$C),  # mean
                               precision = matrix(100,params$N,params$C))   # variance
        list_of_vars$A$mom2 <- list_of_vars$A$mu^2 + 0.01

        # B: T by C, third dimension, normal distribution
        list_of_vars$B <- list(nu = matrix(rnorm(params$T * params$C),params$T,params$C),  # mean
                               precision = matrix(100,params$T,params$C))   # variance
        list_of_vars$B$mom2 <- list_of_vars$B$nu^2 + 0.01

        # Ita: C by 1, variance of B, gamma distribution
        list_of_vars$Ita <- list(a = matrix(params$a,params$C,1),
                                b = matrix(params$b,params$C,1),
                                mom1 = matrix(params$a * params$b,params$C,1))

        # Lamda: N by T, noise variance, gamma distribution
        list_of_vars$Lam <- list(u = matrix(params$u,params$N,params$T),
                                 v = matrix(params$v,params$N,params$T))
        list_of_vars$Lam$mom1 = list_of_vars$Lam$u*list_of_vars$Lam$v

        # Beta: C by 1, variance of W, gamma distribution
        list_of_vars$Beta <- list(e = matrix(params$e,params$C,1),
                                  f = matrix(params$f,params$C,1),
                                  mom1 = matrix(params$e * params$f,params$C,1))

        list_of_vars$Ph <- matrix(0.5,params$C,params$L)
        list_of_vars$Ps <- matrix(0.5,params$C,params$L)
        list_of_vars$Rho <- matrix(rbeta(params$C,params$r,params$z),1,params$C)

        # X: C by L, second dimension, sparse factor matrix, spike and slab prior, normal-bernoulli distribution
        list_of_vars$X <- list(gamma = matrix(0.5,params$C,params$L),  
                                sigma = matrix(c(100),params$C,params$L),
                                m = matrix(list_of_vars$Beta$e *list_of_vars$Beta$f,params$C,params$L),
                                mom1 = matrix(0,params$C,params$L),
                                mom2 = matrix(0,params$C,params$L))
        list_of_vars$X$mom1 = list_of_vars$X$m * list_of_vars$X$gamma   # E[X] = m*gamma
        list_of_vars$X$mom2 = (list_of_vars$X$m**2 + 1/list_of_vars$X$sigma) * list_of_vars$X$gamma   # E[X^2] = (m^2 + 1/sigma)*gamma
        
        return(list_of_vars)
    }

    ############################ local functions ###################################
    ############# Negative Free Energy ##############
    Free_Energy<-function(params){
        #returns negative free energy
        FE = 0

        # E[log(lamda)]
        FE = FE + 0.5 * params$L * sum(data$indicator * (digamma(vars$Lam$u)+log(vars$Lam$v)))

        # E[log(p(y,params))]
        component_1 = (data$data_tensor - outer_product(vars$A$mu, t(vars$X$mom1), vars$B$nu, params$C))^2
        component_2 = outer_product(vars$A$mom2, t(vars$X$mom2), vars$B$mom2,params$C)
        component_3 = outer_product(vars$A$mu^2, t(vars$X$mom1^2), vars$B$nu^2,params$C)
        summation = component_1 + component_2 - component_3     # <(y-sum_c(abx))^2>, N by L by T
        sum_n = matrix(0,params$N, params$T)
        for(t in 1:params$T){
            sum_n[,t] = rowSums(summation[,,t])     # summation over L
        }
        FE = FE - 0.5 * sum(sum_n * vars$Lam$mom1 * data$indicator)   # summation over N, T

        ######## the terms from the prior and approx posteriors
        # FE with respect to A
        FE = FE - 0.5*sum(vars$A$mom2 + log(abs(vars$A$precision))) #+ (params$N * params$C)/2
        
        # FE with respect to B
        tmp = sum((params$a - vars$Ita$a + 0.5 * params$T) * (digamma(vars$Ita$a) + log(abs(vars$Ita$b))) 
                    + vars$Ita$a * (1 - vars$Ita$b / params$b + log(abs(vars$Ita$b))) + lgamma(vars$Ita$a))
        tmp = tmp - 0.5 * sum(log(abs(vars$B$precision)))
        var_b_prime_nu = rbind(array(0, params$C), vars$B$nu[1:(params$T-1),])
        var_b_prime_mom2 = rbind(array(0, params$C), vars$B$mom2[1:(params$T-1),])
        var_b_b_prime = (vars$B$nu - var_b_prime_nu)^2 - vars$B$nu^2 - var_b_prime_nu^2 + vars$B$mom2 + var_b_prime_mom2
        FE = FE + tmp - 0.5 * sum(t(matrix(vars$Ita$mom1, params$C, params$T)) * var_b_b_prime)

        # FE with respect to w
        Wmom2 = vars$X$gamma * ((1/vars$X$sigma) + (vars$X$m)**2) + (1-vars$X$gamma) * matrix(1/vars$Beta$mom1,params$C,params$L)
        FE = FE + 0.5*sum(
            matrix(digamma(vars$Beta$e) + log(vars$Beta$f),params$C,params$L) - matrix(vars$Beta$mom1,params$C,params$L)*Wmom2
        )

        FE = FE + sum(
            -0.5*(
                vars$X$gamma * log(vars$X$sigma) + (1-vars$X$gamma)*log(matrix(vars$Beta$mom1,params$C,params$L))
            )
        )

        # FE with respect to Rho
        FE = FE + sum(
            (params$r - 1)*log(vars$Rho) + (params$z - 1)*log(1-vars$Rho)
        )

        # FE with respect to Psi
        FE = FE + sum(
            (params$g - 1) * log(vars$Ps) + (params$h - 1)*log(1 - vars$Ps)
        )

        # FE with respect to Phi
        FE = FE + sum(
            vars$Ph * matrix(log(vars$Rho),params$C,params$L) + (1-vars$Ph)*matrix(log(1-vars$Rho),params$C,params$L)
        )

        # FE with respect to s
        FE = FE + sum(
            vars$X$gamma*log(vars$Ph*vars$Ps) + (1-vars$X$gamma)*log(1-vars$Ph*vars$Ps)
        )

        Xtmp<- (-(1-vars$X$gamma)*log(1-vars$X$gamma) - vars$X$gamma*log(vars$X$gamma))
        if(any(vars$X$gamma==0 | vars$X$gamma==1)){
            Xtmp[vars$X$gamma==0 | vars$X$gamma==1]=0
        }
        FE = FE + sum(Xtmp)

        # FE with respect to Beta
        FE = FE + sum(
            params$e * log(abs(vars$Beta$f)) + (params$e - vars$Beta$e)*digamma(vars$Beta$e) +
                vars$Beta$e - vars$Beta$mom1/params$f + lgamma(vars$Beta$e)
        )

        # FE with respect to Lamda
        FE = FE + sum(
            params$u*log(vars$Lam$v) + (params$u - vars$Lam$u)*digamma(vars$Lam$u) +
                vars$Lam$u - vars$Lam$u*vars$Lam$v/params$v + lgamma(vars$Lam$u)
        )

        return(FE)
    }

    ################## update functions ####################
    updateA=function(params){
        A=vars$A   # N by C

        # precision
        tmp = (data$indicator * vars$Lam$mom1) %*% vars$B$mom2   # N by C, summation over T  
        A$precision = t(t(tmp) *  rowSums(vars$X$mom2)) + matrix(1,params$N,params$C)   # N by C, summation over T

        # the first term of vars$A$mu
        # mean term 
        mean_term1 = matrix(0, params$N, params$C)
        for(n in 1:params$N){
            tmp = vars$X$mom1 %*% data$data_tensor[n,,]  # C by T, summation over L
            mean_term1[n,] = (data$indicator[n,] * vars$Lam$mom1[n,]) %*% (t(tmp) * vars$B$nu)  # 1 by C, summation over T
        }

        # the second term of vars$A$mu
        for(c in 1:params$C){
            x_product = colSums(vars$X$mom1[c,] * t(vars$X$mom1[-c,]))    # 1 by C-1, summation over L
            b_product = matrix(vars$B$nu[,c],params$T,params$C-1)*vars$B$nu[,-c]  # T by C-1
            tmp = (data$indicator * vars$Lam$mom1) %*% b_product  # N by C-1, summation over T
            A$mu[,c] = (mean_term1[,c] - colSums(t(A$mu[,-c] * tmp) * x_product)) / A$precision[,c]  # summaion over C-1
        }

        A$mom2 = 1/A$precision + A$mu^2

        return(A)
    }

    updateB=function(params){
        B=vars$B   # T by C

        coe = rbind(matrix(2, (params$T -1), params$C), array(1, params$C))
        tmp = t(matrix(vars$Ita$mom1, params$C, params$T)) * coe
        B$precision = t(rowSums(vars$X$mom2) * (t(vars$A$mom2) %*% (data$indicator * vars$Lam$mom1))) + matrix(1, params$T, params$C) + tmp # T by C, summation over L, N

        # mean term 1
        mean_term1 = matrix(0, params$T, params$C)
        for(t in 1:params$T){
            tmp = vars$X$mom1 %*% t(data$data_tensor[,,t]) # C by N, summation over L
            mean_term1[t,] = (data$indicator[,t] * vars$Lam$mom1[,t]) %*% (vars$A$mu * t(tmp))    # 1 by C, summation over N
        }

        for(c in 1:params$C){
            x_product = colSums(vars$X$mom1[c,] * t(vars$X$mom1[-c,]))   # 1 by C-1, summation over L
            a_indicator = t(data$indicator * vars$Lam$mom1) %*% (matrix(vars$A$mu[,c],params$N,params$C-1)*vars$A$mu[,-c])  # T by C-1, summation over N
            tmp = (a_indicator * B$nu[,-c]) %*% x_product  # T by 1, summation over C-1
            var_b_pre_nu = c(0, vars$B$nu[1:(params$T-1),c])
            var_b_after_nu = c(vars$B$nu[2:params$T,c], 0)
            b_time_sum = (var_b_pre_nu + var_b_after_nu) * vars$Ita$mom1[c]
            B$nu[,c] = (mean_term1[,c] - tmp + b_time_sum) / B$precision[,c]
        }
        B$mom2 = 1/B$precision + B$nu^2
        
        return(B)
    }

    updateX=function(params){
        X=vars$X   # C by L

        # sigma term 2
        tmp = colSums(t(data$indicator * vars$Lam$mom1) %*% vars$A$mom2 * vars$B$mom2)  # 1 by C, summation over N, T
        X$sigma = matrix(tmp, params$C, params$L) + matrix(vars$Beta$mom1,params$C,params$L)  # C by L

        # mean term 1
        m_term1 = matrix(0, params$C, params$L)
        for(l in 1:params$L){
            tmp = data$data_tensor[,l,] * data$indicator * vars$Lam$mom1  # N by T
            m_term1[,l] = colSums(vars$A$mu * (tmp %*% vars$B$nu))  # N by C, summation over N, T
        }

        for(c in 1:params$C){
            a_indicator = t(data$indicator * vars$Lam$mom1) %*% (matrix(vars$A$mu[,c],params$N,params$C-1) * vars$A$mu[,-c])  # T by C-1, summation over N
            tmp = colSums(matrix(vars$B$nu[,c],params$T,params$C-1) * vars$B$nu[,-c] * a_indicator)    # 1 by C-1, summation over T
            X$m[c,] = (m_term1[c,] - (tmp %*% X$mom1[-c,])) / X$sigma[c,]  # summation over C-1

            u = -0.5*log(X$sigma[c,]) + 0.5 * (X$m[c,]^2) * X$sigma[c,] 
                + log(vars$Ps[c,]*vars$Ph[c,]) + 0.5*log(vars$Beta$mom1[c]) - log(1-vars$Ps[c,]*vars$Ph[c,])

            X$gamma[c,] = 1/(1+exp(-u))
            X$mom1[c,] = X$m[c,] * X$gamma[c,]
            X$mom2[c,] = (1 / X$sigma[c,] + X$m[c,]^2) * X$gamma[c,]
        } 

        return(X)
    }

    updateLam=function(params){
        Lam=vars$Lam
        Lam$u = params$u + 0.5 * params$L * data$indicator       # N by T

        # vars$Lam$v
        component_1 = (data$data_tensor - outer_product(vars$A$mu, t(vars$X$mom1), vars$B$nu, params$C))^2
        component_2 = outer_product(vars$A$mom2, t(vars$X$mom2), vars$B$mom2,params$C)
        component_3 = outer_product(vars$A$mu^2, t(vars$X$mom1^2), vars$B$nu^2,params$C)
        summation = component_1 + component_2 - component_3

        sum_n = matrix(0,params$N, params$T)
        for(t in 1:params$T){
            sum_n[,t] = rowSums(summation[,,t])     # summation over L
        }
        
        Lam$v = 1.0 / (1.0/params$v + 0.5 * sum_n * data$indicator)  # N by T
        Lam$mom1 = Lam$u * Lam$v

        return(Lam)
    }

    updateBeta=function(params){
        Beta=vars$Beta
        Wmom2 = vars$X$gamma * (1/vars$X$sigma  + vars$X$m^2) +
            (1-vars$X$gamma)*(matrix(1/Beta$mom1,params$C,params$L))

        Beta$e = (params$e + params$L/2)*matrix(1,params$C,1)

        for (c in 1:params$C){
            Beta$f[c] = 1/(1/params$f + 0.5*sum(Wmom2[c,]))
        }
        Beta$mom1 = Beta$e*Beta$f

        return(Beta)
    }

    updateIta = function(params){
        Ita = vars$Ita
        Ita$a = (params$a + 0.5 * params$T) * matrix(1,params$C,1)
        var_b_prime_nu = rbind(array(0, params$C), vars$B$nu[1:(params$T-1),])
        var_b_prime_mom2 = rbind(array(0, params$C), vars$B$mom2[1:(params$T-1),])
        var_b_b_prime = (vars$B$nu - var_b_prime_nu)^2 - vars$B$nu^2 - var_b_prime_nu^2 + vars$B$mom2 + var_b_prime_mom2
        tmp = colSums(var_b_b_prime)  # 1 by C
        Ita$b = 1 / (1/params$b + 0.5 * matrix(tmp, params$C, 1))
        Ita$mom1 = Ita$a * Ita$b

        return(Ita)
    }


    updateRho=function(params){
        Rho=vars$Rho
        for(c in 1:params$C){
            Rho[c] = (params$r - 1 + sum(vars$Ph[c,]))/(params$L + params$r + params$z -2)
        }
        return(Rho)
    }

    updatePhiPsi=function(params){

        Grad=function(x,y,c,l){
            v1 = vars$X$gamma[c,l]/x - (1-vars$X$gamma[c,l])/(1/y -x) + log(vars$Rho[c]) - log(1-vars$Rho[c])
            v2 = vars$X$gamma[c,l]/y - (1-vars$X$gamma[c,l])/(1/x - y) + (params$g - 1)/y - (params$h - 1)/(1-y)
            vec=c(v1,v2)
            return(vec)
        }
        Hess = function(X,c,l){
            mat=matrix(0,2,2)
            mat[1,1] = -vars$X$gamma[c,l]/(X[1]^2) -  (1-vars$X$gamma[c,l])/((1/X[2] - X[1])^2)
            mat[1,2] = -(1-vars$X$gamma[c,l])/(1-X[1]*X[2])^2
            mat[2,1] = mat[1,2]
            mat[2,2] = -vars$X$gamma[c,l]/(X[2]^2) - (1 - vars$X$gamma[c,l])/((1/X[1] - X[2])^2) - (params$g - 1)/(X[2]^2) - (params$h - 1)/((1-X[2])^2)
            return(mat)
        }
        tildeF = function(X,c,l){
            FE =0
            FE = FE + (vars$X$gamma[c,l])*log(X[1]*X[2]) + (1-vars$X$gamma[c,l])*log(1-X[1]*X[2]) +
                (params$g - 1)*log(X[2]) + (params$h -1)*log(1-X[2]) + X[1]*log(vars$Rho[c]) +(1-X[1])*log(1-vars$Rho[c])
            return(FE)
        }
        Phi=vars$Ph
        Psi=vars$Ps
        #use Newton's method for finding Ph, Ps (c,l)
        Xtol = 1e-6
        ftol = 1e-17
        for(c in 1:params$C){
            for (l in 1:params$L){
                tmpY = c(Phi[c,l],Psi[c,l])
                X=tmpY
                g = Grad(tmpY[1],tmpY[2],c,l)
                H=Hess(tmpY,c,l)
                dH = H[1,1]*H[2,2] - H[1,2]*H[2,1]
                Direction = (matrix(c(H[2,2],-H[2,1],-H[1,2],H[1,1]),2,2))%*%g * (1/dH)
                alpha = 0.1
                i=0
                current=tildeF(X,c,l)
                while(alpha^i>1e-10){
                    tmpY=c(X - (alpha^i)*Direction)
                    if(all(tmpY>1e-10) && all(tmpY<1-1e-10)){
                        if(tildeF(tmpY,c,l)>current){
                            Phi[c,l]=tmpY[1]
                            Psi[c,l]=tmpY[2]
                            break
                        }
                    }
                    i=i+1
                }
            }

        }
        return(list(Phi=Phi,Psi=Psi))
    }


    ######################################## main function ############################################

    #### initialise variables ####
    if(params$seed_flag == 'set'){
        set.seed(params$seed)
    }
    
    vars<-initialization(params)
    continue=TRUE
    iteration=1
    
    ##  Iterate until convergence:
    FEcur=Free_Energy(params)
    FEold=FEcur-1
    vars$Neg_FE = c(vars$Neg_FE, FEcur)
    trackingvec=rep(10*track,track) # a vector of length track, listing the changes of modules in the last track iterations

    while(iteration<=maxiter & continue){
        print(paste0("Iteration ",iteration))

        # update Beta
        vars$Beta=updateBeta(params)
        print('Updating Beta')
        if(debugging){
            FEold=FEcur
            FEcur=Free_Energy(params)
            vars$Neg_FE = c(vars$Neg_FE, FEcur)
            print(paste0("Current FE: ", FEcur))
            if(FEcur<FEold){
                vars$Error=1
                print("Error in updating Beta: negative free energy decreased.")
                break
            }
        }

        # update Ita
        vars$Ita=updateIta(params)
        print('Updating Ita')
        if(debugging){
            FEold=FEcur
            FEcur=Free_Energy(params)
            vars$Neg_FE = c(vars$Neg_FE, FEcur)
            print(paste0("Current FE: ", FEcur))
            if(FEcur<FEold){
                vars$Error=1
                print("Error in updating Ita: negative free energy decreased.")
                break
            }
        }

        # update X
        vars$X=updateX(params)
        print('Updating X')
        if(debugging){
            FEold=FEcur
            FEcur=Free_Energy(params)
            print(paste0("Current FE: ", FEcur))
            vars$Neg_FE = c(vars$Neg_FE, FEcur)
            if(FEcur<FEold){
                vars$Error=1
                print("Error in updating X: negative free energy decreased.")
                break
            }
        }

        # update Rho
        vars$Rho=updateRho(params)
        print('Updating Rho')
        if(debugging){
            FEold=FEcur
            FEcur=Free_Energy(params)
            print(paste0("Current FE: ", FEcur))
            vars$Neg_FE = c(vars$Neg_FE, FEcur)
            if(FEcur<FEold){
                vars$Error=1
                print("Error in updating Rho: negative free energy decreased.")
                break
            }
        }

        # update Phi and Psi 
        tmp = updatePhiPsi(params)
        vars$Ph = tmp$Phi
        vars$Ps = tmp$Psi
        print('Updating Phi and Psi')
        if(debugging){
            FEold=FEcur
            FEcur=Free_Energy(params)
            print(paste0("Current FE: ", FEcur))
            vars$Neg_FE = c(vars$Neg_FE, FEcur)
            if(FEcur<FEold){
                vars$Error=1
                print("Error in updating Phi and Psi: negative free energy decreased.")
                break
            }
        }

        # update A
        vars$A=updateA(params)
        print('Updating A')
        if(debugging){
            FEold=FEcur
            FEcur=Free_Energy(params)
            print(paste0("Current FE: ", FEcur))
            vars$Neg_FE = c(vars$Neg_FE, FEcur)
            if(FEcur<FEold){
                vars$Error=1
                print("Error in updating A: negative free energy decreased.")
                break
            }
        }

        # update B
        vars$B=updateB(params)
        print('Updating B')
        if(debugging){
            FEold=FEcur
            FEcur=Free_Energy(params)
            print(paste0("Current FE: ", FEcur))
            vars$Neg_FE = c(vars$Neg_FE, FEcur)
            if(FEcur<FEold){
                vars$Error=1
                print("Error in updating B: negative free energy decreased.")
                break
            }
        }
        
        # update Lamda
        vars$Lam=updateLam(params)
        print('Updating Lamda')
        if(debugging){
            FEold=FEcur
            FEcur=Free_Energy(params)
            print(paste0("Current FE: ", FEcur))
            vars$Neg_FE = c(vars$Neg_FE, FEcur)
            if(FEcur<FEold){
                vars$Error=1
                print("Error in updating Lamda: negative free energy decreased.")
                break
            }
        }

        # evaluate whether to stop (check the changes in the modules)
        PIP = round(vars$X$gamma)
        if(iteration>1){
            indexingvar=iteration%%track
            if(indexingvar==0){indexingvar=track}
            trackingvec[indexingvar]=sum(abs(PIP - PIP_old))
            if(mean(trackingvec)<1){continue=FALSE} 
        }
        PIP_old=PIP
        iteration=iteration+1

        # RMSE of reconstructed tensor
        reconstructed <-  outer_product(vars$A$mu, t(vars$X$mom1), vars$B$nu, params$C)
        tensor_rmse = cal_rmse(data$data_tensor, reconstructed)
        print(paste0('RMSE: ', tensor_rmse))
    }

    vars$maximumiteration = iteration - 1
    vars$data_tensor = data$data_tensor
    vars$last_FE = FEcur
    vars$indicator = data$indicator
    vars$rmse = tensor_rmse
    vars$params = params

    return(vars)
}
