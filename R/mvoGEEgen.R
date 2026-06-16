#------------------------------------------------------------
# function for iterative unweighted GEE
#------------------------------------------------------------

#' @export

mvoGEEgen <- function(mdat, beta0, nbeta, id, Y_var, nresponse, X_mat, corr.str, Z){
  
  match.call()
  beta <- beta0
  bdiff <- 1
  iter <- 0
  
  while(sum(abs(bdiff))>.00000001){
    DWZ <- matrix(0, nrow = nbeta)
    DWD <- matrix(0, ncol = nbeta, nrow = nbeta)
    DWZZWD <- matrix(0, ncol = nbeta, nrow = nbeta)
    
    if (corr.str == "ind")
      R <- diag(nresponse)
    if (corr.str == "unstr")
      R <- estalpha1_unstr_fixed(mdat, id, Z, nresponse)
    if (corr.str == "exch")
      R <- estalpha1_exch_fixed(mdat, id, Z, nresponse)
    
    for (i in id){
      
      csizei <- nlevels(as.factor(mdat[mdat$cluster == i,]$unit))
      wi <- 1/csizei
      unitveci <- unique(mdat[mdat$cluster == i,]$unit)
      
      DWZj <- matrix(0, nrow = nbeta)
      DWDj <- matrix(0, ncol = nbeta, nrow = nbeta)
      
      for(j in unitveci){
        
        y <- as.matrix( Y_var[mdat$cluster==i & mdat$unit == j] )
        x <- as.matrix( X_mat[mdat$cluster==i & mdat$unit == j,] )
        u <- exp( x %*% beta ) # mean function 
        dudb <- exp( x %*% beta ) # derivative
        D <- x[,1] * dudb
        for(p in 2:nbeta) D <- cbind(D, x[,p] * dudb)
        vv <- u 
        V <- matrix(ncol = length(vv), nrow = length(vv), 0)
        diag(V) <- vv
        W <- V ^ (1/2) %*% R %*% V ^ (1/2)
        invW <- solve(W)
        DWZj <- DWZj + t(D) %*% invW %*% ( y - u ) 
        DWDj <- DWDj + t(D) %*% invW %*% D
        
      }
      
      DWZ <- DWZ +DWZj
      DWD <- DWD + DWDj
      DWZZWD <- DWZZWD + ( DWZj ) %*% t( DWZj )
      
    }
    
    invDWD <- solve(DWD)
    bdiff <- invDWD %*% DWZ
    beta <- beta + bdiff
    covbeta <- invDWD %*% DWZZWD %*% invDWD
    U <- exp(X_mat %*% beta)
    var <- U
    Z <- (Y_var - U) / sqrt(var)
    iter <- iter + 1
  }
  
  fit <- list()
  fit$coefficients <- beta
  fit$robust.variance <- covbeta
  fit$robust.se <- sqrt(diag(covbeta))
  fit$wald.chisq <- (beta/sqrt(diag(covbeta)))^2
  fit$p.value <- 1 - pchisq(fit$wald.chisq, df = 1)
  fit$R <- R
  fit$niter <- iter
  fit
  
}
