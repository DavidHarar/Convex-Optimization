# David Harar


rm(list = ls())
gc()

y <- function(s,k,A){
  return(exp((-k)/(10^5)) * pmax(0, 10 - s * (1 + 0.007 * A)))
}

s = 2
N = 100

# i
vmax1 <- function(s = 2, N = 100){
  
  V_list <- list()
  
  A <- seq(-N,N,2) # possible sum xi values
  Y <- y(s, k = N, A)
  
  assign(paste("Y_", N, sep = ""), Y)
  assign(paste("V_", N, sep = ""), Y) # In time N we sell the asset for sure
  assign(paste("A_", N, sep = ""), A)
  
  V_list <- append(V_list, list(get(paste("V_", N, sep = ""))))
  
  for (k in (N-1):0){
    A <- seq(-k,k,2)
    Y <- y(s, k = k, A)
    V <- c()
    
    for (l in (1:length(Y))){  # we always have only one number in n+1 comapred to n
      
      V_prev <- get(paste("V_", k+1, sep = ""))
      
      future_value <- 0.5 * V_prev[l] + 0.5 * V_prev[l+1]
      today_value  <- Y[l]
      decision <- max(today_value, future_value)
      V <- c(V, decision)
    }
    
    assign(paste("Y_", k, sep = ""), Y)
    assign(paste("V_", k, sep = ""), V)
    assign(paste("A_", k, sep = ""), A)
    
    V_list <- append(V_list, list(get(paste("V_", k, sep = ""))))
  }
  
  # V_list[[N+1]] is actually V for t = 0, that is bevause of Rs index method
  return(V_list)
}

V_func <- vmax1(9) # it took t > 0 (that is we went through t = 1)


V_func[[101]]
V_func[[100]]

# ii


