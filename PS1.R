# Computation Methods in Nonlinear Optimization
# PS1
# David Harar


options(scipen=999)
rm(list = ls())

# 1. Find a positive root for the funcion x-sin(x)=0 ----
## (i) Use the bisection method ----
Q1_bisec <- function(a,b,epsilon) {
  "
  This function returns the root for a function with one variable.
  Inputs: 
  - a: The starting point of the interval in which the root is located.
  - b: The ending point of the interval in which the root is located.  
  - epsilon: a small number which defines the extant of accuracy.
  "
  
  q1 <- function(x){
    fx <- x-(2*sin(x))
    return(fx)
  }
  
  i <- 1
  # Starting condition
  fa = q1(a)
  fb = q1(b)
  d = (a+b)/2
  fd = q1(d)
  if (fd == 0) {
    
  }
  while(abs(fd)>epsilon) {
    if (fa*fd < 0) {
      b = d
      d = (a+b)/2
      fd = q1(d)
      
      interval = abs(a-b) # corrent interval
      i <- i + 1
    }
    
    if (fd*fb < 0) {
      a = d
      d = (a+b)/2
      fd = q1(d)
      
      interval = abs(a-b) # corrent interval
      i <- i + 1
    }
  }
  print(paste("The bisection method has converged in ", i, " iterations.",
              "The root is ", d,". and f(d) equals ", fd, sep = ""))
}  
Q1_bisec(a = (-2), b = 3, epsilon = 1/1000000000000)

## (ii) Use the Newton method ----
Q1_Newton <- function(x1, epsilon) {
  "
  This function returns the root for a convex function with one variable.
  Inputs: 
    - x1: a starting point
    - epsilon: a small number which defines the extant of accuracy.
  "
  # Test functions
  # q1 <- function(x){
  #   fx <- (x^2)
  #   return(fx)
  # }
  # q1_tag <- function(x){
  #   fx <- (2*x)
  #   return(fx)
  # }
  
  # Assignment functions
  q1 <- function(x){
    fx <- x-(2*sin(x))
    return(fx)
  }
  q1_tag <- function(x, epsilon){
    fx_tag <- 1-(2*cos(x))
    return(fx_tag)
  }

  
  # Starting conditions
  #x1 <- runif(1,-2,2)
  fx <- q1(x1)
  x_new <- x1
  f_new <- q1(x1)
  fx_tag <- q1_tag(x1)
  
  i <- 1
  
  while (abs(f_new) > epsilon){
    x_old = x_new
    fx_old = q1(x_old)
    fx_tag_old = q1_tag(x_old)
    
    x_new <- x_old - (fx_old/fx_tag_old)
    i <- i + 1
    
    
    # stopping rule
    f_new <- q1(x_new)
  }
  print(paste("The Newton method has converged in ", i, " iterations.",
              "The root is ", x_new, ". f(x) equals to ", f_new, sep = ""))
  return(x_new)
}
Q1_Newton(x1 = runif(1,min = -1,max = 3),epsilon = 1/1000000000000)

## (iii) Which method is faster ----
"
While the bisection method has converged to the value in 42 iterations, for the Newton method
it took only 7 iteration. Therefore the Newton methos is faster.
"

# 2. Find the value of 1/sqrt(2) without using "/" operation. ----
Q2_Newton <- function(x1, epsilon) {
  "
  This function returns the root for a convex function with one variable.
  Inputs: 
  - x1: a starting point
  - epsilon: a small number which defines the extant of accuracy.
  "
  # Assignment functions
  if (x1>=0) {
    print("Input out of range")
  }
  else{
    q1 <- function(x){
      fx <- -1/x - sqrt(2)
      return(fx)
    }
    q1_tag <- function(x, epsilon){
      fx_tag <- 1/(x^2)
      return(fx_tag)
    }
    
    
    # Starting conditions
    #x1 <- runif(1,-2,2)
    fx <- q1(x1)
    x_new <- x1
    f_new <- q1(x1)
    fx_tag <- q1_tag(x1)
    
    i <- 1
    
    while (abs(f_new) > epsilon){
      x_old = x_new
      fx_old = q1(x_old)
      fx_tag_old = q1_tag(x_old)
      
      x_new <- x_old - (fx_old/fx_tag_old)
      i <- i + 1
      
      
      # stopping rule
      f_new <- q1(x_new)
    }
    print(paste("The root of f(x)=x^2-0.5 is ", x_new, ". f(x) equals to ", f_new, sep = ""))
    print(paste("Using R's commands 1/sqrt(2) equals ", 1/sqrt(2), ".", sep = ""))
    return(-x_new)
  }
}
Q2_Newton(x1 = runif(1,1,2),epsilon = 1/1000000)  # Will not work
Q2_Newton(x1 = runif(1,-1,0),epsilon = 1/1000000) # Will work

# Even though we used the "/" operator in f(x), when we plug the values of f(x) and f'(x)
# into the recurrence relation of the Newton method, we are left with the form
# x_{n+1} = 2x_{n} + 2x_{n}^2, a function that defines the relationship between the 
# value of the current iteration (n+1) to the previous one. This function does not include
# any "\" operation.
# Moreover, in order to converge to the root of a funciton, the Newton method requiers that
# the function will be both convex and (increasingly) monotone. Since x<0, choosing 
# f(x) = 1/x - sqrt(2) yields a decreasing monote function. Eventually I chose
# f(x) = -1/x - sqrt(2). We have that f'(x) = 1/x > 0 and that 
# f''(x) = -2/(x^2) > 0 (recall that x<0). Therefore f(x) fulfills the condision.

# 3. Prove that f(x) = e^x - x - 1 has a unique value at x=0 ----
" 
Preliminary step: using function enquiry on the given funtion in order to show it has a 
 unique root at x = 0:
 f(x) = e^x - x - 1
 f'(x) = e^x - 1
 We ask which are the roots of f(x), that is if f(x) = 0 what are the values of x:
 0=e^x - x - 1 -> e^x - x = 1. 
 The equation above fulfiled only under x = 0. 
 first range: x < 0:
  f'(x) < 0 -> f(x) decreasing in that range
 second range (point): x = 0
  f'(x) = 0, it is an extramum. Let's check its nature using f''(x).
  f''(x) = e^x > 0 that the function is convex and the extramum is minimum.
 thired range: x > 0:
  f'(x)>0 then f(x) increases in that range (exponentialy).
 Therefore the function has a unique global minimum at x = 0, which is also a root since
 f(0) = 1 - 0 - 1 = 0.
 Now, lets validate it with the function we have built:
"

## (i) numericly validation for that the function e^x-x-1 has a unique root at 0 ----
Q3_Newton_i <- function(x1, epsilon) {
  "
  This function returns the root for a convex function with one variable.
  Inputs: 
  - x1: a starting point
  - epsilon: a small number which defines the extant of accuracy.
  "
  
  # Assignment functions
  q1 <- function(x){
    fx <- (exp(x) - x - 1)
    return(fx)
  }
  q1_tag <- function(x, epsilon){
    fx_tag <- (exp(x) - 1)
    return(fx_tag)
  }
  
  
  # Starting conditions
  #x1 <- runif(1,-2,2)
  fx <- q1(x1)
  x_new <- x1
  f_new <- q1(x1)
  fx_tag <- q1_tag(x1)
  
  i <- 1
  
  while (abs(f_new) > epsilon){
    x_old = x_new
    fx_old = q1(x_old)
    fx_tag_old = q1_tag(x_old)
    
    x_new <- x_old - (fx_old/fx_tag_old)
    i <- i + 1
    
    
    # stopping rule
    f_new <- q1(x_new)
  }
  print(paste("The Newton method has converged in ", i, " iterations.",
              "The root is ", x_new, ". f(x) equals to ", f_new, sep = ""))
  return(x_new)
}
Q3_Newton_i(x1 = runif(1,-2,2),epsilon = 1/1000000)

## (ii) Mofify the Newton method ----
Q3_Newton_ii <- function(x1, epsilon) {
  "
  This function returns the root for a convex function with one variable.
  Inputs: 
  - x1: a starting point
  - epsilon: a small number which defines the extant of accuracy.
  "
  
  # Assignment functions
  q1 <- function(x){
    fx <- (exp(x) - x - 1)
    return(fx)
  }
  q1_tag <- function(x, epsilon){
    fx_tag <- (exp(x) - 1)
    return(fx_tag)
  }
  
  
  # Starting conditions
  #x1 <- runif(1,-2,2)
  fx <- q1(x1)
  x_new <- x1
  f_new <- q1(x1)
  fx_tag <- q1_tag(x1)
  
  i <- 1
  
  while (abs(f_new) > epsilon){
    x_old = x_new
    fx_old = q1(x_old)
    fx_tag_old = q1_tag(x_old)
    
    x_new <- x_old - 2*(fx_old/fx_tag_old)
    i <- i + 1
    
    
    # stopping rule
    f_new <- q1(x_new)
  }
  print(paste("The Newton method has converged in ", i, " iterations.",
              "The root is ", x_new, ". f(x) equals to ", f_new, sep = ""))
  return(x_new)
}
Q3_Newton_ii(x1 = runif(1,-2,2),epsilon = 1/1000000)

## (iii) ----
x <- runif(1,-3,4)
Q3_Newton_i(x1 = runif(1,-2,2),epsilon = 1/1000000)
Q3_Newton_ii(x1 = runif(1,-2,2),epsilon = 1/1000000)
'
The second algorithm converges to the root faster so in that sense it works better. 
The difference between the two is that the second one updates the x values in bigger steps. 
A possible downside of updating the x values steps which are too big is that the algorithm
can not converge.
'

# 4. Write a code which implement the Newton method for the function f(x) = arctan(x). ----
Q4_Newton <- function(x1, epsilon, quietly) {
  "
  This function returns the root for a convex function with one variable.
  Inputs: 
  - x1: a starting point
  - epsilon: a small number which defines the extant of accuracy.
  - quetly: a binary variable. If quietly == 1 the function produces no output.
  "
  
  # Assignment functions
  q1 <- function(x){
    fx <- (atan(x))
    return(fx)
  }
  q1_tag <- function(x, epsilon){
    fx_tag <- (1/(x^2 + 1))
    return(fx_tag)
  }
  
  
  # Starting conditions
  #x1 <- runif(1,-2,2)
  fx <- q1(x1)
  x_new <- x1
  f_new <- q1(x1)
  fx_tag <- q1_tag(x1)
  
  i <- 1
  
  while (abs(f_new) > epsilon){
    x_old = x_new
    fx_old = q1(x_old)
    fx_tag_old = q1_tag(x_old)
    
    x_new <- x_old - (fx_old/fx_tag_old)
    i <- i + 1
    
    
    # stopping rule
    f_new <- q1(x_new)
    if(is.nan(x_new)){break}
  }
  if (quietly == 0){
    if(!is.nan(x_new)){
      print(paste("The Newton method has converged in ", i, " iterations.",
                  "The root is ", x_new, ". f(x) equals to ", f_new, sep = ""))
    }
    if (is.nan(x_new)){
      print(paste("After ", i, " iterations the Newton method has NOT converged",
                  "x = ", x_new, ". f(x) equals to ", f_new, "also.", sep = ""))
    }
  }
  return(x_new)
}

convergence_summary = data.frame(matrix(NA,nro = 0, ncol = 2))
names(convergence_summary) <- c("x", "convergence")


for (x in seq(-2,2,0.01)){
  x_out <- Q4_Newton(x1 = x ,epsilon = 1/1000000, quietly = 1)
  if (is.nan(x_out)){
    df_temp = as.data.frame(matrix(c(x, "No"), nrow = 1, ncol = 2))
    names(df_temp) <- names(convergence_summary)
    convergence_summary <- rbind(convergence_summary, df_temp)
  }
  if (!is.nan(x_out)){
    df_temp = as.data.frame(matrix(c(x, "Yes"), nrow = 1, ncol = 2))
    names(df_temp) <- names(convergence_summary)
    convergence_summary <- rbind(convergence_summary, df_temp)
  }
}

print(convergence_summary)


"
The function I wrote converges in the interval -1.39 and 1.39.
In any other value we get that x diverge to infinity.
"

