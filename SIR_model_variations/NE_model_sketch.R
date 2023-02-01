library(deSolve)

NE_derv <- function(time, x, parameters) {
  with(as.list(c(x, parameters)), {
    x1 <- -rr*x[2]*x[1] 
    x2 <-  rr*x[3]*x[2]*x[1]*calc_g2(lam,x[1])/calc_g1(lam,x[1])-rr*x[2]*(1-x[2])-x[2]*mm+pp*(x[5]-x[2])
    x3 <-  rr*x[3]*x[2]*(1-x[1]*calc_g2(lam,x[1])/calc_g1(lam,x[1]))+pp*(x[1]*calc_g1(lam,x[1])/calc_g1(lam,1)-x[3])
    x4 <- -rr*x[2]*x[1]*calc_g1(lam,x[1]) 
    x5 <- -mm*x[5]+rr*x[2]*(x[1]**2*calc_g2(lam,x[1])+x[1]*calc_g1(lam,x[1])/calc_g1(lam,1)) 
    x6 <- rr*x[2]*x[1]*calc_g1(lam,x[1])-mm*x[6]
    return(list(c(x1,x2, x3, x4, x5, x6)))
  })
}

calc_g <- function(lamda, x)(0 + 0.5*x + 0*x**2 + 0.5*x**3)

calc_g1 <- function(lamda, x) ( 0.5 + 2*0*x + 3*0.5*x**2)

calc_g2 <- function(lamda, x)(3*0.5*2*x)

lamda = 2.5
r = 0.2
mu = 0.2
ro = 0.2
epsilon = 0.001
time = seq(0,200)

parameters_values <- c(
  lam = lamda,
  rr = r,
  mm = mu,
  pp = ro,
  calc_g = calc_g, 
  calc_g1 = calc_g1, 
  calc_g2 = calc_g2
)

initial_vals <- c(x1 = 1 - epsilon, # proportion susceptible at start
                  x2 = epsilon / (1 - epsilon), #
                  x3 = (1 - 2*epsilon)/(1 - epsilon), #
                  x4 =calc_g(lamda, 1 - epsilon), #
                  x5 = epsilon, #
                  x6 = 1 - calc_g(lamda, 1- epsilon) )

sir_values_1 <- ode(
  y = initial_vals,
  times = time,
  func = NE_derv,
  parms = parameters_values 
)

sir_values_1 <- as.data.frame(sir_values_1)

with(sir_values_1, {
  # plotting the time series of susceptibles:
  plot(time, x4, type = "l", col = "blue",
       xlab = "time (days)", ylab = "proportion of people",
       ylim = c(0,1))
  # adding the time series of infectious:
  lines(time, x6, col = "red")
  # adding the time series of recovered:
  lines(time, 1 - (x4+x6), col = "green")
  
  # adding a legend:
  legend("right", c("susceptibles", "infectious", "recovered"),
         col = c("blue", "red", "green"), lty = 1, bty = "n")
  
})


sir_values_2 <- sir_values_1 * 10000
sir_values_2$time <- sir_values_1$time

with(sir_values_2, {
  # plotting the time series of susceptibles:
  plot(time, x4, type = "l", col = "blue",
       xlab = "time (days)", ylab = "number of people",
       ylim = c(0,10000))
  # adding the time series of infectious:
  lines(time, x6, col = "red")
  # adding the time series of recovered:
  lines(time, 10000 - (x4+x6), col = "green")
  
  # adding a legend:
  legend("right", c("susceptibles", "infectious", "recovered"),
         col = c("blue", "red", "green"), lty = 1, bty = "n")
  
})

