# install.packages('pomp')
library(pomp)
library(ggplot2)

# implementation of a stochastic differential equation

# https://kingaa.github.io/short-course/stochsim/stochsim.html


bsflu <- read.table("http://kingaa.github.io/short-course/stochsim/bsflu_data.txt") 

# bsflu <- as.Date()
head(bsflu)
# B refers to boys confined to bed and C to boys in convalescence.
ggplot(data =bsflu) + 
  geom_line(aes(x = day, y = B, color = 'blue')) +
  geom_line(aes(x = day, y = C)) + 
  labs(x = 'Day',
       y = 'Count',
       title = 'Flu Counts') +
  theme_bw() 

sir_step <- Csnippet("
                     double dN_SI = rbinom(S,1 - exp(-Beta*I/N*dt));
                     double dN_IR = rbinom(I,1 - exp(-gamma*dt) );
                     S -= dN_SI;
                     I += dN_SI - dN_IR;
                     R += dN_IR;
                     ") 

# at day 0 we will assume that I(0) = 1, R() = 0, we dont know how big the school is so we treat N as 
# a param to be estiamted  and let S(0) = N - I

sir_init <- Csnippet("S = N-1;
                     I = 1;
                     R = 0;")
# we fold these csnippets into pompwith the data
sir <- pomp::pomp(bsflu, 
           time = 'day', 
           t0 = 0, 
           rprocess = euler(sir_step, delta.t = 1/6),
           rinit = sir_init,
           paramnames = c('N', "Beta","gamma"),
           statenames = c("S","I","R"))

# now let us add in the H class

sir_step_h <-  Csnippet("
                     double dN_SI = rbinom(S,1 - exp(-Beta*I/N*dt));
                     double dN_IR = rbinom(I,1 - exp(-gamma*dt) );
                     S -= dN_SI;
                     I += dN_SI - dN_IR;
                     R += dN_IR;
                     H += dN_IR;
                     ") 

sir_init_h <- Csnippet("S = N-1;
                     I = 1;
                     R = 0;
                     H = 0;")
sirh <- pomp::pomp(sir, 
                   rprocess=euler(sir_step_h, delta.t = 1.6),
                   rinit  = sir_init_h,
                   paramnames = c("Beta",'gamma',"N"),
                   statenames = c("S","I","R","H"))

sir_h_acc <- pomp(sirh,accumvars="H")

rmeas <- Csnippet("B = rbinom(H, rho);")

sir_measure <- pomp(sir_h_acc, 
                    rmeasure = rmeas, 
                    statenames = "H",
                    paramnames = "rho")

sims <- simulate (sir_measure, params = c(Beta=1.5,gamma=1,rho=0.9,N=2600),
                  nsim = 20, format = "data.frame",
                  include.data = TRUE)
#c loser  geom_line() +
  theme_bw()

