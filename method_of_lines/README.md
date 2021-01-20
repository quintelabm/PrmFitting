# Intra-host viral dynamics

*Goal:* Represent the replication of HCV RNA inside the infected cells

### ODE

- Virus 

  dV/dt = rho( integral(I) . Integral(Rp+Rt) ) - cV

- Target Cells 
    
    dT/dt = s - dT - betaVT

### Advection Equation

- Infected cells 

  dI/dt+dI/da = -delta I

- RNA 
    - translation 

      dRt/dt+dRt/da = theta Rp -(sigma+rho+mu_t)Rt 
      
    - replication positive 

       dRp/dt+dRp/da = 
       alpha Rn + sigma Rt - (theta+rho+mu_c)Rp

    - replication negative 

       dRn/dt+dRn/da = r(1-Rn/Rmax))Rp - (mu_c)Rn


## Solve with Method of Lines

- using finite differences

- Based on example with [Lax Friedrichs](http://people.bu.edu/andasari/courses/numericalpython/Week11Lecture18/laxfriedrichs_periodic.py)