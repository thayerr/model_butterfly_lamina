# Title: thin_films.R
# Author: Rachel C. Thayer, in collaboration with Sam Thayer
# Date: 2019
# Publication: Thayer et al 2020, "Structural color in Junonia butterflies evoves by tuning scale lamina thickness," eLife

# See also: Stavenga 2014 "Thin film and multilayer optics cause structural color in many insects and birds," Mat. Today:Proc

# Purpose: model reflectance off of chitin thin films that compose the lower lamina of Junonia genus butterfly scales.

########################################################################
## DEFINE BASELINE VARIABLES

# wavelengths of light
lambda = seq(300, 700, by=1)

# angles of light that are collected by your microscope objective. 0-30º inclusive is appropriate for an objective lens with Numerical Aperture = 0.5
# granularity here MUST match indexing inside the 'averaging' function later!
angles = seq(0, 30, by=0.03)

# thickness (average and standard deviation, as measured from electron micrographs)
# realistic thicknesses range from about 100 nm for a very thin lower lamina, to about 260 nm for a very thick one
# these examples are for blue dorsal hindwing scales from Junonia hierta
d_mu = 183 #average lamina thickness in nm from replicated measures
d_sigma = 12 #standard deviation in nm among replicated lamina thickness measures
d_list <- seq(d_mu-d_sigma, d_mu+d_sigma, by=0.5) #Use inclusive list of measured thickness values to allow for measurement variation and/or error in the modeled reflectance

# refractive indices
n_chitin = 1.56 + (8.80 * 10^3)/ (lambda^2) # RI of chitin was estimated from glass butterfly scales by Leertouwer 2011
n_air = 1
k = 2*pi/lambda

########################################################################
## DEFINE COMPONENT FUNCTIONS

# this section defines the calculation for the unpolarized reflectance for a single specific angle
# input theta must be given in degrees, not radians
unpolarized_reflectance <- function(theta, d){ 
  theta0 = theta *pi/180 #convert to radians. Input value must be in degrees!
  theta1 = asin((n_air*sin(theta0))/n_chitin)
  theta2 = theta0 #Proof: thetatwo = asin((n_air*sin(thetazero))/n_air)
  phi = k * n_chitin * d * cos(theta1)
  TE_R01_num = n_air * cos(theta0) - n_chitin * cos(theta1)
  TE_R01_den = n_air * cos(theta0) + n_chitin * cos(theta1)
  TE_R01 = TE_R01_num / TE_R01_den
  TE_R12_num = n_chitin * cos(theta1) - n_air * cos(theta2)
  TE_R12_den = n_chitin * cos(theta1) + n_air * cos(theta2)
  TE_R12 = TE_R12_num / TE_R12_den
  
  R_num_TE = TE_R01^2 + TE_R12^2 + 2 * TE_R01 * TE_R12 * cos(2 * phi)
  R_den_TE = 1 + TE_R01^2 * TE_R12^2 + 2 * TE_R01 * TE_R12 * cos(2*phi)
  R_TE = R_num_TE / R_den_TE
  
  TM_R01_num = n_air * cos(theta1) - n_chitin * cos(theta0)
  TM_R01_den = n_air * cos(theta1) + n_chitin * cos(theta0)
  TM_R01 = TM_R01_num / TM_R01_den
  TM_R12_num = n_chitin * cos(theta2) - n_air * cos(theta1)
  TM_R12_den = n_chitin * cos(theta2) + n_air * cos(theta1)
  TM_R12 = TM_R12_num / TM_R12_den
  
  R_num_TM = TM_R01^2 + TM_R12^2 + 2 * TM_R01 * TM_R12 * cos(2 * phi)
  R_den_TM = 1 + TM_R01^2 * TM_R12^2 + 2 * TM_R01 * TM_R12 * cos(2*phi)
  R_TM = R_num_TM / R_den_TM
  
  R_unpolarized = (R_TE + R_TM) /2
  return(R_unpolarized)
}

# this section defines the calculation for averaging reflectance from 0 degrees to the maximal angle, mimicking spectrophotometer readings through a microscope objective lens
# assumes constant light flux, meaning that we assume there is not more total light reaching the probe from 0º than from any other exact angle
averaging <- function(d, maximal_angle){
  raw_r <- t(sapply(angles, unpolarized_reflectance, d=d))
  radiuses <- angles*pi/180 #in radians
  circumferences <- raw_r * 2*pi*radiuses
  inclusive <- circumferences[1:(maximal_angle/.03+1),] #divide by granularity I set at top
  if(is.array(inclusive)){
    sumz <- colSums(inclusive)
    reflectance <- sumz /(2*pi*sum(radiuses[1:(maximal_angle/.03+1)]))
    return(reflectance)
    } else { return(raw_r[1,]) } #theta0 = 0 is a special case. Don't average anything.
}

########################################################################
## Call to apply the component functions

#the output 'average_me' has one column of modeled reflectance per thickness value in d_list
average_me <- as.data.frame(sapply(d_list, averaging, maximal_angle=30))

#useful formatting for graphing the results
model_reflectance <- as.data.frame(lambda)
model_reflectance$type <- "sample_name"
model_reflectance$mean <- apply(average_me, 1, mean)
model_reflectance$min <- apply(average_me, 1, min)
model_reflectance$max <- apply(average_me, 1, max)

#simple plotting
plot(model_reflectance$lambda, model_reflectance$mean)

#plotting reflectance for the range of possible reflectances, accounting for measurement error
library(ggplot2)
ggplot(data=model_reflectance, aes(x=lambda, y=mean, ymin=min, ymax=max)) + 
  geom_line() + 
  geom_ribbon() + 
  xlab("wavelength (nm)") + 
  ylab("reflectance")

########################################################################
### This section calculates the Gaussian-sampled reflectance for a single uneven film. 
# In other words, rather than treating the empirically measured thickness variation as estimate error of a lamina with constant thickness 
  # (calculating reflectance off of films with a uniform thickness of d_mu+ 1 standard deviation, and separately for a film with d_mu - 1 standard deviation), 
  # we assume that one single scale lamina has a variable thickness with average d_mu and d_sigma
# May run slowly
# similar to Siddique 2016 and Stavenga 2014 "Coloration principles of nymphaline butterflies..." figure 4

# first, generate a list of thickness by drawing 400 randomized values from the thickness distribution defined by the empirically measured mu and sigma
d_Gaus <- rnorm(400, d_mu, d_sigma) #n, d_mu, d_sigma

# calculate reflectance for those thicknesses
average_me_Gaus <- sapply(d_Gaus, averaging, maximal_angle=30)
# average the reflectances and plot
plot(lambda, apply(average_me_Gaus, 1, mean), ylim = c(0,0.25), type = "l", col="red", main="Gaussian model")
