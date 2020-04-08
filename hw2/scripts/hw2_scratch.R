
sim_mites_matrix <- function(d, num_experiments, size, p_dose) {
  print(paste("dose: ", d))
  print(num_experiments[d])
  print(size)
  print(p_dose[d])
  
  # Generate a dataframe of simulated data with the same dimension as the mitesA dataframe
  sim_dead <- rbinom(num_experiments[d], size, p_dose[d])
  
  loop_df <- as.matrix(cbind(dose = rep(as.numeric(dose_vec[d])), num_mites=rep(size), num_dead = sim_dead))
}

sim_mites_data <- function(num_experiments, size, p_dose, num_dose_level){
  
  return(data.frame(do.call(rbind, lapply(1:num_dose_level, sim_mites_matrix, num_experiments, 10, p_dose))))

}

test_output<- sim_mites_data(num_experiments, 10, p_dose, num_dose_level)



rbinom_autocorr <- function(n, size, p, lambda){
  q = p
  sample <- c()
  sapply(1:n, )
  for(i in 1:n){
    #print(p*lambda + (1-lambda)*q)
    s <- rbinom(1, size, p*lambda + (1-lambda)*q)
    sample <- c(sample, s)
    q <- s/size
  }
  return(sample)
}