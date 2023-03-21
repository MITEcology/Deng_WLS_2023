# load functions
source("toolbox.R")

# parameters
num_sys <- 100
num_sp_seq <- c(3)
invader <- 1
n_point_3d <- 5000
sd_seq <- c(0.1)

for(num_sp in num_sp_seq){
  print(paste0("num_sp = ", num_sp))
  
  # n_point scales with the surface area of n-dim hypersphere
  sphere_surface <- 2*pi^(num_sp/2)/gamma(num_sp/2)
  sphere_surface_3d <- 4*pi
  n_point <- round(sphere_surface / sphere_surface_3d * n_point_3d)

  for(sd in sd_seq){
    print(paste0("sd = ", sd))
    
    # simultaneous assembly - sphere
    invasion_simu <- c()
    augmentation_simu <- c()
    exclusion_simu <- c()
    # sequential assembly - feasibility region of residents in isolation
    invasion_sequ <- c()
    augmentation_sequ <- c()
    exclusion_sequ <- c()
    
    for (sys in 1:num_sys){
      print(paste0("# sys = ", sys, "/", num_sys))
      # generate interaction matrix
      A <- random_globally_stable_matrix(num_sp, stren = sd)
      
      # simultaneous assembly - sphere
      simultaneous <- simultaneous_sph(A, invader)
      invasion_simu <- c(invasion_simu, simultaneous$survival)
      augmentation_simu <- c(augmentation_simu, simultaneous$persistence)
      exclusion_simu <- c(exclusion_simu, simultaneous$exclusion)
      
      # sequential assembly - feasibility region of residents in isolation
      sequential <- sequential_iso(A, invader, n_point = n_point)
      invasion_sequ <- c(invasion_sequ, sequential$colonization)
      augmentation_sequ <- c(augmentation_sequ, sequential$augmentation)
      exclusion_sequ <- c(exclusion_sequ, sequential$replacement)
    }
    ## output
    output <- data.frame(invasion_simu, augmentation_simu, exclusion_simu, 
                         invasion_sequ, augmentation_sequ, exclusion_sequ)
  }
}