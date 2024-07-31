# inverse PROSPECT
library("Rprospect")

# load the data
reflectance <- read.csv("C://Users/Julian/Desktop/spectra/digitized/malus_domestica.csv", header = T)
#reflectance <- reflectance[grepl("PSME" , reflectance$Species),]

#reflectance_s <- data.frame(x = as.integer(sub("X","",names(reflectance)[grep("X", names(reflectance))])), y = as.numeric(reflectance[1,grep("X", names(reflectance))]))


#reflectance <- rbind(reflectance, c(400, 4.9))

reflectance_s <- spline(x = reflectance[,1], y = reflectance[,2], n = 1000, xout = 360:2400)

plot(reflectance_s)

# calculate N prior
Nprior_R <- prospect::Get_Nprior(lambda = reflectance_s$x, 
                       Refl = reflectance_s$y)

# gap filled trait data
trait_data <- read.csv("C://Users/Julian/Desktop/spectra/try/gap_filled_final_with_names.csv")

# get trait data
trait <- trait_data[trait_data$Species == "Malus domestica",c("chlc_mygcm2","plant_height_m", "LMA_gm2", "EWT_mgcm2", "antc_gcm2", "caroc_gcm2")]
t_names <- names(trait)
#trait <- apply(trait,2, weighted.mean, w= trait$plant_height_m, na.rm = T )
trait <- apply(trait,2, mean, na.rm = T )
names(trait) <- t_names

# calculate spectra
spectra <- Rprospect::prospect4(N = 1.63, trait["chlc_mygcm2"], trait["EWT_mgcm2"]/1000, trait["LMA_gm2"]/10000)

optim_fn <- function(N){
  spectra <- Rprospect::prospect4(N = N, trait["chlc_mygcm2"], trait["EWT_mgcm2"]/1000, trait["LMA_gm2"]/10000)
  dist_df <- merge(spectra, reflectance_s, by.x = "Wavelength", by.y = "x")
  return(sqrt(mean((dist_df$Reflectance - dist_df$y)^2)))
}

optimal_N <- optim(Nprior_R, optim_fn, method = "Brent", lower = 0, upper = 5)$par

spectra <- Rprospect::prospect4(N = optimal_N, trait["chlc_mygcm2"], trait["EWT_mgcm2"]/1000, trait["LMA_gm2"]/10000)

# plot the spectra 
plot(spectra$Wavelength, spectra$Reflectance, type = "l", col = "red", ylim = c(0,0.6))
lines(reflectance_s$x, reflectance_s$y, col = "blue")

# write the spectrum to a file in DART standard
spectra$Wavelength <- spectra$Wavelength / 1000
spectra <- with(spectra, data.frame(wavelength, reflectance, direct_transmittance = 0.0, diffuse_transmittance))
header <- paste("* # Simualted data using Prospect4 model and try information \n* #by Julian Frey WWD University Freiburg\n* # AccSpeciesName Malus_domestica\n* # N", optimal_N, "\n* # CHL", trait["chlc_mygcm2"], "\n* # EWT", trait["EWT_mgcm2"]/1000, "\n* # LMA", trait["LMA_gm2"]/10000, collapse = " ")

writeLines(
  c(
    header,
    paste(colnames(spectra), collapse = "\t"), 
    apply(spectra,1,function(x)paste(format(as.numeric(x), digits = 6, trim = T, scientific = F), collapse  = "\t"))
    ), 
  paste( "c:/Users/Julian/Desktop/spectra/try/spectra/Malus_domestica.txt", sep = "")
  )
