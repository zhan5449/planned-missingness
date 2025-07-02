######################################################################################
# This code generates density plot PDFs comparing imputed and observed data distributions
# How to run:
#   – Put this script and "Population Correlation Matrices.rds" in the same
#     directory and open/run in RStudio *or* call  `Rscript 01_density_plots.R`
#   – The script makes an output folder called  density_plots/and drops PDFs there.
#####################################################################################

set.seed(20250623)  

# Packages
required <- c("MASS", "mice", "tidyverse", "here")
invisible(lapply(required, require, character.only = TRUE))

# File paths
here::i_am("SF vs PM diagnostics_imputedvobserved.R")   # declares project root
data_file  <- here::here("Population Correlation Matrices.rds")
out_dir    <- here::here("density_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load population correlation matrices
P <- readRDS(data_file)

# Two representative conditions: 2-factor with population correlation = 0.4, N = 500
# and 5-factor with cor = 0.6 and N = 300.
conditions <- list(
  list(label = "2factor_rho0.4_n500", P_idx = "P6", n_factors = 2, n = 500),
  list(label = "5factor_rho0.6_n300", P_idx = "P174", n_factors = 5, n = 300)
)

# Loop through each condition x all four missingness levels
for (cond in conditions) {
  for (n_miss in c(2, 4, 6, 8)) { 
    
    # Extract the population matrix 
    pop_mat <- P[[cond$P_idx]]
    n_items <- ncol(pop_mat)
    vars <- paste0("X", 1:n_items)
    
    # Simulate data from multivariate normal population
    dat <- MASS::mvrnorm(n = cond$n, mu = rep(0, n_items), Sigma = pop_mat) |>
      as.data.frame()
    names(dat) <- vars
    
    # Impose planned missingness by sampling N items missing per factor
    dat_pm <- lapply(1:nrow(dat), function(i) {
      x <- dat[i, ]
      idx <- unlist(lapply(0:(cond$n_factors - 1), function(f) {
        sample(1:10 + 10 * f, size = n_miss, replace = FALSE)
      }))
      x[idx] <- NA
      return(x)
    }) |> bind_rows()
    names(dat_pm) <- vars
    
    # Set number of iterations and imputations 
    m <- n_miss*10
    maxit <- 3
    imp <- mice(dat_pm, method = "pmm", m = m, maxit = maxit, 
                print = FALSE, seed = 20250623)
    
    # Set directory to save density plots to PDF
    out_file <- paste0("density_plots/",
                       cond$label, "_miss", n_miss*10, ".pdf")
    # Convert mice object to long format (includes imputed and observed data)
    imp_long <- complete(imp, action = "long", include = TRUE)
    # Open pdf
    pdf(out_file, width = 8, height = 6)
    # Loop over variables and generate plots
    for (var in vars) {
      p <- ggplot(imp_long, aes_string(x = var)) +
        geom_density(data = filter(imp_long, .imp != 0), aes(group = .imp), 
                     alpha = 0.2, color = "red") +  # imputed, semi-transparent
        geom_density(data = filter(imp_long, .imp == 0), 
                     color = "blue", size = 1.2) +  # observed, bold
        labs(title = paste("Density for", var, "-", cond$label),
             subtitle = paste(n_miss * 10, "% Missing"),
             x = var, y = "Density") +
        theme_minimal(base_size = 12)
      
      print(p)
    }

    # Close PDF
    dev.off()
  }
}
