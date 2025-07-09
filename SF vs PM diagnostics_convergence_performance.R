######################################################################################
# This code generates convergence and diagnostic measures, including plots
# How to run:
#  – Put this script and "Population Correlation Matrices.rds" and "Convergence Diagnostic Values.rds" in the same
#     directory
#  – The script makes an output folder /convergence_diagnostic_plots.
#  - The first part of the script creates 'Convergence Diagnostic Values.rds'
#         # This part of the script is computationally very expensive (commented out). 
#         # Running 16 sample conditions across hundreds of iterations takes over 13 hours on an 8 core machine
#         # To skip this part inspect the diagnostic values and generate the figures, skip to line 337, load 'Convergence Diagnostic Values.rds' and run only the plotting code below line 
#####################################################################################

set.seed(20250623)  

# Packages 
required <- c("here", "MASS", "mice", "miceadds", "psych", "tidyverse", "doParallel", "parallel","doRNG","ggplot2", "patchwork")
invisible(lapply(required, require, character.only = TRUE))
# NOTE: mice version 3.18.0 and later introduce a dependency on 'rstan' for convergence()
# This breaks compatibility on some systems without rstan installed.
# To remain platform-agnostic we:
#   #        • compute PSRF manually using the Gelman–Rubin formula
#   #          (identical in closed form to coda::gelman.diag(); a spot-check
#   #          on the P6 | n = 300 | 20 % missing cell showed Δ = 0.01).
#   #        • compute lag-1 autocorrelation (AC) with stats::acf().


# File paths
input_data_file  <- here::here("Population Correlation Matrices.rds")
output_data_file  <- here::here("Convergence Diagnostic Values.rds")
out_dir    <- here::here("convergence_diagnostic_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)



############################ Define sample conditions and run diagnostics #########################################

# # Load population correlation matrices
# P <- read_rds(input_data_file)
# 
# 
# # Define grid of simulation conditions for convergence diagnostics
# # This grid crosses:
# # - 2- and 5-factor population matrices (from P) with known factor structures (true interfactor cor .4)
# # - Sample sizes (n = 300, 1000)
# # - Proportion of items missing per factor (20%, 40%, 60%, and 80%). 
# condition_grid <- expand.grid(
#   P_idx      = c("P6",       # 2-factor: pop_cor ≈ .4
#                  "P126"  # 5-factor
#   ),
#   n          = c(300, 1000),         # Sample sizes
#   missing_pc = c(0.2, .40, 0.6, 0.8),    # Percent missingness per factor
#   stringsAsFactors = FALSE
# ) |>
#   mutate(
#     n_factors = case_when(
#       P_idx == "P6"     ~ 2,
#       P_idx == "P126" ~ 5,
#       TRUE ~ NA_integer_
#     ),
#     
#     # Compute number of items missing per factor (10 items per factor in all cases)
#     n_miss = as.integer(missing_pc * 10), 
#     
#     # Number of imputations set to the missingness rate 
#     m = n_miss*10,
#     
#     # More iterations if missingness is ≥ 60%
#     maxit = ifelse(missing_pc < 0.6, 100, 250)
#   ) |>
#   
#   # Add a column for the actual average population factor correlation
#   rowwise() |>
#   mutate(
#     pop_cor = {
#       mat <- P[[P_idx]]          # Load correlation matrix
#       k <- ncol(mat)             # Total number of observed items
#       f <- n_factors             # Number of latent factors
#       
#       # Create binary keys matrix for computing factor composites
#       # Each factor loads on its own 10 items
#       keys <- matrix(0, nrow = k, ncol = f)
#       for (i in 1:f) {
#         keys[((i - 1) * 10 + 1):(i * 10), i] <- 1
#       }
#       
#       # Use psych::cluster.cor() to estimate factor-level correlation matrix
#       cor_mat <- psych::cluster.cor(keys = keys, r.mat = mat)$cor
#       
#       # Average the off-diagonal elements to get mean population factor correlation
#       mean(cor_mat[lower.tri(cor_mat)])
#     }
#   ) |>
#   ungroup()
# 
# 
# 
# # Simulate and impose missingness
# sim_data_list <- pmap(condition_grid, function(P_idx, n, missing_pc, n_factors, n_miss, m, maxit, pop_cor) {
#   pop_mat <- P[[P_idx]]
#   n_items <- ncol(pop_mat)
#   vars    <- paste0("X", 1:n_items)
#   
#   # Simulate complete data from multivariate normal population
#   # (based on item-level correlation matrix and standard normal marginals)
#   dat_val <- MASS::mvrnorm(n = n, mu = rep(0, n_items), Sigma = pop_mat) |>
#     as.data.frame()
#   names(dat_val) <- vars
#   
#   # Impose planned missingness: for each respondent, drop `n_miss` items per factor
#   # Items assumed to be ordered in blocks of 10 per factor
#   dat_val_pm <- map_dfr(1:n, function(i) {
#     x <- dat_val[i, ]
#     idx <- unlist(map(0:(n_factors - 1), ~ sample(1:10 + 10 * .x, size = n_miss)))
#     x[idx] <- NA
#     return(x)
#   })
#   names(dat_val_pm) <- vars
#   
#   list(
#     P_idx       = P_idx,
#     n           = n,
#     missing_pc  = missing_pc,
#     n_factors   = n_factors,
#     pop_cor     = pop_cor,
#     dat_val     = dat_val,
#     dat_val_pm  = dat_val_pm,
#     m           = m,
#     maxit       = maxit
#   )
# })
# 
# 
# run_diagnostics <- function(condition, simulation_seed = 20250623) {
#   set.seed(simulation_seed)
#   
#   dat_val_pm <- condition$dat_val_pm
#   n_factors  <- condition$n_factors
#   m          <- condition$m
#   maxit      <- condition$maxit
#   n          <- condition$n
#   P_idx      <- condition$P_idx
#   pop_mat    <- P[[P_idx]]
#   
#   # Population correlation from true factor composites
#   k <- ncol(pop_mat)
#   keys <- matrix(0, nrow = k, ncol = n_factors)
#   for (i in 1:n_factors) {
#     keys[((i - 1) * 10 + 1):(i * 10), i] <- 1
#   }
#   pop_cor <- psych::cluster.cor(keys = keys, r.mat = pop_mat)$cor
#   cor_truth <- mean(pop_cor[lower.tri(pop_cor)]) 
#   
#   # Preallocate diagnostics dataframe
#   diagnostic_df <- data.frame(
#     P_idx        = rep(condition$P_idx, maxit),
#     n            = rep(condition$n, maxit),
#     missing_pc   = rep(condition$missing_pc, maxit),
#     n_factors    = rep(n_factors, maxit),
#     iteration_no = 1:maxit,
#     cor_pop      = rep(cor_truth, maxit),
#     cor_pooled   = numeric(maxit),
#     PB           = numeric(maxit),
#     ac           = numeric(maxit),
#     psrf         = numeric(maxit)
#   )
#   
#   # Run MICE iteratively
#   for (i in 1:maxit) {
#     
#     imp <- if (i == 1) {
#       mice(dat_val_pm, method = "pmm",
#            m = m, maxit = 3, print = FALSE, seed = simulation_seed)
#     } else {
#       mice.mids(imp, maxit = 1, print = FALSE)
#     }
#     
#     diagnostic_df$iteration_no[i] <- imp$iteration
#     
#     # Create composite scores from imputed data
#     imp_long <- complete(imp, action = "long", include = TRUE) 
#     for (j in 1:n_factors) {
#       item_range <- ((j - 1) * 10 + 1):(j * 10)
#       imp_long[[paste0("c", j)]] <- rowMeans(imp_long[, item_range]) 
#       
#     }
#     imp_pm <- as.mids(imp_long)
#     comp_vars <- paste0("c", 1:n_factors)
#     cor_mi <- miceadds::micombine.cor(imp_pm, variables = comp_vars)
#     
#     pooled_cor <- mean(cor_mi$r[lower.tri(cor_mi$r)])
#     
#     diagnostic_df$cor_pooled[i] <- pooled_cor
#     diagnostic_df$PB[i]         <- 100 * abs((pooled_cor - cor_truth) / cor_truth)
#     
#     # Manual calculations of convergence diagnostics. Removes mice::convergence() rstan dependence
#     chain_mean <- imp$chainMean                 # dims: m × iter x vars
#     cur_it     <- imp$iteration
#     n_ch       <- dim(chain_mean)[1L]           # m chains
#     n_var      <- dim(chain_mean)[3L]           # 20 or 50 variables 
#     
#     psrf_vals <- ac_vals <- numeric(n_var) 
#     
#     if (cur_it == 1L) { # mimic mice::convergence() behaviour
#       diagnostic_df$psrf[i] <- NA_real_
#       diagnostic_df$ac[i]   <- NA_real_
#       next            # skip to next iteration
#     }
#     
#     for (v in 1:n_var) {
#       cm <- chain_mean[ , 1:cur_it, v, drop = FALSE]  # slice up to cur_it
#       
#       # Skip variable if the entire slice is NA 
#       if (all(is.na(cm))) {
#         psrf_vals[v] <- NA_real_
#         ac_vals[v]   <- NA_real_
#         next
#       }
#       
#       ##  PSRF (R-hat) Same formula as coda::gelman.diag(autoburnin=FALSE). Coda not required here.
#       W <- mean(apply(cm, 1, var, na.rm = TRUE))                    # within-chain variance
#       B <- var(rowMeans(cm, na.rm = TRUE)) * cur_it                 # between-chain variance × n
#       var_hat  <- ((cur_it - 1) / cur_it) * W + B / cur_it
#       psrf_vals[v] <- sqrt(var_hat / W)
#       
#       
#       ##  Lag-1 Autocorrelation
#       ac_chain <- apply(cm, 1, function(series) {
#         if (cur_it > 1)
#           acf(series, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2] #function in stats package
#         else
#           NA_real_
#       })
#       ac_vals[v] <- mean(ac_chain, na.rm = TRUE)
#     }
#     
#     diagnostic_df$psrf[i] <- mean(psrf_vals, na.rm = TRUE)
#     diagnostic_df$ac[i]   <- mean(ac_vals,  na.rm = TRUE)
#     
#   }
#   
#   return(diagnostic_df)
# }
# 
# 
# 
# ## Run the function in parallel 
# 
# # Choose  number of workers
# n_workers <- min(parallel::detectCores() - 1L,    # leave one core free
#                  length(sim_data_list))           
# 
# cl <- parallel::makeCluster(n_workers, type = "PSOCK", outfile = "")
# 
# # One-time setup per worker
# parallel::clusterEvalQ(cl, {
#   library(mice); library(miceadds); library(psych); library(MASS)
#   options(mc.cores = 1)            # prevent nested threading
#   NULL
# })
# 
# parallel::clusterExport(cl, c("run_diagnostics", "P"))
# 
# doParallel::registerDoParallel(cl)
# 
# simulation_seed <- 20250623        # global stream for doRNG
# 
# start_time <- Sys.time()
# 
# diag_imp <- foreach(i = seq_along(sim_data_list),  #seq_along(sim_data_list)                     
#                     .packages    = c("mice","miceadds","psych","MASS"
#                                      ),
#                     .export      = c("run_diagnostics","P"), #ignore warning
#                     .combine     = "list",
#                     .options.RNG = simulation_seed,
#                     .errorhandling = "pass") %dorng% {
#                       run_diagnostics(sim_data_list[[i]],
#                                       simulation_seed = simulation_seed)  
#                     }
# 
# elapsed <- difftime(Sys.time(), start_time, units = "mins")
# print(elapsed)
# 
# 
# parallel::stopCluster(cl)
# 
# 
# # Output is deeply nested list object. Unlist into dataframe before saving.
# ## Helper function to dig deep into each nested list and extract dfs
# get_dfs <- function(x) {
#   if (is.data.frame(x)) {
#     list(x)                    # wrap so we can c() later
#   } else if (is.list(x)) {
#     # recurse into each element and c() results
#     purrr::map(x, get_dfs) |> purrr::flatten()
#   } else {
#     list()                     # ignore non-list, non-df leaves
#   }
# }
# 
# ## Apply helper and bind 
# df_list   <- get_dfs(diag_imp)           # pure list of data frames of diagnostic results for all conditions
# diag_flat <- bind_rows(df_list)
# 
# 
# #Add a flag for each condition where computation of diagnostics failed
# metric_cols  <- c("cor_pooled","PB","ac","psrf")
# 
# fail_flags <- diag_flat %>%
#   group_by(P_idx, n, missing_pc) %>%
#   summarise(
#     across(all_of(metric_cols),
#            ~ all(is.na(.x) | is.nan(.x)),
#            .names = "allNA_{.col}"),
#     .groups = "drop"
#   ) %>%
#   # collapse TRUE metrics into one semicolon-separated string
#   mutate(flag = pmap_chr(
#     select(., starts_with("allNA_")),
#     ~ {
#       bad <- metric_cols[c(...) == TRUE]
#       if (length(bad) == 0) "OK" else paste(bad, collapse = ";")
#     })
#   )
# 
# # join back 
# diag_flat <- diag_flat %>% 
#   left_join(fail_flags, by = c("P_idx","n","missing_pc")) %>% 
#   select(-c(starts_with("allNA_")))
# 
# ### Save diagnostic ouput
# write_rds(diag_flat, output_data_file) #Larde df with all metrics saved for each condition x iteration


################################### Plotting ###################################################

### Read diagnostic output file
diag_flat <- read_rds(output_data_file)

# Plotting function
make_panel <- function(df_cell) {
  
  ## Convert to long format
  metrics_long <- df_cell |>
    pivot_longer(
      cols      = c(ac, psrf, cor_pooled, PB),
      names_to  = "metric",
      values_to = "value"
    ) |>
    mutate(
      metric = factor(metric,
                      levels = c("ac", "psrf", "cor_pooled", "PB"))
    )
  
  ## Flag metrics that are all-NA in this cell (failed)
  na_flags <- metrics_long %>%
    group_by(metric) %>%
    summarise(all_na = all(is.na(value) | is.nan(value)), .groups = "drop")
  
  metrics_long <- left_join(metrics_long, na_flags, by = "metric")
  
  ## Thresholds 
  pop_cor <- unique(df_cell$cor_pop)[1]
  
  meta <- tribble(
    ~metric,        ~ymin, ~ymax, ~yint1,   ~yint2,
    "ac",                 0,     1,    0,        NA,
    "psrf",               1, 3.5, 1.01,   1.10,
    "cor_pooled",       0.2,   0.6, pop_cor,     NA,
    "PB",                 0,  15,    0,        5
  ) |>
    mutate(
      metric = factor(metric, levels(metrics_long$metric))
    )
  
  ## Blank points (one min + one max per facet)
  limits_long <- meta |>
    pivot_longer(cols = c(ymin, ymax),
                        names_to = NULL, values_to = "value")
  
  ## Facet labels 
  base_labels <- c(
    ac         = "'Autocorrelation'",
    psrf       = "hat(R)",
    cor_pooled = "'Average Imputed Factor Correlation'",
    PB         = "'Percent Bias (%)'"
  )
  na_suffix <- ifelse(na_flags$all_na, " *'(NA)*'", "")
  facet_labs <- as_labeller(
    setNames(paste0(base_labels, na_suffix), levels(metrics_long$metric)),
    default = label_parsed
  )

  ## Plot
  ggplot(metrics_long %>%
           filter(!(is.na(value) | is.nan(value))),  # drop NA *and* NaN
         aes(iteration_no, value,
             linetype = factor(missing_pc)))  +
    geom_line() +
    facet_wrap(
      vars(metric),
      nrow = 2, ncol = 2,
      scales = "free_y",
      labeller = labeller(metric = facet_labs)
    ) +
    geom_hline(data = meta,
               aes(yintercept = yint1),
               linetype = "dashed", colour = "grey") +
    geom_hline(data = dplyr::filter(meta, !is.na(yint2)),
               aes(yintercept = yint2),
               linetype = "dashed", colour = "grey") +
    geom_blank(
      data        = limits_long,
      inherit.aes = FALSE,
      aes(x = -Inf, y = value, metric = metric)
    ) +
    ## grey facet strip when metric missing
    theme_bw(base_size = 11) +
    theme(
      legend.position   = "bottom",
      strip.background.x = element_rect(
        fill = ifelse(na_flags$all_na, "grey80", "grey90")
      )
    ) +
    labs(
      linetype = "Missing %",
      x        = "Iteration",
      y        = NULL
    ) 
}



## Loop over 16 conditions and save plots to out_dir
diag_flat %>%
  distinct(P_idx, n, n_factors) %>%   # collapses across missing %
  pwalk(function(P_idx, n, n_factors) {
    
    df_cell <- diag_flat %>%
      filter(P_idx == !!P_idx, n == !!n)
    
    p  <- make_panel(df_cell)
    
    fn <- sprintf("%s_n%d_f%d_panel.png", P_idx, n, n_factors)
    ggsave(
      here::here("convergence_diagnostic_plots", fn),
      p, width = 8, height = 6, dpi = 300
    )
  })





