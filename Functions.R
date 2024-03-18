#' Show Loop Progress
#'
#' Print progress bar, percentage, current iteration, and, optionally, estimation 
#' of remaining time until completion of loop to console.
#'
#' @param i integer. The loop counter. Expected to start at 1.
#' @param n integer. The maximum loop count, i.e. length of the data.
#' @param names A character vector of length n. The names of the dimension
#' being looped over, e.g. row.names. Defaults to iteration number.
#' @param time logical. Indication of whether to print an estimate of time to
#' complete remaining iterations. See Details section for implementation.
#' @param extra_poly integer. Degree of the polynomial used for estimation of 
#' remaining time. Default is 0, i.e. mean duration of loops times 
#' number of remaining loops is used for time estimation. If set to -1, the 
#' default is  up to loop 50 and quadratic estimation for remaining loops. 
#' @details The time estimation is implemented by calculating a time difference 
#' of the current and the last k iterations where k is i-1. Values are stored
#' in a cache for estimation of remaining time. The cache is cleared when the 
#' loop is finished, i.e. when \eqn{i == n}, when the the loop is not finished 
#' (i.e. \eqn{i == n} was not reached) but \eqn{i == 1}, or when the
#' function is called in a different environment. However, the cache is not
#' cleared if the last call did not complete the loop and the current call
#' does not start at \eqn{i == 1}. This may be the case if the loop is
#' continued, but a falsely set \eqn{i = x | x > 1} as start value will not
#' be recognized as an error.
#' @export
#' import crayon

show_progress <- function() {
  #> This function leverages lexical closures. The function returned by this 
  #> higher order function will always have access to the environment it was 
  #> created in. This allows storing objects across calls without manipulating
  #> the global or any other environment. The challenge, however, was to 
  #> determine when objects in this environment should be reset (see @details). 
  
  # convenience function
  pretty_time <- function(x) {
    if (x > 90 * 60 * 24) return(paste(round(x / (60 * 60 * 24)), "days"))
    if (x > 90 * 60) return(paste(round(x / (60 * 60)), "hours"))
    if (x > 90) return(paste(round(x / 60), "minutes"))
    return(paste(round(x), "seconds"))
  }
  # data to enclose by function
  time_cache   <- c()
  time_stamp   <- 0
  loop_count   <- 0
  env_sentinel <- rlang::caller_env()
  # sentinel is set when i == n, i.e. loop completed;
  # if loop terminates in error, cache is cleared on next call
  completed <- FALSE
  # closure to return
  clsr <- function(i, n, names = NA_character_, time = FALSE, extra_poly = 0) {
    if (identical(names,  NA_character_)) names <- as.character(1:n)
    if (!is.numeric(i) || !is.numeric(n) || i <= 0 || n <= 0) {
      stop("Both arguments 'i' and 'n' must be positive integers.")
    }
    if (i > n) {
      stop("Index 'i' must not be greater than length of data 'n'.")
    }
    if (!i == as.integer(i)) {
      warning("Removing decimal places for indexing.")
      i <- as.integer(i)
    }
    if (!is.character(names)) {
      stop("Argument 'names' must be a character vector.")
    }
    # when restarting after error/abort or called in another environment
    if ((completed == FALSE & i == 1) | 
        !identical(env_sentinel, rlang::caller_env())) {
      # cleanup
      time_cache   <<- c()
      time_stamp   <<- 0
      loop_count   <<- 0
      env_sentinel <<- rlang::caller_env()
    }
    loop_count <<- loop_count + 1
    if (loop_count == 1 & i != 1) stop("Index 'i' must start at 1.")
    cat("\014\n")
    extra <- nchar("||100%")
    width <- options()$width
    step  <- round(i / n * (width - extra))
    text  <- sprintf(
      "|%s%s|% 3s%%",
      #crayon::bgMagenta(strrep(" ", step)),
      strrep("=", step),
      strrep(" ", width - step - extra),
      round(i / n * 100)
    )
    if (time == TRUE) {
      if (i == 1) {
        # message for first timediff calculation
        cat(crayon::magenta$bold("Calculating remaining time...\n"))
        time_stamp <<- Sys.time()
      } else {
        now <- Sys.time()
        loop_dur <- difftime(now, time_stamp, units = "secs")
        # smoothing over cached values
        time_cache <<- c(time_cache, loop_dur)
        # switch from linear prediction to polynomial when sufficient data are collected
        if (length(unique(time_cache)) <= extra_poly | 
            (extra_poly == -1 & i < 50) | 
            extra_poly == 0) {
          time_remaining <- mean(time_cache) * ((n - i) + 1)
        } else {
          # adjust degree after switching to quadratic estimation
          if (extra_poly == -1) extra_poly <- 2
          time_data <- data.frame(iterations = 1:length(time_cache),
                                  durations = time_cache)
          pred_data <- data.frame(iterations = (i + 1):n)
          pred_data$predictions <- 
            stats::predict(lm(durations ~ poly(iterations, 
                                               degree = extra_poly, 
                                               raw = TRUE), 
                              data = time_data), 
                           newdata = pred_data)
          time_remaining <- sum(pred_data$predictions)
        }
        cat(crayon::magenta$bold(
          "Approximately",
          pretty_time(time_remaining),
          "remaining.\n"
        ))
        # save for next iteration
        time_stamp <<- Sys.time()
      }
    }
    cat(text)
    cat("\nProcessing ", names[i], "\n")
    if (i == n) {
      # cleanup
      completed    <<- TRUE
      time_cache   <<- c()
      loop_count   <<- 0
      env_sentinel <<- rlang::caller_env()
      cat("\r\n", crayon::green$bold(">>> COMPLETED <<<"), "\n\n")
    } else {
      completed <<- FALSE
    }
  }
  return(clsr)
}


# assign return of higher order function, i.e. function showing progress, to
# desired name; roxygen documents the higher order function (first function
# in file)
show_progress <- show_progress()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot ggplot2 point shapes
point_shapes <- function(fill = "#77cc44", color = "#333333"){
  require("ggplot2")
  ggplot(data.frame(s = c(0:25, 32:126), p = c(0:14, 16:26, 32:126))) +
    scale_shape_identity() +
    geom_point(mapping = aes(x = p %% 16, y = p %/% 16, shape = s), 
               size = 5, fill = fill, color = color) +
    geom_text(mapping = aes(x = p %% 16, y = p %/% 16 + 0.28, label = s), size = 3) +
    geom_point(mapping = aes(x = 12.3, y = 1.1, shape = "\u2615"), size = 13) +
    geom_text(mapping = aes(x = 13.8, y = 1.1, label = "\\u2615"), size = 3) +
    theme_void()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ggplot2 theme 
theme_js <- function(base_size = 12, legend_size = 12, plot_margin = 20, ...) {
  require("ggplot2")
  theme_minimal() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      plot.title       = element_text(size = base_size + 4, hjust = 0),
      plot.subtitle    = element_text(size = base_size, hjust = 0.5),
      plot.caption     = element_text(size = base_size, color = "black"),
      plot.margin      = margin(plot_margin, plot_margin, plot_margin, plot_margin),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", linewidth = 1, fill = NA),
      legend.title     = element_text(size = legend_size, color = "black"),
      legend.text      = element_text(size = legend_size, color = "black"),
      legend.position  = "none",
      axis.line        = element_line(linewidth = 0.1, color = "black"),
      axis.ticks       = element_line(linewidth = 0.5),
      axis.text        = element_text(size = base_size, color = "black"),
      axis.title.y     = element_text(size = base_size, 
                                      margin = margin(t = 0, r = 8, b = 0, l = 0)),
      axis.title.x     = element_text(size = base_size, 
                                      margin = margin(t = 8, r = 0, b = 0, l = 0)),
      strip.text       = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "gray90", linewidth = 1)
    ) +
    theme(...)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
tp <- function(df) {
  tdf <- data.table::transpose(df) # does not always conserve names
  row.names(tdf) <- colnames(df)
  colnames(tdf) <- row.names(df)
  return(tdf)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set_rownames <- function(df, names) {rownames(df) <- names; return(df)}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
bioclim <- function() {
  require(dplyr)
  tibble::tribble(
    ~BIO_var, ~bio_var, ~bio_name,                                                      ~unit,
    "BIO01",  "bio1",   "Annual Mean Temperature",                                      "°C",
    "BIO02",  "bio2",   "Mean Diurnal Range\n(Mean of monthly (max temp - min temp))",  "°C",
    "BIO03",  "bio3",   "Isothermality (bio2/bio7) (×100)",                             NA,
    "BIO04",  "bio4",   "Temperature Seasonality (standard deviation ×100)",            "°C",
    "BIO05",  "bio5",   "Max Temperature of Warmest Month",                             "°C",
    "BIO06",  "bio6",   "Min Temperature of Coldest Month",                             "°C",
    "BIO07",  "bio7",   "Temperature Annual Range (bio5-bio6)",                         "°C",
    "BIO08",  "bio8",   "Mean Temperature of Wettest Quarter",                          "°C",
    "BIO09",  "bio9",   "Mean Temperature of Driest Quarter",                           "°C",
    "BIO10",  "bio10",  "Mean Temperature of Warmest Quarter",                          "°C",
    "BIO11",  "bio11",  "Mean Temperature of Coldest Quarter",                          "°C",
    "BIO12",  "bio12",  "Annual Precipitation",                                         "mm",
    "BIO13",  "bio13",  "Precipitation of Wettest Month",                               "mm",
    "BIO14",  "bio14",  "Precipitation of Driest Month",                                "mm",
    "BIO15",  "bio15",  "Precipitation Seasonality (Coefficient of Variation)",         "mm",
    "BIO16",  "bio16",  "Precipitation of Wettest Quarter",                             "mm",
    "BIO17",  "bio17",  "Precipitation of Driest Quarter",                              "mm",
    "BIO18",  "bio18",  "Precipitation of Warmest Quarter",                             "mm",
    "BIO19",  "bio19",  "Precipitation of Coldest Quarter",                             "mm"
  ) %>% 
    mutate(bio_var = factor(bio_var, levels = .$bio_var)) %>% 
    return()
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
clear <- function(pkgs = FALSE) {
  if (isFALSE(pkgs)) message("clearing environment")
  else message("clearing environment and unloading packages")
  if (isTRUE(pkgs) & !is.null(names(sessionInfo()$otherPkgs))) {
    suppressWarnings( 
      invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)),
                       detach,
                       character.only = TRUE,
                       unload = TRUE,
                       force = TRUE)
      )
    )
  }
  rm(list = ls())
  invisible(gc(verbose = FALSE))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
info <- function(object) {
  message(
    paste0(
      "storage mode:\t", storage.mode(object),
      "\nsize:\t\t", format(object.size(object), unit = "auto"),
      "\nobject type:\t", sloop::otype(object),
      "\nclass:\t\t", paste(class(object), collapse = ", "),
      "\nnamespace:\t", paste(find(deparse(substitute(object))), collapse = ", "),
      "\nattributes:\t", paste(names(attributes(object)), collapse = ", "),
      "\ndimensions:\t", paste(dim(object), collapse = ", ")
    )
  )
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
plot_lm <- function(lm){
  par(mfrow = c(2, 2), mai = c(0.4, 0.8, 0.4, 0.2))
  plot(lm, ask = FALSE, which = 1)
  plot(lm, ask = FALSE, which = 2)
  plot(lm, ask = FALSE, which = 3)
  plot(lm, ask = FALSE, which = 5)
  par(mfrow=c(1, 1), mai = rep(1, 4))
}

# example

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
conyza_population_data <- function() {
  structure(
    list(Range = c("native", "native", "native", "native", 
                   "native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "native", "native", "native", "native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "native", "native", "native", "native", "native", "native", "native", 
                   "native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "native", "native", 
                   "native", "native", "native", "native", "native", "non-native", 
                   "non-native", "native", "native", "native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "native", 
                   "native", "native", "native", "native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "native", 
                   "native", "native", "native", "native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "native", "native", "native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "native", 
                   "native", "native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "native", "native", 
                   "native", "native", "native", "native", "native", "native", "native", 
                   "native", "native", "native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "native", 
                   "native", "native", "native", "native", "native", "native", "native", 
                   "native", "native", "native", "native", "native", "native", "native", 
                   "native", "native", "native", "native", "native", "native", "native", 
                   "native", "native", "non-native", "non-native", "non-native", 
                   "native", "native", "native", "native", "native", "native", "native", 
                   "native", "native", "native", "native", "non-native", "non-native", 
                   "non-native", "native", "native", "native", "native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "native", "native", 
                   "native", "non-native", "non-native", "non-native", "native", 
                   "native", "native", "non-native", "non-native", "non-native", 
                   "native", "native", "native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "native", "native", "native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "native", "native", 
                   "native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "non-native", 
                   "non-native", "non-native", "non-native", "non-native", "native", 
                   "native", "native", "native", "non-native", "non-native", "non-native"
    ), 
    ID_3 = c("ALA-MON", "ALA-TUS", "ALA-VER", "ALB-CZR", "ALB-EDM", 
             "ALP-DAV", "ALP-FIL", "ALP-STA", "AMU-CHU", "AMU-KAN", "AMU-NOV", 
             "AMU-SWO", "ARM-LOR", "ARM-VDZ", "ARM-YER", "ARZ-LEM", "ARZ-NAU", 
             "ARZ-SED", "ARZ-SWE", "BEJ-CHP", "BEJ-FGH", "BEJ-TON", "BEL-BEL", 
             "BEL-KOR", "BEL-NOV", "BIC-ANA", "BIC-CLW", "BIC-CST", "BIC-PAV", 
             "CAL-BOD", "CAL-DAV", "CAL-LIV", "CAL-NAP", "CEU-AAS", "CEU-FOR", 
             "CEU-GIE", "CEU-PLA", "CEU-PRA", "CEU-THO", "CEU-ZED", "COL-AUR", 
             "COL-BOU", "COL-DEN", "COL-ELK", "COL-FAL", "COL-FOR", "COL-PUE", 
             "DEN-AAR", "DEN-MUN", "ECO-EAS", "ECO-MID", "ECO-PHI", "EGY-ALE", 
             "EGY-FAY", "EGY-ISM", "ENG-CAM", "ENG-LIN", "EST-PAR", "EST-RAK", 
             "EST-TAR", "FIN-JOU", "FIN-SAL", "FIN-TUR", "FLO-ANT", "FLO-BIV", 
             "FLO-CIT", "FLO-GAI", "FLO-OCA", "FRA-BAY", "FRA-GIR", "FRA-SAN", 
             "FRU-ANI", "FRU-BUK", "FRU-VLA", "GAN-LON", "GAN-PIN", "GAN-QIN", 
             "GRE-APY", "GRE-KVO", "GRE-PAL", "HOK-DAT", "HOK-FUR", "HOK-TOM", 
             "HUN-BAL", "HUN-HER", "HUN-KAZ", "HUN-MAR", "HUN-RAC", "HUN-ROH", 
             "IDA-BON", "IDA-POC", "ILL-CHI", "ILL-NAP", "IND-GAR", "IRK-BEL", 
             "IRK-CHE", "IRK-MAM", "ITA-LAT", "ITA-POG", "ITA-ROM", "JAP-KYO", 
             "JAP-NAG", "JAP-OKA", "JAP-TOK", "JAP-YAM", "JIL-ANT", "JIL-CHA", 
             "JIL-SON", "JOR-IRB", "JOR-MAL", "JOR-SAR", "JOR-SHA", "KAN-LAW", 
             "KAN-MAN", "KAN-SAL", "KAS-GAN", "KAS-HAZ", "KAS-SON", "KAS-SOP", 
             "KAZ-PAD", "KAZ-PAV", "KEN-CAR", "KEN-GEO", "KEN-STA", "KOS-FER", 
             "KOS-MIR", "KOS-PEJ", "KUR-KER", "KUR-MAR", "KUR-SAN", "KZH-ALM", 
             "KZH-AVA", "KZH-IAK", "LIT-PAB", "LIT-VAS", "MAN-WIN", "MAS-NEW", 
             "MEX-CAB", "MEX-PED", "MIC-ANN", "MIN-FAI", "MIN-MAN", "MIN-REN", 
             "MIN-STA", "MIS-CAH", "MIS-EUR", "MIS-SLO", "MOR-MAR", "MOR-TAN", 
             "MOR-TAZ", "MRU-GOD", "MRU-GOR", "MRU-MOS", "MRU-SUZ", "MRU-ZVE", 
             "NCC-MEI", "NCC-SHU", "NCC-YAN", "NCC-ZHO", "NDA-BRE", "NDA-MEL", 
             "NDA-TRU", "NEB-CAR", "NEB-CRE", "NEB-HAV", "NEB-LIN", "NEB-NOR", 
             "NEB-OMA", "NEV-REN", "NEW-GEN", "NEW-WNY", "NOS-COR", "NOS-KEA", 
             "NOS-MAI", "NOS-MUS", "NWU-ALB", "NWU-BOZ", "NWU-BUT", "NWU-SMA", 
             "OHI-TRB", "OKL-ELR", "OKL-NOR", "OKL-SHA", "OMS-NOV", "OMS-OKT", 
             "OMS-TRO", "ONT-FIN", "ONT-RID", "ONT-SAR", "ORE-POM", "ORE-SAL", 
             "PEN-BRO", "PEN-CHA", "PEN-CLI", "PEN-LAT", "PEN-NKE", "PEN-PIT", 
             "POL-TOR", "POR-COI", "POR-PED", "QUE-EDM", "QUE-LAB", "QUE-LAP", 
             "QUE-SOT", "RAL-BER", "RAL-GEN", "RAL-GRE", "ROM-HAL", "ROM-HER", 
             "ROM-IAS", "ROM-MUR", "ROM-PIE", "ROM-TIF", "SAB-LOM", "SAB-OXN", 
             "SAB-SAB", "SAR-ARB", "SAR-BAD", "SAR-PLA", "SAS-MEL", "SAS-PAL", 
             "SAS-REG", "SAU-ALB", "SAU-BAL", "SAU-FAI", "SCA-ABV", "SCA-CLE", 
             "SCA-GRE", "SCC-AOQ", "SCC-CJC", "SCC-QUA", "SCC-XIS", "SCC-ZHU", 
             "SDA-BRO", "SIB-ACZ", "SIB-BAR", "SIB-NOV", "SPA-CAS", "SPA-ELE", 
             "SPA-FRO", "SPA-TUI", "SRU-EGI", "SRU-PYA", "SRU-VED", "SWE-ALM", 
             "SWE-HAS", "SWE-LAH", "TUN-AGB", "TUN-BIR", "TUN-SUD", "TUN-TUN", 
             "TUR-ARS", "TUR-BAL", "TUR-CAP", "TUR-GIR", "TUR-ORD", "TUR-SOG", 
             "TUR-TRA", "UFA-BEZ", "UFA-GUR", "UKR-HOL", "UKR-KAM", "UKR-KYI", 
             "UKR-LVI", "UKR-RUB", "UKR-STA", "UTA-MOR", "UTA-SPV", "UTA-TRM", 
             "UZB-ASU", "UZB-CHU", "UZB-ESK", "UZB-GAL", "UZB-SAM", "UZB-TAS", 
             "VIR-BLB", "VIR-BLS", "VIR-PAI", "VOL-ARD", "VOL-ELA", "VOL-KAM", 
             "VOL-KVA", "VOL-MAM", "VOL-MIP", "VOL-NIZ", "VOL-RYB", "VOL-SAR", 
             "VOL-SYZ", "VOL-TOG", "VOL-TSI", "VOL-YOS", "WAS-OLY", "WAS-PEN", 
             "WAS-STE", "WAS-VAN", "XIN-HUT", "XIN-SAE", "XIN-XIN"), 
    Lat = c(32.36572, 
            32.47231, 32.74125, 52.37603, 53.52295, 46.81089, 46.67447, 46.77651, 
            50.15381, 50.19381, 50.49731, 51.34515, 40.83083, 39.74111, 40.21056, 
            32.45093, 35.19203, 34.99198, 32.27927, 40.15444, 39.75333, 39.89389, 
            50.53, 50.84, 50.8, 49.01254, 51.6466, 50.22348, 49.30789, 38.30685, 
            38.539, 37.65673, 38.36847, 54.39191, 53.399, 50.60756, 51.325, 
            49.999, 48.96683, 49.69, 39.609, 40.034, 39.77641, 39.565, 38.9403, 
            40.563, 38.3335, 56.13197, 56.30106, 38.74536, 39.44136, 39.97358, 
            31.00021, 30.26767, 30.61928, 52.16882, 53.22781, 58.38883, 59.35556, 
            58.39334, 61.14109, 60.388, 60.44396, 29.26286, 29.62799, 29.41433, 
            29.61545, 29.20794, 43.4948, 45.03165, 46.57151, 43.17621, 42.60636, 
            43.19717, 34.97184, 35.45889, 34.84417, 40.59278, 40.18862, 40.78326, 
            42.489, 43.221, 42.581, 46.793, 48.25559, 48.25559, 45.879, 47.026, 
            48.01679, 43.55293, 42.86699, 41.8575, 41.77889, 41.62139, 52.81243, 
            53.00929, 52.28872, 41.45007, 42.24659, 41.81532, 35.0174, 36.222, 
            34.67192, 35.652, 35.736, 42.40917, 43.86667, 44.57972, 32.48397, 
            32.64902, 32.15913, 32.05578, 38.95565, 39.10215, 38.95798, 34.21635, 
            34.11861, 34.26994, 34.53793, 52.27969, 52.29057, 38.3248, 38.101, 
            37.5423, 42.27806, 42.84944, 42.70306, 34.30575, 35.54839, 35.28087, 
            43.22107, 43.39573, 43.25189, 55.03046, 54.35054, 49.8942409016656, 
            42.33716, 18.91148, 19.32115, 42.29348, 43.60883, 44.18739, 44.46002, 
            33.4446, 38.66, 38.526, 38.684, 35.62583, 35.7775, 34.26528, 
            56.687, 56.20712, 55.829, 56.43199, 55.771, 34.3069, 34.2412, 
            34.2356, 34.4505, 47.91848, 47.9043, 46.98232, 40.307, 42.43767, 
            40.86142, 40.893, 41.04694, 41.296, 39.55389, 42.86663, 42.33693, 
            44.4533, 44.69524, 44.65107, 44.77603, 47.01729, 45.41537, 45.5918, 
            47.31305, 41.20968, 35.52821, 35.20016, 35.35112, 54.83762, 55.00102, 
            54.87825, 42.64875, 42.44972, 42.9966, 45.45683, 44.92506, 41.80474, 
            40.64086, 41.03391, 40.32004, 40.58962, 40.45745, 52.99919, 40.19351, 
            39.9103, 47.37528, 46.91111, 47.37639, 46.03611, 46.86684, 46.18571, 
            45.19673, 47.59032, 44.86069, 47.18042, 46.2312, 45.71114, 45.8677, 
            34.695583, 34.163911, 34.403982, 39.77635, 40.96068, 40.82007, 
            52.9487, 53.22182, 50.418, 19.998, 19.8493, 17.246, 34.20592, 
            34.66839, 34.76849, 25.9351, 26.0462, 29.4452, 29.1155, 29.2986, 
            44.31361, 52.159, 53.341, 54.837, 37.39122, 37.49992, 37.20081, 
            42.07609, 42.87694, 44.06222, 42.69778, 56.55255, 56.16989, 56.44689, 
            35.76113, 36.42033, 35.79858, 36.80799, 37.78639, 37.86167, 38.621, 
            40.96911, 40.9777, 37.93528, 41.04662, 54.58328, 54.95461, 49.43553, 
            49.29308, 50.45171, 49.79258, 49.02048, 49.28128, 41.07572, 40.14583, 
            41.72004, 41.36809, 41.07043, 40.92088, 41.54085, 39.64854, 41.34104, 
            37.23213, 37.083, 37.58866, 54.7561, 55.79202, 55.21017, 52.48927, 
            55.1951, 55.09053, 55.71169, 55.45039, 54.15021, 53.17507, 53.47375, 
            55.89047, 56.48471, 47.05305, 46.53591, 45.61701, 45.70073, 44.178, 
            43.471, 43.937), 
    Lon = c(-86.1478, -85.6908, -86.5274, -110.781, 
            -113.508, 9.843248, 9.679603, 10.25234, 128.2016, 127.6304, 127.5298, 
            127.95, 44.55944, 45.34944, 44.55806, -110.783, -111.654, -111.725, 
            -111.024, 116.251389, 116.010556, 116.702778, 36.65, 37.23, 37.87, 
            -119.354, -120.018, -119.24, -124.323, -123.058, -121.779, -121.727, 
            -122.34, 12.43947, 12.801, 8.653955, 12.323, 14.406, 8.504937, 
            10.941, -104.705, -105.29, -105.157, -104.588, -104.578, -105.037, 
            -104.626, 10.14901, 9.65887, -76.0703, -75.7431, -75.117, 29.72934, 
            32.36937, 32.27338, 0.110284, -0.54918, 24.50741, 26.34594, 26.71286, 
            28.48892, 23.09934, 22.21261, -82.057, -82.3532, -82.1482, -82.34, 
            -82.1511, -1.49136, 0.47762, 0.433536, 132.7981, 131.2206, 131.9207, 
            104.6886, 106.905, 105.6556, 22.97112, 23.80595, 22.253, 140.834, 
            142.4361, 141.427, 17.712, 20.64861, 20.64861, 18.249, 18.918, 
            22.14199, -112.033, -112.429, -87.6808, -88.1789, -87.2247, 103.5328, 
            103.0435, 104.1339, 12.8768, 12.6354, 12.45426, 135.5866, 138.484, 
            133.9577, 139.457, 138.714, 128.1169, 125.35, 123.5583, 35.98485, 
            35.72239, 35.99811, 35.92976, -95.2716, -96.6058, -97.5575, 74.77194, 
            74.835, 75.15406, 74.72583, 76.97253, 76.93459, -84.087, -84.5124, 
            -84.6557, 21.20722, 20.90833, 20.31722, 47.03805, 46.14753, 46.99572, 
            76.92122, 77.26108, 77.47357, 25.69882, 24.57465, -97.1233035571736, 
            -71.2093, -98.4493, -99.1945, -83.8121, -94.4814, -94.0077, -94.7933, 
            -88.7956, -90.076, -90.56, -90.36, -5.27444, -5.79833, -3.96639, 
            43.424, 42.67468, 37.701, 40.44289, 36.806, 107.7341, 108.4012, 
            108.0501, 107.5975, -97.11617, -98.11443, -96.88205, -97.677, 
            -97.8947, -96.5942, -96.571, -100.746, -96.196, -119.834, -76.98068, 
            -77.3331, -64.4182, -63.7057, -63.6521, -63.1776, -114.514, -111.031, 
            -112.321, -116.56, -80.5268, -98.1883, -97.4611, -97.0719, 73.37272, 
            73.36855, 73.33226, -81.3557, -81.8812, -82.3425, -118.965, -123.072, 
            -77.0838, -77.9535, -77.5164, -79.3812, -79.7576, -79.9092, 18.60659, 
            -8.41631, -8.95237, -68.3014, -71.3222, -69.9683, -73.1283, 7.364357, 
            6.153596, 5.81512, 23.03184, 22.40317, 27.54991, 24.97564, 27.12129, 
            27.09742, -120.365412, -119.060907, -119.698138, 8.576433, 8.8782, 
            8.516353, -105.01, -105.715, -104.616, 41.464, 41.5838, 43.088, 
            -82.4124, -82.8945, -82.3572, 112.5633, 112.432, 112.762, 117.9545, 
            112.6914, -96.7289, 79.711, 83.776, 83.104, -6.32642, -4.21792, 
            -6.91517, -8.62302, 44.93028, 43.05528, 45.57361, 14.16545, 13.81425, 
            12.92146, 10.82714, 10.56823, 10.66047, 10.19361, 27.3075, 27.78583, 
            34.805, 38.62107, 37.97725, 28.37556, 39.24485, 55.85843, 55.68667, 
            27.45858, 25.60344, 30.57651, 23.97841, 38.3738, 38.90665, -111.722, 
            -111.611, -112.141, 69.33908, 71.35236, 71.14165, 69.85774, 66.96551, 
            69.3175, -80.4358, -77.971, -75.8208, 46.27097, 52.20709, 49.24563, 
            48.11027, 45.92193, 45.97367, 48.95233, 50.2998, 45.16878, 48.54478, 
            49.56385, 47.43543, 48.04161, -122.912, -117.777, -122.043, -122.731, 
            86.963, 87.278, 87.432)), 
    row.names = 1:298, 
    class = "data.frame")
}
