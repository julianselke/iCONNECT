
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Functions used throughout the analysis scripts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




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
#' remaining time. Default is 0, i.e. mean duration of loops \eqn{\times} 
#' number of remaining loops is used for time estimation. If set to -1, the 
#' default is used up to loop 50 and quadratic estimation for remaining loops. 
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
    if (x > 90 * 60 * 24 * 30 * 12) return(paste(round(x / (60 * 60 * 24 * 30 * 12)), "years"))
    if (x > 90 * 60 * 24 * 30) return(paste(round(x / (60 * 60 * 24 * 30)), "months"))
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
      crayon::bgMagenta(strrep(" ", step)),
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
# desired name; NOTE: roxygen documents the higher order function (first function
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

theme_js <- function(base_size = 12, legend_size = 12, ...) {
  require("ggplot2")
  theme_minimal() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      plot.title       = element_text(size = base_size + 4, hjust = 0),
      plot.subtitle    = element_text(size = base_size, hjust = 0.5),
      plot.caption     = element_text(size = base_size, color = "black"),
      plot.margin      = margin(20, 20, 20, 20),
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

plot_lm <- function(lm){
  par(mfrow = c(2, 2), mai = c(0.4, 0.8, 0.4, 0.2))
  plot(lm, ask = FALSE, which = 1)
  plot(lm, ask = FALSE, which = 2)
  plot(lm, ask = FALSE, which = 3)
  plot(lm, ask = FALSE, which = 5)
  par(mfrow=c(1, 1), mai = rep(1, 4))
}

# example


physeqList <- function(...) {
  psl <- rlang::dots_list(.named = TRUE, ...)
  class(psl) <- c(class(psl), "psl")
  return(psl)
}

rare_fun <- function(ps, depth = NA_integer_) {
  ps %>% 
    rarefy_even_depth(sample.size = depth, 
                      rngseed = 2023, 
                      replace = FALSE, 
                      trimOTUs = TRUE, 
                      verbose = TRUE) %>% 
    prune_taxa(taxa_sums(.) > 0, .) %>% 
    prune_samples(sample_sums(.) > 0, .) 
}

ps_filter <- function(physeq, 
                      min.reads.across.samples = 0, 
                      min.confidence = 0, 
                      remove.unidentified = NULL, 
                      min.occurrences.across.samples = 0,
                      min.relative.abundance.per.samples = 0,
                      min.relative.abundance.across.samples = 0) {
  if (!"phyloseq" %in% class(physeq)) {
    stop(paste0("object '", deparse(substitute(physeq)), "' not of class 'phyloseq'"))
  }
  require("phyloseq")
  
  if (isFALSE(physeq@otu_table@taxa_are_rows)) physeq <- t(physeq)
  
  if (min.reads.across.samples > 0) {
    physeq <- prune_taxa(rowSums(data.frame(otu_table(physeq))) >= min.reads.across.samples, physeq)  
  }
  
  if (min.confidence > 0) {
    physeq <- prune_taxa(as.numeric(data.frame(tax_table(physeq))$confidence) >= confidence, physeq)
  }
  
  if (!is.null(remove.unidentified)) {
    physeq <- prune_taxa(data.frame(tax_table(physeq))[[remove.unidentified]] != "unidentified", 
                         physeq)
  }
  
  if (min.occurrences.across.samples > 0) {
    physeq <- prune_taxa(
      rowSums(
        ifelse(data.frame(otu_table(physeq, taxa_are_rows = TRUE)) > 0, 1, 0)
      ) >= min.occurrences.across.samples, 
      physeq
    )
  }
  
  if (min.relative.abundance.across.samples > 0) {
    total <- physeq %>% otu_table() %>% data.frame() %>% sum()
    physeq_names_rel <- physeq %>% 
      otu_table() %>% 
      data.frame() %>% 
      rowSums() %>% 
      data.frame(ab = ., ft_id = names(.)) %>%
      mutate(rel_ab = ab / total) %>% 
      filter(rel_ab > min.relative.abundance.across.samples) %>% 
      pull(ft_id)
    physeq_names <- physeq %>% tax_table() %>% data.frame() %>% row.names()
    physeq <- prune_taxa(x = physeq, physeq_names %in% physeq_names_rel)
  }
  
  if (min.relative.abundance.per.samples > 0) {
    physeq_names_rel <- physeq %>% 
      transform_sample_counts(function(x) x / sum(x) ) %>% 
      filter_taxa(function(x) sum(x) > min.relative.abundance.per.samples, prune = TRUE) %>% 
      tax_table() %>% 
      row.names()
    physeq_names <- physeq %>% tax_table() %>% row.names()
    physeq <- prune_taxa(x = physeq, physeq_names %in% physeq_names_rel)
  }
  
  return(physeq %>% prune_taxa(taxa_sums(.) > 0, .) %>% prune_samples(sample_sums(.) > 0, .) ) 
}


# get helpful information about an R object
info <- function(object) UseMethod("info")


info.default <- function(object) {
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


info.psl <- function(psl) {
  data.frame(
    "n.samples" = unlist(lapply(psl, nsamples)),
    "n.taxa" = unlist(lapply(psl, ntaxa))
  ) %>% 
    print()
}


info.phyloseq <- function(object) {
  
  suppressPackageStartupMessages(require(phyloseq))
  suppressPackageStartupMessages(require(tidyverse))
  suppressPackageStartupMessages(require(rlang))
  suppressPackageStartupMessages(require(crayon))
  
  ft <- phyloseq::otu_table(object, taxa_are_rows = TRUE)
  ft_bool <- ifelse(ft > 0, 1, 0)
  
  sparcity = sprintf("%1.2f%%", 100 * (1 - sum(ft_bool)/(ncol(ft_bool) * nrow(ft_bool))))
  
  tab <- function(tax, p.object = object) {
    p.tax <- rlang::as_name(enquo(tax))
    taxa_abund(physeq = p.object, tax = p.tax, add_meta = FALSE) %>% 
      mutate(n.taxa = nrow(.)) %>% 
      pivot_longer(-c({{tax}}, n.taxa)) %>%
      #pivot_longer(-c(class, n.taxa)) %>%
      
      group_by(name) %>% 
      mutate(sam_size = sum(value)) %>% 
      ungroup() %>% 
      group_by({{tax}}) %>% 
      #group_by(class) %>% 
      mutate(prevalence = (value/sam_size) * 100) %>%
      mutate(prevalence_sd = paste0("%  (± ", round(sd(prevalence), 2), "%)")) %>%
      mutate(prevalence = paste0(round(mean(value/sam_size * 100), 2), prevalence_sd)) %>%
      select(-prevalence_sd) %>% 
      ungroup() %>% 
      
      group_by({{tax}}) %>% 
      #group_by(class) %>% 
      mutate(detection = sum(value > 0)) %>% 
      reframe(#taxon = unique({{tax}}),
        count = sum(value), 
              n.taxa = unique(n.taxa),
              detection = unique(detection),
              prevalence = unique(prevalence)) %>% 
      mutate(count = count/sum(count)) %>% 
      dplyr::arrange(desc(count)) %>% 
      dplyr::slice(1) %>% 
      mutate(count = sprintf("%1.2f%%", 100 * count))
  }
  
  taxon_info <- data.frame(Rank = rank_names(object)[1:7], 
                           Taxon = NA, 
                           Rel.Abund. = NA,
                           n.Taxa = NA,
                           Detection = NA,
                           Prevalence = NA)
  
  for (i in 1:length(taxon_info$Rank)) taxon_info[i, 2:6] <- tab(!!sym(taxon_info$Rank[i]))
  
  taxon_info <- taxon_info %>% 
    .[, c(1,2,4,3,5,6)] %>% 
    mutate(Taxon = str_replace(Taxon, "_", ".")) %>% 
    {capture.output(print(.))} %>% 
    {gsub("^\\d", " ", .)} 
  
  orange <- make_style("#FF9900")
  rose   <- make_style("deeppink")
  leaf   <- make_style("#44cc44")
  sky    <- make_style("#3366ff")
  navy   <- make_style("#001122", bg = TRUE)
  
  f <- function(x, l = 12, s = "left") str_pad(format(x, unit = "auto", big.mark = ","), l, s, " ")
  
  f2 <- function(x, perc = FALSE, ll = 12) {
    str_split(str_replace_all(str_trim(x), "\\s+", " "), " ") %>%
      unlist() %>% 
      {ifelse(!is.na(suppressWarnings(as.numeric(.))), as.numeric(.), .)} %>%  
      {if (isTRUE(perc)) scales::label_percent()(./100) else .} %>%  
      sapply(f, l = ll, s = "left") %>% 
      paste0(collapse = "")
  }
  
  f3 <- function(x, perc = FALSE, ll = 12) {
    tmp <- str_split(str_replace_all(str_trim(x), "\\s+", " "), " ") %>% unlist()
    for (i in 1:length(ll)) tmp[i] <- f(tmp[i], l = ll[i], s = "left")

      paste0(tmp, collapse = " ")
  }
  
  filler <- function(n = 1, m = 0) paste0(rep("\t", n), collapse = "")
  
  info_message <- paste0(
    
    "\n", 
    
    # line 1
    bold(f2(capture.output(print(quantile(colSums(ft))))[1])),
    filler(2),
    paste0(rose$bold("object name:\t"), deparse(substitute(object))),
    
    # line 2
    leaf$bold("\nreads per sample:\t"), 
    filler(6),
    #navy(str_pad("", 26)),
    "\n", 
    
    # line 3
    navy(f2(capture.output(print(quantile(colSums(ft))))[2])),
    filler(2), 
    navy(paste0(orange$bold("size:\t\t"), f(object.size(object)))),
    
    # line 4
    leaf$bold("\nreads per feature:"), 
    filler(7),
    navy(str_pad("", 28)),
    "\n",
    
    #line 5
    navy(f2(capture.output(print(quantile(rowSums(ft))))[2])),
    filler(2),
    navy(paste0(orange$bold("samples:\t"), f(phyloseq::nsamples(object)))),
    
    # line 6
    leaf$bold("\nfeatures per samples:\t"), 
    filler(6),
    navy(paste0(orange$bold("features:\t"), f(phyloseq::ntaxa(object)))),
    "\n",
    
    # line 7
    navy(f2(capture.output(print(quantile(colSums(ft_bool))))[2])),
    filler(2),
    navy(str_pad("", 28)),
    
    # line 8
    leaf$bold("\nsamples per feature:\t"),
    filler(6),
    navy(paste0(orange$bold("sparcity:\t      "), sparcity)),
    "\n",
    
    #line 9
    navy(f2(capture.output(print(quantile(rowSums(ft_bool))))[2])),
    filler(2),
    navy(str_pad("", 28)),
    
    # line 10
    leaf$bold("\nassignment confience:\t"),
    filler(6),
    navy(paste0(orange$bold("taxa.are.rows:\t"), f(object@otu_table@taxa_are_rows))),
    "\n", 
    navy(f2(capture.output(print(quantile(as.numeric(
      data.frame(tax_table(object))$confidence) * 100)))[2], 
            TRUE)),
    filler(2),
    navy(paste0(orange$bold("has.phytree:\t"), f(!is.null(object@phy_tree)))),
    "\n", 
    sky$bold("\ntop taxa:"),
    "\n",
    sky(navy(f3(taxon_info[1], f3, ll = c(10, 31, 10, 12, 10, 22)))),
    "\n",
    sapply(taxon_info[2], f3, ll = c(10, 31, 10, 12, 10, 11)), "\n",
    navy(sapply(taxon_info[3], f3, ll = c(10, 31, 10, 12, 10, 11))), "\n",
    sapply(taxon_info[4], f3, ll = c(10, 31, 10, 12, 10, 11)), "\n",
    navy(sapply(taxon_info[5], f3, ll = c(10, 31, 10, 12, 10, 11))), "\n",
    sapply(taxon_info[6], f3, ll = c(10, 31, 10, 12, 10, 11)), "\n",
    navy(sapply(taxon_info[7], f3, ll = c(10, 31, 10, 12, 10, 11))), "\n",
    sapply(taxon_info[8], f3, ll = c(10, 31, 10, 12, 10, 11)), "\n",
    "\n"
    )

    cat(info_message)
  
}

diff <- function(...)  UseMethod("diff")

diff.phyloseq <- function(object1, object2) {
  
  suppressPackageStartupMessages(require(phyloseq))
  suppressPackageStartupMessages(require(tidyverse))
  suppressPackageStartupMessages(require(rlang))
  suppressPackageStartupMessages(require(crayon))
  
  
  f <- function(x, y = NULL, l = 12, is.size = FALSE, simple = TRUE) {
    
    f.to_num <- function(x) str_extract(x, pattern = "\\d+\\.?\\d*") %>% as.numeric()
    
    f.format <- function(x, keep = FALSE) {
      if (keep) return(format(x, unit = "auto", big.mark = ","))
      return(format(round(f.to_num(x), 2), unit = "auto", big.mark = ","))
    }
    
    f.sub <- function(x, is_size = FALSE, keep = FALSE) {
      if(isTRUE(is_size)) return(format(x, unit = "auto", big.mark = ","))
      is_numeric_or_percent <- str_detect(x, "[0-9][^a-zA-Z]?%?") && str_detect(x, "[0-9]+")
      is_percent <- str_detect(x, ".*%$")
      x <- x %>% 
        str_replace_all("%", "") %>% 
        str_trim() %>% 
        str_replace_all("\\s+", " ") %>% 
        str_split(" ") %>%
        unlist()
      if (is_percent && keep) return(paste0(sapply(x, f.format, is_size), "%"))
      if (is_numeric_or_percent && keep) return(sapply(x, f.format))
      if (is_numeric_or_percent) return(as.numeric(sapply(x, f.to_num)))
      return(x)
    }
    
    f.len_out <- function(x, l.out = l) str_pad(x, width = l.out, side = "left")
    
    f.compare <- function(x, y, l.out = l, is_size = is.size) {
      is_percent <- str_detect(x, ".*%$")
      x_num <- f.sub(x)
      y_num <- f.sub(y)
      if (is.numeric(x_num) && is.numeric(y_num)){
        if (x_num - y_num > 0) {
          if (is_percent) return(magenta(f.len_out(paste0(sapply(y, f.format, is_size), "%"), l.out)))
          return(magenta(f.len_out(f.format(y, is_size), l.out)))
        } else if (x_num - y_num < 0) {
          if (is_percent) return(green(f.len_out(paste0(sapply(y, f.format, is_size), "%"), l.out)))
          return(green(f.len_out(f.format(y, is_size), l.out)))
        } else {
          if (is_percent) return(silver(f.len_out(paste0(sapply(y, f.format, is_size), "%"), l.out)))
          return(silver(f.len_out(f.format(y, is_size), l.out)))
        }
      } else {
        if (y == "___________") return(black(y))
        if (x == y) {
          return(silver(f.len_out(y, l.out)))
        } else {
          return(magenta(f.len_out(y, l.out)))
        }
      }
    }
    
    # format only
    if (is.null(y)) {
      if (isTRUE(is.size)) {
        return(format(x, unit = "auto", big.mark = ",") %>% f.len_out())
      } else {
        sapply(x, f.sub, is.size, keep = TRUE) %>% 
          f.len_out() %>% 
          paste0(collapse = "") %>% 
          return()
      }
    } else {
      if (isTRUE(is.size)) {
        # pipe breaks format
        return(f.compare(x, y, l.out = l, is_size = is.size))
      }
      if (simple == 0) {
        tmp <- map2_chr(x, y, f.compare, l = 1, is_size = is.size)
        ll <- rep(l, length(tmp))
        for (i in 1:length(tmp)) tmp[i] <- f.len_out(x = tmp[i], l.out = ll[i])
        return(paste0(tmp, collapse = ""))
      }
      map2_chr(x, y, f.compare, l = 12, is_size = is.size) %>% f.len_out(l = 12) %>% paste0(collapse = "") %>% return()
    }
    
  }
  
  ft1 <- phyloseq::otu_table(object1, taxa_are_rows = TRUE)
  ft2 <- phyloseq::otu_table(object2, taxa_are_rows = TRUE)
  ft_bool1 <- ifelse(ft1 > 0, 1, 0)
  ft_bool2 <- ifelse(ft2 > 0, 1, 0)
  
  sparcity1 = sprintf("%1.2f%%", 100 * (1 - sum(ft_bool1)/(ncol(ft_bool1) * nrow(ft_bool1))))
  sparcity2 = sprintf("%1.2f%%", 100 * (1 - sum(ft_bool2)/(ncol(ft_bool2) * nrow(ft_bool2))))
  
  sparcity_diff <- sprintf("%1.2f%%", 
                           as.numeric(gsub("%", "", sparcity2)) - as.numeric(gsub("%", "", sparcity1))) %>% 
    {ifelse(str_detect(., "-"), magenta(f(sparcity2)), green(f(sparcity2)))}
    
  
  obj1.nsamples <- phyloseq::nsamples(object1)
  obj2.nsamples <- phyloseq::nsamples(object2)
  obj1.ntaxa <- phyloseq::ntaxa(object1)
  obj2.ntaxa <- phyloseq::ntaxa(object2)
  
  tab <- function(tax, p.object = object) {
    p.tax <- rlang::as_name(enquo(tax))
    taxa_abund(physeq = p.object, tax = p.tax, add_meta = FALSE) %>% 
      mutate(n.taxa = nrow(.)) %>% 
      pivot_longer(-c({{tax}}, n.taxa)) %>%
      group_by(name) %>% 
      mutate(sam_size = sum(value)) %>% 
      ungroup() %>% 
      group_by({{tax}}) %>% 
      mutate(prevalence = (value/sam_size) * 100) %>%
      mutate(prevalence_sd = paste0("%__(±", round(sd(prevalence), 2), "%)")) %>%
      mutate(prevalence = paste0(round(mean(value/sam_size * 100), 2), prevalence_sd)) %>%
      select(-prevalence_sd) %>% 
      ungroup() %>% 
      group_by({{tax}}) %>% 
      mutate(detection = sum(value > 0)) %>% 
      reframe(count = sum(value), 
              n.taxa = unique(n.taxa),
              detection = unique(detection),
              prevalence = unique(prevalence)) %>% 
      mutate(count = count/sum(count)) %>% 
      dplyr::arrange(desc(count)) %>% 
      dplyr::slice(1) %>% 
      mutate(count = sprintf("%1.2f%%", 100 * count))
  }
  
  taxon_info1 <- data.frame(Rank = rank_names(object1)[1:7], 
                           Taxon = NA, 
                           Rel.Abund. = NA,
                           n.Taxa = NA,
                           Detection = NA,
                           Prevalence = NA)
  
  for (i in 1:length(taxon_info1$Rank)) taxon_info1[i, 2:6] <- tab(!!sym(taxon_info1$Rank[i]), object1)
  
  taxon_info1 <- taxon_info1 %>% 
    mutate(Taxon = str_replace(Taxon, "_", ".")) %>% 
    separate(Prevalence, c("Prevalence", "SD"), "__")
  
  taxSD1 <- taxon_info1$SD
  taxon_info1 <- taxon_info1 %>% select(-SD)
  
  taxon_info2 <- data.frame(Rank = rank_names(object2)[1:7], 
                           Taxon = NA, 
                           Rel.Abund. = NA,
                           n.Taxa = NA,
                           Detection = NA,
                           Prevalence = NA)
  
  for (i in 1:length(taxon_info2$Rank)) taxon_info2[i, 2:6] <- tab(!!sym(taxon_info2$Rank[i]), object2)
  
  taxon_info2 <- taxon_info2 %>% 
    mutate(Rank = "___________") %>%  
    mutate(Taxon = str_replace(Taxon, "_", ".")) %>% 
    separate(Prevalence, c("Prevalence", "SD"), "__")

  taxSD2 <- taxon_info2$SD
  taxon_info2 <- taxon_info2 %>% select(-SD)
  
  
  gold <- make_style("#FF9900")
  sky  <- make_style("#3366ff")
  navy <- make_style("#001122", bg = TRUE)
  
  filler <- function(n = 2, m = 0) paste0(rep("\t", n), collapse = "")
  
  spacer <- c(10, 31, 15, 10, 13, 11)
  spacer2 <- c(10, 31+7, 15+8, 10+8, 13+8, 10+9)
  
  info_message <- paste0(
    
    "\n",
    paste0("\t", 
           gold$bold(deparse(substitute(object1))), 
           silver$bold(" vs. "), 
           gold$bold(deparse(substitute(object2)))),
    "\n",
    "\n",
    "\n",
    
    # line 1
    bold(f(c("0%", "25%", "50%", "75%", "100%"))),

    # line 2
    sky$bold("\nsequences per sample:"), 
    filler(7), 
    navy(paste0(sky$bold("size:\t\t"), f(object.size(object1), is.size = TRUE))),
    "\n", 
    
    # line 3
    navy(f(capture.output(print(quantile(colSums(ft1))))[2])),
    filler(),
    navy(paste0("\t\t", f(object.size(object1), object.size(object2), is.size = TRUE))),
    # line 4
    "\n",
    navy(f(quantile(colSums(ft1)), quantile(colSums(ft2)))),
    filler(), 
    navy(str_pad("", 28)),
    
    # line 5
    sky$bold("\nsequences per feature:"), 
    filler(7),
    navy(paste0(sky$bold("samples:\t"), f(nsamples(object1)))),
    "\n",
    
    # line 6
    navy(f(capture.output(print(quantile(rowSums(ft1))))[2])),
    filler(),
    navy(paste0("\t\t", f(nsamples(object1), nsamples(object2)))),
    "\n",
    # line 7
    navy(f(quantile(rowSums(ft1)), quantile(rowSums(ft2)))),
    filler(),
    navy(paste0(sky$bold("features:\t"), f(ntaxa(object1)))),
     
    # line 8
    sky$bold("\nfeatures per samples:\t"),
    filler(6),
    navy(paste0("\t\t", f(ntaxa(object1), ntaxa(object2)))),
    "\n",

    # line 9
    navy(f(capture.output(print(quantile(colSums(ft_bool1))))[2])),
    filler(),
    navy(str_pad("", 28)),
    "\n",
    # line 10
    navy(paste0(f(quantile(colSums(ft_bool1)), quantile(colSums(ft_bool2))))),
    filler(),
    navy(paste0(sky$bold("sparcity:\t      "), sparcity1)),

    # line 11
    sky$bold("\nsamples per feature:\t"),
    filler(6),
    navy(paste0("\t\t", sparcity_diff)),
    "\n",
    #line 12
    navy(f(capture.output(print(quantile(rowSums(ft_bool1))))[2])),
    filler(),
    navy(str_pad("", 28)),
    # line 13
    "\n",
    navy(f(quantile(rowSums(ft_bool1)), quantile(rowSums(ft_bool2)))),
    filler(),
    navy(paste0(sky$bold("taxa.are.rows:\t"), f(object1@otu_table@taxa_are_rows))),

    # line 14
    sky$bold("\nassignment confience:\t"),
    filler(6),
    navy(paste0("\t\t", f(object1@otu_table@taxa_are_rows, object2@otu_table@taxa_are_rows))),
    "\n",
    # line 15
    navy(f(paste0(quantile(as.numeric(data.frame(tax_table(object1))$confidence) * 100), "%"))),
    filler(),
    navy(paste0(sky$bold("has.phytree:\t"), f(!is.null(object1@phy_tree)))),
    "\n",
    # line 16
    navy(f(
      paste0(quantile(as.numeric(data.frame(tax_table(object1))$confidence) * 100), "%"),
      paste0(quantile(as.numeric(data.frame(tax_table(object2))$confidence) * 100), "%"))),
    filler(),
    navy(paste0(sky$bold("\t\t"), f(!is.null(object1@phy_tree), !is.null(object2@phy_tree)))
    ),
    "\n", 
    sky$bold("\ntop taxa:"),
    "\n",
    sky(navy(f(colnames(taxon_info1), l = c(10, 31, 15, 10, 13, 21)))),
    "\n",
    paste(f(taxon_info1[1,], l = spacer), silver(taxSD1[1])), "\n",
    paste(f(taxon_info1[1,], taxon_info2[1,], l = spacer2, simple = FALSE), silver(taxSD2[1])), 
    "\n",
    navy(paste(f(taxon_info1[2,], l = spacer)), silver(taxSD1[2])), "\n",
    navy(paste(f(taxon_info1[2,], taxon_info2[2,], l = spacer2, simple = FALSE)), silver(taxSD2[2])), 
    "\n",
    paste(f(taxon_info1[3,], l = spacer), silver(taxSD1[3])), "\n",
    paste(f(taxon_info1[3,], taxon_info2[3,], l = spacer2, simple = FALSE), silver(taxSD2[3])), 
    "\n",
    navy(paste(f(taxon_info1[4,], l = spacer)), silver(taxSD1[4])), "\n",
    navy(paste(f(taxon_info1[4,], taxon_info2[4,], l = spacer2, simple = FALSE)), silver(taxSD2[4])), 
    "\n",
    paste(f(taxon_info1[5,], l = spacer), silver(taxSD1[5])), "\n",
    paste(f(taxon_info1[5,], taxon_info2[5,], l = spacer2, simple = FALSE), silver(taxSD2[5])), 
    "\n",
    navy(paste(f(taxon_info1[6,], l = spacer)), silver(taxSD1[6])), "\n",
    navy(paste(f(taxon_info1[6,], taxon_info2[6,], l = spacer2, simple = FALSE)), silver(taxSD2[6])), "\n",
    paste(f(taxon_info1[7,], l = spacer), silver(taxSD1[7])), "\n",
    paste(f(taxon_info1[7,], taxon_info2[7,], l = spacer2, simple = FALSE), silver(taxSD2[7])), "\n"
  )
  
  cat(info_message)
  
}



geom_boxvio <- function(mapping = NULL, 
                        box.notch = TRUE,
                        out.dot.alpha = 0.65,
                        out.dot.size = 0.9,
                        line.width = 0.4) {
  call_env <- rlang::caller_env()
  if (!is.null(mapping)) {
    mapping <- ggplot2:::validate_mapping(mapping, call_env)
  }
  geoms <- list(
    geom_violin(mapping = mapping, 
                linewidth = line.width, 
                alpha = 0.3, 
                width = 0.9
    ),
    geom_boxplot(mapping = mapping, 
                 linewidth = line.width,
                 outlier.alpha = out.dot.alpha, 
                 outlier.size= out.dot.size,
                 width = 0.35, 
                 notch = box.notch
    )
  )
  return(geoms)
}

taxa_plot_facet <- function(physeq, 
                            tax, 
                            variable, 
                            facet_var = NULL, 
                            add_meta = TRUE,
                            exclude_nan = FALSE,
                            palette) {
  require(reshape2)
  data <- taxa_abund(physeq = physeq, tax = tax, add_meta = add_meta)
  d <- data[, c(
    1:(ncol(data) - ncol(sample_data(physeq))),
    which(colnames(data) == variable),
    which(colnames(data) == facet_var)
  )]
  d <- reshape2::melt(d, id.vars = c(variable, facet_var))
  f <- function(...) as.list(substitute(...()))
  arglist <- f(data = d, formula = paste(variable, "~ ."), fun.aggregate = sum)
  arglist$formula <- eval(arglist$formula)
  arglist$id.var <- c("Water_code", "Range", "variable")
  x <- do.call("recast", args = arglist)
  colnames(x) <- c(variable, "total")
  z <- merge(d, x)

  if (is.null(facet_var)) { stop("no facet provided")  }
  
  if (!is.null(facet_var)) {
    colnames(z) <- c("data", "facet", "variable", "value", "total")
    z <- z %>%
      group_by(facet, data) %>%
      dplyr::mutate(facet_total = sum(value)) %>% 
      ungroup()
    z <- z %>%
      group_by(facet) %>%
      dplyr::mutate(plot_width = length(unique(data))) %>% 
      ungroup()
    plot_widths <- z %>%
      group_by(facet) %>%
      dplyr::reframe(plot_width)
    plot_widths <- unique(as.data.frame(plot_widths))[, 2]
    plot_widths <- plot_widths / sum(plot_widths)
    z <- z %>%
      group_by(facet, data, variable) %>%
      dplyr::mutate(rel_ab = sum(value / facet_total * 100)) %>% 
      ungroup()
    if (isTRUE(exclude_nan)) z <- z %>% filter(!is.nan(rel_ab))
    p <- ggplot(z, aes(
      x = data, y = value / facet_total,
      fill = variable, color = variable,
      text = paste(
        data, "\n",
        "Rank:                    ", variable, "\n",
        "Rel. Abundance:   ", rel_ab, "%"
      ))) +
      geom_col() +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette) +
      theme_bw() +
      facet_grid(. ~ facet, scales = "free_x", space = "free_x")
  }
  p <- p + scale_y_continuous(limits = c(0,1), expand = c(0,0))
  return(p)
}

parse_taxa <- function(taxa) {
  tx_names <- taxa$Feature.ID
  # split taxonomy
  taxa <- tidyr::separate(
    data = taxa,
    col  = Taxon,
    into = c(
      "kingdom",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    ),
    sep = ";"
  )
  taxa <- taxa[,-1]
  taxa[is.na(taxa)] <- "unidentified"
  # replace empty fields
  taxa <- as.data.frame(lapply(taxa, 
                               function(x) gsub("[a-z]__$", "unidentified", x)))
  taxa <- as.data.frame(lapply(taxa, 
                               function(x) gsub("[a-z]__", "", x)))
  taxa <- as.data.frame(sapply(taxa, 
                               function(x) {gsub("^ ", "", x)}))
  rownames(taxa) <- tx_names
  return(taxa)
}

df_to_tax <- function(tax_table) {
  if(!inherits(tax_table, "data.frame")){
    stop(paste(substitute(tax_table), "is not of class 'data.frame'"))
  }
  if(ncol(tax_table) != 8){
    stop(paste("expected 8 columns;", 
               substitute(tax_table), 
               "has", 
               ncol(tax_table), 
               "columns"))
  }
  taxa <- as.matrix(tax_table)
  colnames(taxa) <- c(
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "confidence"
  )
  rownames(taxa) <- rownames(tax_table)
  return(tax_table(taxa))
}

taxa_abund <- function(physeq, tax, add_meta = TRUE) {
  if (!physeq@otu_table@taxa_are_rows) {
    x <- t(as(otu_table(physeq), "matrix"))
  } else {
    x <- as(otu_table(physeq), "matrix")
  }
  x <- aggregate(x, data.frame(tax_table(physeq)[, tax, drop = FALSE]), sum)
  if (add_meta == TRUE) {
    x <- t(x) %>%
      as.data.frame() %>%
      `colnames<-`(.[1, ]) %>%
      `<-`(., .[-1, ]) %>%
      lapply(., as.numeric)
    x <- cbind(x, sample_data(physeq))
  }
  return(x)
}


make_physeq <- function(meta_path, ftbl_path, taxa_path, conf_filter = 0.9, min_abund = 10) {
  md <- read.delim(meta_path) %>% 
    `rownames<-`(., .$sample_name) %>% 
    mutate(sample_name_alt = gsub("-", ".", .$sample_name)) %>% 
    mutate(Comp_code = case_when(Comp_code == "1:0" ~ "alone",
                                 Comp_code == "5:0" ~ "intra",
                                 Comp_code == "1:4" ~ "inter")) %>% 
    filter(Sample_type == "Field" | 
             (Sample_type == "GH" & Fine_root_dry_weight < 0.015))
  ft <- qiime2R::read_qza(ftbl_path)$data %>% as.data.frame()
  ft  <- ft[ , colnames(ft) %in% md$sample_name]  
  tx <- qiime2R::read_qza(taxa_path)$data %>% 
    data.frame() %>%
    parse_taxa() %>% 
    mutate(Confidence = as.numeric(.$Confidence)) %>% 
    mutate_at(vars(!matches("Confidence")), 
              function(l) ifelse(.$Confidence < conf_filter, "unidentified", l))
  ps <- phyloseq(otu_table(ft, taxa_are_rows = TRUE), 
                 df_to_tax(tx),
                 sample_data(md)) %>% 
    subset_samples(Sample_type == "GH" & Comp_code == "alone") %>% 
    prune_taxa(taxa_sums(.) > min_abund, .) %>% 
    prune_samples(sample_sums(.) > 0, .)
  return(ps)
}


myMarginal <- function (p, 
                        data, 
                        x, 
                        y, 
                        type = c("density", 
                                 "histogram", 
                                 "boxplot", 
                                 "violin", 
                                 "densigram"), 
                        margins = c("both", "x", "y"), 
                        size = 5, 
                        ..., 
                        xparams = list(), 
                        yparams = list(), 
                        groupColour = FALSE, 
                        groupFill = FALSE, 
                        intercepts = c(1000, 2000)) 
{
  type <- match.arg(type)
  margins <- match.arg(margins)
  prmL <- ggExtra:::toParamList(list(...), xparams, yparams)
  prmL <- ggExtra:::reconcileColParamApply(prmL)
  scatP <- ggExtra:::reconcileScatPlot(p, data, x, y) + 
    ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0.25, 0.25), "cm"))
  scatPbuilt <- ggplot2::ggplot_build(scatP)
  labels <- scatPbuilt$plot$labels
  hasTitle <- (!is.null(labels$title) || !is.null(labels$subtitle))
  if (hasTitle) {
    titleGrobs <- ggExtra:::getTitleGrobs(p)
    scatP$labels$title <- NULL
    scatP$labels$subtitle <- NULL
  }
  if (margins != "y") {
    plt <- ggExtra:::MarginalPlot$new("x", type, scatPbuilt, prmL, 
                                      groupColour, groupFill)
    # START MANIPULATION
    top <- plt$build() +
      geom_vline(xintercept = intercepts, linetype = "dotted", linewidth = 1)
    # END MANIPULATION
  }
  if (margins != "x") {
    plt <- ggExtra:::MarginalPlot$new("y", type, scatPbuilt, prmL, 
                                      groupColour, groupFill)
    right <- plt$build()
  }
  pGrob <- ggplot2::ggplotGrob(scatP)
  withCallingHandlers({
    suppressMessages({
      if (margins == "both") {
        ggxtraTmp <- ggExtra:::addTopMargPlot(pGrob, top, size)
        ggxtraNoTtl <- ggExtra:::addRightMargPlot(ggxtraTmp, right, 
                                                  size) 
        
      }
      else if (margins == "x") {
        ggxtraTmp <- 
          gtable::gtable_add_padding(pGrob, 
                                     grid::unit(c(0, 0.5, 0, 0), "lines"))
        ggxtraNoTtl <- ggExtra:::addTopMargPlot(ggxtraTmp, top, 
                                                size)
      }
      else if (margins == "y") {
        ggxtraTmp <- 
          gtable::gtable_add_padding(pGrob, 
                                     grid::unit(c(0.5, 0, 0, 0), "lines"))
        ggxtraNoTtl <- ggExtra:::addRightMargPlot(ggxtraTmp, right, 
                                                  size)
      }
    })
  }, warning = function(w) {
    if (grepl("did you forget aes", w, ignore.case = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
  if (hasTitle) {
    ggExtraPlot <- ggExtra:::addTitleGrobs(ggxtraNoTtl, titleGrobs)
  }
  else {
    ggExtraPlot <- ggxtraNoTtl
  }
  class(ggExtraPlot) <- c("ggExtraPlot", class(ggExtraPlot))
  ggExtraPlot
}

rarefaction_plot <- function(pseq, palette, title  = "", xlim = 18E3) {
  md <- sample_data(pseq)
  x <- otu_table(pseq, taxa_are_rows = TRUE) %>% data.frame() %>% nrow()
  palette <- unname(palette)
  this_pal <- rep(palette, ceiling(x/length(palette)))
  otu_table(pseq) %>% 
    data.frame(check.names = FALSE) %>% 
    t() %>% 
    data.frame(check.names = FALSE) %>% 
    rarecurve(step = 50, tidy = TRUE) %>%
    merge(., md, by.x = "Site", by.y = "sample_name", all.x = TRUE) %>%
    group_by(Site) %>%
    dplyr::mutate(max_y = max(Species)) %>%
    dplyr::mutate(max_x = max(Sample)) %>%
    ungroup() %>% 
    group_by(Sample_type, Sample) %>%
    dplyr::mutate(means = mean(Species)) %>% 
    {ggplot(.) + 
        geom_point(aes(x = max_x, y = max_y, color = Sample_type), 
                   size = 0.7, alpha = 0) +
        geom_line(aes(x = Sample, y = Species, group = Site, color = Site), 
                  alpha = 0.7, linewidth = 0.8) + 
        geom_vline(xintercept = c(2000), linetype = "dotted", linewidth = 1) +
        scale_x_continuous(breaks = seq(0, xlim, 2000),
                           limits = c(0, xlim)) +
                           #sec.axis = dup_axis(name = NULL, labels = NULL)) +
        #scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
        scale_color_manual(values = this_pal) +
        xlab("Rarefaction Depth") +
        ylab("Feature Count") +
        ggtitle(title) +
        theme_js()} %>% 
    #myMarginal(p, fill = "#33335555", margins = "x", intercepts = 2000)
    #ggMarginal(p, fill = "#3366ff55", margins = "x")
    #ggMarginal(p, fill = "#33333355", margins = "x")
    ggMarginal(p, fill = "#33335555", margins = "x", linewidth = 1)
}

my_dist <- function(physeq,
                    dis = "bray",
                    part.weights = "uniform",
                    ilr.weights = "uniform",
                    pseudocount = 1e-12,
                    normalized = TRUE,
                    parallel = FALSE,
                    ...) {
  # Check orientation of feature table in phyloseq object since there is
  # no internal control for assumption "samples are rows" in philr.
  # method dispatch:
  # philr::philr --> philr.phyloseq --> philr.data.frame
  if (isTRUE(physeq@otu_table@taxa_are_rows)) physeq <- phyloseq::t(physeq)
  ft <- otu_table(physeq)
  ft_bool <- ifelse(ft > 0, 1, 0)
  x <- rowSums(ft_bool) %>% data.frame()
  if (any(x == 0)) {
    warning("One or more samples contain zero features and will be removed.")
  }
  ft <- ft[x > 0, ]
  if (dis %in% c("aitchison", "robust.aitchison")) {
    dist <- vegan::vegdist(ft, method = dis, pseudocount = pseudocount, ...)
  } else if (!dis %in% c("u.unifrac", "w.unifrac", "philr")) {
    dist <- vegan::vegdist(ft, method = dis, ...)
  } else {
    if (dis != "philr") {
      dist <- phyloseq::UniFrac(physeq,
                                weighted = (dis == "w.unifrac"),
                                normalized = normalized,
                                parallel = parallel
      )
    } else {
      otu_table(physeq) <- otu_table(physeq) + pseudocount
      ph.ft <- philr::philr(otu_table(physeq),
                            phy_tree(physeq),
                            part.weights = part.weights,
                            ilr.weights = ilr.weights,
                            ...
      )
      dist <- vegan::vegdist(ph.ft, method = "euclidian")
    }
  }
  return(dist)
}


pcoa_dat_fun <- function(ps, Range, Comp){
  pcoa <- my_dist(ps %>% get_sample_subsets("Comp_code", {{Comp}}) %>% 
                    get_sample_subsets("Range", {{Range}}), "bray") %>% 
    ape::pcoa(.)
  xlab <- (paste0("PCo ", 1, " - ", 
                  round(pcoa$values$Relative_eig[[1]] * 100, 2),
                  "%"))
  ylab <- (paste0("\n\nPCo ", 2, " - ", 
                  round(pcoa$values$Relative_eig[[2]] * 100, 2),
                  "%"))
  data <- data.frame(
    X = pcoa$vectors[, 1],
    Y = pcoa$vectors[, 2],
    Z = pcoa$vectors[, 3]
  ) %>%
    mutate(Sample = rownames(.)) %>%
    `colnames<-`(c("X", "Y", "Z", "sampleid")) %>%
    merge(., sample_data(ps), by = "row.names")
  return(list(data, xlab, ylab))
}

#> convenience function for plotting
#> 
get_sample_subsets <- function(ps, col, val){
  # https://github.com/joey711/phyloseq/issues/487#issuecomment-441099026
  sam_dat <- sample_data(ps) %>% data.frame()
  sample_subset <- 
    sam_dat[which(sam_dat[[col]] == val),]
  phy_subset <- phyloseq(otu_table(ps), sample_data(sample_subset))
  return(phy_subset)
}


pcoa_plot_fun <- function(dat, title = "", subtitle = ""){
  ggplot(dat[[1]], aes(x = X, y = Y, color = Watering_treatment)) + 
    geom_point(size = 2) + 
    stat_ellipse(type = "norm") +
    theme_bw() + 
    scale_color_manual(values = c("#edae49", "#66a182")) +
    xlab(dat[[2]]) + 
    ylab(dat[[3]]) +
    labs(title = title,
         subtitle = subtitle) +
    jstheme()
}

do_procrustes <- function(physeq_1, physeq_2, water_code, range){
  
  gh_ids <- do.call("subset_samples", 
                    list(quote(physeq_2), 
                         substitute(Watering_treatment == water_code & 
                                      Range == range))) %>% 
    sample_data() %>%
    pull(sample_name)
  
  dat_1 <- do.call("subset_samples", 
                    list(quote(physeq_1), 
                         substitute(sample_name %in% gh_ids))) %>% 
    my_dist("bray") %>%
    ape::pcoa()
  
  pop_ids <- do.call("subset_samples", 
                     list(quote(physeq_1), 
                          substitute(sample_name %in% gh_ids))) %>% 
    sample_data() %>%
    pull(sample_name)
  
  dat_2 <- do.call("subset_samples", 
                    list(quote(physeq_2), 
                         substitute(Watering_treatment == water_code & 
                                      Range == range &
                                      sample_name %in% pop_ids))) %>% 
    my_dist("bray") %>% 
    ape::pcoa()
  
  res <- list(dat_1, dat_2, length(pop_ids), pop_ids)
  
  return(res)
}


#' Generate QC Summaries
#'
#' Summarize quality information of one or multiple fastq files.
#'
#' @md
#' @param file_paths A character vector of relative or absolute paths of
#' fastq files.
#' @param qs A numeric vector of quantiles of PHRED scores to compute.
#' @param meta_data A dataframe storing metadata.
#' @param id_col Sample identifier column in metadata.
#' @param extract_pattern A string specifying the sample identifier in the
#' file name. See details.
#' @param rm_pattern A string specifying the remainder to remove. See details.
#' @return A List. Elements are:
#' * "read.counts": A dataframe of read counts per file and respective
#' metadata.
#' * "read.means": A dataframe of mean PHRED scores per cycle.
#' * "read.quants": A (molten) dataframe of PHRED scores per cycle and
#' quantile.
#' * "quality.tiles": A (molten) dataframe of counts per PHRED score and cycle.
#' * "max.read.length": Integer.
#' @details Sample identifiers are assumed to be a substring of the respective
#' filenames. **extract_pattern** is used to match and extract this pattern
#' from file names. Since it might be necessary to match field
#' delimiters as well, **rm_pattern** can be specified to remove a pattern from
#' sample identifiers. Both patterns accept regular expressions.
#' @keywords MBT
#' @export
#' @examples readQC(c("NGSXX_041_lib123_1_.fastq",
#' "NGSXX_042_lib123_1_.fastq"), meta.dat, "SampleID", "_\\d\\d\\d_", "_")
#' @importFrom ShortRead qa
#' @importFrom reshape2 melt
#' @importFrom dplyr %>%
#' @importFrom stringr str_extract str_pad


readQC <- function(file_paths,
                   qs = c(0.09, 0.25, 0.5, 0.75, 0.91),
                   meta_data = NULL,
                   id_col,
                   extract_pattern,
                   rm_pattern = "",
                   show_progress = TRUE) {
  n <- length(file_paths)
  # init results
  res <- vector(mode = "list", length = n)
  # sequences are assumed to be Illumina reads of variable length but shorter
  # than 1000 bases; is trimmed to max.read.length later
  # max supported PHRED value is 40
  quality.tiles <- rep(list(rep(0, 40)), 1000)
  max.read.length <- 0
  # store only relevant info per file to reduce footprint
  for (i in 1:n) {
    tmp <- ShortRead::qa(file_paths[i], sample = FALSE)
    res[[i]][["read.count"]] <- c(tmp[["readCounts"]]$read, file_paths[i])
    pcq <- tmp[["perCycle"]]$quality
    res[[i]][["quality.means"]] <- data.frame(
      meanQ = rowsum(pcq$Score * pcq$Count, pcq$Cycle) / rowsum(pcq$Count, pcq$Cycle)
      )
    res[[i]][["quantiles"]] <-
      # apply per quantile and cycle
      lapply(qs, function(q) {
        by(pcq, pcq$Cycle, function(f) {
          .get_quant(f$Score, f$Count, q)
        }, simplify = TRUE) %>%
          unlist() %>%
          as.vector()
      })
    names(res[[i]][["quantiles"]]) <- as.character(qs)
    max.read.length <- max(
      max.read.length,
      max(tmp@.srlist$perCycle$quality$Cycle)
    )
    # add to counts in each iteration
    .phred <- tmp@.srlist$perCycle$quality[, c(1, 3, 4)] %>% data.frame()
    for (j in 1:nrow(.phred)) {
      .cycle <- .phred[j, 1]
      .score <- .phred[j, 2]
      .count <- .phred[j, 3]
      quality.tiles[[.cycle]][[.score]] <-
        quality.tiles[[.cycle]][[.score]] + .count
    }
    if (isTRUE(show_progress)) show_progress(i, n, file_paths, time = TRUE)
  }
  
  cat("\n")
  # aggregate qc results
  read.counts <- lapply(res, `[[`, 1) %>%
    data.frame() %>%
    t() %>%
    data.frame() %>%
    dplyr::mutate(X1 = as.numeric(X1)) %>%
    `names<-`(c("Read.Counts", "SampleID"))
  # add metadata
  if (!is.null(meta_data)){
    read.counts$Sample <- stringr::str_extract(read.counts$SampleID,
                                               extract_pattern) %>%
      gsub(rm_pattern, "", .)
    read.counts <- merge(read.counts, meta_data, by.x = "Sample", by.y = id_col)
  }
  read.means <- lapply(res, `[[`, 2) %>%
    lapply(., `[[`, 1) %>%
    lapply(., function(l) c(l, rep(NA, max.read.length-length(l)))) %>%
    data.frame() %>%
    t() %>%
    melt()
  read.quants <- lapply(res, `[[`, 3) %>%
    lapply(., as.data.frame, check.names = F) %>%
    mapply(function(x) "[<-"(x, "Cycle", value = 1:nrow(x)),
           .,
           SIMPLIFY = FALSE
    ) %>%
    mapply(melt, ., id.vars = "Cycle", SIMPLIFY = FALSE) %>%
    mapply(cbind, ., "SampleID" = seq_len(n), SIMPLIFY = FALSE) %>%
    bind_rows()
  quality.tiles <- quality.tiles[seq_len(max.read.length)]
  quality.tiles <- quality.tiles %>%
    as.data.frame() %>%
    t() %>%
    `rownames<-`(1:max.read.length) %>%
    melt() %>%
    `names<-`(c("Cycle", "Score", "Count"))
  quality.tiles$Count[quality.tiles$Count == 0] <- NA
  res <- list(
    read.counts,
    read.means,
    read.quants,
    quality.tiles,
    max.read.length
  )
  names(res) <- c(
    "read.counts",
    "read.means",
    "read.quants",
    "quality.tiles",
    "max.read.length"
  )
  return(res)
}

# helper function
.get_quant <- function(score, count, q) {
  score[which(cumsum(count) / sum(count) >= q)][[1]]
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#' Plot QC Summaries
#'
#' Various plotting options for results of readQC.
#'
#' @md
#' @param data A list returned by readQC().
#' @param option One of:
#' * "qp" | 1: Quantile plot.
#' * "qr" | 2: Quantile ribbons.
#' * "mq" | 3: Mean quality per cycle.
#' * "rc" | 4: Read count per sample.
#' * "bp" | 5: Quantile boxplot (similar to quality plot in qzv produced
#' by **qiime demux summarize**)
#' * "hm" | 6: Quality heatmap
#' @param meta_data_col Character specifying the the metadata column to color
#' read counts by.
#' @param fill_log A boolean. Scale fill of option "hm".
#' @param add_mean A boolean. Plot means on top of heatmap ("hm")
#' @param q Integer. One of 1:length(qs)
#' @param qs Numeric vector. Specifies quantiles calculated in readQC().
#' @param ... Passed to \link[viridis]{scale_fill_viridis_c}
#' @return A ggplot object.
#' @export
#' @keywords MBT
#' @examples plotQC(data)
#' @import dplyr
#' @import forcats
#' @import ggplot2
#' @import ggsci
#' @import viridis

plotQC <- function(data,
                   option,
                   meta_data_col,
                   fill_log = FALSE,
                   add_mean = TRUE,
                   q = 3,
                   qs = c(0.09, 0.25, 0.5, 0.75, 0.91),
                   alpha = 0.1,
                   ...) {
  
  if (option == "qr" | option == 1) p <- plot_qr(data = data, ...)
  if (option == "qp" | option == 2) p <- plot_qp(data = data,
                                                 q = q,
                                                 qs = qs, ...)
  if (option == "mq" | option == 3) p <- plot_mq(data = data,
                                                 alpha = alpha,
                                                 ...)
  if (option == "rc" | option == 4) p <- plot_rc(data = data,
                                                 meta_data_col = meta_data_col,
                                                 ...)
  if (option == "bp" | option == 5) p <- plot_bp(data = data, ...)
  if (option == "hm" | option == 6) p <- plot_hm(data = data,
                                                 fill_log = fill_log,
                                                 add_mean = add_mean,
                                                 ...)
  return(p)
}

plot_qr <- function(data = data, ...){
  .dat <- data[["read.quants"]] %>%
    .[order(.[, 1], .[, 2]), ] %>%
    group_by(Cycle, variable) %>%
    dplyr::mutate(value = mean(value))
  p <- ggplot(data = .dat[.dat$variable == "0.5",],
              aes(x = Cycle, y = value)) +
    geom_ribbon(aes(ymin = .dat$value[.dat$variable == "0.09"],
                    ymax = .dat$value[.dat$variable == "0.5"]),
                fill = "darkblue", alpha = 0.1) +
    geom_ribbon(aes(ymin = .dat$value[.dat$variable == "0.25"],
                    ymax = .dat$value[.dat$variable == "0.5"]),
                fill = "darkblue", alpha = 0.1) +
    geom_ribbon(aes(ymin = .dat$value[.dat$variable == "0.5"],
                    ymax = .dat$value[.dat$variable == "0.91"]),
                fill = "#4682b4", alpha = 0.3) +
    geom_ribbon(aes(ymin = .dat$value[.dat$variable == "0.5"],
                    ymax = .dat$value[.dat$variable == "0.75"]),
                fill = "#4682b4", alpha = 0.3) +
    geom_line(data = .dat[.dat$variable == "0.5",] %>%
                group_by(Cycle) %>%
                mutate(mean = mean(value)),
              aes(x = Cycle, y = mean, group = SampleID),
              alpha = 0.1) +
    scale_x_continuous(breaks = seq(20, max(.dat$Cycle), 20)) +
    scale_y_continuous(breaks = seq(0, 40, 5), limits = c(0, 40)) +
    coord_cartesian(expand = FALSE) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  return(p)
}

plot_qp <- function(data = data,
                    q = 3,
                    qs = c(0.09, 0.25, 0.5, 0.75, 0.91),
                    alpha = 0.1,
                    ...){
  .dat <- data[["read.quants"]] %>%
    .[.$variable == as.character(qs[q]), ] %>%
    .[order(.[, 1], .[, 2]), ]
  p <- ggplot(data = .dat) +
    theme_bw() +
    geom_line(aes(x = Cycle, y = value, group = SampleID),
              alpha = alpha) +
    scale_x_continuous(breaks = seq(20, max(.dat$Cycle), 20)) +
    scale_y_continuous(breaks = seq(0, 40, 5), limits = c(0, 40)) +
    coord_cartesian(expand = FALSE) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      legend.position = "none"
    )
  return(p)
}

plot_mq <- function(data = data, alpha = 0.1, ...){
  p <- ggplot(data[["read.means"]]) +
    geom_line(aes(group = Var1, x = Var2, y = value), alpha = alpha) +
    theme_bw() +
    scale_x_continuous(breaks = seq(20, max(data[[2]]$Var2), 20)) +
    scale_y_continuous(breaks = seq(0, 40, 5), limits = c(0, 40)) +
    coord_cartesian(expand = FALSE) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      legend.position = "none"
    )
  return(p)
}

plot_rc <- function(data = data, meta_data_col, ...){
  if (!meta_data_col %in% colnames(data[["read.counts"]])) {
    did_you_mean(colnames(data[["read.counts"]]), meta_data_col)
  }
  p <- data[["read.counts"]] %>%
    dplyr::mutate(name = fct_reorder(Sample, Read.Counts, .desc = TRUE)) %>%
    ggplot() +
    geom_col(aes(
      x = name, y = Read.Counts,
      fill = get(meta_data_col)
    )) +
    theme_bw() +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.1)),
      breaks = seq(0, 1000000, 10000)
    ) +
    theme(
      axis.text.x = element_text(angle = -90),
      legend.position = "bottom"
    ) +
    ylab("Read Count") +
    xlab("Sample") +
    {
      if (!is.numeric(data[["read.counts"]][[eval(meta_data_col)]])) {
        scale_fill_d3(name = eval(meta_data_col))
      } else {
        scale_fill_viridis_c(name = eval(meta_data_col), ...)
      }
    }
  return(p)
}

plot_bp <- function(data = data, ...){
  .dat <- data[["read.quants"]]
  p <- ggplot() +
    geom_boxplot(
      aes(
        x = .dat[.dat$variable == "0.09", ]$Cycle %>%
          unique() %>%
          sort(),
        group = .dat[.dat$variable == "0.09", ]$Cycle %>%
          unique() %>%
          sort(),
        ymin = .dat[.dat$variable == "0.09", ] %>%
          group_by(Cycle) %>%
          dplyr::summarize(Mean = mean(value, na.rm = TRUE)) %>%
          .$Mean,
        lower = .dat[.dat$variable == "0.25", ] %>%
          group_by(Cycle) %>%
          dplyr::summarize(Mean = mean(value, na.rm = TRUE)) %>%
          .$Mean,
        middle = .dat[.dat$variable == "0.5", ] %>%
          group_by(Cycle) %>%
          dplyr::summarize(Mean = mean(value, na.rm = TRUE)) %>%
          .$Mean,
        upper = .dat[.dat$variable == "0.75", ] %>%
          group_by(Cycle) %>%
          dplyr::summarize(Mean = mean(value, na.rm = TRUE)) %>%
          .$Mean,
        ymax = .dat[.dat$variable == "0.91", ] %>%
          group_by(Cycle) %>%
          dplyr::summarize(Mean = mean(value, na.rm = TRUE)) %>%
          .$Mean
      ),
      stat = "identity",
      fill = "#8aa5b8",
      linetype = 1,
      color = "gray20",
      lwd = 0.2
    ) +
    theme_bw() +
    xlab("Cycle") +
    scale_x_continuous(breaks = seq(20, 1000, 20)) +
    scale_y_continuous(
      breaks = seq(0, 40, 5),
      limits = c(0, 40),
      minor_breaks = seq(1, 40, 1)
    ) +
    coord_cartesian(expand = FALSE) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "none"
    )
  return(p)
}

plot_hm <- function(data = data,
                    fill_log = FALSE,
                    add_mean = TRUE,
                    alpha = 0.1,
                    ...){
  p <- ggplot() +
    {
      if (isTRUE(fill_log)) {
        geom_tile(
          data = data[["quality.tiles"]],
          aes(x = Cycle, y = Score, fill = log(Count))
        )
      }
    } +
    {
      if (isFALSE(fill_log)) {
        geom_tile(
          data = data[["quality.tiles"]],
          aes(x = Cycle, y = Score, fill = Count)
        )
      }
    } +
    {
      if (isTRUE(add_mean)) {
        geom_line(
          data = data[["read.means"]],
          aes(group = Var1, x = Var2, y = value),
          color = "red", alpha = alpha
        )
      }
    } +
    scale_fill_viridis_c(
      direction = -1,
      option = "mako",
      na.value = "white"
    ) +
    theme_bw() +
    scale_x_continuous(breaks = seq(20, 1000, 20)) +
    scale_y_continuous(breaks = seq(0, 40, 5), limits = c(0, 40)) +
    coord_cartesian(expand = FALSE) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  return(p)
}





sample_effort <- function(physeq, ..., seed = 42, n = 1) {
  
  set.seed(seed)
  
  .f <- function(physeq, ...) {
    groups = enquos(...)
    
    physeq %>% 
      otu_table() %>% 
      data.frame() %>% 
      {ifelse(. > 0, 1, 0)} %>% 
      data.frame() %>% 
      rownames_to_column("feature") %>% 
      pivot_longer(cols = 2:ncol(.), names_to = "sample", values_to = "ftcount") %>% 
      filter(ftcount > 0) %>% 
      left_join(physeq %>% 
                  sample_data() %>% 
                  data.frame() %>%
                  # column name!
                  select(!!!groups, sample_name_alt),
                by = join_by("sample" == "sample_name_alt")) %>% 
      group_by(sample, !!!groups) %>% 
      mutate(features = list(unique(feature))) %>% 
      select(-feature, -ftcount) %>% 
      unique() %>% 
      ungroup() %>% 
      group_by(!!!groups) %>% 
      mutate(nsamples = sample(n(), n())) %>% 
      dplyr::arrange(nsamples) %>%
      # inverse of reduce => combine feature lists of samples  
      mutate(nfeatures = accumulate(features, ~ c(.x, .y))) %>% 
      group_by(nsamples, !!!groups) %>% 
      mutate(nfeatures = length(unique(unlist(nfeatures)))) %>% 
      ungroup() %>% 
      return()
  }
  
  if (n == 1) return(.f(physeq, ...))
  
  res <- .f(physeq, ...) %>% mutate(iter = 1)
  
  for (i in seq_len(n-1)) {
    res <- rbind(res, .f(physeq, ...) %>% mutate(iter = i))
  }
  
  groups = enquos(...)
  
  res <- res %>% 
    group_by(!!!groups, nsamples) %>% 
    reframe(nft_mean = mean(nfeatures),
            nft_n = n(),
            nft_sd = sd(nfeatures),
            nft_se = nft_sd/sqrt(nft_n),
            alpha = 0.05,
            nft_df = nft_n - 1,
            nft_t = qt(p = alpha / 2, df = nft_df, lower.tail = F),
            nft_margin = nft_t * nft_se,
            nft_lower = nft_mean - nft_margin,
            nft_upper = nft_mean + nft_margin) %>% 
    select(-c(nft_n, nft_sd, nft_se, alpha, nft_df, nft_t, nft_margin))
  
  return(res)
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

e_rich <- function(x) {
  estimate_richness(x, 
                    measures = c("Observed", 
                                 "InvSimpson", 
                                 "Shannon",
                                 "Evenness")) %>% 
    dplyr::mutate(sample_name = gsub("\\.", "-", row.names(.)),
                  Evenness = Shannon / log(Observed)) %>% 
    reshape2::melt() %>% 
    merge(sample_data(x) %>% data.frame(), 
          by = "sample_name", 
          all = FALSE)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#> stepwise selection for linear mixed-effects models, i.e. outputting the 
#> model without any non-significant fixed effects  

lmer_select <- function(formula, dat){
  formula <- deparse(substitute(formula))
  lm0 <- lmerTest::lmer(formula = as.formula(paste(formula, 
                                                             collapse = " ")), 
                                  data = dat)
  repeat{
    lm_test <- lmerTest::step(lm0, reduce.random = FALSE) 
    if(sum(lm_test$fixed$Eliminated) == 0) break
    if(sum(lm_test$fixed$Eliminated != 0) == nrow(lm_test$fixed)){
      warning("all terms eliminated")
      return(NA)
    }
    new_formula <- lm_test %>% 
      attr("model") %>% 
      `@`(call) %>% 
      as.character() %>% 
      .[2] 
    lm0 <- lmerTest::lmer(formula = as.formula(new_formula), data = dat)
  }
  return(lm0)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# wrapper for differential abundance (DA) testing with ANCOMBC and DESeq2;
# groups by range affiliation and tests DA between watering treatment levels

diff_abund <- function(physeq, taxon){
  
  ps_tmp <- physeq %>% 
    subset_samples(Range == "native", drop = TRUE)
  ps_tmp <- prune_taxa(taxa_sums(ps_tmp) > 0, ps_tmp) 
  ps_tmp <- microbiome::aggregate_taxa(ps_tmp, taxon)
  
  out <- ancombc2(
    data = ps_tmp,
    tax_level = taxon,
    fix_formula = "Watering_treatment", 
    p_adj_method = "fdr", 
    prv_cut = 0.05, # prev filtering has been done above already
    lib_cut = 0, 
    struc_zero = FALSE, # TRUE causes error if only "native"
    neg_lb = TRUE,
    iter_control = list(tol = 1e-2, max_iter = 20, verbose = FALSE),
    em_control = list(tol = 1e-5, max_iter = 200), # use max_iter >= 100 on real data 
    alpha = 0.05, 
    global = FALSE, # multi group comparison will be deactivated automatically 
    n_cl = 4,
    verbose = TRUE
  )
  # store the results in res 
  ancombc_res_native <- out$res
  ancombc_res_native$l2fc <- 
    log2(exp(ancombc_res_native$lfc_Watering_treatmentmesic))
  ancombc_res_native$l2fc_se <- 
    log2(exp(ancombc_res_native$se_Watering_treatmentmesic))
  
  ds <-  phyloseq_to_deseq2(ps_tmp, ~ Watering_treatment)
  # https://support.bioconductor.org/p/63229/#63230
  # https://support.bioconductor.org/p/100201/#100202
  # alternative: add pseudo count of 1
  ds <-  DESeq(ds, test = "Wald", fitType = "parametric", sfType = 'poscounts')
  
  res <- results(ds, alpha = 0.05)
  res <- res[order(res$padj, na.last = NA), ]
  
  deseq2_res_native <- res@listData %>% as.data.frame()
  
  taxa_sig <- rownames(res) %>% 
    strsplit("_") %>% 
    lapply(tail, 4) %>% 
    lapply(paste, sep = "", collapse = "_") %>% 
    unlist()
  
  deseq2_res_native$species <- taxa_sig
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  ps_tmp <- physeq %>% 
    subset_samples(Range == "non-native", drop = TRUE)
  ps_tmp <- prune_taxa(taxa_sums(ps_tmp) > 0, ps_tmp) 
  ps_tmp <- microbiome::aggregate_taxa(ps_tmp, taxon)
  
  out <- ancombc2(
    data = ps_tmp,
    tax_level = taxon,
    fix_formula = "Watering_treatment", 
    p_adj_method = "fdr", 
    prv_cut = 0.05, # prev filtering has been done above already
    lib_cut = 0, 
    struc_zero = FALSE, # TRUE causes error if only "native"
    neg_lb = TRUE,
    iter_control = list(tol = 1e-2, max_iter = 20, verbose = FALSE),
    em_control = list(tol = 1e-5, max_iter = 200), # use max_iter >= 100 on real data 
    alpha = 0.05, 
    global = FALSE, # multi group comparison will be deactivated automatically 
    n_cl = 4,
    verbose = TRUE
  )
  # store the results in res 
  ancombc_res_non_native <- out$res
  ancombc_res_non_native$l2fc <- 
    log2(exp(ancombc_res_non_native$lfc_Watering_treatmentmesic))
  ancombc_res_non_native$l2fc_se <- 
    log2(exp(ancombc_res_non_native$se_Watering_treatmentmesic))
  
  
  ds <- phyloseq_to_deseq2(ps_tmp, ~ Watering_treatment)
  # https://support.bioconductor.org/p/63229/#63230
  # https://support.bioconductor.org/p/100201/#100202
  # alternative: add pseudo count of 1
  ds <- DESeq(ds, test = "Wald", fitType = "parametric", sfType = 'poscounts')
  
  res <- results(ds, alpha = 0.05)
  res <- res[order(res$padj, na.last = NA), ]
  
  deseq2_res_non_native <- res@listData %>% as.data.frame()
  
  taxa_sig <- rownames(res) %>% 
    strsplit("_") %>% 
    lapply(tail, 4) %>% 
    lapply(paste, sep = "", collapse = "_") %>% 
    unlist()
  
  deseq2_res_non_native$species <- taxa_sig
  
  
  daa_res <- rbind(
    ancombc_res_native %>% 
      filter(p_Watering_treatmentmesic < 0.05) %>% 
      dplyr::select(taxon, l2fc, l2fc_se) %>% 
      dplyr::mutate(method = "ANCOMBC", Range = "native"),
    ancombc_res_non_native %>% 
      filter(p_Watering_treatmentmesic < 0.05) %>% 
      dplyr::select(taxon, l2fc, l2fc_se) %>% 
      dplyr::mutate(method = "ANCOMBC", Range = "non-native"),
    deseq2_res_native %>% 
      filter(padj < 0.05) %>% 
      dplyr::select(species, log2FoldChange, lfcSE) %>% 
      dplyr::rename(taxon = species, l2fc = log2FoldChange, l2fc_se = lfcSE) %>% 
      dplyr::mutate(method = "DESeq2", Range = "native"),
    deseq2_res_non_native %>% 
      filter(padj < 0.05) %>% 
      dplyr::select(species, log2FoldChange, lfcSE) %>% 
      dplyr::rename(taxon = species, l2fc = log2FoldChange, l2fc_se = lfcSE) %>% 
      dplyr::mutate(method = "DESeq2", Range = "non-native")
  ) %>% 
    dplyr::mutate(direction = ifelse(.$l2fc > 0, "pos", "neg"))
  
  return(daa_res)
}


