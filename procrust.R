

# margin=1 if rows are samples; margin=2 otherwise
hellinger_transform <- function(x, margin=1) return(sqrt(x/apply(x,margin,sum)))

library(vegan)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(qs2) # save R objects like you do with saveRDS but very fast

# rownames pattern:
# in words: 3-digit population code -- watering code -- seed family -- range
# as regex: ([A-Z]{3}-[A-Z]{3})_(dry|wet)_(\\d)_([native|non-native])

# load matrices
matrix_fungi <- qs_read("matrix_fungi.qs2")
matrix_exudates <- qs2::qs_read("matrix_exudates.qs2")
matrix_root <- qs2::qs_read("matrix_root.qs2")
# remove all samples with missing data
matrix_root <- matrix_root[apply(matrix_root, MARGIN = 1, FUN = \(x) !any(is.na(x))), ]
# subset to samples present across all three matrices
samples_to_keep <- Reduce(base::intersect, lapply(list(matrix_fungi, matrix_exudates, matrix_root), rownames))
matrix_fungi <- matrix_fungi[samples_to_keep, ]
matrix_exudates <- matrix_exudates[samples_to_keep, ]
matrix_root <- matrix_root[samples_to_keep, ]

# assert same sample order
all(rownames(matrix_exudates) == rownames(matrix_fungi)) & all(rownames(matrix_exudates) == rownames(matrix_root))


id <- data.frame(id = rownames(matrix_root)) %>% tidyr::separate_wider_delim(id, "_", names = c("pop_code", "water", "sf", "range"), cols_remove = FALSE)

variable_data <- read.delim("data.csv", sep = ",") %>%
  dplyr::select(CWD, Admixture_total, sample_ID) %>%
  mutate(pop_code = str_extract(sample_ID, "[A-Z]{3}-[A-Z]{3}"),
         water = ifelse(str_detect(sample_ID, "_D_"), "dry", "wet")) %>%
  right_join(id, by = join_by(pop_code, water))


ScreePlot <- function(data) {
  ord_x <- data$ordination_x
  ord_y <- data$ordination_y
  if (!inherits(ord_x, "pcoa")) x <- summary(ord_x)$importance[3,] else x <- ord_x$values$Cum_corr_eig
  if (!inherits(ord_y, "pcoa")) y <- summary(ord_y)$importance[3,] else y <- ord_y$values$Cum_corr_eig
  scree_data <- data.frame(x = seq_len(length(x)), y = c(x, y), dataset = rep(unique(data$points$dataset), each = length(x)))
  ggplot(scree_data) +
    aes(x, y, color = dataset, group = dataset) +
    geom_line(linewidth = 0.4) +
    geom_point(shape = 21) +
    labs(title = "Screeplot", x = "Axis", y = "Cumul. Expl. Var.")
}




#' @param x_data matrix, e.g., matrix_fungi
#' @param x_dist "bray" or "eucl"
#' @param x_ord "pca" or "pcoa"
#' @param y_data matrix, e.g., matrix_fungi
#' @param y_dist "bray" or "eucl"
#' @param y_ord "pca" or "pcoa"
#' @param water_subset_pattern "both", "dry", or "wet"
#' @param range_subset_pattern "both", "native", or "non-native"
#'
#' @return
#' @export
#'
#' @examples
ProcrustesByGroup <- function(x_data,
                              x_dist,
                              x_ord,
                              y_data,
                              y_dist,
                              y_ord,
                              water = "both",
                              range = "both"
) {

  if (!x_dist %in% c("bray", "eucl") | !y_dist %in% c("bray", "eucl")) {
    stop(sprintf("(x,y)_dist must be one of 'bray' or 'eucl' for Bray-Curtis Dissimilarity and Euclidian Distance, respectively."))
  }
  if (!x_ord %in% c("pca", "pcoa") | !y_ord %in% c("pca", "pcoa")) {
    stop(sprintf("(x,y)_ord must be one of 'pca' or 'pcoa' for Principal Component Analysis and Principal Coordinate Analysis, respectively."))
  }
  if (!water %in% c("both", "dry", "wet")) {
    stop(sprintf("water_subset_pattern must be one of 'both', 'dry', 'wet'"))
  }
  if (!range %in% c('both', 'native', 'non-native')) {
    stop(sprintf("range_subset_pattern must be one of 'both', 'native', 'non-native'"))
  }

  range_subset_pattern <- switch(range,
                                 "both" = "(native)|(non-native)",
                                 "native" = "(?<!-)native",
                                 "non-native" = "non-native")
  water_subset_pattern <- switch(water,
                                 "both" = "(dry)|(wet)",
                                 "wet" = "wet",
                                 "dry" = "dry")

  require(stringr)
  require(dplyr)
  require(tidyr)

  library(ape) # fail if not available
  pcoa <- ape::pcoa
  pca <- stats::prcomp # bind to arg name

  # compute distances and ordinations
  dist_x <- vegan::vegdist(x_data, na.rm = TRUE, method = x_dist)
  ord_x <- do.call(x_ord, args = list(dist_x))
  dist_y <- vegan::vegdist(y_data, na.rm = TRUE, method = y_dist)
  ord_y <- do.call(y_ord, args = list(dist_y))

  #mantel test
  x_mantel_data <- x_data[str_detect(rownames(x_data), water_subset_pattern) & str_detect(rownames(x_data), range_subset_pattern), ]
  y_mantel_data <- y_data[str_detect(rownames(y_data), water_subset_pattern) & str_detect(rownames(y_data), range_subset_pattern), ]
  dist_mantel_x <- vegan::vegdist(x_mantel_data, na.rm = TRUE, method = x_dist)
  dist_mantel_y <- vegan::vegdist(y_mantel_data, na.rm = TRUE, method = y_dist)
  mantel_result <- ade4::mantel.randtest(dist_mantel_x, dist_mantel_y)

  # extract and tidy up ordination data for plotting
  f <- function(x, ord) {
    if (ord == "pca") {
      re <- data.frame(
        PC_1 = x$rotation[,1],
        PC_2 = x$rotation[,2]
      )
    }
    if (ord == "pcoa"){
      re <- data.frame(
        PCoA_1 = x$vectors[,1],
        PCoA_2 = x$vectors[,2]
      )
    }
    re
  }

  x <- list(f(ord_x, x_ord), f(ord_y, y_ord))
  names(x) <- c(as.character(deparse(substitute(x_data))), as.character(deparse(substitute(y_data))))

  point_data <- data.frame(
    Axis_1 = c(x[[1]][,1], x[[2]][,1]),
    Axis_2 = c(x[[1]][,2], x[[2]][,2]),
    dataset = c(rep(names(x)[1], length(x[[1]][,1])), rep(names(x)[2], length(x[[2]][,1])))
  )
  segment_data <- data.frame(
    from_x = x[[1]][,1],
    to_x = x[[2]][,1],
    from_y = x[[1]][,2],
    to_y = x[[2]][,2]
  )
  id <- data.frame(id = rownames(x[[1]])) %>% tidyr::separate_wider_delim(id, "_", names = c("pop_code", "water", "sf", "range"))

  point_data <- cbind(point_data, id)
  segment_data <- cbind(segment_data, id)


  stats_data <- data.frame(water = NA, range = NA, mantel_pvalue = NA, mantel_expl_var = NA, protest_pvalue = NA)
  # subset to matches of range and water patterns
  protest_results <- procrustes_results <- mantel_results <- list()
  set.seed(0)
  #mantel test
  x_mantel_data <- x_data[str_detect(rownames(x_data), water_subset_pattern) & str_detect(rownames(x_data), range_subset_pattern), ]
  y_mantel_data <- y_data[str_detect(rownames(y_data), water_subset_pattern) & str_detect(rownames(y_data), range_subset_pattern), ]
  dist_mantel_x <- vegan::vegdist(x_mantel_data, na.rm = TRUE, method = x_dist)
  dist_mantel_y <- vegan::vegdist(y_mantel_data, na.rm = TRUE, method = y_dist)
  mantel_result <- ade4::mantel.randtest(dist_mantel_x, dist_mantel_y)

  set.seed(0)
  #procrust test
  x_procrust_data = if (x_ord == "pca") ord_x[["rotation"]] else ord_x[["vectors"]]
  y_procrust_data = if (y_ord == "pca") ord_y[["rotation"]] else ord_y[["vectors"]]
  x_procrust_data <- x_procrust_data[str_detect(rownames(x_procrust_data), water_subset_pattern) & str_detect(rownames(x_procrust_data), range_subset_pattern), ]
  y_procrust_data <- y_procrust_data[str_detect(rownames(y_procrust_data), water_subset_pattern) & str_detect(rownames(y_procrust_data), range_subset_pattern), ]
  suppressWarnings(procrustes_result <- procrustes(X = x_procrust_data, Y = y_procrust_data))
  suppressWarnings(protest_result <- protest(X = x_procrust_data, Y = y_procrust_data))

  protest_results[[paste0(water_subset_pattern, "_", str_remove(range_subset_pattern, "\\(.+\\)"))]] <- protest_result
  procrustes_results[[paste0(water_subset_pattern, "_", str_remove(range_subset_pattern, "\\(.+\\)"))]] <- procrustes_result
  mantel_results[[paste0(water_subset_pattern, "_", str_remove(range_subset_pattern, "\\(.+\\)"))]] <- mantel_result

  stats_data$water <- water_subset_pattern
  stats_data$range <- str_remove(range_subset_pattern, "\\(.+\\)")
  stats_data$mantel_pvalue <- mantel_results[[1]]$pvalue
  stats_data$mantel_expl_var <- mantel_results[[1]]$expvar[3]
  stats_data$protest_pvalue <- protest_results[[1]]$signif
  stats_data$min_x <- min(point_data$Axis_1)
  stats_data$min_y <- min(point_data$Axis_2)
  stats_data$max_x <- max(point_data$Axis_1)
  stats_data$max_y <- max(point_data$Axis_2)

  return(list(stats = stats_data,
              points = point_data,
              segments = segment_data,
              procrustes = procrustes_result,
              protest = protest_result,
              mantel = mantel_result,
              ordination_x = ord_x,
              ordination_y = ord_y))
}


x <- ProcrustesByGroup(hellinger_transform(matrix_fungi), "bray", "pcoa",
                       # scale(matrix_root), "eucl", "pca",
                       hellinger_transform(matrix_exudates), "bray", "pcoa",
                       water = "wet",
                       range = "native")

x$stats
x$mantel
x$procrustes
x$protest

ScreePlot(x)

ggplot() +
  geom_segment(data = x$segments, aes(x = from_x, xend = to_x, y = from_y, yend = to_y), color = "#999999", linewidth = 0.2) +
  geom_point(data = x$points, aes(Axis_1, Axis_2, color = dataset), size = 1.5) +
  stat_ellipse(geom = "polygon", data = x$points, aes(Axis_1, Axis_2, color = dataset), size = 1, fill = NA) +
  facet_grid(water~range) +
  labs(x = "Axis 1", y = "Axis 2") +
  theme_linedraw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#ffffff"), strip.text = element_text(color = "#000000")) +
  geom_label(data = x$stats, aes(x =0, y = min_y*1.2, label = paste0("protest: p=",protest_pvalue,"\tmantel: p=", mantel_pvalue)))




d <- vegdist(hellinger_transform(matrix_fungi), method = "bray")
#d <- vegdist(matrix_fungi, method = "bray")
vegan::dbrda(d ~ range * water + sqrt(CWD) * water, data = variable_data) %>% anova(by = "term")
vegan::betadisper(d, paste0(variable_data$range, "_", variable_data$water)) %>% vegan::permutest(permutations = how(nperm = 999))
vegan::adonis2(d ~ range * water + sqrt(CWD) * water, variable_data, permutations = 999, by = "terms")

d <- vegdist(hellinger_transform(matrix_exudates), method = "bray")
# d <- vegdist(matrix_exudates, method = "bray")
vegan::dbrda(d ~ range * water + sqrt(CWD) * water, data = variable_data) %>% anova(by = "term")
vegan::betadisper(d, paste0(variable_data$range, "_", variable_data$water)) %>% vegan::permutest(permutations = how(nperm = 999))
vegan::adonis2(d ~ range * water + sqrt(CWD) * water, variable_data, permutations = 999, by = "terms")

d <- vegdist(matrix_root, method = "eucl")
vegan::dbrda(d ~ range * water + sqrt(CWD) * water, data = variable_data) %>% anova(by = "term")
vegan::betadisper(d, paste0(variable_data$range, "_", variable_data$water)) %>% vegan::permutest(permutations = how(nperm = 999))
vegan::adonis2(d ~ range * water + sqrt(CWD) * water, variable_data, permutations = 999, by = "terms")
