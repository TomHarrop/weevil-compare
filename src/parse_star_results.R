# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


library(data.table)

#############
# FUNCTIONS #
#############

FindSpeciesDirs <- function(x){
    no_idx_dirs[grepl(x, no_idx_dirs)]
}

ReadStarLog <- function(x){
    my_dt <- fread(x,
                   sep = "|",
                   fill = TRUE,
                   header = FALSE,
                   col.names = c("variable", "value"))
    my_dt[, value := gsub("[[:space:]]", "", value)]
    my_dt[, value := as.numeric(gsub("%", "", value))]
    return(my_dt[!is.na(value)])
}

###########
# GLOBALS #
###########

af_file <- snakemake@input[["assembly_filenames"]]
rnaseq_sample_key_file <- snakemake@input[["sample_key"]]
star_dir <- snakemake@params[["star_dir"]]

full_results_file <- snakemake@output[[1]]

# dev
# af_file <- "data/assembly_filenames.txt"
# rnaseq_sample_key_file <- "data/full_sample_key.csv"
# star_dir <- "output/030_map"

########
# MAIN #
########

# search the directory for log files
all_dirs <- list.dirs(star_dir, recursive = FALSE)
no_idx_dirs <- all_dirs[!grepl("^star-index", basename(all_dirs))]


sample_key <- fread(rnaseq_sample_key_file)
af <- fread(af_file)

# get the star results per-species
spec_results <- af[, .(spec_result_dir = FindSpeciesDirs(species_name)),
                   by = species_name]
spec_results[, sample_id := gsub(paste0("_", species_name),
                                 "",
                                 basename(spec_result_dir)),
             by = .(spec_result_dir, species_name)]

spec_results[, star_result_file := list.files(path = spec_result_dir,
                                              pattern = "Log.final.out",
                                              recursive = FALSE,
                                              full.names = TRUE),
             by = spec_result_dir]

# merge the sample key
spec_results[, ogbf_id := gsub("^[^-]*-([[:digit:]]+-[[:digit:]]+).*",
                               "\\1",
                               sample_id)]
spec_results_with_tissue <- merge(spec_results,
                                  sample_key[, .(ogbf_id = OGF_sample_ID,
                                                 tissue = Tissue)],
                                  by = "ogbf_id",
                                  all.x = TRUE,
                                  all.y = FALSE)

# generate a list of files to read
file_list <- spec_results_with_tissue[!is.na(star_result_file),
                                      structure(star_result_file,
                                                names = star_result_file)]

# read the STAR results
star_results_list <- lapply(file_list, ReadStarLog)
star_results <- rbindlist(star_results_list, idcol = "star_result_file")

# merge star logs with metadata
full_results <- merge(spec_results_with_tissue,
                      star_results, by = "star_result_file")


fwrite(full_results, full_results_file)

sessionInfo()

quit("no", status = 0)

###########
# NOT RUN #
###########

# generate plot data
vars_to_plot <- c("Uniquely mapped reads %" = "Uniquely mapped",
  "% of reads mapped to multiple loci" = "Multimapped", 
  "% of reads unmapped: too short" = "Unmapped")

pd <- full_results[variable %in% names(vars_to_plot)]

pd[, grp := factor(plyr::revalue(variable, vars_to_plot),
                   levels = vars_to_plot)]

ggplot(pd, aes(x = species_name, y = value, colour = tissue)) +
    coord_flip() +
    facet_wrap(~ grp, ncol = 1) +
    xlab(NULL) + ylab("%") +
    scale_colour_viridis_d(guide = guide_legend(title = NULL)) +
    geom_point(position = position_jitter(width = 0.2))
