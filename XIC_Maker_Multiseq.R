# Packages ------------------------------------------------------------------------------------

library(here)
library(rawDiag)
library(MSnbase)
library(enviPat)
library(protViz)
library(Peptides)
library(OrgMassSpecR)
library(magrittr)
library(tictoc)
library(glue)
library(fs)
library(tibble)
library(gridExtra)
library(glue)
library(RSQLite)
library(purrr)
library(furrr)
library(readr)
library(openxlsx)
library(fs)
library(ggplot2)
library(stringr)
library(dplyr)
library(tictoc)
library(timeR)

# Initialize Parameters ---------------------------------------------------

# Raw File Directory

rawFileDir <- 
  "raw/"

#Raw File Name

rawFileName <- 
  "ZP_09032019_Norm_F3F4_A_01.raw"

# AA sequence to use for calculation of isotopic pattern

target_seq_file <- 
  "input/201909_EcoliBL21_GF_Norm_F01-12_2run_nomod_F3F4_tightAbsMass_proteoforms.csv"

# How many pforms to search (arranged by Q value)

top_n_pforms <- 4

# Use table with depleted isotope ratios?

use_depleted_isotopes <- TRUE

# Charge to use for calculating isotopic distribution

target_charges <- 
  c(1:30)

# m/z range in raw file

mz_range <- 
  c(600, 2000)

# XIC mass tolerance in ppm

XIC_tol <- 
  25

# Use iodoacetylation of all cysteines?

use_IAA <- 
  FALSE

# Abundance cutoff for isotopic distribution (percentage)

abund_cutoff <- 
  5

# Size of window for zoomed mass spectra

mz_window <- 4

# Make zoomed mass spectra with multiple scans?

make_multiscan_MS <- FALSE

# How many scans to average for zoomed masss spectra

nscans <- 3

# Output file size (width,height in inches)

outputWidth <- 16

# Output file DPI

outputDPI <- 150

# Functions ---------------------------------------------------------------

mem_change_selection <- function() {
  
  context <-  
    rstudioapi::getActiveDocumentContext()
  
  out <- 
    paste0("pryr::mem_change(~{
           ",
      context$selection[[1]]$text,
      "
      })")
  
  rstudioapi::modifyRange(context$selection[[1]]$range, out, id = context$id)

}

kickout <- function(list) {
  
  # This function removes any element from the list of input files
  # (from root/input) which does not have one of the allowed
  # extensions or which has "deprecated"
  
  allowed_ext <- c("raw")
  
  for (i in rev(seq_along(list))) {
    
    if (!(tools::file_ext(list[[i]]) %in% allowed_ext)) {
      
      list[[i]] <- NULL 
      
    } else if (str_detect(list[[i]], fixed("deprecated", TRUE)) == TRUE) {
      
      list[[i]] <- NULL 
      
    }
  }
  
  return(list)
}

kickoutXIC <- function(list) {
  
  # This function removes any element from the list of input files
  # (from root/input) which does not have one of the allowed
  # extensions or which has "deprecated"
  
  for (i in rev(seq_along(list))) {
    
    if (is.null(list[[i]]$intensities) == TRUE) {
      
      list[[i]] <- NULL 
      
    }
    
  }
  
  return(list)
}

stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
}

make_XIC_plot = function(df, x, y, seq_name,  theme = NULL) {
  
  df %>% 
    ggplot(aes({{x}}, {{y}})) +
    geom_line() +
    labs(
      title = glue::glue("{seq_name}") 
    ) +
    theme
  
}

make_spectrum_top1 = function(df, x, y, accession = "NA", scan_num = 0, charge = 0, xrange, theme  = NULL) {
  
  {
    xmin <- 
      df %>% 
      filter({{x}} == min({{x}})) %>%
      pull({{x}})
    
    xmax <- 
      df %>% 
      filter({{x}} == max({{x}})) %>%
      pull({{x}})
    
    if (xrange[[1]] < xmin) {
      
      xrange[[1]] <- xmin
      
    }
    
    if (xrange[[2]] > xmax) {
      
      xrange[[2]] <- xmax
      
    }
    
    ymax <-
      df %>% 
      filter({{x}} >= xrange[[1]] & {{x}} <= xrange[[2]]) %>%
      filter({{y}} == max({{y}})) %>%
      pull({{y}}) %>% 
      extract(1)
    
    if (length(ymax) == 0) {
      
      ymax <- 100
      
    } 
    
    if (all.equal(0, ymax) == TRUE) {
      
      ymax <- 100
      
    }
    
    scan_cap <- 
      df %>% 
      pull({{scan_num}}) %>% 
      .[[1]]
    
  }
  
  df %>%
    ggplot(aes({{x}}, {{y}})) +
    geom_line() +
    geom_vline(
      aes(xintercept = xrange[[1]]+((xrange[[2]]-xrange[[1]])/2), color = "red", alpha = 0.5)
    ) +
    annotate(
      "text",
      x = xrange[[2]],
      y = ymax,
      label = glue("{accession}\n Scan #{scan_cap}\n Charge +{charge}"),
      vjust="inward",
      hjust="inward",
      size = 3,
      alpha = 0.5
    ) +
    lims(
      x = xrange,
      y = c(0, ymax)
    ) +
    guides(
      color = "none",
      size = "none",
      alpha = "none"
    ) +
    labs(
      x = "m/z",
      y = "Intensity"
    ) +
    theme
  
}

make_spectrum_topN = function(df, x, y, accession = "NA", charge = 0, xrange, theme  = NULL) {
  
  {
    xmin <- 
      df %>% 
      filter({{x}} == min({{x}})) %>%
      pull({{x}})
    
    xmax <- 
      df %>% 
      filter({{x}} == max({{x}})) %>%
      pull({{x}})
    
    if (xrange[[1]] < xmin) {
      
      xrange[[1]] <- xmin
      
    }
    
    if (xrange[[2]] > xmax) {
      
      xrange[[2]] <- xmax
      
    }
    
    ymax <-
      df %>% 
      filter({{x}} >= xrange[[1]] & {{x}} <= xrange[[2]]) %>%
      filter({{y}} == max({{y}})) %>%
      pull({{y}}) %>% 
      extract(1)
    
    if (length(ymax) == 0) {
      
      ymax <- 100
      
    } 
    
    if (all.equal(0, ymax) == TRUE) {
      
      ymax <- 100
      
    }
    
  }
  
  df %>%
    ggplot(aes({{x}}, {{y}})) +
    geom_line() +
    geom_vline(
      aes(xintercept = xrange[[1]]+((xrange[[2]]-xrange[[1]])/2), color = "red", alpha = 0.5)
    ) +
    annotate(
      "text",
      x = xrange[[2]],
      y = ymax,
      label = glue("{accession}\n Charge +{charge}"),
      vjust="inward",
      hjust="inward",
      size = 3,
      alpha = 0.5
    ) +
    
    coord_cartesian(xlim = xrange, ylim = c(0, ymax), expand = TRUE) +
    guides(
      color = "none",
      size = "none",
      alpha = "none"
    ) +
    labs(
      x = "m/z",
      y = "Intensity"
    ) +
    theme
  
}


# ggplot themes -----------------------------------------------------------


XICtheme <- 
  list(
    geom_text(
      data = ~filter(.x, int_sum == max(int_sum)),
      aes(
        x = times,
        y = int_sum,
        label = glue(
          "Max TIC: {format(int_sum, scientific = TRUE, nsmall = 4, digits = 4)}\n RT@Max: {format(times, scientific = FALSE, nsmall = 4, digits = 3)}"),
        alpha = 0.5,
        size = 3
      ),
      vjust="inward",
      hjust="inward",
      nudge_x = 5
    ),
    ggthemes::theme_clean(base_size = 10),
    # theme(
    #   text = element_text(size = 6),
    #   title = element_text(size = 8),
    #   axis.text = element_text(size = 8),
    #   axis.title = element_text(size = 8)
    # ),
    labs(
      x = "Retention Time (min)",
      y = "Total Ion Current"
    ),
    guides(
      alpha = "none",
      size = "none"
    )
  )

MStheme01 <- 
  list(
    ggthemes::theme_clean(base_size = 12),
    labs(
      x = "m/z",
      y = "Intensity"
    ),
    guides(
      alpha = "none",
      size = "none"
    )
  )


# Make future workers and timer -------------------------------------------

timer <- 
  createTimer(verbose = FALSE)

timer$start("Make future workers")

if (top_n_pforms < 10) plan(multisession(workers = as.integer(top_n_pforms), gc = TRUE, persistent = FALSE))

if (top_n_pforms >= 10) plan(multisession(workers = 10L, gc = TRUE, persistent = FALSE))

timer$stop("Make future workers")

# Check filename and path -------------------------------------------------

rawFilesInDir <- 
  fs::dir_ls(
    rawFileDir, 
    recurse = TRUE,
    type = "file",
    regexp = c("[.]raw$")
  )

if (length(rawFilesInDir) == 0) {
  stop("No .raw files found in raw file directory")
}

if (rawFilesInDir %>% 
    str_detect(rawFileName) %>% 
    any() == FALSE) {
  stop("Input raw file not found in raw file directory")
}

rawFile <- 
  rawFilesInDir %>% 
  str_subset(rawFileName)

# Process target sequences ------------------------------------------------

target_seqs <- 
  target_seq_file %>% 
  readr::read_csv() %>% 
  group_by(UNIPROTKB) %>% 
  filter(GlobalQvalue == max(GlobalQvalue)) %>%
  ungroup() %>% 
  top_n(top_n_pforms, desc(GlobalQvalue)) %>% 
  select(UNIPROTKB, ProteoformSequence) %>% 
  deframe() %>% 
  as.list()

# Get elemental composition -----------------------------------------------

timer$start("Chemical formulas and isotopic distributions")

if (use_IAA == TRUE) {
  
  chemform_temp <- 
    map(
      target_seqs,
      ConvertPeptide,
      IAA = TRUE
    )
  
} else {
  
  chemform_temp <- 
    map(
      target_seqs,
      ConvertPeptide,
      IAA = FALSE
    )
  
}

chemform_names <- 
  chemform_temp %>% 
  map(unlist) %>% 
  map(names)

chemform <- 
  chemform_temp %>% 
  map(unlist) %>%
  map2(
    chemform_names,
    ~map2_chr(.y, .x, paste0)
  ) %>% 
  map(
    paste0,
    collapse = ""
  )

H_to_merge <- 
  target_charges %>% 
  map_chr(~paste0("H", .x)) %>% 
  list() %>% 
  rep(length(target_seqs))

chemform_to_merge <- 
  chemform %>% 
  map(
    ~rep(.x, length(target_charges))
  )

chemform_withH <- 
  map2(
    chemform_to_merge,
    H_to_merge,
    ~map2_chr(.x, .y, ~mergeform(.x, .y))
  )

for (i in seq_along(chemform_withH)) {
  
  if (
    str_detect(chemform_withH[[i]], fixed("S0"))
  ) {
    
    chemform_withH[[i]] <- 
      str_remove_all(chemform_withH[[i]], "S0")
    
  }
  
}

rm(chemform_temp)
rm(chemform_names)
rm(chemform_to_merge)
rm(H_to_merge)

# Calculate isotopic dist -------------------------------------------------

if (use_depleted_isotopes == FALSE) {
  
  data(isotopes)
  
  isotopes_to_use <- isotopes
  
} else {
  
  isotopes_to_use <-
    readr::read_csv("input/depleted_isotopes.csv") %>% 
    as.data.frame()
  
}

check_chemform(isotopes = isotopes_to_use, chemforms = chemform) %>% 
  pull(warning) %>% 
  any() %>% 
  `if`(., stop("Problem with chemical formula"))

iso_dist <- 
  future_map2(
    chemform_withH,
    list(target_charges) %>% rep(length(target_seqs)),
    ~map2(
      .x,
      .y,
      ~isopattern(
        isotopes_to_use,
        chemform=.x,
        threshold=0.1,
        plotit=FALSE,
        charge=.y,
        emass=0.00054858,
        algo=1
      ) 
    )
  )

target_charges_list <- 
  list(target_charges) %>%
  rep(length(target_seqs)) %>% 
  map(as.list)

iso_dist_list_union <-
  iso_dist %>% 
  modify_depth(3, as_tibble) %>% 
  map2(
    .y = target_charges_list,
    ~map2(
      .x = .x,
      .y = .y,
      ~map2(
        .x = .x,
        .y = .y,
        ~mutate(.x, charge = .y)
      )
    )
  ) %>%
  map(reduce, union_all) %>%
  map(reduce, union_all) %>%
  map(
    ~ filter(
      .x,
      abundance > abund_cutoff & `m/z` > mz_range[[1]] & `m/z` < mz_range[[2]]
    )
  )

timer$stop("Chemical formulas and isotopic distributions")

# Calculate XIC -----------------------------------------------------------

timer$start("Calculate and plot XICs")


XIC_target_mz <- 
  iso_dist_list_union %>% 
  map(~pull(.x, `m/z`)) %>% 
  map(unique)

XIC <- 
  future_map(
    XIC_target_mz,
    ~readXICs(
      rawfile = rawFile,
      masses = .x,
      tol = XIC_tol
    )
  )


XIC_nonull <- 
  XIC %>%
  map(kickoutXIC)

sumXIC1 <- 
  XIC_nonull %>% 
  modify_depth(2, as_tibble) %>% 
  modify_depth(2, ~select(.x, times, intensities)) %>%
  future_map(reduce, full_join, .progress = TRUE)

sumXIC2 <-  
  sumXIC1 %>% 
  future_map(~group_by(.x, times)) %>% 
  future_map(~summarize(.x, int_sum = sum(intensities)))

sumXIC_summary <- 
  sumXIC2 %>% 
  set_names(chemform) %>% 
  map(
    ~summarize(.x, max_TIC = max(int_sum))
  ) %>% 
  map2(
    names(.),
    ~mutate(.x, chem_form = .y)
  ) %>% 
  reduce(union_all) %>% 
  mutate(Sequence = unlist(target_seqs)) %>% 
  mutate(seq_name = names(target_seqs)) %>% 
  select(seq_name, Sequence, chem_form, max_TIC)

# Get RT corresponding to maximum TIC for each sequence for later use in
# making spectra

RT_of_maxTIC <- 
  sumXIC2 %>% 
  map(
    ~filter(.x, int_sum == max(int_sum))
  ) %>% 
  map(
    ~pull(.x, times) 
  )

RT_of_maxTIC_topN <- 
  sumXIC2 %>% 
  map(
    ~top_n(.x, nscans, int_sum)
  ) %>% 
  map(
    ~pull(.x, times) 
  )

# Plot XICs ---------------------------------------------------------------

if (dir.exists("output/")== FALSE) dir.create("output/")
if (dir.exists("output/XIC_summary/")== FALSE) dir.create("output/XIC_summary/")

systime <- format(Sys.time(), "%Y%m%d_%H%M")

XIC_plots <- 
  sumXIC2 %>% 
  imap(
    ~make_XIC_plot(
      .x,
      times,
      int_sum,
      .y,
      XICtheme
    )
  )

# MultiArrange XIC grobs -------------------------------------------------------

XIC_groblist <- 
  marrangeGrob(
    grobs = XIC_plots,
    ncol = 5,
    nrow = 3,
    top = rawFileName
  )

timer$stop("Calculate and plot XICs")

# Save XIC results --------------------------------------------------------

timer$start("Save XIC results")

if (dir_exists(paste0("output/XIC_summary/", systime, "/")) == FALSE) dir_create(paste0("output/XIC_summary/", systime, "/"))

# Save spreadsheet data

readr::write_csv(
  sumXIC_summary,
  path = 
    paste(
      "output/XIC_summary/",
      systime,
      "_",
      fs::path_ext_remove(rawFileName),
      "_XIC_summary.csv",
      sep = ""
    )
)

writexl::write_xlsx(
  iso_dist_list_union,
  path = 
    paste(
      "output/XIC_summary/",
      systime,
      "_",
      fs::path_ext_remove(rawFileName),
      "_isotopic_distributions.xlsx",
      sep = ""
    )
)

# Save chromatograms, all together

XIC_groblist_filename <- 
  paste(
    "output/XIC_summary/",
    systime,
    "/",
    fs::path_ext_remove(rawFileName),
    "_XICs",
    sep = ""
  )

ggsave(
  filename = paste0(XIC_groblist_filename, ".pdf"),
  plot = XIC_groblist,
  width = 20,
  height = 12,
  limitsize = FALSE
)

timer$stop("Save XIC results")

# Make MS, Top 1 ---------------------------------------------------

timer$start("Make MS, Top 1 most intense PART 1")

scanNumber_and_RT <- 
  read.raw(
    rawFile,
    rawDiag = FALSE
  ) %>% 
  as_tibble() %>% 
  select(scanNumber, StartTime)

scanNumsToRead <- 
  map(
    RT_of_maxTIC,
    ~filter(scanNumber_and_RT, StartTime == .x)
  ) %>% 
  map(
    ~pull(.x, scanNumber) 
  )

scansToPlot <- 
  future_map(
    scanNumsToRead,
    ~readScans(
      rawFile,
      scans = .x
    )
  )

spectra_highestTIC <- 
  imap(
    scansToPlot,
    ~tibble(
      UNIPROTKB = .y,
      scan = .x[[1]]$scan,
      mz = .x[[1]]$mZ,
      intensity = .x[[1]]$intensity)
  )

mz_max_abund <- 
  iso_dist_list_union %>% 
  map(
    ~filter(.x, abundance == 100) 
  ) %>% 
  map(
    ~pull(.x, `m/z`) 
  ) %>% 
  map(as.list)

mz_max_abund_charge <- 
  iso_dist_list_union %>% 
  map(
    ~filter(.x, abundance == 100) 
  ) %>% 
  map(
    ~pull(.x, charge) 
  ) %>% 
  map(as.list)

spectra_highestTIC_list <- 
  spectra_highestTIC %>% 
  map(list) %>% 
  map2(
    mz_max_abund,
    ~rep(.x, length(.y))
  )

timer$stop("Make MS, Top 1 most intense PART 1")

timer$start("Make MS, Top 1 most intense, PART 2")

spectra_highestTIC_plots <- 
  future_pmap(
    list(
      spectra_highestTIC_list,
      names(spectra_highestTIC_list) %>% as.list(),
      mz_max_abund,
      mz_max_abund_charge
    ),
    ~pmap(
      list(
        ..1,
        ..2,
        ..3,
        ..4
      ),
      ~make_spectrum_top1(
        df = ..1,
        x = mz,
        y = intensity,
        accession = ..2,
        scan_num = scan,
        charge = ..4,
        xrange = c(..3 - (mz_window/2), ..3 + (mz_window/2)),
        theme = NULL
      )
    )
  )

timer$stop("Make MS, Top 1 most intense, PART 2")

# Arrange MS Grobs, Top 1 ----------------------------------------------

timer$start("Arrange MS grobs, Top 1 most intense")

tablegrob_list_top1 <-
  future_map(
    spectra_highestTIC_plots,
    ~arrangeGrob(
      grobs = .x,
      ncol = 4,
      top = rawFileName
    )
  )

tablegrob_list_top1_arranged <-
  future_map(
    tablegrob_list_top1,
    ~gridExtra::grid.arrange(.x),
    .progress = TRUE
  )

timer$stop("Arrange MS grobs, Top 1 most intense")

# Save arranged MS, Top 1 -------------------------------------------------

timer$start("Save MS, Top 1")

if (dir_exists(paste0("output/XIC_summary/", systime, "/mass_spec/")) == FALSE) {
  
  dir_create(paste0("output/XIC_summary/", systime, "/mass_spec/"))
  
}
  
tablegrob_filenames <-
  names(target_seqs) %>%
  as.list() %>%
  map(
    ~paste(
      "output/XIC_summary/",
      systime,
      "/mass_spec/",
      fs::path_ext_remove(rawFileName),
      "_",
      .x,
      "_spectrum_zooms.png",
      sep = ""
    )
  )

tablegrob_heights <-
  tablegrob_list_top1 %>%
  map(use_series, "heights") %>%
  map(length)

future_pwalk(
  list(
    tablegrob_filenames,
    tablegrob_list_top1_arranged,
    tablegrob_heights
  ),
  ~ggsave(
    filename = ..1 ,
    plot = ..2,
    width = outputWidth,
    height = ..3 * 3,
    dpi = outputDPI
  )
)

message("\n\n Done with top 1!")

timer$stop("Save MS, Top 1")
