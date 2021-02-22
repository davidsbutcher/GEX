#' make_every_XIC_MS2
#'
#' @param rawFileDir
#' @param rawFileName
#' @param targetSeqData
#' @param outputDir
#' @param target_col_name
#' @param target_sequence_col_name
#' @param PTMname_col_name
#' @param PTMformula_col_name1
#' @param PTMformula_col_name2
#' @param isoAbund
#' @param fragment_mass_range
#' @param fragment_charges
#' @param fragment_mz_range
#' @param fragment_pos_cutoff
#' @param XIC_tol_MS2
#' @param use_IAA
#' @param save_output
#' @param abund_cutoff
#' @param hClust_height
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples

make_every_XIC_MS2 <-
   function(
      rawFileDir = NULL,
      rawFileName = NULL,
      targetSeqData = NULL,
      outputDir = getwd(),
      target_col_name = c("UNIPROTKB"),
      target_sequence_col_name = c("ProteoformSequence"),
      PTMname_col_name = c("PTMname"),
      PTMformula_col_name1 = c("FormulaToAdd"),
      PTMformula_col_name2 = c("FormulaToSubtract"),
      isoAbund = c("12C" = 0.9893, "14N" = 0.99636),
      fragment_mass_range = c(0,100000),
      fragment_charges = c(1:5),
      fragment_mz_range = c(600,2000),
      fragment_pos_cutoff = 20,
      XIC_tol_MS2 = 10,
      XIC_TIC_cutoff = 25000,
      resPowerMS2 = 150000,
      mz_window = 3,
      use_IAA = FALSE,
      save_output = TRUE,
      abund_cutoff = 5,
      hClust_height = 0.005
   ) {

      # Load rawrr package

      library(rawrr)

      options(dplyr.summarise.inform = FALSE)

      # Assertions --------------------------------------------------------------

      assertthat::assert_that(
         assertthat::is.dir(rawFileDir),
         msg = "rawFileDir is not a recognized path"
      )


      assertthat::assert_that(
         assertthat::is.string(rawFileName),
         msg = "rawFileName is not a string"
      )

      assertthat::assert_that(
         assertthat::has_extension(targetSeqData, "csv") |
            is.data.frame(targetSeqData),
         msg = "targetSeqData is not a CSV file or dataframe"
      )

      if (
         assertthat::see_if(assertthat::has_extension(targetSeqData, "csv"))
         == TRUE
      ) {

         assertthat::assert_that(
            assertthat::is.readable(targetSeqData),
            msg = "targetSeqData is not a readable file"
         )

      }

      assertthat::assert_that(
         assertthat::is.dir(outputDir),
         msg = "outputDir is not a recognized directory"
      )

      assertthat::assert_that(
         assertthat::is.flag(use_IAA),
         msg = "use_IAA should be TRUE or FALSE"
      )

      assertthat::assert_that(
         assertthat::is.flag(save_output),
         msg = "save_output should be TRUE or FALSE"
      )


      # Create necessary quosures -----------------------------------------------

      target_col_name <- rlang::enquo(target_col_name)
      target_sequence_col_name_sym <- rlang::sym(target_sequence_col_name)
      target_sequence_col_name <- rlang::enquo(target_sequence_col_name)
      PTMformula_col_name1 <- rlang::enquo(PTMformula_col_name1)
      PTMformula_col_name2 <- rlang::enquo(PTMformula_col_name2)
      PTMname_col_name <- rlang::enquo(PTMname_col_name)

      # Check filename and path -------------------------------------------------

      ## purrr::as_mapper(purrr::reduce)

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
          stringr::str_detect(rawFileName) %>%
          any() == FALSE) {
         stop("Input raw file not found in raw file directory")
      }

      rawFile <-
         rawFilesInDir %>%
         stringr::str_subset(rawFileName)

      # ggplot themes -----------------------------------------------------------

      XICtheme <-
         list(
            ggplot2::geom_text(
               data = ~dplyr::filter(.x, int_sum == max(int_sum)),
               ggplot2::aes(
                  x = times,
                  y = int_sum,
                  label = glue::glue(
                     "Max TIC: {format(int_sum, scientific = TRUE, nsmall = 4, digits = 4)}\n RT@Max: {format(times, scientific = FALSE, nsmall = 4, digits = 3)}\n PTM:{!!PTMname}"),
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
            ggplot2::labs(
               x = "Retention Time (min)",
               y = "Total Ion Current"
            ),
            ggplot2::guides(
               alpha = "none",
               size = "none"
            )
         )

      MStheme01 <-
         list(
            ggthemes::theme_clean(base_size = 12),
            ggplot2::labs(
               x = "m/z",
               y = "Intensity"
            ),
            ggplot2::guides(
               alpha = "none",
               size = "none"
            )
         )


      # Load isotopes table -----------------------------------------------------

      data(isotopes, package = "enviPat")

      isotopes_to_use <- isotopes

      # Get positions of specific isotopes in dataframe

      index_12C <-
         which(isotopes_to_use$isotope == "12C") %>%
         .[[1]]

      index_13C <-
         which(isotopes_to_use$isotope == "13C") %>%
         .[[1]]

      index_14N <-
         which(isotopes_to_use$isotope == "14N") %>%
         .[[1]]

      index_15N <-
         which(isotopes_to_use$isotope == "15N") %>%
         .[[1]]

      # Replace values in dataframe with supplied values

      isotopes_to_use$abundance[index_12C] <-
         isoAbund[which(names(isoAbund) == "12C")]

      isotopes_to_use$abundance[index_13C] <-
         1 - isoAbund[which(names(isoAbund) == "12C")]

      isotopes_to_use$abundance[index_14N] <-
         isoAbund[which(names(isoAbund) == "14N")]

      isotopes_to_use$abundance[index_15N] <-
         1 - isoAbund[which(names(isoAbund) == "14N")]

      # Process target sequences ------------------------------------------------

      target_seqs_df <-
         targetSeqData %>%
         {if (!is.data.frame(.)) readr::read_csv(.) else .} %>%
         dplyr::mutate(
            MonoisoMass =
               Peptides::mw(!!target_sequence_col_name_sym, monoisotopic = TRUE)
         )

      target_seqs <-
         target_seqs_df %>%
         dplyr::select(
            !!target_col_name, !!target_sequence_col_name
         ) %>%
         tibble::deframe() %>%
         as.list()

      # Get elemental compositions ----------------------------------------------

      # Generate all possible fragments

      fragments <-
         purrr::map2(
            target_seqs,
            rep(list(fragment_charges), length(target_seqs)),
            ~purrr::map2(
               .x,
               .y,
               ~MSnbase::calculateFragments(
                  .x,
                  z = .y,
                  type=c("b", "y"),
                  neutralLoss = NULL
               ) %>%
                  tibble::as_tibble() %>%
                  dplyr::filter(
                     mz > fragment_mz_range[[1]] & mz < fragment_mz_range[[2]],
                     pos >= fragment_pos_cutoff[[1]] & pos <= fragment_pos_cutoff[[2]]
                  )
            ) %>%
               purrr::reduce(dplyr::union_all)
         )

      # Generate chemical formulas

      chemformulas_noH <-
         purrr::map(
            fragments,
            ~dplyr::pull(.x, seq) %>%
               as.list() %>%
               purrr::map(
                  ~OrgMassSpecR::ConvertPeptide(.x) %>%
                     unlist() %>%
                     purrr::imap_chr(
                        ~paste0(.y, .x, collapse = '')
                     ) %>%
                     paste0(collapse = '')
               ) %>%
               purrr::map(
                  ~stringr::str_remove_all(.x, "C0|H0|N0|O0|P0|S0")
               ) %>%
               unlist()
         )

      # Add one H per charge to chemical formulas

      charges <-
         purrr::map(
            fragments,
            ~dplyr::pull(.x, z)
         )

      chemformulas_withH <-
         purrr::map2(
            chemformulas_noH,
            charges,
            ~purrr::map2_chr(
               .x,
               .y,
               ~enviPat::mergeform(
                  .x,
                  paste0(
                     "H",
                     .y
                  )
               )
            )
         )

      # Add chemical formulas to fragments data

      fragments2 <-
         purrr::map2(
            fragments,
            chemformulas_withH,
            ~dplyr::mutate(
               .x,
               chemform = .y
            )
         )

      # Generate isotopic distributions ---------


      iso_dist <-
         purrr::map_depth(
            fragments2,
            1,
            ~enviPat::isopattern(
               isotopes_to_use,
               chemform=.x$chemform,
               threshold=0.1,
               plotit=FALSE,
               charge=.x$z,
               emass=0.00054858,
               algo=1,
               verbose = F
            ) %>%
               {
                  purrr::pmap(
                     list(
                        .,
                        .x$z,
                        .x$ion
                     ),
                     ~tibble::as_tibble(..1) %>%
                        dplyr::select(`m/z`, abundance) %>%
                        dplyr::filter(abundance > abund_cutoff) %>%
                        dplyr::mutate(
                           charge = ..2,
                           ion = ..3
                        )
                  ) %>%
                     purrr::reduce(dplyr::union_all)
               }
         )

      # Cluster isotopic distributions

      iso_dist_cluster <-
         purrr::imap(
            iso_dist,
            ~dplyr::mutate(
               .x,
               cluster =
                  cutree(
                     hclust(
                        dist(
                           `m/z`, method = "maximum"), method = "centroid"
                     ),
                     h = hClust_height
                  )
            ) %>%
               dplyr::group_by(cluster) %>%
               dplyr::summarise(
                  `m/z` = mean(`m/z`),
                  abundance = sum(abundance), # All clustered peak abundances are summed
                  charge = mean(charge),
                  ion = dplyr::first(ion)
               ) %>%
               dplyr::filter(charge %% 1 == 0) %>% # Remove all incorrectly clustered peaks
               dplyr::ungroup()
         )

      # Split all isotopic distributions by ion name

      fragment_groups <-
         purrr::map(
            iso_dist_cluster,
            # ~dplyr::group_by(.x, ion, charge) %>%
            ~dplyr::group_by(.x, ion) %>%
               dplyr::group_split()
         )

      # Read XICs --------

      XIC_MS2 <-
         purrr::map(
            fragment_groups,
            ~purrr::map(
               .x,
               ~rawrr::readChromatogram(
                  rawfile = rawFile,
                  mass = .x$`m/z`,
                  filter = "ms2",
                  tol = XIC_tol_MS2
               )
            )
         )

      # Make blank tibble to replace NULL tibbles

      blank_tibble <-
         tibble::tibble(times = 0, intensities = 0)

      # Process XICs by replacing NULL tibbles and adding together all XICs
      # for same ion and charge

      # charges_to_delete_A <-
      #    XIC_MS2 %>%
      #    purrr::map_depth(
      #       3,
      #       ~tibble::tibble(
      #          times = .x$times,
      #          intensities = .x$intensities
      #       )
      #    )
      #
      # charges_to_delete_B <-
      #    purrr::map_depth(
      #       charges_to_delete_A,
      #       3,
      #       ~which(nrow(.x) == 0)
      #    )

      XIC_MS2_sum <-
         XIC_MS2 %>%
         purrr::map_depth(
            3,
            ~tibble::tibble(
               times = .x$times,
               intensities = .x$intensities
            )
         ) %>%
         purrr::map_depth(
            3,
            ~{if (is.null(dplyr::pull(.x, times)) == TRUE) blank_tibble else .x}
         ) %>%
         purrr::map_depth(
            2,
            dplyr::bind_rows
         ) %>%
         purrr::map_depth(
            2,
            ~dplyr::group_by(.x, times) %>%
               dplyr::summarize(int_sum = sum(intensities))
         )

      # XIC_MS2 is large and no longer needed, delete it

      rm(XIC_MS2)

      # Remove fragments with TIC < XIC_TIC_cutoff from
      # fragment_groups and XIC_MS2_sum

      fragments_to_retain <-
         XIC_MS2_sum %>%
         purrr::map_depth(
            2,
            ~dplyr::filter(.x, int_sum == max(int_sum)) %>%
               dplyr::pull(int_sum)
         ) %>%
         purrr::map(
            unlist
         ) %>%
         purrr::map(
            ~which(.x > XIC_TIC_cutoff)
         )

      fragment_groups_trunc <-
         purrr::map2(
            fragment_groups,
            fragments_to_retain,
            ~.x[.y]
         )

      XIC_MS2_sum_trunc <-
         purrr::map2(
            XIC_MS2_sum,
            fragments_to_retain,
            ~.x[.y]
         )

      # Get ion names, charges, and max TIC from XIC for plotting

      ion_names <-
         fragment_groups_trunc %>%
         purrr::map_depth(
            2,
            ~dplyr::pull(.x, ion) %>%
               dplyr::first()
         )

      ion_charges <-
         fragment_groups_trunc %>%
         purrr::map_depth(
            2,
            ~dplyr::pull(.x, charge) %>%
               dplyr::first()
         )

      XIC_maxTIC_MS2 <-
         XIC_MS2_sum_trunc %>%
         purrr::map_depth(
            2,
            ~dplyr::filter(.x, int_sum == max(int_sum)) %>%
               dplyr::pull(int_sum)
         )


      # Plot XICs -------------------------------------------------------------

      XIC_plots_MS2 <-
         purrr::pmap(
            list(
               XIC_MS2_sum_trunc,
               as.list(names(fragment_groups_trunc)),
               ion_names,
               ion_charges,
               XIC_maxTIC_MS2
            ),
            ~purrr::pmap(
               list(
                  ..1,
                  ..2,
                  ..3,
                  ..4,
                  ..5
               ),
               ~make_XIC_plot_MS2(
                  ..1,
                  times,
                  int_sum,
                  ..2,
                  ..3,
                  ..4,
                  ..5
               )
            )
         )


      XIC_plots_MS2_marrange <-
         gridExtra::marrangeGrob(
            grobs = purrr::flatten(XIC_plots_MS2),
            ncol = 5,
            nrow = 3,
            top = rawFileName
         )

      # Save XIC results ------

      ## Create save directory ---------------------------------------------------

      systime <- format(Sys.time(), "%Y%m%d")

      if (all(isoAbund == c("12C" = 0.9893, "14N" = 0.99636))) {

         saveDir <-
            fs::path(
               outputDir,
               paste0(
                  fs::path_ext_remove(rawFileName),
                  "_",
                  length(target_seqs),
                  "seqs_IsoNorm"
               )
            ) %>%
            stringr::str_trunc(246, "right", ellipsis = "")

      } else if (isoAbund[[1]] < 0.9893 | isoAbund[[2]] < 0.99636) {

         saveDir <-
            fs::path(
               outputDir,
               paste0(
                  fs::path_ext_remove(rawFileName),
                  "_",
                  length(target_seqs),
                  "seqs_IsoDep"
               )
            ) %>%
            stringr::str_trunc(246, "right", ellipsis = "")

      } else {

         saveDir <-
            fs::path(
               outputDir,
               paste0(
                  fs::path_ext_remove(rawFileName),
                  "_",
                  length(target_seqs),
                  "seqs"
               )
            ) %>%
            stringr::str_trunc(246, "right", ellipsis = "")

      }

      if (dir.exists(saveDir) == FALSE) dir.create(saveDir)

      ## Save target seqs -----

      message(
         paste0(
            "Saving target sequences to ",
            saveDir,
            "/",
            fs::path_ext_remove(
               stringr::str_trunc(
                  rawFileName, 90, "right", ellipsis = ""
               )
            ),
            '_seqs.csv'
         )
      )

      readr::write_csv(
         target_seqs_df,
         path =
            fs::path(
               saveDir,
               paste0(
                  fs::path_ext_remove(
                     stringr::str_trunc(
                        rawFileName, 90, "right", ellipsis = ""
                     )
                  ),
                  '_seqs_MS2.csv'
               )
            )
      )

      ## Save theo iso distributions -----

      writexl::write_xlsx(
         iso_dist_cluster,
         path =
            fs::path(
               saveDir,
               paste0(
                  fs::path_ext_remove(
                     stringr::str_trunc(
                        rawFileName, 90, "right", ellipsis = ""
                     )
                  ),
                  '_isodist_MS2.xlsx'
               )
            )
      )


      ## Save XIC plots -----

      ggplot2::ggsave(
         filename =
            fs::path(
               saveDir,
               paste0(
                  fs::path_ext_remove(
                     stringr::str_trunc(
                        rawFileName, 90, "right", ellipsis = ""
                     )
                  ),
                  '_MS2_XICs.pdf'
               )
            ),
         plot = XIC_plots_MS2_marrange,
         width = 20,
         height = 12,
      )


      # Extract scans and make spectra ------------------------------------------

      rawFileMetadata <-
         rawrr::readIndex(
            rawFile
         ) %>%
         tibble::as_tibble() %>%
         dplyr::mutate(RT_min = round(rtinseconds/60, digits = 5))

      scan_nums <-
         rawFileMetadata %>%
         dplyr::pull(scan)

      RT_of_maxTIC_MS2 <-
         XIC_MS2_sum_trunc %>%
         purrr::map_depth(
            2,
            ~dplyr::filter(.x, int_sum == max(int_sum)) %>%
               dplyr::rename("RT_min" = times)
         )

      fragments_maxTIC <-
         purrr::map_depth(
            RT_of_maxTIC_MS2,
            2,
            ~dplyr::mutate(
               .x,
               scan =
                  scan_nums[MALDIquant::match.closest(.x$RT_min, rawFileMetadata$RT_min)]
            )
         ) %>%
         purrr::map2(
            fragment_groups_trunc,
            ~purrr::map2(
               .x,
               .y,
               ~dplyr::bind_cols(
                  .x,
                  dplyr::filter(.y, abundance == max(abundance))
               )
            )
         ) %>%
         purrr::map_depth(
            2,
            ~dplyr::group_by(.x, charge) %>%
               dplyr::group_split()
         )


      mz_to_read <-
         purrr::map_depth(
            fragments_maxTIC,
            3,
            ~dplyr::pull(.x, `m/z`)
         )

      scan_to_read <-
         purrr::map_depth(
            fragments_maxTIC,
            3,
            ~dplyr::pull(.x, scan) %>%
               dplyr::first()
         )


      scans_to_plot <-
         purrr::pmap(
            list(
               scan_to_read,
               mz_to_read
            ),
            ~purrr::pmap(
               list(
                  ..1,
                  ..2
               ),
               ~purrr::pmap(
                  list(
                     ..1,
                     ..2
                  ),
                  ~rawrr::readSpectrum(
                     rawfile = rawFile,
                     scan = ..1
                  ) %>%
                     purrr::flatten() %>%
                     {
                        tibble::tibble(
                           mz = .$mZ,
                           intensity = .$intensity
                        ) %>%
                           dplyr::filter(
                              mz >= ..2 - mz_window/2,
                              mz <= ..2 + mz_window/2
                           )
                     }
               )
            )
         )

      # Make list of charges in fragments_maxTIC and filter iso_dist_cluster_trunc
      # to remove unneeded charges

      charges_to_retain <-
         fragments_maxTIC %>%
         purrr::map_depth(
            3,
            ~dplyr::pull(.x, charge)
         ) %>%
         purrr::map_depth(
            2,
            ~unlist(.x)
         )

      iso_dist_cluster_trunc <-
         purrr::map_depth(
            iso_dist_cluster,
            1,
            ~dplyr::group_by(.x, ion) %>%
               dplyr::group_split()
         ) %>%
         purrr::map_depth(
            2,
            ~dplyr::group_by(.x, charge) %>%
               dplyr::group_split()
         ) %>%
         purrr::map2(
            fragments_to_retain,
            ~.x[.y]
         ) %>%
         purrr::map_depth(
            2,
            ~dplyr::bind_rows(.x)
         ) %>%
         purrr::map2(
            charges_to_retain,
            ~purrr::map2(
               .x,
               .y,
               ~dplyr::filter(.x, charge %in% .y)
            )
         ) %>%
         purrr::map_depth(
            2,
            ~dplyr::group_by(.x, charge) %>%
               dplyr::group_split()
         )




      ion_names_MS2 <-
         iso_dist_cluster_trunc %>%
         purrr::map_depth(
            3,
            ~dplyr::pull(.x, ion) %>%
               dplyr::first()
         )

      ion_charges_MS2 <-
         iso_dist_cluster_trunc %>%
         purrr::map_depth(
            3,
            ~dplyr::pull(.x, charge) %>%
               dplyr::first()
         )


      mz_all_abund_MS2 <-
         purrr::map_depth(
            iso_dist_cluster_trunc,
            3,
            ~dplyr::pull(.x, `m/z`)
         )

      spectra_MS2 <-
         purrr::pmap(
            list(
               scans_to_plot,
               as.list(names(fragment_groups_trunc)),
               mz_to_read,
               scan_to_read,
               ion_names_MS2,
               ion_charges_MS2,
               mz_all_abund_MS2,
               mz_to_read
            ),
            ~purrr::pmap(
               list(
                  ..1,
                  ..2,
                  ..3,
                  ..4,
                  ..5,
                  ..6,
                  ..7,
                  ..8
               ),
               ~purrr::pmap(
                  list(
                     ..1,
                     ..2,
                     ..3,
                     ..4,
                     ..5,
                     ..6,
                     ..7,
                     ..8
                  ),
                  ~make_spectrum_MS2(
                     ..1,
                     mz,
                     intensity,
                     name = ..2,
                     max_mz = ..3,
                     scan = ..4,
                     ion = ..5,
                     charge = ..6,
                     mz_all_abund = ..7,
                     isotopologue_window = ..8/resPowerMS2,
                     theme = MStheme01
                  )
               )
            )
         ) %>%
         purrr::map_depth(
            2,
            ~ggplot_list_checker(.x)
         )

      message("Making tablegrob list for PDF output")

      spectra_MS2_marrange <-
         unlist(spectra_MS2, recursive = F) %>%
         purrr::flatten() %>%
         gridExtra::marrangeGrob(
            grobs = .,
            ncol = 5,
            nrow = 3,
            top = paste0(rawFileName)
         )

      ggplot2::ggsave(
         filename =
            fs::path(
               saveDir,
               paste0(
                  fs::path_ext_remove(rawFileName),
                  "_specZoom_MS2.pdf"
               )
            ),
         plot = spectra_MS2_marrange,
         width = 20,
         height = 12,
         limitsize = FALSE
      )


   }
