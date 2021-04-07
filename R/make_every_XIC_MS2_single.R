#' make_every_XIC_MS2_single
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
#' @param fragment_charges
#' @param fragment_mz_range
#' @param fragment_pos_cutoff
#' @param XIC_tol_MS2
#' @param use_IAA
#' @param abund_cutoff
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples

make_every_XIC_MS2_single <-
   function(
      rawFileDir = NULL,
      rawFileName = NULL,
      targetSeqData = NULL,
      outputDir = getwd(),
      target_col_name = "UNIPROTKB",
      target_sequence_col_name = "ProteoformSequence",
      PTMname_col_name = "PTMname",
      PTMformula_col_name1 = "FormulaToAdd",
      PTMformula_col_name2 = "FormulaToSubtract",
      isoAbund = c("12C" = 0.9893, "14N" = 0.99636),
      fragment_charges = c(1:50),
      fragment_types = c("b", "y"),
      fragment_mz_range = c(300,2000),
      fragment_pos_cutoff = c(1, 50),
      XIC_tol_MS2 = 10,
      XIC_cutoff = 0.0005,
      scoreMFAcutoff = 0.3,
      cosinesimcutoff = 0.99,
      SNcutoff = 20,
      resPowerMS2 = 150000,
      isotopologue_window_multiplier = 6,
      mz_window = 5,
      use_IAA = FALSE,
      abund_cutoff = 5,
      rawrrTemp = tempdir()
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

      # Fix this later lmao

      # assertthat::assert_that(
      #    assertthat::has_extension(targetSeqData, "csv") |
      #       is.data.frame(targetSeqData),
      #    msg = "targetSeqData is not a CSV file or dataframe"
      # )

      # if (
      #    assertthat::see_if(assertthat::has_extension(targetSeqData, "csv"))
      #    == TRUE
      # ) {
      #
      #    assertthat::assert_that(
      #       assertthat::is.readable(targetSeqData),
      #       msg = "targetSeqData is not a readable file"
      #    )
      #
      # }

      assertthat::assert_that(
         assertthat::is.dir(outputDir),
         msg = "outputDir is not a recognized directory"
      )

      assertthat::assert_that(
         assertthat::is.flag(use_IAA),
         msg = "use_IAA should be TRUE or FALSE"
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

      # Get positions of specific isotopes in dataframe

      indices <-
         list(
            "12C" = which(isotopes$isotope == "12C") %>% .[[1]],
            "13C" = which(isotopes$isotope == "13C") %>% .[[1]],
            "14N" = which(isotopes$isotope == "14N") %>% .[[1]],
            "15N" = which(isotopes$isotope == "15N") %>% .[[1]]
         )

      # Replace values in dataframe with supplied values

      isotopes$abundance[indices[["12C"]]] <-
         isoAbund[which(names(isoAbund) == "12C")]

      isotopes$abundance[indices[["13C"]]] <-
         1 - isoAbund[which(names(isoAbund) == "12C")]

      isotopes$abundance[indices[["14N"]]] <-
         isoAbund[which(names(isoAbund) == "14N")]

      isotopes$abundance[indices[["15N"]]] <-
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

      # Create save directory ---------------------------------------------------

      systime <- format(Sys.time(), "%Y%m%d_%H%M")
      systime2 <- Sys.time()
      rand_string <- stringi::stri_rand_strings(1, 3, '[A-Z]')

      saveDir <-
         fs::path(
            outputDir,
            paste0(
               systime,
               "_",
               fs::path_ext_remove(rawFileName),
               "_",
               length(target_seqs),
               "seqs_MS2"
            )
         ) %>%
         stringr::str_trunc(246, "right", ellipsis = "")


      if (dir.exists(saveDir) == FALSE) dir.create(saveDir)

      ## Save target seqs

      message(
         paste0(
            "\nSaving target sequences to ",
            saveDir,
            "/",
            fs::path_ext_remove(
               stringr::str_trunc(
                  rawFileName, 90, "right", ellipsis = ""
               )
            ),
            '_target_seqs.csv'
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
                  '_target_seqs_MS2_',
                  rand_string
               ),
               ext = "csv"
            )
      )

      ## Prepare summary datasheet

      spectra_MS2_dataframe_union <-
         tibble::tibble()


      # Get elemental compositions ----------------------------------------------

      for (i in seq_along(target_seqs)) {

         # Generate all possible fragments

         fragments <-
            purrr::map2(
               target_seqs[[i]],
               list(fragment_charges),
               ~purrr::map2(
                  .x,
                  .y,
                  ~MSnbase::calculateFragments(
                     .x,
                     z = .y,
                     type = fragment_types,
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
                     ~OrgMassSpecR::ConvertPeptide(
                        .x,
                        IAA = use_IAA
                     ) %>%
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
                  chemform = .y,
                  ion_chemform = paste(ion, chemform, sep = "_")
               )
            )

         # Prepare fragments data for combination with iso_dist_MS2 later

         # fragments3_old <-
         #    fragments2 %>%
         #    purrr::map(
         #       ~dplyr::group_by(.x, type, pos, mz) %>%
         #          dplyr::group_split()
         #    )

         fragments3 <-
            fragments2 %>%
            purrr::map(
               ~dplyr::group_by(.x, type, pos, mz) %>%
                  dplyr::group_split()
            )

         chemform_names <-
            purrr::map_depth(
               fragments3,
               2,
               ~paste(dplyr::pull(.x, ion), dplyr::pull(.x, chemform), sep = "_")
            )

         # Need to reorder fragments3 to match the order of fragments2/iso_dist_MS2

         fragments3 <-
            purrr::map2(
               fragments3,
               chemform_names,
               ~purrr::set_names(.x, .y)
            ) %>%
            purrr::map2(
               fragments2,
               ~.x[.y$ion_chemform]
            )

         # Generate isotopic distributions ---------

         iso_dist_MS2 <-
            purrr::map(
               fragments2,
               ~enviPat::isopattern(
                  isotopes,
                  chemform=.x$chemform,
                  threshold=0.1,
                  plotit=FALSE,
                  charge=.x$z,
                  emass=0.00054858,
                  algo=1,
                  verbose = F
               ) %>%
                  enviPat::envelope(
                     dmz = "get",
                     resolution = resPowerMS2,
                     verbose = F
                  )
            ) %>%
            purrr::modify_depth(
               2,
               ~new(
                  "Spectrum1",
                  mz = .x[,1],
                  intensity = .x[,2],
                  centroided = FALSE
               ) %>%
                  MSnbase::pickPeaks(
                     SNR = 1,
                     method = "MAD",
                     refineMz = "kNeighbors",
                     k = 2
                  )
            ) %>%
            purrr::modify_depth(
               2,
               ~tibble::tibble(
                  `m/z` = MSnbase::mz(.x),
                  abundance = MSnbase::intensity(.x)
               ) %>%
                  dplyr::filter(abundance > abund_cutoff)
            )


         iso_dist_MS2_cluster <-
            purrr::map2(
               iso_dist_MS2,
               fragments3,
               ~purrr::map2(
                  .x = .x,
                  .y = .y,
                  ~dplyr::mutate(
                     .x,
                     charge = .y$z,
                     ion = .y$ion
                  )
               )
            ) %>%
            purrr::map(
               dplyr::bind_rows
            )



         # Cluster isotopic distributions

         # fragment_groups_old <-
         #    purrr::map(
         #       iso_dist_MS2_cluster,
         #       ~dplyr::group_by(.x, ion) %>%
         #          dplyr::group_split()
         #    )

         fragment_groups <-
            purrr::map(
               iso_dist_MS2_cluster,
               ~dplyr::group_by(.x, ion) %>%
                  dplyr::group_split() %>%
                  purrr::map(
                     ~dplyr::group_by(.x, charge) %>%
                        dplyr::filter(abundance == max(abundance))
                  )
            )

         # Read XICs --------

         # Get max TIC for MS2s from rawfile

         max_TIC_MS2 <-
            rawrr::readChromatogram(
               rawfile = rawFile,
               filter = "ms2",
               type = "tic",
               tol = XIC_tol_MS2
            ) %>%
            .$intensities %>%
            max()

         XIC_TIC_cutoff <-
            max_TIC_MS2 * XIC_cutoff

         XIC_MS2 <-
            purrr::map(
               fragment_groups,
               ~purrr::map(
                  .x,
                  ~rawrr::readChromatogram(
                     rawfile = rawFile,
                     mass = .x$`m/z`,
                     filter = "ms2",
                     type = "xic",
                     tol = XIC_tol_MS2
                  )
               )
            )

         # Make blank tibble to replace NULL tibbles

         blank_tibble <-
            tibble::tibble(times = 0, intensities = 0)

         # Process XICs by replacing NULL tibbles and adding together all XICs
         # for same ion and charge

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

         # fragments_to_retain <-
         #    XIC_MS2_sum %>%
         #    purrr::map_depth(
         #       2,
         #       ~dplyr::filter(.x, int_sum == max(int_sum)) %>%
         #          dplyr::pull(int_sum)
         #    ) %>%
         #    purrr::map(
         #       unlist
         #    ) %>%
         #    purrr::map(
         #       ~which(.x > XIC_TIC_cutoff)
         #    )

         fragments_to_retain <-
            XIC_MS2_sum %>%
            purrr::map_depth(
               2,
               ~dplyr::filter(.x, int_sum == max(int_sum)) %>%
                  dplyr::pull(int_sum) %>%
                  dplyr::first()
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
                  unique()
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
                  as.list(names(target_seqs)[[i]]),
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

         message(glue::glue("\nSaving XICs for {names(target_seqs)[[i]]}"))
         message(paste("Memory used:", pryr::mem_used()/1000000, "MB"))


         XIC_plots_MS2_marrange <-
            gridExtra::marrangeGrob(
               grobs = purrr::flatten(XIC_plots_MS2),
               ncol = 5,
               nrow = 3,
               top = rawFileName
            )

         # Save XIC results ------

         saveDirSub <-
            fs::path(
               saveDir,
               "XIC_MS2"
            )

         if (fs::dir_exists(saveDirSub) == FALSE) fs::dir_create(saveDirSub)

         ## Save theo iso distributions -----

         writexl::write_xlsx(
            iso_dist_MS2_cluster,
            path =
               fs::path(
                  saveDirSub,
                  paste0(
                     stringr::str_trunc(
                        names(target_seqs)[[i]], 90, "right", ellipsis = ""),
                     '_isodist_MS2.xlsx'
                  )
               )
         )

         ## Save XIC plots -----

         ggplot2::ggsave(
            filename =
               fs::path(
                  saveDirSub,
                  paste0(
                     stringr::str_trunc(
                        names(target_seqs)[[i]], 90, "right", ellipsis = ""
                     ),
                     '_MS2_XICs.pdf'
                  )
               ),
            plot = XIC_plots_MS2_marrange,
            width = 20,
            height = 12,
         )

         # XIC_plots_MS2_marrange is large and no longer needed

         rm(XIC_plots_MS2_marrange)




         # Extract scans and make spectra ------------------------------------------

         rawFileMetadata <-
            rawrr::readIndex(
               rawFile,
               tmpdir = rawrrTemp
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
                     dplyr::group_by(.y, charge) %>%
                        dplyr::filter(abundance == max(abundance))
                  )
               )
            ) %>%
            purrr::map_depth(
               2,
               ~dplyr::group_by(.x, charge) %>%
                  dplyr::group_split()
            )

         # mz_to_read is really the m/z of the theoretical most abundant
         # fragment ion peak

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
            ) %>%
            purrr::map_depth(
               2,
               ~magrittr::extract2(.x, 1)
            )

         # Read the scan with max TIC for each fragment and subset it to the
         # desired range


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
                  ~rawrr::readSpectrum(
                     rawfile = rawFile,
                     scan = ..1,
                     tmpdir = rawrrTemp
                  ) %>%
                     purrr::flatten() %>%
                     {
                        tibble::tibble(
                           mz = .$mZ,
                           intensity = .$intensity
                        )
                     } %>%
                     list() %>%
                     purrr::map2(
                        .y = ..2,
                        ~dplyr::filter(
                           .x,
                           mz >= .y - mz_window/2,
                           mz <= .y + mz_window/2
                        )
                     )
               )
            )

         # Get noise for every scan being plotted to use for S/N calculation

         noise_obs_MS2 <-
            purrr::pmap(
               list(
                  scan_to_read
               ),
               ~purrr::pmap(
                  list(
                     ..1
                  ),
                  ~purrr::pmap(
                     list(
                        ..1
                     ),
                     ~rawrr::readSpectrum(
                        rawfile = rawFile,
                        scan = ..1,
                        tmpdir = rawrrTemp
                     ) %>%
                        purrr::flatten() %>%
                        {
                           MALDIquant::createMassSpectrum(
                              mass = .[["mZ"]],
                              intensity = .[["intensity"]]
                           ) %>%
                              MALDIquant::estimateNoise(method = "MAD") %>%
                              .[,2] %>%
                              .[[1]]
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
               iso_dist_MS2_cluster,
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


         # Extract maxY and maxX from within isotopologue windows to use for
         # cosine sim and scoreMFA scoring


         mz_obs_MS2 <-
            purrr::pmap(
               list(
                  scans_to_plot,
                  iso_dist_cluster_trunc
               ),
               ~purrr::pmap(
                  list(
                     ..1 = ..1,
                     ..2 = ..2
                  ),
                  ~purrr::pmap(
                     list(
                        ..1 = ..1,
                        ..2 = ..2
                     ),
                     ~get_maxX_in_Xrange_vector(
                        df = ..1,
                        x = mz,
                        y = intensity,
                        mz = ..2[["m/z"]],
                        res_power = resPowerMS2,
                        isotopologueWinMultiplier = isotopologue_window_multiplier
                     )
                  )
               )
            )


         intensity_obs_MS2 <-
            purrr::pmap(
               list(
                  scans_to_plot,
                  iso_dist_cluster_trunc
               ),
               ~purrr::pmap(
                  list(
                     ..1 = ..1,
                     ..2 = ..2
                  ),
                  ~purrr::pmap(
                     list(
                        ..1 = ..1,
                        ..2 = ..2
                     ),
                     ~get_maxY_in_Xrange_vector(
                        df = ..1,
                        x = mz,
                        y = intensity,
                        mz = ..2[["m/z"]],
                        res_power = resPowerMS2,
                        isotopologueWinMultiplier = isotopologue_window_multiplier
                     )
                  )
               )
            )


         mz_theo_MS2 <-
            purrr::map_depth(
               iso_dist_cluster_trunc,
               3,
               ~dplyr::pull(.x, `m/z`)
            )

         intensity_theo_MS2 <-
            purrr::map_depth(
               iso_dist_cluster_trunc,
               3,
               ~dplyr::pull(.x, abundance)
            )

         # Calculate scaling factor for intensity_theo_MS2

         intensity_theo_scaling_factors <-
            purrr::map2(
               intensity_obs_MS2,
               intensity_theo_MS2,
               ~purrr::map2(
                  .x = .x,
                  .y = .y,
                  ~purrr::map2(
                     .x = .x,
                     .y = .y,
                     ~max(.x)/max(.y)
                  )
               )
            )

         intensity_theo_MS2_scaled <-
            purrr::map2(
               intensity_theo_MS2,
               intensity_theo_scaling_factors,
               ~purrr::map2(
                  .x = .x,
                  .y = .y,
                  ~purrr::map2(
                     .x = .x,
                     .y = .y,
                     ~.x*.y
                  )
               )
            )

         # Calculate S/N for observed and theoretical intensities


         SN_estimate_obs_MS2 <-
            purrr::map2(
               intensity_obs_MS2,
               noise_obs_MS2,
               ~purrr::map2(
                  .x,
                  .y,
                  ~purrr::map2(
                     .x,
                     .y,
                     ~(.x/.y) %>%
                        {if (length(.) == 0 | is.nan(.) | is.null(.) | is.infinite(.)) 0 else .}
                  )
               )
            ) %>%
            purrr::map(
               ~purrr::map_if(
                  .x,
                  ~length(.x) == 0,
                  ~0
               )
            )


         SN_estimate_theo_MS2 <-
            purrr::map2(
               intensity_theo_MS2_scaled,
               noise_obs_MS2,
               ~purrr::map2(
                  .x,
                  .y,
                  ~purrr::map2(
                     .x,
                     .y,
                     ~(.x/.y) %>%
                        {if (length(.) == 0 | is.nan(.) | is.null(.) | is.infinite(.)) 0 else .}
                  )
               )
            ) %>%
            purrr::map(
               ~purrr::map_if(
                  .x,
                  ~length(.x) == 0,
                  ~0
               )
            )



         # Attempt to determine RP for each isotopologue peak

         rp_obs_MS2 <-
            purrr::pmap(
               list(
                  scans_to_plot,
                  iso_dist_cluster_trunc
               ),
               ~purrr::pmap(
                  list(
                     ..1 = ..1,
                     ..2 = ..2
                  ),
                  ~purrr::pmap(
                     list(
                        ..1 = ..1,
                        ..2 = ..2
                     ),
                     ~get_rp_in_Xrange_vector(
                        df = ..1,
                        x = mz,
                        y = intensity,
                        mz = ..2[["m/z"]],
                        res_power = resPowerMS2,
                        isotopologueWinMultiplier = isotopologue_window_multiplier
                     )
                  )
               )
            )

         # Convert all resolving powers into a data frame for linear modeling

         rp_obs_MS2_df <-
            purrr::map2(
               rp_obs_MS2,
               mz_to_read,
               ~purrr::map2(
                  .x,
                  .y,
                  ~tibble::tibble(
                     rp = .x,
                     mz = .y
                  ) %>%
                     tidyr::unnest(cols = c(rp, mz)) %>%
                     dplyr::filter(rp != 0)
               )
            ) %>%
            dplyr::bind_rows() %>%
            dplyr::arrange(mz)


         # Calculate linear model for RP based on determined non-zero values

         x <- rp_obs_MS2_df$mz
         y <- rp_obs_MS2_df$rp

         resolving_power_model_MS2 <-
            fit <- lm(y ~ I(1/x))

         # Calculate theoretical resolving powers for every isotopologue peak

         rp_theo_MS2 <-
            purrr::map_depth(
               mz_theo_MS2,
               3,
               ~predict(
                  resolving_power_model_MS2,
                  newdata =
                     data.frame(
                        x = .x
                     )
               )
            )

         # Replace all zeros in rp_obs_MS1 with theoretical values from the fit

         rp_obs_theo_hybrid_MS2 <-
            purrr::map2(
               rp_obs_MS2,
               rp_theo_MS2,
               ~purrr::map2(
                  .x,
                  .y,
                  ~purrr::map2(
                     .x,
                     .y,
                     ~purrr::map2_dbl(
                        .x,
                        .y,
                        ~{if (.x == 0) .y else .x}
                     )
                  )
               )
            )



         # Calculate ScoreMFA

         score_MFA_MS2 <-
            purrr::pmap(
               list(
                  mz_obs_MS2,
                  mz_theo_MS2,
                  rp_obs_theo_hybrid_MS2,
                  SN_estimate_obs_MS2,
                  SN_estimate_theo_MS2
               ),
               ~purrr::pmap(
                  list(
                     ..1,
                     ..2,
                     ..3,
                     ..4,
                     ..5
                  ),
                  ~purrr::pmap(
                     list(
                        ..1,
                        ..2,
                        ..3,
                        ..4,
                        ..5
                     ),
                     ~ScoreMFA(
                        ..1,
                        ..2,
                        ..3,
                        ..4,
                        ..5,
                        rp_mult = 2
                     ) %>%
                        {if (is.nan(.) | is.null(.) | length(.) == 0) 0 else .}
                  )
               )
            )

         # Calculate cosine similarity

         cosine_sim_MS2 <-
            purrr::pmap(
               list(
                  intensity_obs_MS2,
                  intensity_theo_MS2_scaled
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
                     ~calculate_cosine_similarity(
                        ..1,
                        ..2
                     ) %>%
                        {if (is.nan(.) | is.null(.) | length(.) == 0) 0 else .}
                  )
               )
            )


         # Prepare data for plotting

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

         scan_to_read_special <-
            purrr::map2(
               scan_to_read,
               ion_names_MS2,
               ~purrr::map2(
                  .x,
                  .y,
                  ~purrr::map2(
                     .x,
                     .y,
                     ~rep(.x, times = length(.y))
                  )
               )
            )

         spectra_MS2 <-
            purrr::pmap(
               list(
                  scans_to_plot,
                  as.list(names(target_seqs)[[i]]),
                  mz_to_read,
                  scan_to_read_special,
                  ion_names_MS2,
                  ion_charges_MS2,
                  mz_to_read,
                  score_MFA_MS2,
                  mz_theo_MS2,
                  intensity_theo_MS2_scaled,
                  cosine_sim_MS2,
                  SN_estimate_obs_MS2
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
                     ..8,
                     ..9,
                     ..10,
                     ..11,
                     ..12
                  ),
                  ~purrr::pmap(
                     list(
                        ..1[..8 > scoreMFAcutoff & ..11 > cosinesimcutoff],                   # By subsetting this way, spectra
                        rep(..2, length(..3))[..8 > scoreMFAcutoff & ..11 > cosinesimcutoff], # with low scoreMFAs are removed
                        ..3[..8 > scoreMFAcutoff & ..11 > cosinesimcutoff],
                        ..4[..8 > scoreMFAcutoff & ..11 > cosinesimcutoff],
                        ..5[..8 > scoreMFAcutoff & ..11 > cosinesimcutoff],
                        ..6[..8 > scoreMFAcutoff & ..11 > cosinesimcutoff],
                        ..7[..8 > scoreMFAcutoff & ..11 > cosinesimcutoff],
                        ..8[..8 > scoreMFAcutoff & ..11 > cosinesimcutoff],
                        ..9[..8 > scoreMFAcutoff & ..11 > cosinesimcutoff],
                        ..10[..8 > scoreMFAcutoff & ..11 > cosinesimcutoff],
                        ..11[..8 > scoreMFAcutoff & ..11 > cosinesimcutoff],
                        ..12[..8 > scoreMFAcutoff & ..11 > cosinesimcutoff]
                     ),
                     ~make_spectrum_MS2(
                        ..1,
                        name = ..2,
                        max_mz = ..3,
                        scan = ..4,
                        ion = ..5,
                        charge = ..6,
                        isotopologue_window = isotopologue_window_multiplier*(..7/resPowerMS2),
                        ScoreMFA = ..8,
                        cosine_sim = ..11,
                        xrange = c(..3 - (mz_window/2), ..3 + (mz_window/2)),
                        mz_theo = ..9,
                        int_theo = ..10,
                        sn_estimate = ..12,
                        theme = MStheme01
                     )
                  )
               )
            ) %>%
            purrr::map_depth(
               2,
               ~ggplot_list_checker(.x)
            ) %>%
            purrr::map(
               purrr::compact
            )


         # Put all of the info together into a data frame

         spectra_MS2_dataframe <-
            purrr::pmap(
               list(
                  scan_to_read_special,
                  mz_obs_MS2,
                  mz_theo_MS2,
                  ion_names_MS2,
                  ion_charges_MS2,
                  score_MFA_MS2,
                  cosine_sim_MS2,
                  SN_estimate_obs_MS2
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
                        ..1[..6 > scoreMFAcutoff & ..7 > cosinesimcutoff], # By subsetting this way, spectra
                        ..2[..6 > scoreMFAcutoff & ..7 > cosinesimcutoff], # with low scoreMFAs are removed
                        ..3[..6 > scoreMFAcutoff & ..7 > cosinesimcutoff],
                        ..4[..6 > scoreMFAcutoff & ..7 > cosinesimcutoff],
                        ..5[..6 > scoreMFAcutoff & ..7 > cosinesimcutoff],
                        ..6[..6 > scoreMFAcutoff & ..7 > cosinesimcutoff],
                        ..7[..6 > scoreMFAcutoff & ..7 > cosinesimcutoff],
                        ..8[..6 > scoreMFAcutoff & ..7 > cosinesimcutoff]
                     ),
                     ~tibble::tibble(
                        raw_filename = rawFileName,
                        sequence_name = names(target_seqs)[[i]],
                        scan = ..1,
                        max_obs_mz = max(..2),
                        max_theo_mz = max(..3),
                        ppm_error = calculate_mma_ppm(max(..2), max(..3)),
                        ion = ..4,
                        charge = ..5,
                        scoreMFA = ..6,
                        cosineSim = ..7,
                        mean_SN_estimate = mean(..8)
                     )
                  )
               )
            ) %>%
            purrr::map(
               purrr::compact
            ) %>%
            purrr::map(
               dplyr::bind_rows
            ) %>%
            dplyr::bind_rows()


         # Add data frame to summary data frame of assigned fragments

         spectra_MS2_dataframe_union <-
            dplyr::bind_rows(
               spectra_MS2_dataframe_union,
               spectra_MS2_dataframe
            )

         message(glue::glue("\nSaving spectra for {names(target_seqs)[[i]]}"))
         message(paste("Memory used:", pryr::mem_used()/1000000, "MB"))

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
                  saveDirSub,
                  paste0(
                     stringr::str_trunc(
                        names(target_seqs)[[i]], 90, "right", ellipsis = ""
                     ),
                     '_specZoom_MS2.pdf'
                  )
               ),
            plot = spectra_MS2_marrange,
            width = 20,
            height = 12,
            limitsize = FALSE
         )

      }

      # Save data frame of fragments which beat ScoreMFA cutoff

      message("\nWriting all 'assigned' fragment ions to spreadsheet")

      sum_path <-
         fs::path(
            saveDir,
            paste0(
               fs::path_ext_remove(
                  stringr::str_trunc(
                     rawFileName, 90, "right", ellipsis = ""
                  )
               ),
               '_assigned_fragments_MS2.xlsx'
            )
         )

      if (fs::file_exists(sum_path)) {

         sum_path <-
            fs::path(
               saveDir,
               paste0(
                  fs::path_ext_remove(
                     stringr::str_trunc(
                        rawFileName, 90, "right", ellipsis = ""
                     )
                  ),
                  '_assigned_fragments_MS2_',
                  rand_string
               ),
               ext = "xlsx"
            )

      }

      writexl::write_xlsx(
         spectra_MS2_dataframe_union,
         path = sum_path
      )

      # Write params to text file

      params_path <-
         fs::path(
            saveDir,
            'GEX_params.txt'
         )

      if (!fs::file_exists(params_path)) {

         readr::write_lines(
            glue::glue(
               "
               ----------------
               {systime2}
               ----------------
               rawFileDir = {toString(rawFileDir)}
               rawFileName = {toString(rawFileName)}
               targetSeqData = {toString(targetSeqData)}
               outputDir = {toString(outputDir)}
               target_col_name = {rlang::as_name(target_col_name)}
               target_sequence_col_name = {rlang::as_name(target_sequence_col_name)}
               PTMname_col_name = {rlang::as_name(PTMname_col_name)}
               PTMformula_col_name1 = {rlang::as_name(PTMformula_col_name1)}
               PTMformula_col_name2 = {rlang::as_name(PTMformula_col_name2)}
               isoNames = {toString(names(isoAbund))}
               isoAbund = {toString(isoAbund)}
               fragment_charges = {toString(fragment_charges)}
               fragment_types = {toString(fragment_types)}
               fragment_mz_range = {toString(fragment_mz_range)}
               fragment_pos_cutoff = {toString(fragment_pos_cutoff)}
               XIC_tol_MS2 = {XIC_tol_MS2}
               XIC_cutoff = {XIC_cutoff}
               scoreMFAcutoff = {scoreMFAcutoff}
               cosinesimcutoff = {cosinesimcutoff}
               resPowerMS2 = {resPowerMS2}
               isotopologue_window_multiplier = {isotopologue_window_multiplier}
               mz_window = {mz_window}
               use_IAA = {use_IAA}
               abund_cutoff = {abund_cutoff}
               "
            ),
            file = params_path
         )

      }

      # Trying to deal with problem with futures not freeing up
      # memory correctly!

      rm(list = ls(all.names = TRUE))
      gc()


   }


