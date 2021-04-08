#' make_every_XIC_MS1
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
#' @param sample_n_pforms
#' @param mass_range
#' @param target_charges
#' @param mz_range
#' @param XIC_tol
#' @param use_IAA
#' @param save_output
#' @param abund_cutoff
#'
#' @return
#' @export
#'
#' @import MSnbase
#' @importFrom magrittr %>%
#'
#' @examples

make_every_XIC_MS1 <-
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
      target_charges = c(1:50),
      mass_range = c(0,100000),
      mz_range = c(600,2000),
      abund_cutoff = 15,
      sample_n_pforms = NULL,
      XIC_tol = 2,
      use_IAA = FALSE,
      save_output = TRUE,
      scoreMFAcutoff = 0.3,
      cosinesimcutoff = 0.8,
      SN_cutoff = 10,
      resPowerMS1 = 300000,
      isotopologue_window_multiplier = 6,
      mz_window = 3,
      return_timers = TRUE,
      save_spec_object = FALSE,
      rawrrTemp = tempdir()
   ) {

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

      if (is.null(sample_n_pforms) == FALSE && sample_n_pforms > 0) {

         assertthat::assert_that(
            assertthat::is.count(sample_n_pforms),
            msg = "sample_n_pforms should be a positive integer"
         )

      }

      # Create necessary quosures -----------------------------------------------

      target_col_name <- rlang::enquo(target_col_name)
      target_sequence_col_name_sym <- rlang::sym(target_sequence_col_name)
      target_sequence_col_name <- rlang::enquo(target_sequence_col_name)
      PTMformula_col_name1 <- rlang::enquo(PTMformula_col_name1)
      PTMformula_col_name2 <- rlang::enquo(PTMformula_col_name2)
      PTMname_col_name <- rlang::enquo(PTMname_col_name)

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
            ggplot2::labs(
               x = "Retention Time (min)",
               y = "Total Ion Current"
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

      indices <-
         list(
            "12C" = which(isotopes$isotope == "12C") %>% .[[1]],
            "13C" = which(isotopes$isotope == "13C") %>% .[[1]],
            "14N" = which(isotopes$isotope == "14N") %>% .[[1]],
            "15N" = which(isotopes$isotope == "15N") %>% .[[1]]
         )

      # Replace values in dataframe with supplied values

      isotopes_to_use$abundance[indices[["12C"]]] <-
         isoAbund[which(names(isoAbund) == "12C")]

      isotopes_to_use$abundance[indices[["13C"]]] <-
         1 - isoAbund[which(names(isoAbund) == "12C")]

      isotopes_to_use$abundance[indices[["14N"]]] <-
         isoAbund[which(names(isoAbund) == "14N")]

      isotopes_to_use$abundance[indices[["15N"]]] <-
         1 - isoAbund[which(names(isoAbund) == "14N")]

      # Process target sequences ------------------------------------------------

      target_seqs_df <-
         targetSeqData %>%
         {if (!is.data.frame(.)) readr::read_csv(.) else .} %>%
         dplyr::mutate(
            MonoisoMass =
               Peptides::mw(!!target_sequence_col_name_sym, monoisotopic = TRUE)
         ) %>%
         dplyr::filter(
            MonoisoMass >= mass_range[[1]],
            MonoisoMass <= mass_range[[2]]
         ) %>%
         {if (!is.null(sample_n_pforms)) dplyr::sample_n(., size = sample_n_pforms) else .}

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

      saveDir <-
         fs::path(
            outputDir,
            paste0(
               systime,
               "_",
               fs::path_ext_remove(rawFileName),
               "_",
               length(target_seqs),
               "seqs"
            )
         ) %>%
         stringr::str_trunc(246, "right", ellipsis = "")


      if (fs::dir_exists(saveDir) == FALSE) fs::dir_create(saveDir)

      timer <-
         timeR::createTimer(verbose = FALSE)


      # Get elemental composition -----------------------------------------------

      timer$start("Chemical formulas and isotopic distributions")

      chemform_temp <-
         purrr::map(
            target_seqs,
            OrgMassSpecR::ConvertPeptide,
            IAA = use_IAA
         )

      chemform_names <-
         chemform_temp %>%
         purrr::map(unlist) %>%
         purrr::map(names)

      chemform <-
         chemform_temp %>%
         purrr::map(unlist) %>%
         purrr::map2(
            chemform_names,
            ~purrr::map2_chr(.y, .x, paste0)
         ) %>%
         purrr::map(
            paste0,
            collapse = ""
         )

      ## Extract PTM names and make list for later use

      PTM_names_list <-
         target_seqs_df %>%
         dplyr::select(!!target_col_name, !!PTMname_col_name) %>%
         tibble::deframe() %>%
         as.list()

      PTM_names_list[is.na(PTM_names_list) == TRUE] <- "None"

      ## Make list of PTM formulas to ADD

      PTM_form_to_add <-
         target_seqs_df %>%
         dplyr::pull(!!PTMformula_col_name1) %>%
         as.list()

      PTM_form_to_add[is.na(PTM_form_to_add) == TRUE] <- "C0"

      PTM_form_to_add <-
         purrr::map(
            PTM_form_to_add,
            ~enviPat::check_chemform(isotopes = isotopes_to_use, chemforms = .x) %>%
               dplyr::pull(new_formula)
         )

      ## Make list of PTM formulas to SUBTRACT

      PTM_form_to_sub <-
         target_seqs_df %>%
         dplyr::pull(!!PTMformula_col_name2) %>%
         as.list()

      PTM_form_to_sub[is.na(PTM_form_to_sub) == TRUE] <- "C0"

      PTM_form_to_sub <-
         purrr::map(
            PTM_form_to_sub,
            ~enviPat::check_chemform(isotopes = isotopes_to_use, chemforms = .x) %>%
               dplyr::pull(new_formula)
         )

      ## Add and subtract PTM formulas as needed

      chemform <-
         purrr::map2(
            chemform,
            PTM_form_to_add,
            ~enviPat::mergeform(.x, .y)
         ) %>%
         purrr::map2(
            PTM_form_to_sub,
            ~enviPat::subform(.x, .y)
         )

      # Add one H per positive charge

      H_to_merge <-
         target_charges %>%
         purrr::map_chr(~paste0("H", .x)) %>%
         list() %>%
         rep(length(target_seqs))

      chemform_to_merge <-
         chemform %>%
         purrr::map(
            ~rep(.x, length(target_charges))
         )

      chemform_withH <-
         purrr::map2(
            chemform_to_merge,
            H_to_merge,
            ~purrr::map2_chr(.x, .y, ~enviPat::mergeform(.x, .y))
         )

      # Remove any elements with 0 count, they cause enviPat::isopattern to fail

      for (i in seq_along(chemform_withH)) {

         if (
            any(stringr::str_detect(chemform_withH[[i]], "C0|H0|N0|O0|P0|S0")) == TRUE
         ) {
            chemform_withH[[i]] <-
               stringr::str_remove_all(chemform_withH[[i]], "C0|H0|N0|O0|P0|S0")
         }
      }

      # Get rid of some unneeded junk

      rm(chemform_temp)
      rm(chemform_names)
      rm(chemform_to_merge)
      rm(H_to_merge)

      # Calculate isotopic dist -------------------------------------------------

      purrr::map(
         chemform,
         ~enviPat::check_chemform(isotopes = isotopes_to_use, chemforms = .x) %>%
            dplyr::pull(warning) %>%
            any() %>%
            `if`(., stop("Problem with chemical formula"))
      )

      iso_dist <-
         purrr::map2(
            chemform_withH,
            list(target_charges) %>% rep(length(target_seqs)),
            ~purrr::map2(
               .x,
               .y,
               ~enviPat::isopattern(
                  isotopes_to_use,
                  chemform=.x,
                  threshold=0.1,
                  plotit=FALSE,
                  charge=.y,
                  emass=0.00054858,
                  algo=1,
                  verbose = F
               ) %>%
                  enviPat::envelope(
                     dmz = "get",
                     resolution = resPowerMS1,
                     verbose = F
                  )
            )
         ) %>%
         purrr::modify_depth(
            3,
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
         )


      target_charges_list <-
         list(target_charges) %>%
         rep(length(target_seqs))

      iso_dist_list_union <-
         iso_dist %>%
         purrr::modify_depth(
            3,
            ~tibble::tibble(
               `m/z` = MSnbase::mz(.x),
               abundance = MSnbase::intensity(.x)
            )
         ) %>%
         purrr::map2(
            .y = target_charges_list,
            ~purrr::map2(
               .x = .x,
               .y = .y,
               ~purrr::map2(
                  .x = .x,
                  .y = .y,
                  ~dplyr::mutate(.x, charge = .y) %>%
                     dplyr::filter(
                        `m/z` > mz_range[[1]] & `m/z` < mz_range[[2]],
                        abundance > abund_cutoff
                     )
               )
            )
         )

      # Remove iso_dist to save space hopefully

      rm(iso_dist)

      # Keep the next three chunks separate, doesn't work if they are condensed

      iso_dist_list_union <-
         purrr::map(
            iso_dist_list_union,
            purrr::reduce,
            dplyr::union_all
         )

      # Save this for later use in making spectra

      iso_dist_vlines <-
         iso_dist_list_union


      iso_dist_list_union <-
         purrr::map(
            iso_dist_list_union,
            ~tibble_list_checker(.x)
         )


      iso_dist_list_union <-
         purrr::map(
            iso_dist_list_union,
            purrr::reduce,
            dplyr::union_all
         )

      timer$stop("Chemical formulas and isotopic distributions")

      # Calculate XIC -----------------------------------------------------------

      timer$start("Calculate and plot XICs")

      # Not sure whether all isotopologue peaks should be included in generating
      # this XIC. Use most abundant peak for every charge state for now

      XIC_target_mz <-
         iso_dist_list_union %>%
         purrr::map(~dplyr::group_by(.x, charge) %>%
                       dplyr::filter(abundance == max(abundance)) %>%
                       dplyr::pull(`m/z`)) %>%
         purrr::map(unique)

      XIC <-
         purrr::map(
            XIC_target_mz,
            ~rawrr::readChromatogram(
               rawfile = rawFile,
               mass = .x,
               filter = "ms",
               type = "xic",
               tol = XIC_tol
            )
         )

      XIC_nonull <-
         XIC %>%
         purrr::map(kickoutXIC)

      # SAVE SPACE BY DELETING XICS!

      rm(XIC)

      sumXIC1 <-
         XIC_nonull %>%
         purrr::modify_depth(
            2,
            ~tibble::tibble(
               times = .x$times,
               intensities = .x$intensities
            )
         ) %>%
         purrr::modify_depth(2, ~dplyr::select(.x, times, intensities)) %>%
         purrr::map(purrr::reduce, dplyr::full_join)

      sumXIC2 <-
         sumXIC1 %>%
         purrr::map(~dplyr::group_by(.x, times)) %>%
         purrr::map(~dplyr::summarize(.x, int_sum = sum(intensities)))

      # Get RT corresponding to maximum TIC for each sequence for later use in
      # making spectra

      RT_of_maxTIC <-
         sumXIC2 %>%
         purrr::map(
            ~dplyr::filter(.x, int_sum == max(int_sum))
         ) %>%
         purrr::map(
            ~dplyr::pull(.x, times)
         )

      rawFileMetadata <-
         rawrr::readIndex(
            rawFile,
            tmpdir = rawrrTemp
         ) %>%
         tibble::as_tibble() %>%
         dplyr::mutate(RT_min = rtinseconds/60)

      scanNumber_and_RT <-
         rawFileMetadata %>%
         dplyr::select(scan, RT_min)

      scanNumber_and_RT_vec <-
         scanNumber_and_RT$RT_min %>%
         purrr::set_names(scanNumber_and_RT$scan)


      scanNumsToRead <-
         purrr::map(
            RT_of_maxTIC,
            ~which(
               abs(.x - scanNumber_and_RT_vec) ==
                  min(abs(.x - scanNumber_and_RT_vec))
            )
         )

      # Make summary of XIC data

      sumXIC_summary <-
         sumXIC2 %>%
         rlang::set_names(chemform) %>%
         purrr::map(
            ~dplyr::summarize(.x, max_TIC = max(int_sum))
         ) %>%
         purrr::map2(
            names(.),
            ~dplyr::mutate(.x, chem_form = .y)
         ) %>%
         purrr::reduce(dplyr::union_all) %>%
         dplyr::mutate(Sequence = unlist(target_seqs)) %>%
         dplyr::mutate(seq_name = names(target_seqs)) %>%
         dplyr::mutate(RT_of_maxTIC = unlist(RT_of_maxTIC)) %>%
         dplyr::mutate(scan_of_maxTIC = unlist(scanNumsToRead)) %>%
         dplyr::select(seq_name, Sequence, chem_form, max_TIC, tidyr::everything()) %>%
         dplyr::arrange(dplyr::desc(max_TIC))

      # Arrange XICs to go in order of decreasing intensity

      sumXIC2 <-
         sumXIC2 %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)]


      PTM_names_list <-
         PTM_names_list %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)]

      iso_dist_list_union <-
         iso_dist_list_union %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)]

      scanNumsToRead <-
         scanNumsToRead %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)]

      # Plot XICs ---------------------------------------------------------------

      XIC_plots <-
         purrr::pmap(
            list(
               sumXIC2,
               PTM_names_list,
               names(sumXIC2)
            ),
            ~make_XIC_plot(
               ..1,
               times,
               int_sum,
               ..2,
               ..3
            )
         )

      ########## EVERYTHING AFTER THIS IS CONTROLLED BY save_output!

      if (save_output == TRUE) {

         # MultiArrange XIC grobs -------------------------------------------------------

         XIC_groblist <-
            gridExtra::marrangeGrob(
               grobs = XIC_plots,
               ncol = 5,
               nrow = 3,
               top = rawFileName
            )

         timer$stop("Calculate and plot XICs")

         # Save XIC results --------------------------------------------------------

         timer$start("Save XIC results")

         # Save target sequences to saveDir

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
            file =
               fs::path(
                  saveDir,
                  paste0(
                     fs::path_ext_remove(
                        stringr::str_trunc(
                           rawFileName, 90, "right", ellipsis = ""
                        )
                     ),
                     '_seqs.csv'
                  )
               )
         )

         # Save XIC data to saveDir

         readr::write_csv(
            sumXIC_summary,
            path =
               fs::path(
                  saveDir,
                  paste0(
                     fs::path_ext_remove(
                        stringr::str_trunc(
                           rawFileName, 90, "right", ellipsis = ""
                        )
                     ),
                     '_XICsum.csv'
                  )
               )
         )

         # Save theo isotope distributions to saveDir

         purrr::imap(
            iso_dist_list_union,
            ~dplyr::mutate(.x, accession = .y)
         ) %>%
            purrr::reduce(dplyr::union_all) %>%
            dplyr::select(accession, tidyr::everything()) %>%
            readr::write_csv(
               path =
                  fs::path(
                     saveDir,
                     paste0(
                        fs::path_ext_remove(
                           stringr::str_trunc(
                              rawFileName, 90, "right", ellipsis = ""
                           )
                        ),
                        '_isodist.csv'
                     )
                  )
            )

         # Save chromatograms, all together

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
                     '_XICs.pdf'
                  )
               ),
            plot = XIC_groblist,
            width = 20,
            height = 12,
         )

         timer$stop("Save XIC results")

      }


      # make_every_spectrum -----------------------------------------------------

      ## ggplot themes -----------------------------------------------------------


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



      ## Analyze data ------------------------------------------------------------


      timer$start("Make MS, Top 1 most intense PART 1")

      scansToPlot <-
         purrr::map(
            scanNumsToRead,
            ~rawrr::readSpectrum(
               rawfile = rawFile,
               scan = .x,
               tmpdir = rawrrTemp
            )
         )

      spectra_highestTIC <-
         purrr::imap(
            scansToPlot,
            ~tibble::tibble(
               UNIPROTKB = .y,
               scan = .x[[1]]$scan,
               mz = .x[[1]]$mZ,
               intensity = .x[[1]]$intensity)
         ) %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)]

      # Estimate noise for all spectra using MALDIquant::estimateNoise

      spectra_noiseEstimates <-
         spectra_highestTIC %>%
         purrr::map(
            ~MALDIquant::createMassSpectrum(
               mass = .x$mz,
               intensity = .x$intensity,
               metaData = list(name = .x$UNIPROTKB[[1]])
            ) %>%
               MALDIquant::estimateNoise(method = "MAD") %>%
               .[,2] %>%
               .[[1]]
         )

      # Add noise estimates to spectra

      spectra_highestTIC <-
         purrr::map2(
            spectra_highestTIC,
            spectra_noiseEstimates,
            ~dplyr::mutate(.x, noise = .y)
         )

      # m/z value of highest abundance theoretical isotopologue peaks.
      # These are the standard for comparison of list lengths!

      message("Getting highest intensity isotopologues")

      mz_max_abund <-
         iso_dist_list_union %>%
         purrr::map(
            ~dplyr::group_by(.x, charge) %>%
               dplyr::filter(abundance == max(abundance)) %>%
               dplyr::slice(1)
         ) %>%
         purrr::map(
            ~dplyr::pull(.x, `m/z`)
         ) %>%
         purrr::map(as.list) %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)]

      standard_lengths <-
         purrr::map(mz_max_abund, length)

      message("Getting charges of highest intensity isotopologues")

      mz_max_abund_charge <-
         iso_dist_list_union %>%
         purrr::map(
            ~dplyr::group_by(.x, charge) %>%
               dplyr::filter(abundance == max(abundance)) %>%
               dplyr::slice(1)
         ) %>%
         purrr::map(
            ~dplyr::pull(.x, charge)
         ) %>%
         purrr::map(as.list) %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)] %>%
         purrr::map2(
            standard_lengths,
            ~fix_list_length(.x, .y)
         )

      # Get m/zs for all theoretical isotopologues by charge state
      # for use in getting observed intensities for each isotopologue

      message("Getting relative abundances of isotopologues")

      mz_theo_MS1 <-
         iso_dist_list_union %>%
         purrr::map(
            ~dplyr::group_by(.x, charge) %>%
               dplyr::group_split() %>%
               purrr::map(
                  ~dplyr::select(.x, `m/z`) %>%
                     tibble::deframe()
               )
         ) %>%
         purrr::map2(
            standard_lengths,
            ~fix_list_length(.x, .y)
         )

      # Get abundances for all theoretical isotopologues by charge state
      # for use in calculating cosine similarities

      message("Getting isotopologue abundance per charge state")

      intensity_theo_MS1 <-
         iso_dist_list_union %>%
         purrr::map(
            ~dplyr::group_by(.x, charge) %>%
               dplyr::group_split() %>%
               purrr::map(
                  ~dplyr::select(.x, abundance) %>%
                     tibble::deframe()
               )
         ) %>%
         purrr::map2(
            standard_lengths,
            ~fix_list_length(.x, .y)
         )


      # Get first noise value from spectrum to use in results_chargestateTICs2
      # Theyre all the same anyway

      message("Getting noise estimate for highest intensity isotopologue")

      mz_max_abund_noise <-
         spectra_highestTIC %>%
         purrr::map(
            ~dplyr::pull(.x, noise) %>%
               .[[1]]
         ) %>%
         purrr::map(as.list) %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)]

      # Get values to use for vertical lines on plots

      message("Get values to use for vertical lines on plots")

      iso_dist_vlines2 <-
         iso_dist_vlines %>%
         purrr::map(
            tibble_list_checker
         ) %>%
         purrr::modify_depth(
            2,
            ~dplyr::filter(
               .x,
               abundance != 100
            ) %>%
               dplyr::top_n(10, abundance)
         ) %>%
         purrr::modify_depth(
            2,
            ~dplyr::pull(.x, `m/z`)
         ) %>%
         purrr::map(as.list) %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)] %>%
         purrr::map2(
            standard_lengths,
            ~fix_list_length(.x, .y)
         )

      spectra_highestTIC_names <-
         spectra_highestTIC %>%
         names() %>%
         as.list()

      timer$stop("Make MS, Top 1 most intense PART 1")

      ## Extract highest TICs from MS made in the previous step from within a
      ## narrow window of the theoretical most abundant peak (mz_window * mz_window_scaling),
      ## add all charge states together

      timer$start("Make MS, Top 1 most intense, PART 2")

      message("Getting intensity for theoretical highest abundance isotopologue")

      results_chargestateTICs <-
         purrr::map2(
            spectra_highestTIC %>% purrr::map(list),
            mz_max_abund,
            ~purrr::map2(
               .x,
               .y,
               ~get_maxY_in_Xrange(
                  df = .x,
                  x = mz,
                  y = intensity,
                  xrange =
                     c(
                        .y - (.y/resPowerMS1), .y + (.y/resPowerMS1)
                     )
               )
            )
         ) %>%
         purrr::map2(
            standard_lengths,
            ~fix_list_length(.x, .y)
         )

      ## Extract highest TICs from MS made in the previous step from within a
      ## narrow window of the all isotopologue peaks (mz_window * mz_window_scaling)

      message("Getting intensities for all isotopologues")

      intensity_obs_MS1 <-
         purrr::map2(
            spectra_highestTIC %>% purrr::map(list),
            mz_theo_MS1,
            ~purrr::map2(
               .x,
               .y,
               ~get_maxY_in_Xrange_vector(
                  df = .x,
                  x = mz,
                  y = intensity,
                  mz = .y,
                  res_power = resPowerMS1,
                  isotopologueWinMultiplier = isotopologue_window_multiplier
               )
            )
         ) %>%
         purrr::map2(
            standard_lengths,
            ~fix_list_length(.x, .y)
         )

      ## Extract highest TICs from MS made in the previous step from within a
      ## narrow window of the all isotopologue peaks (mz_window * mz_window_scaling)

      message("Getting m/z values at max intensity for all isotopologues")

      mz_obs_MS1 <-
         purrr::map2(
            spectra_highestTIC %>% purrr::map(list),
            mz_theo_MS1,
            ~purrr::map2(
               .x,
               .y,
               ~get_maxX_in_Xrange_vector(
                  df = .x,
                  x = mz,
                  y = intensity,
                  mz = .y,
                  res_power = resPowerMS1,
                  isotopologueWinMultiplier = isotopologue_window_multiplier
               )
            )
         ) %>%
         purrr::map2(
            standard_lengths,
            ~fix_list_length(.x, .y)
         )


      # Get abundances for theoretical isotopologues by charge state
      # for use in plotting by scaling theoretical isotopologue
      # intensity distribution by observed values

      message("Scaling theoretical isotopic abundances")


      intensity_theo_scaling_factors <-
         purrr::map2(
            intensity_obs_MS1,
            intensity_theo_MS1,
            ~purrr::map2(
               .x = .x,
               .y = .y,
               ~max(.x)/max(.y)
            )
         )


      intensity_theo_MS1_scaled <-
         purrr::pmap(
            list(
               intensity_theo_MS1,
               intensity_theo_scaling_factors,
               intensity_obs_MS1
            ),
            ~purrr::pmap(
               list(
                  ..1,
                  ..2,
                  ..3
               ),
               ~as.numeric(.x)*as.numeric(.y) %>%
                  {if (length(.) == 0) rep(0, length(..3)) else .}
            )
         ) %>%
         purrr::map2(
            standard_lengths,
            ~fix_list_length(.x, .y)
         )

      ## Try to determine peak resolving power using linear interpolation

      rp_obs_MS1 <-
         purrr::pmap(
            list(
               spectra_highestTIC,
               mz_theo_MS1
            ),
            ~purrr::pmap(
               list(
                  ..1 = list(..1),
                  ..2 = ..2
               ),
               ~purrr::pmap(
                  list(
                     ..1 = list(..1),
                     ..2 = ..2
                  ),
                  ~get_rp_in_Xrange_vector(
                     df = ..1,
                     x = mz,
                     y = intensity,
                     mz = ..2,
                     res_power = resPowerMS1,
                     isotopologueWinMultiplier = isotopologue_window_multiplier
                  )
               ) %>%
                  unlist()
            )
         )

      # Convert all resolving powers into a data frame for linear modeling

      rp_obs_MS1_df <-
         purrr::map2(
            rp_obs_MS1,
            mz_theo_MS1,
            ~purrr::map2(
               .x,
               .y,
               ~tibble::tibble(
                  rp = .x,
                  mz = .y
               ) %>%
                  tidyr::unnest() %>%
                  dplyr::filter(rp != 0)
            )
         ) %>%
         dplyr::bind_rows() %>%
         dplyr::arrange(mz)

      # Calculate linear model for RP based on determined non-zero values

      x <- rp_obs_MS1_df$mz
      y <- rp_obs_MS1_df$rp

      resolving_power_model <-
         fit <- lm(y ~ I(1/x))

      # Calculate theoretical resolving powers for every isotopologue peak

      rp_theo_MS1 <-
         purrr::map_depth(
            mz_theo_MS1,
            2,
            ~predict(
               resolving_power_model,
               newdata =
                  data.frame(
                     x = .x
                  )
            )
         )

      # Replace all zeros in rp_obs_MS1 with theoretical values from the fit

      rp_obs_theo_hybrid_MS1 <-
         purrr::map2(
            rp_obs_MS1,
            rp_theo_MS1,
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



      ## Calculate cosine similarities

      message("Calculating cosine similarities")

      cosine_sims <-
         purrr::map2(
            intensity_obs_MS1,
            intensity_theo_MS1,
            ~purrr::map2(
               .x,
               .y,
               ~calculate_cosine_similarity(
                  .x,
                  .y
               ) %>%
                  {if (length(.) == 0 | is.nan(.) | is.null(.)) 0 else .}
            )
         ) %>%
         purrr::map2(
            standard_lengths,
            ~fix_list_length(.x, .y)
         )

      ## Calculate S/N estimates

      SN_estimate_obs_MS1 <-
         purrr::map2(
            intensity_obs_MS1,
            mz_max_abund_noise,
            ~purrr::map2(
               .x,
               .y,
               ~(.x/.y) %>%
                  {if (length(.) == 0 | is.nan(.) | is.null(.)) 0 else .}
            )
         ) %>%
         purrr::map(
            ~purrr::map_if(
               .x,
               ~length(.x) == 0,
               ~0
            )
         ) %>%
         purrr::map2(
            standard_lengths,
            ~fix_list_length(.x, .y)
         )


      SN_estimate_theo_MS1 <-
         purrr::pmap(
            list(
               intensity_theo_MS1_scaled,
               mz_max_abund_noise,
               SN_estimate_obs_MS1
            ),
            ~purrr::pmap(
               list(
                  ..1,
                  ..2,
                  ..3
               ),
               ~(..1/..2) %>%
                  {if (length(.) == 0) rep(0, length(..3)) else .}
               # {if (is.nan(.) | is.null(.) | is.infinite(.)) rep(0, length(.x)) else .}
            )
         ) %>%
         purrr::map2(
            standard_lengths,
            ~fix_list_length(.x, .y)
         )

      ## Calculate Mel's ScoreMFA

      score_MFA <-
         purrr::pmap(
            list(
               mz_obs_MS1,
               mz_theo_MS1,
               rp_obs_theo_hybrid_MS1,
               SN_estimate_obs_MS1,
               SN_estimate_theo_MS1
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
         ) %>%
         purrr::map2(
            standard_lengths,
            ~fix_list_length(.x, .y)
         )

      ## Save score_MFA data for further analysis, TESTING PURPOSES ONLY

      # saveRDS(score_MFA, paste0(saveDir, "/score_MFA.rds"))
      #
      # saveRDS(mz_obs_MS1, paste0(saveDir, "/mz_obs_MS1.rds"))
      # saveRDS(mz_theo_MS1, paste0(saveDir, "/mz_theo_MS1.rds"))
      #
      # saveRDS(SN_estimate_obs_MS1, paste0(saveDir, "/SN_estimate_obs_MS1.rds"))
      # saveRDS(SN_estimate_theo_MS1, paste0(saveDir, "/SN_estimate_theo_MS1.rds"))

      ## Make table with all charge state TICs

      message("Making table with all charge state TICs")

      results_chargestateTICs2 <-
         purrr::pmap(
            list(
               spectra_highestTIC_names,
               mz_max_abund,
               mz_max_abund_charge,
               PTM_names_list,
               results_chargestateTICs,
               mz_max_abund_noise,
               cosine_sims,
               score_MFA
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
               ~tibble::tibble(
                  name = ..1,
                  mz_max_abund = ..2,
                  mz_max_abund_charge = ..3,
                  PTM_name = ..4,
                  maxTIC = ..5,
                  noise_estimate = ..6,
                  `estimated_S/N` = round(maxTIC/noise_estimate, digits = 0),
                  cosine_sim = ..7,
                  score_MFA = ..8
               )
            )
         ) %>%
         purrr::reduce(dplyr::union_all) %>%
         purrr::reduce(dplyr::union_all) %>%
         dplyr::distinct()

      results_chargestateTICs_summary <-
         results_chargestateTICs2 %>%
         dplyr::group_by(name) %>%
         dplyr::filter(
            `estimated_S/N` > SN_cutoff,
            cosine_sim != 0,
            score_MFA != 0
         ) %>%
         dplyr::summarize(
            charge_states = paste(unique(mz_max_abund_charge), collapse = ", "),
            PTM_name = PTM_name[[1]],
            maxTICsum = sum(maxTIC),
            `highest_S/N` = max(`estimated_S/N`),
            `lowest_S/N` = min(`estimated_S/N`),
            `mean_S/N` = mean(`estimated_S/N`),
            mean_cosine_sim = mean(cosine_sim),
            mean_score_mfa = mean(score_MFA)
         ) %>%
         dplyr::arrange(mean_cosine_sim)

      results_chargestateTICs_summary_filtered <-
         results_chargestateTICs_summary %>%
         dplyr::filter(
            `highest_S/N` > SN_cutoff,
            mean_score_mfa > scoreMFAcutoff,
            mean_cosine_sim > cosinesimcutoff
         )

      ## Make all MS based on mz of max abundance for each charge state

      message("Making plots")

      spectra_highestTIC_plots <-
         purrr::pmap(
            list(
               spectra_highestTIC %>% purrr::map(list),
               spectra_highestTIC_names,
               mz_max_abund,
               mz_max_abund_charge,
               PTM_names_list,
               iso_dist_vlines2,
               results_chargestateTICs,
               cosine_sims,
               mz_theo_MS1,
               intensity_theo_MS1_scaled,
               score_MFA
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
                  ..11
               ),
               ~{
                  if (..11 > scoreMFAcutoff & ..8 > cosinesimcutoff & mean(dplyr::pull(..1, noise)) > SN_cutoff) {
                     make_spectrum_top1(
                        df = ..1,
                        x = mz,
                        y = intensity,
                        noise = noise,
                        accession = ..2,
                        scan_num = scan,
                        charge = ..4,
                        xrange = c(..3 - (mz_window/2), ..3 + (mz_window/2)),
                        ..5,
                        vlines = ..6,
                        chargestateTIC = ..7,
                        cosine_sims = ..8,
                        mz_theo = ..9,
                        intensity_theo = ..10,
                        isotopologue_window = isotopologue_window_multiplier*(..9/resPowerMS1),
                        theme = MStheme01,
                        score_mfa = ..11
                     )
                  }
               }
            ) %>%
               purrr::set_names(..4)
         ) %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)] %>%
         purrr::map(
            ~ggplot_list_checker(.x)
         )

      # SAVE PLOTS OBJECT, FOR DEV PURPOSES

      if (save_spec_object == TRUE) {

         saveRDS(spectra_highestTIC_plots, paste0(saveDir, "/spectra_highestTIC_plots.rds"))
         saveRDS(target_seqs, paste0(saveDir, "/target_seqs.rds"))

      }

      timer$stop("Make MS, Top 1 most intense, PART 2")

      ## Arrange MS Grobs, Top 1 ----------------------------------------------

      timer$start("Arrange MS grobs, Top 1 most intense")

      message("Making tablegrob list for PDF output")

      tablegrob_list_multi <-
         spectra_highestTIC_plots %>%
         unlist(recursive = FALSE) %>%
         gridExtra::marrangeGrob(
            grobs = .,
            ncol = 5,
            nrow = 3,
            top = paste0(rawFileName)
         )

      timer$stop("Arrange MS grobs, Top 1 most intense")

      ## Save arranged MS, Top 1 -------------------------------------------------

      # Save max TIC for charge states

      writexl::write_xlsx(
         list(
            "TIC per CS" = results_chargestateTICs2,
            "TIC summary" = results_chargestateTICs_summary,
            "TIC summary filtered" = results_chargestateTICs_summary_filtered
         ),
         fs::path(
            saveDir,
            paste0(
               fs::path_ext_remove(rawFileName),
               "_maxTIC.xlsx"
            )
         ),
         format_headers = TRUE
      )

      # Multi-arranged

      timer$start("Save MS, Top 1, PDF")

      ggplot2::ggsave(
         filename =
            fs::path(
               saveDir,
               paste0(
                  fs::path_ext_remove(rawFileName),
                  "_specZoom.pdf"
               )
            ),
         plot = tablegrob_list_multi,
         width = 20,
         height = 12,
         limitsize = FALSE
      )


      # Write params to text file

      params_path <-
         fs::path(
            saveDir,
            'GEX_params_MS1.txt'
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
               target_charges = {toString(target_charges)}
               mass_range = {toString(mass_range)}
               mz_range = {toString(mz_range)}
               abund_cutoff = {abund_cutoff}
               XIC_tol = {XIC_tol}
               use_IAA = {use_IAA}
               save_output = {save_output}
               scoreMFAcutoff = {scoreMFAcutoff}
               cosinesimcutoff = {cosinesimcutoff}
               SN_cutoff = {SN_cutoff}
               resPowerMS1 = {resPowerMS1}
               isotopologue_window_multiplier = {isotopologue_window_multiplier}
               mz_window = {mz_window}
               return_timers = {return_timers}
               abund_cutoff = {abund_cutoff}
               rawrrTemp = {rawrrTemp}
               "
            ),
            file = params_path
         )

      }

      message("\n\n make_every_spectrum done")

      timer$stop("Save MS, Top 1, PDF")


      if (return_timers == TRUE) {
         return(
            timeR::getTimer(timer)
         )
      } else {
         return(
            spectra_highestTIC_plots
         )
      }


   }
