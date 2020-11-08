#' make_every_spectrum
#'
#' @param make_every_XIC_output Direct output from make_every_XIC.
#' @param mz_window Width of the m/z window to be used to make mini spectra.
#' @param mz_window_scaling Scaling factor for window for determination of intensity
#' for individual isotopologue peaks. Multiplied by mz_window value.
#' @param MSoutputWidth Width of png output. Only used is makePNG is TRUE.
#' @param MSoutputDPI DPI of png output. Only used is makePNG is TRUE.
#' @param makePNG Whether to make png output. Probably deprecated soon.
#' @param return_timers Controls whether timers or spectra are returned.
#'
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples

make_every_spectrum <-
   function(
      make_every_XIC_output,
      mz_window = 3,
      mz_window_scaling = 0.002,
      MSoutputWidth = 18,
      MSoutputDPI = 200,
      makePNG = FALSE,
      return_timers = TRUE
   ) {

      library(rawDiag)

      scanNumsToRead <-
         make_every_XIC_output[[1]]

      iso_dist_list_union <-
         make_every_XIC_output[[2]]

      timer <-
         make_every_XIC_output[[3]]

      rawFile <-
         make_every_XIC_output[[4]]

      rawFileMetadata <-
         make_every_XIC_output[[5]]

      saveDir <-
         make_every_XIC_output[[6]]

      rawFileName <-
         make_every_XIC_output[[7]]

      target_seqs <-
         make_every_XIC_output[[8]]

      sample_n_pforms <-
         make_every_XIC_output[[9]]

      PTM_names_list <-
         make_every_XIC_output[[10]]

      sumXIC_summary <-
         make_every_XIC_output[[11]]

      iso_dist_vlines <-
         make_every_XIC_output[[12]]


      # Make future workers -----------------------------------------------------

      if (makePNG == TRUE) {

         timer$start("Make future workers")

         if (is.null(sample_n_pforms) == FALSE) {

            if (sample_n_pforms < 10) {

               future::plan(
                  future::multisession(
                     workers = as.integer(sample_n_pforms),
                     gc = TRUE,
                     persistent = FALSE
                  )
               )

            }

            if (sample_n_pforms >= 10) {

               future::plan(
                  future::multisession(
                     workers = 10L,
                     gc = TRUE,
                     persistent = FALSE
                  )
               )
            }

         } else {

            future::plan(
               future::multisession(
                  workers = 5L,
                  gc = TRUE,
                  persistent = FALSE
               )
            )
         }

         timer$stop("Make future workers")

      }

      # ggplot themes -----------------------------------------------------------


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



      # Analyze data ------------------------------------------------------------


      timer$start("Make MS, Top 1 most intense PART 1")

      scansToPlot <-
         purrr::map(
            scanNumsToRead,
            ~rawDiag::readScans(
               rawFile,
               scans = .x
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
               .[,2]
         )

      # Add noise estimates to spectra

      spectra_highestTIC <-
         purrr::map2(
            spectra_highestTIC,
            spectra_noiseEstimates,
            ~dplyr::mutate(.x, noise = .y)
         )

      # Delete noise estimates to save memory

      # mz_max_abund <-
      #    iso_dist_list_union %>%
      #    purrr::map(
      #          ~dplyr::filter(.x, abundance == 100)
      #    ) %>%
      #    purrr::map(
      #       ~dplyr::pull(.x, `m/z`)
      #    ) %>%
      #    purrr::map(as.list) %>%
      #    .[sumXIC_summary %>% dplyr::pull(seq_name)]

      mz_max_abund <-
         iso_dist_list_union %>%
         purrr::map(
            ~dplyr::group_by(.x, charge) %>%
               dplyr::filter(abundance == 100) %>%
               dplyr::slice(1)
         ) %>%
         purrr::map(
            ~dplyr::pull(.x, `m/z`)
         ) %>%
         purrr::map(as.list) %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)]


      # mz_max_abund_charge <-
      #    iso_dist_list_union %>%
      #    purrr::map(
      #       ~dplyr::filter(.x, abundance == 100)
      #    ) %>%
      #    purrr::map(
      #       ~dplyr::pull(.x, charge)
      #    ) %>%
      #    purrr::map(as.list) %>%
      #    .[sumXIC_summary %>% dplyr::pull(seq_name)]

      mz_max_abund_charge <-
         iso_dist_list_union %>%
         purrr::map(
            ~dplyr::group_by(.x, charge) %>%
               dplyr::filter(abundance == 100) %>%
               dplyr::slice(1)
         ) %>%
         purrr::map(
            ~dplyr::pull(.x, charge)
         ) %>%
         purrr::map(as.list) %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)]


      # Get m/zs for all theoretical isotopologues by charge state
      # for use in getting observed intensities for each isotopologue

      mz_all_abund <-
         iso_dist_list_union %>%
         purrr::map(
            ~dplyr::group_by(.x, charge) %>%
               dplyr::group_split() %>%
               purrr::map(
                  ~dplyr::select(.x, `m/z`) %>%
                     tibble::deframe()
               )
         )
      # %>%
      #    purrr::map(
      #       rev
      #    )
      # %>%
      #    purrr::map2(
      #       mz_max_abund_charge,
      #       ~purrr::set_names(.x, .y)
      #    )

      # Get abundances for all theoretical isotopologues by charge state
      # for use in calculating cosine similarities

      iso_abund_per_charge_state <-
         iso_dist_list_union %>%
         purrr::map(
            ~dplyr::group_by(.x, charge) %>%
               dplyr::group_split() %>%
               purrr::map(
                  ~dplyr::select(.x, abundance) %>%
                     tibble::deframe()
               )
         )
      # %>%
      #    purrr::map2(
      #       mz_max_abund_charge,
      #       ~purrr::set_names(.x, .y)
      #    )


      # Get first noise value from spectrum to use in results_chargestateTICs2
      # Theyre all the same anyway

      mz_max_abund_noise <-
         spectra_highestTIC %>%
         purrr::map(
            ~dplyr::pull(.x, noise) %>%
               .[[1]]
         ) %>%
         purrr::map(as.list) %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)]

      # Get values to use for vertical lines on plots

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
               # dplyr::mutate(`m/z_round` = round(`m/z`, digits = 1)) %>%
               # dplyr::distinct(`m/z_round`, .keep_all = TRUE) %>%
               dplyr::top_n(10, abundance)
         ) %>%
         purrr::modify_depth(
            2,
            ~dplyr::pull(.x, `m/z`)
         ) %>%
         purrr::map(as.list) %>%
         # purrr::map(rev) %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)]


      ## ADDITIONAL PROCESSING OF ISO_DIST_VLINES2
      ## Lengths of list elements in iso_dist_vlines2 occasionally don't match other arguments to the pmap
      ## call that generates spectra_highestTIC_plots. This makes iso_dist_vlines2 match mz_max_abund

      for (i in seq_along(iso_dist_vlines2)) {

         for (j in rev(seq_along(iso_dist_vlines2[[i]]))) {

            if (length(iso_dist_vlines2[[i]][[j]]) == 0) {

               iso_dist_vlines2[[i]][[j]] <- NULL

            }

         }

      }

      for (i in seq_along(iso_dist_vlines2)) {

         if (length(iso_dist_vlines2[[i]]) > length(mz_max_abund[[i]])) {

            iso_dist_vlines2[[i]] <-
               iso_dist_vlines2[1:length(mz_max_abund[[i]])]

         } else if (length(iso_dist_vlines2[[i]]) < length(mz_max_abund[[i]])) {

            while (length(iso_dist_vlines2[[i]]) != length(mz_max_abund[[i]])) {

               iso_dist_vlines2[[i]] <-
                  append(iso_dist_vlines2[[i]], 0)

            }

         }

      }

      spectra_highestTIC_names <-
         spectra_highestTIC %>%
         names() %>%
         as.list()

      # Find potential MS2 scans

      rawFileMetadataMS2 <-
         rawFileMetadata %>%
         dplyr::select(scanNumber, StartTime, MSOrder, MS2IsolationWidth, PrecursorMass, TIC) %>%
         dplyr::filter(MSOrder == "Ms2") %>%
         dplyr::mutate(isoWindowLow = PrecursorMass - MS2IsolationWidth/2) %>%
         dplyr::mutate(isoWindowHigh = PrecursorMass + MS2IsolationWidth/2)

      potential_MS2 <-
         mz_max_abund %>%
         purrr::map(
            ~purrr::map(
               .x,
               ~{
                  dplyr::filter(
                     rawFileMetadataMS2,
                     .x > isoWindowLow & .x < isoWindowHigh
                  ) %>%
                     dplyr::arrange(dplyr::desc(TIC)) %>%
                     dplyr::pull(scanNumber) %>%
                     as.character() %>%
                     paste(collapse = ", ")
               }
            )
         ) %>%
         purrr::map2(
            mz_max_abund_charge,
            ~rlang::set_names(.x, .y)
         )

      potential_MS2_output <-
         purrr::map(
            potential_MS2,
            tibble::enframe,
            name = "charge",
            value = "potential_MS2_scans"
         ) %>%
         purrr::map(
            tidyr::unnest,
            cols = c(potential_MS2_scans)
         ) %>%
         purrr::imap(
            ~dplyr::mutate(.x, name = .y)
         ) %>%
         purrr::reduce(dplyr::union_all) %>%
         dplyr::filter(potential_MS2_scans != "") %>%
         dplyr::select(name, charge, potential_MS2_scans)

      timer$stop("Make MS, Top 1 most intense PART 1")

      ## Extract highest TICs from MS made in the previous step from within a
      ## narrow window of the theoretical most abundant peak (mz_window * mz_window_scaling),
      ## add all charge states together

      timer$start("Make MS, Top 1 most intense, PART 2")

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
                        .y - (mz_window*mz_window_scaling), .y + (mz_window*mz_window_scaling)
                     )
               )
            )
         )

      ## Extract highest TICs from MS made in the previous step from within a
      ## narrow window of the all isotopologue peaks (mz_window * mz_window_scaling)

      results_isotopologueTICs <-
         purrr::map2(
            spectra_highestTIC %>% purrr::map(list),
            mz_all_abund,
            ~purrr::map2(
               .x,
               .y,
               ~get_maxY_in_Xrange_vector(
                  df = .x,
                  x = mz,
                  y = intensity,
                  mz = .y,
                  mz_window = mz_window,
                  mz_window_scaling = mz_window_scaling
               )
            )
         )

      ## Extract highest TICs from MS made in the previous step from within a
      ## narrow window of the all isotopologue peaks (mz_window * mz_window_scaling)

      results_isotopologueMZ <-
         purrr::map2(
            spectra_highestTIC %>% purrr::map(list),
            mz_all_abund,
            ~purrr::map2(
               .x,
               .y,
               ~get_maxX_in_Xrange_vector(
                  df = .x,
                  x = mz,
                  y = intensity,
                  mz = .y,
                  mz_window = mz_window,
                  mz_window_scaling = mz_window_scaling
               )
            )
         )


      # Get abundances for theoretical isotopologues by charge state
      # for use in plotting by scaling theoretical isotopologue
      # intensity distribution by observed values

      iso_abund_theoretical_scaled <-
         purrr::map2(
            iso_abund_per_charge_state,
            results_isotopologueTICs,
            ~purrr::map2(
               .x,
               .y,
               ~scales::rescale(.x, c(min(.y), max(.y)))
            )
         )


      ## Calculate cosine similarities

      cosine_sims <-
         purrr::map2(
            results_isotopologueTICs,
            iso_abund_per_charge_state,
            ~purrr::map2(
               .x,
               .y,
               ~coop::cosine(
                  .x,
                  .y,
               ) %>%
                  {if (is.nan(.) | is.null(.) | length(.) == 0) 0 else .}
            )
         )

      ## Calculate spectral contrast angle

      spectral_contrast_angles <-
         purrr::pmap(
            list(
               results_isotopologueMZ,
               results_isotopologueTICs,
               mz_all_abund,
               iso_abund_per_charge_state
            ),
            ~purrr::pmap(
               list(
                  ..1,
                  ..2,
                  ..3,
                  ..4
               ),
               ~MicroRaman::SCA(
                  ..2 %>%
                     purrr::set_names(..1),
                  ..3 %>%
                     purrr::set_names(..4)
               ) %>%
                  {if (is.nan(.) | is.null(.) | length(.) == 0) 1 else .}
            )
         )

      ## Make table with all charge state TICs

      results_chargestateTICs2 <-
         purrr::pmap(
            list(
               spectra_highestTIC_names,
               mz_max_abund,
               mz_max_abund_charge,
               PTM_names_list,
               results_chargestateTICs,
               mz_max_abund_noise
            ),
            ~purrr::pmap(
               list(
                  ..1,
                  ..2,
                  ..3,
                  ..4,
                  ..5,
                  ..6
               ),
               ~tibble::tibble(
                  name = ..1,
                  mz_max_abund = ..2,
                  mz_max_abund_charge = ..3,
                  PTM_name = ..4,
                  maxTIC = ..5,
                  noise_estimate = ..6,
                  `estimated_S/N` = round(maxTIC/noise_estimate, digits = 0)
               )
            )
         ) %>%
         purrr::reduce(dplyr::union_all) %>%
         purrr::reduce(dplyr::union_all) %>%
         dplyr::distinct()

      results_chargestateTICs_summary <-
         results_chargestateTICs2 %>%
         dplyr::group_by(name) %>%
         dplyr::summarize(
            charge_states = paste(unique(mz_max_abund_charge), collapse = ", "),
            PTM_name = PTM_name[[1]],
            maxTICsum = sum(maxTIC),
            `highest_S/N` = max(`estimated_S/N`),
            `lowest_S/N` = min(`estimated_S/N`)
         ) %>%
         dplyr::arrange(desc(maxTICsum))

      ## Make all MS based on mz of max abundance for each charge state

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
               mz_all_abund,
               iso_abund_theoretical_scaled,
               spectral_contrast_angles
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
               ~make_spectrum_top1(
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
                  mz_all_abund = ..9,
                  iso_abund_theoretical_scaled = ..10,
                  spectral_contrast_angles = ..11,
                  theme = MStheme01
               )
            ) %>%
               purrr::set_names(..4)
         ) %>%
         .[sumXIC_summary %>% dplyr::pull(seq_name)] %>%
         purrr::map(
            ~ggplot_list_checker(.x)
         )

      timer$stop("Make MS, Top 1 most intense, PART 2")

      # Arrange MS Grobs, Top 1 ----------------------------------------------

      timer$start("Arrange MS grobs, Top 1 most intense")

      if (makePNG == TRUE) {

         tablegrob_list_top1 <-
            purrr::map(
               spectra_highestTIC_plots,
               ~gridExtra::arrangeGrob(
                  grobs = .x,
                  ncol = 4,
                  top = rawFileName
               )
            )

         tablegrob_list_top1_arranged <-
            purrr::map(
               tablegrob_list_top1,
               ~gridExtra::grid.arrange(.x),
               .progress = TRUE
            )

      }

      # Multi-arranged

      # spectra_highestTIC_plots <<- spectra_highestTIC_plots
      #
      # rawFileName <<- rawFileName

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

      # Save arranged MS, Top 1 -------------------------------------------------

      if (makePNG == TRUE) {

         tablegrob_filenames <-
            names(target_seqs) %>%
            as.list() %>%
            purrr::map(
               ~paste(
                  saveDir,
                  "mass_spec/",
                  fs::path_ext_remove(rawFileName),
                  "_",
                  .x,
                  "_spectrum_zooms.png",
                  sep = ""
               )
            )

         tablegrob_heights <-
            tablegrob_list_top1 %>%
            purrr::map(use_series, "heights") %>%
            purrr::map(length)

      }

      if (makePNG == TRUE) {

         if (dir.exists(paste0(saveDir, "mass_spec/")) == FALSE) {

            dir.create(paste0(saveDir, "mass_spec/"))

         }

         timer$start("Save MS, Top 1, PNG")

         furrr::future_pwalk(
            list(
               tablegrob_filenames,
               tablegrob_list_top1_arranged,
               tablegrob_heights
            ),
            ~ggplot2::ggsave(
               filename = ..1 ,
               plot = ..2,
               width = MSoutputWidth,
               height = ..3 * 3,
               dpi = MSoutputDPI
            )
         )

         timer$stop("Save MS, Top 1, PNG")

      }

      # Save potential MS2 scans info

      readr::write_csv(
         potential_MS2_output,
         paste0(
            saveDir,
            fs::path_ext_remove(rawFileName),
            "_MS2.csv"
         )
      )

      # Save max TIC for charge states

      writexl::write_xlsx(
         list(
            "TIC per CS" = results_chargestateTICs2,
            "TIC summary" = results_chargestateTICs_summary
         ),
         paste0(
            saveDir,
            fs::path_ext_remove(rawFileName),
            "_maxTIC.xlsx"
         ),
         format_headers = TRUE
      )

      # Multi-arranged

      timer$start("Save MS, Top 1, PDF")

      ggplot2::ggsave(
         filename = paste0(
            saveDir,
            fs::path_ext_remove(rawFileName),
            "_specZoom.pdf"
         ),
         plot = tablegrob_list_multi,
         width = 20,
         height = 12,
         limitsize = FALSE
      )

      message("\n\n make_every_spectrum done")

      timer$stop("Save MS, Top 1, PDF")


      # readr::write_lines(
      #    timer %>%  tibble::as_tibble(),
      #    paste0(
      #       saveDir,
      #       fs::path_ext_remove(rawFileName),
      #       "_GEX_timer.txt"
      #    )
      # )

      if (return_timers == TRUE) {
         return(
            timeR::getTimer(timer)
         )
      } else {
         return(
            list(
               spectra_highestTIC_plots
            )
         )
      }
   }
