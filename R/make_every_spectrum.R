#' make_every_spectrum
#'
#' @param make_every_XIC_output Direct output from make_every_XIC.
#' @param mz_window Width of the m/z window to be used to make mini spectra.
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
      SN_cutoff = 10,
      mean_cosine_sim_cutoff = 0.8,
      mean_score_mfa_cutoff = 0.3,
      resPowerMS1 = 300000,
      isotopologue_window_multiplier = 6,
      MSoutputWidth = 18,
      MSoutputDPI = 200,
      makePNG = FALSE,
      return_timers = TRUE,
      save_spec_object = FALSE
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

      # timer$start("Make future workers")
      #
      # future::plan("future::multisession", workers = future_workers)
      #
      # timer$stop("Make future workers")

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

      mz_all_abund <-
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

      iso_abund_per_charge_state <-
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

      # Use new method to get noise estimate based on a window centered on
      # the m/z of max abundance with width 10*mz_window

      # THIS IS DEPRECATED FOR NOW. With reduced profile mode data, the dearth
      # of points in the vicinity of some peaks makes local noise estimation
      # problematic

      # spectra_noiseEstimates_window <-
      #    purrr::map2(
      #       spectra_highestTIC,
      #       mz_max_abund,
      #       ~purrr::map2(
      #          list(.x),
      #          .y,
      #          ~dplyr::filter(
      #             .x,
      #             mz > .y - (mz_window*100)/2,
      #             mz < .y + (mz_window*100)/2
      #          ) %>%
      #             {MALDIquant::createMassSpectrum(
      #                mass = .$mz,
      #                intensity = .$intensity
      #             )} %>%
      #             MALDIquant::estimateNoise(method = "MAD") %>%
      #             {if (length(.) == 1) c(0) else .[,2]} %>%
      #             .[[1]] %>%
      #             {if (. == 0) 1E9 else .}
      #       )
      #    )

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

      # Find potential MS2 scans

      message("Trying to find potential MS2 scans")

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
            results_isotopologueTICs,
            iso_abund_per_charge_state,
            ~purrr::map2(
               .x = .x,
               .y = .y,
               ~max(.x)/max(.y)
            )
         )


      iso_abund_theoretical_scaled <-
         purrr::pmap(
            list(
               iso_abund_per_charge_state,
               intensity_theo_scaling_factors,
               results_isotopologueTICs
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
               mz_all_abund
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
            mz_all_abund,
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
            mz_all_abund,
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
            results_isotopologueTICs,
            iso_abund_per_charge_state,
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

      # SN_estimate_obs_OLD <-
      #    purrr::map2(
      #       results_isotopologueTICs,
      #       mz_max_abund_noise,
      #       ~purrr::map2(
      #          .x,
      #          .y,
      #          ~(.x/.y) %>%
      #             {if (length(.) == 0 | is.nan(.) | is.null(.)) 0 else .}
      #       )
      #    ) %>%
      #    purrr::map(
      #       ~purrr::map_if(
      #          .x,
      #          ~length(.x) == 0,
      #          ~0
      #       )
      #    ) %>%
      #    purrr::map2(
      #       standard_lengths,
      #       ~fix_list_length(.x, .y)
      #    )

      SN_estimate_obs <-
         purrr::map2(
            results_isotopologueTICs,
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


      SN_estimate_theo <-
         purrr::pmap(
            list(
               iso_abund_theoretical_scaled,
               mz_max_abund_noise,
               SN_estimate_obs
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
         # purrr::map2(
         #    SN_estimate_obs,
         #    ~purrr::map_if(
         #       .x,
         #       ~length(.x) == 0,
         #       ~rep(0, length(.y))
         #    )
         # ) %>%
         # purrr::map(
         #    ~purrr::map_if(
         #       .x,
      #       ~any(is.nan(.x) | is.null(.) | is.infinite(.)),
      #       ~rep(0, length(.x))
      #    )
      # ) %>%
      purrr::map2(
         standard_lengths,
         ~fix_list_length(.x, .y)
      )


      ## Calculate Mel's ScoreMFA

      score_MFA <-
         purrr::pmap(
            list(
               results_isotopologueMZ,
               mz_all_abund,
               rp_obs_theo_hybrid_MS1,
               SN_estimate_obs,
               SN_estimate_theo
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
      # saveRDS(results_isotopologueMZ, paste0(saveDir, "/results_isotopologueMZ.rds"))
      # saveRDS(mz_all_abund, paste0(saveDir, "/mz_all_abund.rds"))
      #
      # saveRDS(SN_estimate_obs, paste0(saveDir, "/SN_estimate_obs.rds"))
      # saveRDS(SN_estimate_theo, paste0(saveDir, "/SN_estimate_theo.rds"))

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
            # `highest_S/N` > SN_cutoff,
            mean_cosine_sim > mean_cosine_sim_cutoff,
            mean_score_mfa > mean_score_mfa_cutoff
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
               mz_all_abund,
               iso_abund_theoretical_scaled,
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
                  isotopologue_window = isotopologue_window_multiplier*(..9/resPowerMS1),
                  theme = MStheme01,
                  score_mfa = ..11
               )
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

         if (dir.exists(paste0(saveDir, "/mass_spec/")) == FALSE) {

            dir.create(paste0(saveDir, "/mass_spec/"))

         }

         timer$start("Save MS, Top 1, PNG")

         purrr::pwalk(
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
         fs::path(
            saveDir,
            paste0(
               fs::path_ext_remove(rawFileName),
               "_MS2.csv"
            )
         )
      )

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
            spectra_highestTIC_plots
         )
      }
   }
