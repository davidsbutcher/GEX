#' make_every_spectrum
#'
#' @param make_every_XIC_output
#' @param mz_window
#' @param MSoutputWidth
#' @param MSoutputDPI
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
      mz_window = 4,
      MSoutputWidth = 18,
      MSoutputDPI = 200,
      makePNG = FALSE
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

      top_n_pforms <-
         make_every_XIC_output[[9]]


      # Make future workers -----------------------------------------------------

      if (makePNG == TRUE) {

         timer$start("Make future workers")

         if (is.null(top_n_pforms) == FALSE) {

            if (top_n_pforms < 10) {

               future::plan(
                  future::multisession(
                     workers = as.integer(top_n_pforms),
                     gc = TRUE,
                     persistent = FALSE
                  )
               )

            }

            if (top_n_pforms >= 10) {

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
         )

      mz_max_abund <-
         iso_dist_list_union %>%
         purrr::map(
            ~dplyr::filter(.x, abundance == 100)
         ) %>%
         purrr::map(
            ~dplyr::pull(.x, `m/z`)
         ) %>%
         purrr::map(as.list)

      mz_max_abund_charge <-
         iso_dist_list_union %>%
         purrr::map(
            ~dplyr::filter(.x, abundance == 100)
         ) %>%
         purrr::map(
            ~dplyr::pull(.x, charge)
         ) %>%
         purrr::map(as.list)

      spectra_highestTIC_names <-
         spectra_highestTIC %>%
         names() %>%
         purrr::map(list) %>%
         purrr::map2(
            mz_max_abund,
            ~rep(.x, length(.y))
         )

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
                     isoWindowLow < .x & isoWindowHigh > .x
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

      timer$start("Make MS, Top 1 most intense, PART 2")

      spectra_highestTIC_plots <-
         purrr::pmap(
            list(
               spectra_highestTIC %>% map(list),
               spectra_highestTIC_names,
               mz_max_abund,
               mz_max_abund_charge
            ),
            ~purrr::pmap(
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
                  theme = MStheme01
               )
            )
         )

      timer$stop("Make MS, Top 1 most intense, PART 2")

      # Arrange MS Grobs, Top 1 ----------------------------------------------

      timer$start("Arrange MS grobs, Top 1 most intense")

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

      # Multi-arranged

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

      {

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
            "_potential_MS2_scans.csv"
         )
      )

      # Multi-arranged

      timer$start("Save MS, Top 1, PDF")

      ggplot2::ggsave(
         filename = paste0(
            saveDir,
            fs::path_ext_remove(rawFileName),
            "_spectrum_zooms.pdf"
         ),
         plot = tablegrob_list_multi,
         width = 20,
         height = 12,
         limitsize = FALSE
      )

      message("\n\n make_every_spectrum done")

      timer$stop("Save MS, Top 1, PDF")

      return(
         timeR::getTimer(timer)
      )

   }
