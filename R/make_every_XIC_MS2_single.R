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
#' @param fragment_types
#' @param fragment_mz_range
#' @param fragment_pos_cutoff
#' @param abund_cutoff
#' @param XIC_tol_MS2
#' @param XIC_cutoff
#' @param use_IAA
#' @param include_PTMs
#' @param scoreMFAcutoff
#' @param cosinesimcutoff
#' @param SN_cutoff
#' @param resPowerMS2
#' @param isotopologue_window_multiplier
#' @param mz_window
#' @param rawrrTemp
#' @param save_spec_object
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
      isoAbund = c("12C" = 0.9889, "14N" = 0.99636),
      fragment_charges = c(1:50),
      fragment_types = c("b", "y"),
      fragment_mz_range = c(300,2000),
      fragment_pos_cutoff = c(1, 50),
      abund_cutoff = 5,
      XIC_tol_MS2 = 10,
      XIC_cutoff = 0.000001,
      use_IAA = FALSE,
      include_PTMs = TRUE,
      scoreMFAcutoff = 0.3,
      cosinesimcutoff = 0.99,
      SN_cutoff = 10,
      resPowerMS2 = 150000,
      isotopologue_window_multiplier = 6,
      mz_window = 5,
      rawrrTemp = tempdir(),
      save_spec_object = FALSE
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
      # tidyr::separate(
      #    !!PTMname_col_name,
      #    c("PTM", "PTMpos"),
      #    sep = "@",
      #    remove = FALSE
      # ) %>%
      # dplyr::mutate(PTMpos = as.numeric(PTMpos))


      target_seqs <-
         target_seqs_df %>%
         dplyr::select(
            !!target_col_name, !!target_sequence_col_name
         ) %>%
         tibble::deframe() %>%
         as.list()

      if (include_PTMs == TRUE) {

         target_PTM_split <-
            as.list(target_seqs_df$PTMname) %>%
            purrr::map(
               ~stringr::str_split(.x, "; ")
            )

         target_PTM <-
            target_PTM_split %>%
            purrr::map(
               ~stringr::str_extract(.x[[1]], "[^@]+")
            )

         target_PTM_pos <-
            target_PTM_split %>%
            purrr::map(
               ~stringr::str_extract(.x[[1]], "(?<=@).*") %>%
                  as.numeric()
            )

         target_PTM_chemform <-
            target_seqs_df %>%
            dplyr::pull("FormulaToAdd") %>%
            as.list() %>%
            purrr::map(
               ~stringr::str_split(.x, "; ")[[1]]
            )

         rm(target_PTM_split)

      }

      # target_PTM <-
      #    target_seqs_df %>%
      #    dplyr::pull("PTM")
      #
      # target_PTM_pos <-
      #    target_seqs_df %>%
      #    dplyr::pull("PTMpos")
      #
      # target_PTM_chemform <-
      #    target_seqs_df %>%
      #    dplyr::pull("FormulaToAdd") %>%
      #    tidyr::replace_na("C0")


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


      if (fs::dir_exists(saveDir) == FALSE) fs::dir_create(saveDir)

      ## Save target seqs ----

      seqs_path  <-
         fs::path(
            saveDir,
            paste0(
               fs::path_ext_remove(
                  stringr::str_trunc(
                     rawFileName, 90, "right", ellipsis = ""
                  )
               ),
               '_target_seqs_MS2'
            ),
            ext = "csv"
         )

      if (fs::file_exists(seqs_path)) {

         seqs_path <-
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

      }

      message(
         paste0(
            "\nSaving target sequences to ",
            seqs_path
         )
      )

      readr::write_csv(
         target_seqs_df,
         file = seqs_path
      )

      ## Write params to text file ----

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
               SN_cutoff = {SN_cutoff}
               resPowerMS2 = {resPowerMS2}
               isotopologue_window_multiplier = {isotopologue_window_multiplier}
               mz_window = {mz_window}
               use_IAA = {use_IAA}
               include_PTMs = {include_PTMs}
               abund_cutoff = {abund_cutoff}
               "
            ),
            file = params_path
         )

      }


      ## Prepare summary datasheet ----

      spectra_MS2_dataframe_union <-
         tibble::tibble()


      # Create save path for summary of assigned fragments

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


      ## Wrap functions ----------------------------------------------------------

      possibly_write_xlsx <-
         purrr::possibly(
            writexl::write_xlsx,
            otherwise = NULL
         )


      ## Prepare timer -----------------------------------------------------------

      timer <-
         timeR::createTimer(verbose = FALSE)

      timer_df <-
         tibble::tibble()

      # Create save path for timer summary

      timer_path <-
         fs::path(
            saveDir,
            'GEX_timers.csv'
         )

      # Get elemental compositions ----------------------------------------------

      for (i in seq_along(target_seqs)) {

         timer$start("1. Fragment library generation and isodist simulation")

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

         # Add PTM chemical formulas to appropriate fragments

         if (include_PTMs == TRUE) {

            fragments2[[1]] <-
               fragments2[[1]] %>%
               dplyr::mutate(
                  chemform_noPTM = chemform,
                  PTM = "",
                  PTM_chemform = ""
               )

            for (j in seq_along(target_PTM[[i]])) {

               fragments2_temp <-
                  fragments2[[1]] %>%
                  dplyr::mutate(
                     PTM_chemform_temp =
                        dplyr::case_when(
                           type == "b" & pos >= target_PTM_pos[[i]][j] ~ target_PTM_chemform[[i]][j],
                           type == "y" & pos > nchar(target_seqs[[i]]) - target_PTM_pos[[i]][j] ~ target_PTM_chemform[[i]][j],
                           TRUE ~ "C0"
                        ),
                     PTM =
                        dplyr::case_when(
                           type == "b" & pos >= target_PTM_pos[[i]][j] ~ paste0(PTM, target_PTM[[i]][j], sep = "; "),
                           type == "y" & pos > nchar(target_seqs[[i]]) - target_PTM_pos[[i]][j] ~ paste0(PTM, target_PTM[[i]][j], sep = "; "),
                           TRUE ~ ""
                        ),
                     PTM_chemform =
                        dplyr::case_when(
                           type == "b" & pos >= target_PTM_pos[[i]][j] ~ paste0(PTM_chemform, target_PTM_chemform[[i]][j], sep = "; "),
                           type == "y" & pos > nchar(target_seqs[[i]]) - target_PTM_pos[[i]][j] ~ paste0(PTM_chemform, target_PTM_chemform[[i]][j], sep = "; "),
                           TRUE ~ ""
                        )
                  )


               new_chemform <-
                  purrr::map2_chr(
                     fragments2_temp$chemform,
                     fragments2_temp$PTM_chemform_temp,
                     ~enviPat::mergeform(.x, .y)
                  )

               fragments2[[1]] <-
                  fragments2_temp %>%
                  dplyr::mutate(
                     chemform = new_chemform,
                     ion_chemform = paste(ion, chemform, sep = "_")
                  )

            }

            fragments2[[1]] <-
               fragments2[[1]] %>%
               dplyr::mutate(
                  ion_chemform = paste(ion, chemform, sep = "_")
               ) %>%
               dplyr::select(-PTM_chemform_temp)

         }


         # Prepare fragments data for combination with iso_dist_MS2 later

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

         timer$stop("1. Fragment library generation and isodist simulation")

         # Read XICs --------

         timer$start("2. Read XICs")

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

         timer$stop("2. Read XICs")

         timer$start("3. Process XICs")

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

         timer$stop("3. Process XICs")

         # Plot XICs -------------------------------------------------------------

         timer$start("4. Plot XICs and save XICs/isodists")

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

         message(glue::glue("\n\nSaving XICs for {names(target_seqs)[[i]]}"))
         message(paste("Memory used:", pryr::mem_used()/1000000, "MB"))


         XIC_plots_MS2_marrange <-
            gridExtra::marrangeGrob(
               grobs = purrr::flatten(XIC_plots_MS2),
               ncol = 5,
               nrow = 3,
               top = rawFileName
            )

         # Save XIC results ------

         saveDirXIC <-
            fs::path(
               saveDir,
               "XIC_MS2"
            )

         saveDirIsoDist <-
            fs::path(
               saveDir,
               "IsoDist_MS2"
            )

         saveDirFragLib <-
            fs::path(
               saveDir,
               "fragment_library_MS2"
            )

         if (fs::dir_exists(saveDirXIC) == FALSE) fs::dir_create(saveDirXIC)
         if (fs::dir_exists(saveDirIsoDist) == FALSE) fs::dir_create(saveDirIsoDist)
         if (fs::dir_exists(saveDirFragLib) == FALSE) fs::dir_create(saveDirFragLib)

         ## Save theo iso distributions -----

         writexl::write_xlsx(
            iso_dist_MS2_cluster,
            path =
               fs::path(
                  saveDirIsoDist,
                  paste0(
                     stringr::str_trunc(
                        names(target_seqs)[[i]], 90, "right", ellipsis = ""),
                     '_isodist_MS2.xlsx'
                  )
               )
         )

         writexl::write_xlsx(
            fragments2,
            path =
               fs::path(
                  saveDirFragLib,
                  paste0(
                     stringr::str_trunc(
                        names(target_seqs)[[i]], 90, "right", ellipsis = ""),
                     '_fragment_library_MS2.xlsx'
                  )
               )
         )

         ## Save XIC plots -----

         ggplot2::ggsave(
            filename =
               fs::path(
                  saveDirXIC,
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

         timer$stop("4. Plot XICs and save XICs/isodists")

         # Extract scans and make spectra ------------------------------------------

         timer$start("5. Extract MS data")

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

         timer$stop("5. Extract MS data")

         timer$start("6. Process MS data")


         # Process MS data ----

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

         # Calculate mean noise for all MS2 scans

         noise_obs_MS2_mean <-
            noise_obs_MS2 %>%
            unlist() %>%
            .[. != 0] %>%
            mean

         # Replace all zeros with the mean noise

         noise_obs_MS2 <-
            purrr::map_depth(
               noise_obs_MS2,
               3,
               ~{if (.x == 0) noise_obs_MS2_mean else .}
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


         ## Calculate spectral scores -----------------------------------------------

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
                        rp_mult = 2,
                        rm_zero = TRUE
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
                        ..2,
                        rm_zero = TRUE
                     ) %>%
                        {if (is.nan(.) | is.null(.) | length(.) == 0) 0 else .}
                  )
               )
            )

         # Calculate MMA for every peak

         MMA_MS2 <-
            purrr::pmap(
               list(
                  mz_obs_MS2,
                  mz_theo_MS2
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
                     ~calculate_mma_ppm(
                        ..1,
                        ..2
                     )
                        # {if (is.nan(.) | is.null(.) | length(.) == 0) NA else .}
                  )
               )
            )

         timer$stop("6. Process MS data")


         # Plot MS -----------------------------------------------------------------


         # Prepare data for plotting

         timer$start("7. Plot MS data, save plots and spreadsheets")

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
                  SN_estimate_obs_MS2,
                  mz_obs_MS2,
                  intensity_obs_MS2,
                  MMA_MS2
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
                     ..12,
                     ..13,
                     ..14,
                     ..15
                  ),
                  ~purrr::pmap(
                     list(
                        ..1,
                        rep(..2, length(..3)),
                        ..3,
                        ..4,
                        ..5,
                        ..6,
                        ..7,
                        ..8,
                        ..9,
                        ..10,
                        ..11,
                        ..12,
                        ..13,
                        ..14,
                        ..15
                     ),
                     ~{
                        if (..8 > scoreMFAcutoff & ..11 > cosinesimcutoff & max(..12) > SN_cutoff) {
                           make_spectrum_MS2(
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
                              mma = ..15,
                              mz_obs = ..13,
                              int_obs = ..14,
                              theme = MStheme01
                           )
                        } else {
                           NULL
                        }
                     }
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
                        ..1,
                        ..2,
                        ..3,
                        ..4,
                        ..5,
                        ..6,
                        ..7,
                        ..8
                     ),
                     ~{
                        if (..6 > scoreMFAcutoff & ..7 > cosinesimcutoff & max(..8) > SN_cutoff) {
                           tibble::tibble(
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
                              max_SN_estimate = max(..8)
                           )
                        } else {
                           NULL
                        }
                     }
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


         # Save MS results ---------------------------------------------------------

         # Save fragments dataframe for this sequence to saveDirFrag

         saveDirFrag <-
            fs::path(
               saveDir,
               "assigned_fragments_MS2"
            )

         if (fs::dir_exists(saveDirFrag) == FALSE) fs::dir_create(saveDirFrag)

         writexl::write_xlsx(
            spectra_MS2_dataframe,
            path =
               fs::path(
                  saveDirFrag,
                  paste0(
                     stringr::str_trunc(
                        names(target_seqs)[[i]], 90, "right", ellipsis = ""),
                     '_assigned_fragments_MS2.xlsx'
                  )
               )
         )

         # Add data frame to summary data frame of assigned fragments

         spectra_MS2_dataframe_union <-
            dplyr::bind_rows(
               spectra_MS2_dataframe_union,
               spectra_MS2_dataframe
            )

         # Save data frame of fragments which beat ScoreMFA cutoff

         message("\nWriting all 'assigned' fragment ions to spreadsheets")

         # writexl::write_xlsx(
         #    spectra_MS2_dataframe_union,
         #    path = sum_path
         # )

         possibly_write_xlsx(
            spectra_MS2_dataframe_union,
            path = sum_path
         )

         # Check memory usage


         message(glue::glue("\n\nSaving spectra for {names(target_seqs)[[i]]}"))
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

         saveDirSpec <-
            fs::path(
               saveDir,
               "specZoom_MS2"
            )

         if (fs::dir_exists(saveDirSpec) == FALSE) fs::dir_create(saveDirSpec)

         ggplot2::ggsave(
            filename =
               fs::path(
                  saveDirSpec,
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

         timer$stop("7. Plot MS data, save plots and spreadsheets")

         # Process timers

         timer_df_temp <-
            timeR::getTimer(timer) %>%
            tibble::as_tibble()

         timer_df <-
            dplyr::bind_rows(
               list(
                  timer_df,
                  timer_df_temp
               )
            )

         # SAVE PLOTS OBJECT, FOR DEV PURPOSES

         if (save_spec_object == TRUE) {

            message(glue::glue("\n\nSaving spectra object for {names(target_seqs)[[i]]}"))

            saveDirSpecObject <-
               fs::path(
                  saveDirSpec,
                  paste0(
                     stringr::str_trunc(
                        names(target_seqs)[[i]], 90, "right", ellipsis = ""
                     ),
                     '_specObject_MS2.rds'
                  )
               )

            saveRDS(spectra_MS2, saveDirSpecObject)

         }

      }

      timer_df_sum <-
         timer_df %>%
         dplyr::group_by(event) %>%
         dplyr::summarize(
            time_elapsed_sec = sum(timeElapsed)
         ) %>%
         dplyr::mutate(
            time_elapsed_min = time_elapsed_sec/60
         )

      readr::write_csv(
         timer_df_sum,
         file = timer_path
      )

      # Trying to deal with problem with futures not freeing up
      # memory correctly!

      rm(list = ls(all.names = TRUE))
      gc()


   }


