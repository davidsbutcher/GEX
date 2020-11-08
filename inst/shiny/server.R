
options(repos = BiocManager::repositories())

options(shiny.maxRequestSize = 1000*1024^2)

# Functions ---------------------------------------------------------------

find_newest_file <-
   function(
      path
   ) {
      purrr::map_chr(
         path,
         ~fs::dir_info(.x) %>%
            dplyr::filter(modification_time == max(modification_time)) %>%
            dplyr::pull("path")
      )
   }

kickout <-
   function(
      list,
      allowed_ext = c("tdReport", "csv", "xlsx")
   ) {

      # This function removes any element from the list of input files
      # which does not have one of the allowed
      # extensions or which has "deprecated" in its name

      for (i in rev(seq_along(list))) {

         if (!(tools::file_ext(list[[i]]) %in% allowed_ext)) {

            list[[i]] <- NULL

         } else if (any(grep("deprecated", list[[i]], fixed = TRUE)) == TRUE) {

            list[[i]] <- NULL

         }
      }

      return(list)
   }

get_data_path <-
   function(
      filedir,
      filename,
      extension
   ) {

      filesindir <-
         fs::dir_ls(
            filedir, recurse = TRUE, type = "file",
            regexp = paste0("[.]", extension, "$")
         ) %>%
         purrr::as_vector()

      if (purrr::map(
         filename,
         ~stringr::str_detect(filesindir, .x)
      ) %>%
      purrr::map(any) %>%
      purrr::as_vector() %>%
      all() == FALSE) stop("One or more input files not found")

      if (filename %>%
          as.list() %>%
          purrr::map(~stringr::str_subset(filesindir, .x)) %>%
          purrr::map(~any(length(.x) > 1)) %>%
          unlist() %>%
          any() == TRUE) stop("One or more input files found in multiple locations")

      filelist <-
         filename %>%
         purrr::map_chr(~stringr::str_subset(filesindir, .x)) %>%
         as.list() %>%
         kickout()

      names(filelist) <- seq(1, length(filelist))

      return(filelist)

   }


# Server ------------------------------------------------------------------

shinyServer(
   function(input, output, session) {

      # Initial Setup -----------------------------------------------------------

      # Establish params to use for shinyFiles input (local only)

      volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())

      shinyFileChoose(
         input = input,
         "rawfile",
         session = session,
         roots = volumes
      )

      shinyFileChoose(
         input = input,
         "targetseqs",
         session = session,
         roots = volumes
      )

      disable("downloadReport")
      disable("updateplot")

      # Listeners ---------------------------------------------------------------

      listener_selectraw <-
         reactive(
            {
               list(
                  input$rawfile
               )
            }
         )

      listener_selectcsv <-
         reactive(
            {
               list(
                  input$targetseqs
               )
            }
         )

      # Reactives ---------------------------------------------------------------

      is_rawfile_valid <-
         reactive(
            {
               file_exists(rawfile_path()) &&
                  length(rawfile_path() > 0) &&
                  has_extension(rawfile_path(), "raw")
            }
         )

      is_csv_valid <-
         reactive(
            {
               file_exists(targetseqs_path()) &&
                  length(targetseqs_path() > 0) &&
                  has_extension(targetseqs_path(), "csv")
            }
         )

      rawfile_path <-
         reactive(
            {
               parseFilePaths(volumes, input$rawfile)$datapath
            }
         )

      rawfile_name <-
         reactive(
            {
               parseFilePaths(volumes, input$rawfile)$name
            }
         )

      targetseqs_path <-
         reactive(
            {
               parseFilePaths(volumes, input$targetseqs)$datapath
            }
         )

      targetseqs_name <-
         reactive(
            {
               parseFilePaths(volumes, input$targetseqs)$name
            }
         )

      activate_startbutton <-
         reactive(
            {
               if (is_rawfile_valid() && is_csv_valid()) {
                  enable("GEXstart")
               } else {
                  disable("GEXstart")
               }
            }
         )

      spectrum_names <-
         reactive(
            {
               validate(
                  need(input$plotchoice != "Generate XICs first", "")
               )
               names(Spectra()[[input$plotchoice]])
            }
         )


      output$spectrumchoice <-
         renderUI(
            {
               validate(
                  need(input$plotchoice != "Generate XICs first", "")
               )
               selectInput(
                  "spectrumchoice",
                  "Choose charge state",
                  choices = spectrum_names(),
                  selectize = FALSE,
                  size = 5
               )
            }
         )

      # Observers ---------------------------------------------------------------

      observeEvent(
         listener_selectraw(),
         {
            if(
               length(rawfile_name()) > 0 &&
               !is_rawfile_valid()
            ) {
               showModal(
                  modalDialog(
                     title = "Invalid file",
                     "You must upload a Thermo .raw file"
                  )
               )
            }

         }
      )

      observeEvent(
         listener_selectcsv(),
         {
            if(
               length(targetseqs_name()) > 0 &&
               !is_csv_valid()
            ) {
               showModal(
                  modalDialog(
                     title = "Invalid target sequences file",
                     "You must upload a .csv file"
                  )
               )
            }

         }
      )

      XICs <- reactiveVal()
      Spectra <- reactiveVal()
      MEXoutput <- reactiveVal()

      observeEvent(
         input$GEXstart,
         {
            disable("GEXstart")
            disable("rawfile")
            disable("targetseqs")
            disable("downloadReport")
            disable("updateplot")

            withProgress(
               message = "Analyzing raw file, please wait...",
               value = 0.33,
               style = getShinyOption("progress.style", default = "notification"),
               {
                  MEXout <-
                     make_every_XIC(
                        rawFileDir = fs::path_dir(rawfile_path()),
                        rawFileName = rawfile_name(),
                        targetSeqData = targetseqs_path(),
                        outputDir = tempdir(),
                        target_col_name = input$targetcolname,
                        target_sequence_col_name = input$targetseqname,
                        PTMname_col_name = input$ptmcolname,
                        PTMformula_col_name1 = input$ptmaddname,
                        PTMformula_col_name2 = input$ptmsubname,
                        use_depleted_isotopes = as.logical(input$usedepleted),
                        sample_n_pforms = input$sample_pforms,
                        mass_range = input$massrange,
                        target_charges =
                           c(input$charges[[1]]:input$charges[[2]]),
                        mz_range = input$mzrange,
                        abund_cutoff = input$abundcutoff
                     )

                  setProgress(value = 0.66)

                  MESout <-
                     make_every_spectrum(
                        MEXout,
                        mz_window = 3,
                        MSoutputWidth = 18,
                        MSoutputDPI = 200,
                        makePNG = FALSE,
                        return_timers = FALSE
                     )

               }
            )

            XICs(MEXout[[13]])
            MEXoutput(MEXout)
            Spectra(MESout[[1]])

            # MIXoutput[[13]] is the XIC plots object

            updateSelectInput(
               session,
               "plotchoice",
               choices = names(MEXout[[1]])
            )

            enable("GEXstart")
            enable("rawfile")
            enable("targetseqs")
            enable("downloadReport")
            enable("updateplot")

            setProgress(value = 1)

         }
      )

      output$outputPlot1 <-
         renderPlot(
            {
               validate(
                  need(input$plotchoice != "Generate XICs first", "")
               )
               XICs()[[input$plotchoice]] +
                  theme(text = element_text(size = 16))
            },
            width = 600,
            height= 400
         )

      output$outputPlot2 <-
         renderPlot(
            {
               validate(
                  need(input$plotchoice != "Generate XICs first", ""),
                  need(input$spectrumchoice != "Select an XIC", "")
               )
               Spectra()[[input$plotchoice]][[input$spectrumchoice]] +
                  theme(text = element_text(size = 16))
            },
            width = 600,
            height= 400
         )


      # Testers -----------------------------------------------------------------

      output$confirm <-
         renderText(
            {
               c(
                  is_rawfile_valid(),
                  is_csv_valid(),
                  input$usedepleted,
                  rawfile_path(),
                  targetseqs_path(),
                  spectrum_names()
               )
            }
         )



   }
)
