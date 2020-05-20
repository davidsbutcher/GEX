kickoutXIC <- function(list) {

   # This function removes any element from the list of input files
   # (from root/input) which does not have one of the allowed
   # extensions or which has "deprecated"

   for (i in rev(seq_along(list))) {

      if (is.null(magrittr::use_series(list[[i]], intensities)) == TRUE) {

         list[[i]] <- NULL

      }

   }

   return(list)

}

make_XIC_plot = function(df, x, y, PTMname, seq_name) {

   # x_var <-  rlang::enquo(x)
   # y_var <-  rlang::enquo(y)

   ggplot2::ggplot(df, ggplot2::aes({{x}}, {{y}})) +
      ggplot2::geom_line() +
      ggplot2::labs(
         title = glue::glue("{seq_name}")
      ) +
      ggplot2::geom_text(
         data = ~dplyr::filter(.x, int_sum == max(int_sum)),
         ggplot2::aes(
            x = times,
            y = int_sum,
            label = glue::glue(
               "Max TIC: {format(int_sum, scientific = TRUE, nsmall = 4, digits = 4)}\n RT@Max: {format(times, scientific = FALSE, nsmall = 4, digits = 3)}\n PTM: {stringr::str_wrap(PTMname, width = 10)}"),
            alpha = 0.5,
            size = 3
         ),
         vjust="inward",
         hjust="inward",
         nudge_x = 5
      ) +
      ggthemes::theme_clean(base_size = 10) +
      ggplot2::labs(
         x = "Retention Time (min)",
         y = "Total Ion Current"
      ) +
      ggplot2::guides(
         alpha = "none",
         size = "none"
      )

}


make_spectrum_top1 =
   function(
      df, x,
      y,
      accession = "NA",
      scan_num = 0,
      charge = 0,
      xrange = c(0,0),
      PTMname = NULL,
      vlines = NULL,
      theme  = NULL
   ) {
      {
         xmin <-
            df %>%
            dplyr::filter({{x}} == min({{x}})) %>%
            dplyr::pull({{x}})

         xmax <-
            df %>%
            dplyr::filter({{x}} == max({{x}})) %>%
            dplyr::pull({{x}})

         if (xrange[[1]] < xmin) {

            xrange[[1]] <- xmin

         }

         if (xrange[[2]] > xmax) {

            xrange[[2]] <- xmax

         }

         ymax <-
            df %>%
            dplyr::filter({{x}} >= xrange[[1]] & {{x}} <= xrange[[2]]) %>%
            dplyr::filter({{y}} == max({{y}})) %>%
            dplyr::pull({{y}}) %>%
            magrittr::extract(1)

         if (length(ymax) == 0) {

            ymax <- 10

         }

         if (all.equal(0, ymax) == TRUE) {

            ymax <- 10

         }

         scan_cap <-
            df %>%
            dplyr::pull({{scan_num}}) %>%
            .[[1]]

      }

      df %>%
         dplyr::filter({{x}} >= xrange[[1]] & {{x}} <= xrange[[2]]) %>%
         ggplot2::ggplot(ggplot2::aes({{x}}, {{y}})) +
         ggplot2::geom_line(
            ggplot2::aes(size = 0.5)
         ) +
         ggplot2::geom_vline(
            ggplot2::aes(
               xintercept = xrange[[1]]+((xrange[[2]]-xrange[[1]])/2),
               color = "red",
               alpha = 0.5,
               size = 0.5
            )
         ) +
         # ggplot2::geom_vline(
         #    data = vlines,
         #    ggplot2::aes(
         #       xintercept = `m/z`,
         #       color = "red",
         #       alpha = 0.5,
         #       linetype = "longdash"
         #    )
         # ) +
         ggplot2::annotate(
            "segment",
            x = vlines,
            xend = vlines,
            y = 0,
            yend = ymax,
            color = "red",
            alpha = 0.5,
            size = 0.5,
            linetype = "longdash"
         ) +
         # ggplot2::annotate(
         #    "vline",
         #    x = vlines,
         #    xintercept = vlines,
         #    color = "red",
         #    alpha = 0.5,
         #    linetype = "longdash"
         # ) +
         ggplot2::annotate(
            "text",
            x = xrange[[2]],
            y = ymax,
            label = glue::glue("{accession}\n Scan #{scan_cap}\n Charge +{charge}\n Max TIC: {format(ymax, scientific = TRUE, nsmall = 4, digits = 4)}\n PTM: {PTMname}"),
            vjust="inward",
            hjust="inward",
            size = 3,
            alpha = 0.5
         ) +
         ggplot2::lims(
            x = xrange,
            y = c(0, ymax)
         ) +
         ggplot2::guides(
            color = "none",
            size = "none",
            alpha = "none"
         ) +
         ggplot2::labs(
            x = "m/z",
            y = "Intensity"
         ) +
         ggplot2::scale_size_identity() +
         theme

   }

get_maxY_in_Xrange <-
   function(
      df, x, y, xrange = c(0,0)
   ) {

      df %>%
         dplyr::filter({{x}} >= xrange[[1]] & {{x}} <= xrange[[2]]) %>%
         dplyr::filter({{y}} == max({{y}})) %>%
         dplyr::pull({{y}})

   }

ggplot_list_checker <-
   function(
      ggplot_list
   ) {

      # Checks a list of ggplot objects for "blank" plots and removes them

      for (i in rev(seq_along(ggplot_list))) {

         if ((ggplot_list[[i]] %>% .$data %>% nrow()) == 0) {

            ggplot_list[[i]] <- NULL

         }

      }

      return(ggplot_list)

   }

tibble_list_checker <-
   function(
      tibble_list
   ) {

      # Checks a list of ggplot objects for "blank" plots and removes them

      for (i in rev(seq_along(tibble_list))) {

         if (tibble_list[[i]] %>% dplyr::pull(`m/z`) %>% length() == 0) {

            tibble_list[[i]] <- NULL

         }

      }

      return(tibble_list)

   }
