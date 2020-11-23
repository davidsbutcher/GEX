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
            size = 2
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
      df,
      x,
      y,
      noise,
      accession = "NA",
      scan_num = 0,
      charge = 0,
      xrange = c(0,0),
      PTMname = NULL,
      vlines = NULL,
      chargestateTIC = NULL,
      cosine_sims = NULL,
      mz_all_abund = NULL,
      iso_abund_theoretical_scaled = NULL,
      isotopologue_window = NULL,
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

         # ymax <-
         #    chargestateTIC

         if (length(ymax) > 1) {

            ymax <- ymax[[1]]

         }

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

         mean_noise <-
            df %>%
            dplyr::filter({{x}} >= xrange[[1]] & {{x}} <= xrange[[2]]) %>%
            dplyr::pull({{noise}}) %>%
            mean()

      }

      theo_points <-
         iso_abund_theoretical_scaled %>%
         purrr::set_names(mz_all_abund) %>%
         tibble::enframe(name = "mz", value = "intensity") %>%
         dplyr::mutate(mz = as.double(mz))

      df %>%
         dplyr::filter({{x}} >= xrange[[1]] & {{x}} <= xrange[[2]]) %>%
         ggplot2::ggplot(ggplot2::aes({{x}}, {{y}})) +
         ggplot2::geom_line(
            ggplot2::aes(size = 0.25)
         ) +
         ggplot2::geom_point(
            data = theo_points,
            mapping = ggplot2::aes(x = mz, y = intensity),
            color = "red",
            shape = 1,
            size = 2.5,
            alpha = 0.5
         ) +
         # ggplot2::annotate(
         #    "tile",
         #    x = mean(xrange),
         #    y = 0,
         #    width = isotopologue_window*10,
         #    height = ymax,
         #    fill = "red4",
         #    alpha = 0.5,
         #    linetype = "blank"
         # ) +
         ggplot2::annotate(
            "rect",
            xmin = mean(xrange)-(isotopologue_window/2),
            xmax = mean(xrange)+(isotopologue_window/2),
            ymin = 0,
            ymax = ymax,
            fill = "blue",
            alpha = 0.5,
            linetype = "blank"
         ) +
         ggplot2::annotate(
            "rect",
            xmin = vlines-(isotopologue_window/2),
            xmax = vlines+(isotopologue_window/2),
            ymin = 0,
            ymax = ymax,
            fill = "red",
            alpha = 0.5,
            linetype = "blank"
         ) +
         # ggplot2::annotate(
         #    "tile",
         #    x = vlines,
         #    y = rep(0, length(vlines)),
         #    width = isotopologue_window*10,
         #    height = ymax,
         #    fill = "red",
         #    alpha = 0.5,
         #    linetype = "blank"
         # ) +
         ggplot2::geom_hline(
            ggplot2::aes(
               yintercept = noise,
               alpha = 0.5,
               size = 0.25
            ),
            color = "blue"
         ) +
         # ggplot2::annotate(
         #    "segment",
         #    x = vlines,
         #    xend = vlines,
         #    y = 0,
         #    yend = ymax,
         #    color = "red",
         #    alpha = 0.35,
         #    size = 0.25,
         #    linetype = "longdash"
         # ) +
      ggplot2::annotate(
         "text",
         x = xrange[[2]],
         y = ymax,
         label = glue::glue("{accession}\n Scan #{scan_cap}\n Charge +{charge}\n Theo. Max TIC: {format(chargestateTIC, scientific = TRUE, nsmall = 3, digits = 3)}\n Est. S/N: {round((chargestateTIC/mean_noise), digits = 0)}\n Cos. Sim.: {round((cosine_sims), digits = 3)}\nPTM: {PTMname}"),
         vjust="inward",
         hjust="inward",
         size = 2,
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

get_maxY_in_Xrange_vector <-
   function(
      df, x, y, mz = 0, mz_window, mz_window_scaling
   ) {

      out <- vector(mode = "numeric", length = length(mz))

      for (i in seq_along(mz)) {

         xrange <-
            c(
               mz[[i]] - (mz_window*mz_window_scaling),
               mz[[i]] + (mz_window*mz_window_scaling)
            )

         out[[i]] <-
            df %>%
            dplyr::filter({{x}} >= xrange[[1]] & {{x}} <= xrange[[2]]) %>%
            dplyr::filter({{y}} == max({{y}})) %>%
            dplyr::pull({{y}}) %>%
            {if (length(.) == 0) 0 else .} %>%
            {if (length(.) > 1) 0 else .}

      }

      return(out)

   }

get_maxX_in_Xrange_vector <-
   function(
      df, x, y, mz = 0, mz_window, mz_window_scaling
   ) {

      out <- vector(mode = "numeric", length = length(mz))

      for (i in seq_along(mz)) {

         xrange <-
            c(
               mz[[i]] - (mz_window*mz_window_scaling),
               mz[[i]] + (mz_window*mz_window_scaling)
            )

         out[[i]] <-
            df %>%
            dplyr::filter({{x}} >= xrange[[1]] & {{x}} <= xrange[[2]]) %>%
            dplyr::filter({{y}} == max({{y}})) %>%
            dplyr::pull({{x}}) %>%
            {if (length(.) == 0) 0 else .} %>%
            {if (length(.) > 1) 0 else .}

      }

      return(out)

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

fix_list_length <-
   function(
      list, desired_length, fill = c(0)
   ) {

      if (length(list) < desired_length) {

         while (length(list) < desired_length) {
            list[[length(list)+1]] <- fill
         }

      } else if (length(list) > desired_length) {

         while (length(list) > desired_length) {
            list[[length(list)]] <- NULL
         }

      }

      return(list)

   }

#' calculate_cosine_similarity
#'
#' Calculate cosine similarity of two vectors. Function adapted from the
#' MicroRaman package.
#'
#' @param a
#' @param b
#'
#' @importFrom magrittr %>%

calculate_cosine_similarity <-
   function(
      a,
      b
   ) {

      numerator <- sum(a * b)

      quada <-
         sapply(
            a,
            FUN = function(x) x ^ 2
         ) %>%
         sum()

      quadb <-
         sapply(
            b,
            FUN = function(x) x ^ 2
         ) %>%
         sum()

      denominator <- sqrt(quada * quadb)

      costh <- numerator / denominator

      return(costh)
   }
