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
               "Max TIC: {format(int_sum, scientific = TRUE, nsmall = 4, digits = 4)}\n RT@Max: {format(times, scientific = FALSE, nsmall = 4, digits = 3)}\n PTM: {stringr::str_wrap(PTMname, width = 30)}"),
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
      ) +
      ggplot2::scale_size_identity()


}

make_XIC_plot_MS2 =
   function(
      df,
      x,
      y,
      seq_name,
      ion_name,
      ion_charge,
      RT_of_maxTIC
   ) {

      ggplot2::ggplot(df, ggplot2::aes({{x}}, {{y}})) +
         ggplot2::geom_line() +
         ggplot2::labs(
            title = glue::glue("{seq_name}")
         ) +
         ggplot2::geom_text(
            data = ~dplyr::filter(.x, {{y}} == max({{y}})),
            ggplot2::aes(
               x = {{x}},
               y = {{y}},
               label = glue::glue(
                  "{ion_name}\nCharge {stringr::str_wrap(toString(paste0(ion_charge,'+')), width = 40)}\nMax TIC: {format(int_sum, scientific = TRUE, nsmall = 4, digits = 4)}\n RT@Max: {format(times, scientific = FALSE, nsmall = 4, digits = 3)}"
               ),
               alpha = 0.35,
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
         ) +
         ggplot2::scale_size_identity()

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
      score_mfa = NULL,
      mz_theo = NULL,
      int_theo = NULL,
      mz_obs = NULL,
      int_obs = NULL,
      mma = NULL,
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

         # ymax <-
         #    df %>%
         #    dplyr::filter({{x}} >= xrange[[1]] & {{x}} <= xrange[[2]]) %>%
         #    dplyr::filter({{y}} == max({{y}})) %>%
         #    dplyr::pull({{y}}) %>%
         #    magrittr::extract(1)

         ymax <-
            max(int_obs) + 0.05*max(int_obs)

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
         int_theo %>%
         purrr::set_names(mz_theo) %>%
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
         ggplot2::annotate(
            "rect",
            xmin = mean(xrange)-(isotopologue_window/2),
            xmax = mean(xrange)+(isotopologue_window/2),
            ymin = 0,
            ymax = max(theo_points$intensity),
            fill = "yellow",
            alpha = 0.05,
            linetype = "blank"
         ) +
         ggplot2::annotate(
            "rect",
            xmin = mz_theo-(isotopologue_window/2),
            xmax = mz_theo+(isotopologue_window/2),
            ymin = 0,
            ymax = theo_points$intensity,
            fill = "red",
            alpha = 0.25,
            linetype = "blank"
         ) +
         ggplot2::annotate(
            "text",
            x = mz_theo,
            y = int_theo+(int_theo*0.05),
            label = signif(as.numeric(mma), digits = 2),
            vjust="inward",
            color = "red",
            size = 1.5,
            alpha = 0.65
         ) +
         ggplot2::annotate(
            "segment",
            x = mz_obs[1:length(mz_obs)-1],
            xend = mz_obs[2:length(mz_obs)],
            y = int_obs[1:length(int_obs)-1]/2,
            yend = int_obs[1:length(int_obs)-1]/2,
            color = "blue",
            size = 0.75,
            alpha = 0.5
         ) +
         ggplot2::annotate(
            "text",
            x = mz_obs[1:length(mz_obs)-1]+(diff(mz_obs)/2),
            y = (int_obs[1:length(int_obs)-1]/2)+(int_obs[1:length(int_obs)-1]*0.03),
            label = round(diff(mz_obs)*1000),
            vjust="inward",
            size = 1.25,
            alpha = 0.5,
            color = "blue"
         ) +
         ggplot2::geom_hline(
            ggplot2::aes(
               yintercept = noise,
               alpha = 0.5,
               size = 0.25
            ),
            color = "blue"
         ) +
         ggplot2::annotate(
            "text",
            x = xrange[[2]],
            y = ymax,
            label = glue::glue("{accession}\n Scan #{scan_cap}\n Charge +{charge}\n THAI TIC: {format(chargestateTIC, scientific = TRUE, nsmall = 3, digits = 3)}\n Est. S/N: {round((chargestateTIC/mean_noise), digits = 0)}\n Cos. Sim.: {round((cosine_sims), digits = 3)}\n ScoreMFA: {round((score_mfa), digits = 3)}\nPTM: {stringr::str_wrap(PTMname, width = 30)}"),
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

make_spectrum_MS2 =
   function(
      df,
      name = "",
      max_mz = NULL,
      xrange = c(0,0),
      scan = NULL,
      ion = NULL,
      charge = NULL,
      isotopologue_window = NULL,
      ScoreMFA = NULL,
      cosine_sim = NULL,
      mz_theo = NULL,
      int_theo = NULL,
      sn_estimate = NULL,
      mma = NULL,
      mz_obs = NULL,
      int_obs = NULL,
      noise = NULL,
      theme  = NULL
   ) {

      if (is.null(max_mz)) {

         empty_plot <- ggplot2::ggplot() + ggplot2::theme_void()

         return(empty_plot)

      }

      {

         if (is.null(xrange) == FALSE & length(xrange) > 1) {

            xmax <- xrange[[2]]

         } else if (length(xrange) < 2) {

            xmax <- max_mz

         } else {

            xmax <- max_mz

         }

         # Calculate Ymax

         ymax <-
            max(int_theo)+(0.1*max(int_theo))


         if (length(ymax) > 1) {

            ymax <- ymax[[1]]

         } else if (length(ymax) == 0) {

            ymax <- 10

         } else if (all.equal(0, ymax) == TRUE) {

            ymax <- 10

         }

         if (is.null(ScoreMFA)) {

            ScoreMFA <- 0

         }

         if (is.null(cosine_sim)) {

            cosine_sim <- 0

         }

      }

      ggplot2::ggplot(df, ggplot2::aes(mz, intensity)) +
         ggplot2::geom_line(
            ggplot2::aes(size = 0.25)
         ) +
         ggplot2::annotate(
            "text",
            x = xmax,
            y = ymax,
            label = glue::glue("{name}\n{ion}\nScan #{scan}\nCharge {charge}+\nS/N estimate: {signif(max(sn_estimate), digits = 2)}\nCosine Sim: {signif(cosine_sim, digits = 3)}\nScoreMFA: {signif(ScoreMFA, digits = 3)}"),
            vjust="inward",
            hjust="inward",
            size = 2,
            alpha = 0.5
         ) +
         ggplot2::annotate(
            "point",
            x = mz_theo,
            y = int_theo,
            color = "red",
            shape = 1,
            size = 2.5,
            alpha = 0.5
         ) +
         ggplot2::annotate(
            "rect",
            xmin = mz_theo - (isotopologue_window/2),
            xmax = mz_theo + (isotopologue_window/2),
            ymin = 0,
            ymax = int_theo,
            fill = "red",
            alpha = 0.25,
            linetype = "blank"
         ) +
         ggplot2::annotate(
            "text",
            x = mz_theo,
            y = int_theo+(int_theo*0.05),
            label = signif(mma, digits = 2),
            vjust="inward",
            color = "red",
            size = 1.5,
            alpha = 0.65
         ) +
         ggplot2::annotate(
            "segment",
            x = mz_obs[1:length(mz_obs)-1],
            xend = mz_obs[2:length(mz_obs)],
            y = int_obs[1:length(int_obs)-1]/2,
            yend = int_obs[1:length(int_obs)-1]/2,
            color = "blue",
            size = 0.75,
            alpha = 0.5
         ) +
         ggplot2::annotate(
            "text",
            x = mz_obs[1:length(mz_obs)-1]+(diff(mz_obs)/2),
            y = (int_obs[1:length(int_obs)-1]/2)+(int_obs[1:length(int_obs)-1]*0.03),
            label = round(diff(mz_obs)*1000),
            vjust="inward",
            size = 1.25,
            alpha = 0.5,
            color = "blue"
         ) +
         ggplot2::geom_hline(
            ggplot2::aes(
               yintercept = noise,
               alpha = 0.5,
               size = 0.25
            ),
            color = "blue"
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
      df, x, y, mz = 0, res_power, isotopologueWinMultiplier = 2
   ) {

      out <- vector(mode = "numeric", length = length(mz))

      for (i in seq_along(mz)) {

         xrange <-
            c(
               mz[[i]] - ((mz[[i]]/res_power)*isotopologueWinMultiplier)/2,
               mz[[i]] + ((mz[[i]]/res_power)*isotopologueWinMultiplier)/2
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
      df, x, y, mz = 0, res_power, isotopologueWinMultiplier = 2
   ) {

      out <- vector(mode = "numeric", length = length(mz))

      for (i in seq_along(mz)) {

         xrange <-
            c(
               mz[[i]] - ((mz[[i]]/res_power)*isotopologueWinMultiplier)/2,
               mz[[i]] + ((mz[[i]]/res_power)*isotopologueWinMultiplier)/2
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

get_rp_in_Xrange_vector <-
   function(
      df, x, y, mz = 0, res_power, isotopologueWinMultiplier = 2
   ) {

      out <- vector(mode = "numeric", length = length(mz))

      for (i in seq_along(mz)) {

         xrange <-
            c(
               mz[[i]] - ((mz[[i]]/res_power)*isotopologueWinMultiplier)/2,
               mz[[i]] + ((mz[[i]]/res_power)*isotopologueWinMultiplier)/2
            )

         df_trunc <-
            df %>%
            dplyr::filter({{x}} >= xrange[[1]] & {{x}} <= xrange[[2]])

         # Get halfmax for the single isotopologue peak

         halfmax <-
            df_trunc %>%
            dplyr::filter({{y}} == max({{y}})) %>%
            dplyr::pull({{y}}) %>%
            {if (length(.) == 0) 0 else .} %>%
            {if (length(.) > 1) 0 else .} %>%
            {./2}

         maxX <-
            df_trunc %>%
            dplyr::filter({{y}} == max({{y}})) %>%
            dplyr::pull({{x}})

         # Split the spectrum into two halves, left and right of the max,
         # both including max

         df_trunc_left <-
            df_trunc %>%
            dplyr::filter({{x}} <= maxX)

         df_trunc_right <-
            df_trunc %>%
            dplyr::filter({{x}} >= maxX)

         # Find the values closest to half max, then perform linear interpolation
         # to approximate the m/z value at the halfmax intensity

         df_trunc_left2 <-
            df_trunc_left %>%
            dplyr::mutate(delta_max = abs(intensity - halfmax)) %>%
            dplyr::arrange(delta_max) %>%
            dplyr::slice_head(n = 2)


         df_trunc_right2 <-
            df_trunc_right %>%
            dplyr::mutate(delta_max = abs(intensity - halfmax)) %>%
            dplyr::arrange(delta_max) %>%
            dplyr::slice_head(n = 2)

         # Perform linear interpolation safely

         possibly_approx <-
            purrr::possibly(approx, otherwise = list(y = Inf))

         deltaX_left <-
            possibly_approx(
               y = df_trunc_left2$mz,
               x = df_trunc_left2$intensity,
               xout = halfmax
            ) %>%
            .$y

         deltaX_right <-
            possibly_approx(
               y = df_trunc_right2$mz,
               x = df_trunc_right2$intensity,
               xout = halfmax
            ) %>%
            .$y

         deltaX <-  deltaX_right - deltaX_left

         # Calculate estimated resolving power for the isotopologue peak

         result <-
            maxX/deltaX

         if (length(result) == 0) {
            out[[i]] <- 0
         } else {
            out[[i]] <- result[[1]]
         }

      }

      out[is.na(out)] <- 0
      out[is.nan(out)] <- 0
      out[is.infinite(out)] <- 0

      return(out)

   }

ggplot_list_checker <-
   function(
      ggplot_list
   ) {

      # Checks a list of ggplot objects for "blank" plots and removes them

      for (i in rev(seq_along(ggplot_list))) {

         if (length(ggplot_list[[i]]$data) == 0) {

            ggplot_list[[i]] <- NULL

         } else if (nrow(ggplot_list[[i]]$data) == 0) {

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

calculate_cosine_similarity <-
   function(
      a,
      b,
      rm_zero = TRUE
   ) {

      if (rm_zero == TRUE) {

         length_a_untrunc <- length(a)

         # If any values are 0, remove them from calculation

         b <- b[a !=0]
         a <- a[a !=0]

         # Check to see if there are any values left

         if (length(a) == 0) return(0)

         # Check to see if there is only one value left

         if (length(a) == 1) return(0)

         # Check to see if the number of zero values is less than half of the total

         if (length(a) < length_a_untrunc/2) return(0)

      }

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

ovln <- function(theo_mz, exp_mz, rp, rp_mult = 6){
   sd <- exp_mz / (rp_mult * rp)
   Zn <- -abs(theo_mz - exp_mz) / (2 * sd)
   return(2*pnorm(Zn, mean=0, sd=1))
}

ovlhn <- function(theo_mz, exp_mz, theo_sn, exp_sn, rp, rp_mult = 6){

   # Added by DSB to handle bug with missing values in exp_sn

   if (is.na(theo_sn) | is.na(exp_sn)) return(0)
   if (length(theo_sn) != length(exp_sn)) return(0)
   if (is.null(theo_sn) | is.null(exp_sn)) return(0)

   if(theo_sn==exp_sn){
      return(ovln(theo_mz,exp_mz, rp, rp_mult = rp_mult))
   } else {
      sd <- exp_mz / (rp_mult * rp)

      mu1t <- theo_mz - trunc(exp_mz)
      mu2t <- exp_mz - trunc(exp_mz)

      Z <- (mu1t + mu2t) / 2 + ((sd * sd) * log(theo_sn / exp_sn)) / (mu2t - mu1t)
      #Note log=ln

      muy <- min(mu1t, mu2t)
      mux <- max(mu1t, mu2t)

      if (muy == mu1t) {
         sn_muy <- theo_sn
         sn_mux <- exp_sn
      } else{
         sn_muy <- exp_sn
         sn_mux <- theo_sn
      }

      #Overlap of Normals without h
      ovln_sem_sn = pnorm(Z, mux, sd) + 1 - pnorm(Z, muy, sd)

      Iarea <- sn_mux*pnorm(Z, mux, sd) + sn_muy*(1 - pnorm(Z, muy, sd))

      Tarea <- sn_mux*(1 - pnorm(Z, mux, sd)) +  sn_muy*pnorm(Z, muy, sd)

      return(Iarea/Tarea)
   }

}

ScoreMFA <- function(vexp_mz, vtheo_mz, vrp, vexp_sn, vtheo_sn, rp_mult = 6, rm_zero = TRUE) {

   # If experimental and theoretical S/N vectors are not the same length,
   # return a score of 0

   if (length(vexp_sn) != length(vtheo_sn)) return(0)

   if (rm_zero == TRUE) {

      # If any values for experimental m/z are 0, remove it from calculation

      vtheo_mz <- vtheo_mz[vexp_mz !=0]
      vrp <- vrp[vexp_mz !=0]
      vexp_sn <- vexp_sn[vexp_mz !=0]
      vtheo_sn <- vtheo_sn[vexp_mz !=0]
      vexp_mz <- vexp_mz[vexp_mz !=0]

      # Check to see if there are any values left

      if (length(vexp_mz) == 0) return(0)

   }

   vsd <- vexp_mz / (rp_mult*vrp)
   ntheo <- length(vtheo_mz)
   nexp <- length(vexp_mz)

   if(nexp < ntheo){
      vtheo_mz <- vtheo_mz[1:nexp]
      vtheo_sn <- vtheo_sn[1:nexp]
   }
   vovlhn <- numeric(nexp)
   for(i in 1:ntheo){
      vovlhn[i] <- ovlhn(theo_mz=vtheo_mz[i], exp_mz=vexp_mz[i], theo_sn=vtheo_sn[i], exp_sn=vexp_sn[i], rp=vrp[i], rp_mult = rp_mult)
   }
   Score <- sum(vtheo_sn*vovlhn)/sum(vtheo_sn)
   return(Score)

}

calculate_mma_ppm <-
   function(
      obs_mass,
      theo_mass
   ) {

      if (all(obs_mass == 0) | all(theo_mass == 0)) return(rep(NA, length(obs_mass)))

      mma <- ((obs_mass - theo_mass)/theo_mass) * 1E6

      # mma[mma > 100 | mma < -100] <- NA

      purrr::map_dbl(
         mma,
         ~{if (.x > 100 | .x < -100) NA else .}
      )

   }
