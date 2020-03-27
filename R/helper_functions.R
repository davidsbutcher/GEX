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
   
   make_XIC_plot = function(df, x, y, seq_name,  theme = NULL) {
      
      df %>% 
         ggplot2::ggplot(aes({{x}}, {{y}})) +
         ggplot2::geom_line() +
         ggplot2::labs(
            title = glue::glue("{seq_name}") 
         ) +
         theme
      
   }
   
   
   make_spectrum_top1 = 
      function(df, x, y, accession = "NA", scan_num = 0, charge = 0, xrange = c(0,0), theme  = NULL) {
         
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
               
               ymax <- 100
               
            } 
            
            if (all.equal(0, ymax) == TRUE) {
               
               ymax <- 100
               
            }
            
            scan_cap <- 
               df %>% 
               dplyr::pull({{scan_num}}) %>% 
               .[[1]]
            
         }
         
         df %>%
            dplyr::filter({{x}} >= xrange[[1]] & {{x}} <= xrange[[2]]) %>%
            ggplot2::ggplot(aes({{x}}, {{y}})) +
            ggplot2::geom_line() +
            ggplot2::geom_vline(
               ggplot2::aes(
                  xintercept = xrange[[1]]+((xrange[[2]]-xrange[[1]])/2),
                  color = "red", alpha = 0.5
               )
            ) +
            ggplot2::annotate(
               "text",
               x = xrange[[2]],
               y = ymax,
               label = glue::glue("{accession}\n Scan #{scan_cap}\n Charge +{charge}"),
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
            theme
         
      }