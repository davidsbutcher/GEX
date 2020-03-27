rawFileDir = "D:/Z Drive Backup 20191213/21T data"
rawFileName = "DSB_201900906_M9-J-20190404_PEPPI_F03_01.raw"
target_seq_file = "C:/Users/ranar/Documents/meXICan_output/meXICan_spectrum_testset.csv"
outputDir = "C:/Users/ranar/Documents/meXICan_output"
target_col_name = c("UNIPROTKB", "PFR")
target_sequence_col_name = c("ProteoformSequence")
use_depleted_isotopes = FALSE
top_n_pforms = NULL
target_charges = c(1:30)
mz_range = c(600,2000)
XIC_tol = 25
use_IAA = FALSE
abund_cutoff = 5

make_every_XIC(
   rawFileDir = "D:/Z Drive Backup 20191213/21T data",
   rawFileName = "DSB_201900906_M9-J-20190404_PEPPI_F03_01.raw",
   target_seq_file = "C:/Users/ranar/Documents/meXICan_output/meXICan_spectrum_testset.csv",
   outputDir = "C:/Users/ranar/Documents/meXICan_output",
   target_col_name = c("UNIPROTKB", "PFR"),
   target_sequence_col_name = c("ProteoformSequence"),
   use_depleted_isotopes = FALSE,
   top_n_pforms = NULL,
   target_charges = c(1:30),
   mz_range = c(600,2000),
   XIC_tol = 25,
   use_IAA = FALSE,
   abund_cutoff = 5
   ) 
   