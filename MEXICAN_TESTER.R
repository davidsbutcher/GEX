# Function arguments for testing functions

rawFileDir = "D:/Z Drive Backup 20191213/21T data"
rawFileName = "DSB_201900906_M9-J-20190404_PEPPI_F03_01.raw"
targetSeqData = "C:/Users/ranar/Documents/meXICan_output/meXICan_spectrum_testset_short-01.csv"
outputDir = "C:/Users/ranar/Documents/meXICan_output"
target_col_name = c("UNIPROTKB", "PFR")
target_sequence_col_name = c("ProteoformSequence")
PTMname_col_name = c("PTMname")
PTMformula_col_name1 = c("FormulaToAdd")
PTMformula_col_name2 = c("FormulaToSubtract")
use_depleted_isotopes = FALSE
top_n_pforms = NULL
mass_range = c(10000,15000)
target_charges = c(1:30)
mz_range = c(600,2000)
XIC_tol = 25
use_IAA = FALSE
abund_cutoff = 5

# Running a test on a single raw file

MEXIC_TEST <-
   make_every_XIC(
      rawFileDir = "C:/Users/ranar/Documents/zdrive_local",
      rawFileName = "ZP_09032019_Norm_F3F4_A_01.raw",
      targetSeqData = "C:/Users/ranar/Documents/meXICan_output/meXICan_spectrum_testset_short-02.csv",
      outputDir = "C:/Users/ranar/Documents/meXICan_output",
      target_col_name = c("UNIPROTKB", "PFR"),
      target_sequence_col_name = c("ProteoformSequence"),
      PTMname_col_name = c("PTMname"),
      PTMformula_col_name1 = c("FormulaToAdd"),
      PTMformula_col_name2 = c("FormulaToSubtract"),
      use_depleted_isotopes = FALSE,
      top_n_pforms = NULL,
      mass_range = c(10000,20000),
      target_charges = c(1:40),
      mz_range = c(600,2000),
      XIC_tol = 10,
      use_IAA = FALSE,
      abund_cutoff = 5
   )

timer <-
   make_every_spectrum(
      MEXIC_TEST,
      mz_window = 3,
      makePNG = FALSE
   )

# Running a test on a single raw file 2

MEXIC_TEST <-
   make_every_XIC(
      rawFileDir = "C:/Users/ranar/Documents/zdrive_local",
      rawFileName = "ZP_09032019_Norm_F3F4_A_01.raw",
      targetSeqData = "C:/Users/ranar/Documents/meXICan_output/proteome625_sequences_only_trunc.csv",
      outputDir = "C:/Users/ranar/Documents/meXICan_output",
      target_col_name = c("UNIPROTKB"),
      target_sequence_col_name = c("ProteoformSequence"),
      PTMname_col_name = c("PTMname"),
      PTMformula_col_name1 = c("FormulaToAdd"),
      PTMformula_col_name2 = c("FormulaToSubtract"),
      use_depleted_isotopes = FALSE,
      top_n_pforms = NULL,
      mass_range = c(10000,20000),
      target_charges = c(1:40),
      mz_range = c(600,2000),
      XIC_tol = 10,
      use_IAA = FALSE,
      abund_cutoff = 5
   )

timer <-
   make_every_spectrum(
      MEXIC_TEST,
      mz_window = 3,
      makePNG = FALSE
   )

# Running a test on a single raw file on Anderson-T5810

MEXIC_TEST <-
   make_every_XIC(
      rawFileDir = "Z:/ICR/Zeljka Popovic/Zubarev Project Fall 2019 Iso Dep/Data Collected AugSept 2019/raw files",
      rawFileName = "ZP_09032019_Norm_F3F4_A_01.raw",
      targetSeqData = "ignore/proteome625_sequences_only.csv",
      outputDir = "C:/Users/dbutcher/Documents/meXICan_output",
      target_col_name = c("UNIPROTKB"),
      target_sequence_col_name = c("ProteoformSequence"),
      PTMname_col_name = c("PTMname"),
      PTMformula_col_name1 = c("FormulaToAdd"),
      PTMformula_col_name2 = c("FormulaToSubtract"),
      use_depleted_isotopes = FALSE,
      top_n_pforms = NULL,
      mass_range = c(10000,20000),
      target_charges = c(1:40),
      mz_range = c(600,2000),
      XIC_tol = 10,
      use_IAA = FALSE,
      abund_cutoff = 5
   )

timer <-
   make_every_spectrum(
      MEXIC_TEST,
      mz_window = 3,
      makePNG = FALSE
   )

RPushbullet::pbPost(type = "note", body = "mXIC Done")

# Running a test on a single ID raw file on Anderson-T5810

MEXIC_TEST <-
   make_every_XIC(
      rawFileDir = "Z:/ICR/Zeljka Popovic/Zubarev Project Fall 2019 Iso Dep/Data Collected AugSept 2019/raw files",
      rawFileName = "ZP_09042019_ID_F3F4_A_01.raw",
      targetSeqData = "ignore/proteome625_sequences_only.csv",
      outputDir = "C:/Users/dbutcher/Documents/meXICan_output",
      target_col_name = c("UNIPROTKB"),
      target_sequence_col_name = c("ProteoformSequence"),
      PTMname_col_name = c("PTMname"),
      PTMformula_col_name1 = c("FormulaToAdd"),
      PTMformula_col_name2 = c("FormulaToSubtract"),
      use_depleted_isotopes = TRUE,
      top_n_pforms = NULL,
      mass_range = c(10000,20000),
      target_charges = c(1:40),
      mz_range = c(600,2000),
      XIC_tol = 10,
      use_IAA = FALSE,
      abund_cutoff = 5
   )

timer <-
   make_every_spectrum(
      MEXIC_TEST,
      mz_window = 3,
      makePNG = FALSE
   )

RPushbullet::pbPost(type = "note", body = "mXIC ID Done")

# Running a test on a single raw file 4/27/20

rawFileDir = "C:/Users/ranar/Documents/zdrive_local"
rawFileName = "ZP_09032019_Norm_F3F4_A_01.raw"
targetSeqData = "ignore/proteome625_sequences_only_trunc.csv"
outputDir = "C:/Users/ranar/Documents/meXICan_output"
target_col_name = c("UNIPROTKB")
target_sequence_col_name = c("ProteoformSequence")
PTMname_col_name = c("PTMname")
PTMformula_col_name1 = c("FormulaToAdd")
PTMformula_col_name2 = c("FormulaToSubtract")
use_depleted_isotopes = FALSE
top_n_pforms = NULL
mass_range = c(10000,20000)
target_charges = c(1:40)
mz_range = c(600,2000)
XIC_tol = 10
use_IAA = FALSE
abund_cutoff = 5

MEXIC_TEST <-
   make_every_XIC(
      rawFileDir = "C:/Users/ranar/Documents/zdrive_local",
      rawFileName = "ZP_09032019_Norm_F3F4_A_01.raw",
      targetSeqData = "ignore/proteome625_sequences_only_trunc.csv",
      outputDir = "C:/Users/ranar/Documents/meXICan_output",
      target_col_name = c("UNIPROTKB"),
      target_sequence_col_name = c("ProteoformSequence"),
      PTMname_col_name = c("PTMname"),
      PTMformula_col_name1 = c("FormulaToAdd"),
      PTMformula_col_name2 = c("FormulaToSubtract"),
      use_depleted_isotopes = FALSE,
      top_n_pforms = NULL,
      mass_range = c(10000,20000),
      target_charges = c(1:40),
      mz_range = c(600,2000),
      XIC_tol = 10,
      use_IAA = FALSE,
      abund_cutoff = 5
   )

timer <-
   make_every_spectrum(
      MEXIC_TEST,
      mz_window = 3,
      makePNG = FALSE
   )
