#' @name eval_mism
#' @title Evaluate HLA Mismatches
#' @description Calculates the number of mismatches between donor and recipient HLA alleles.
#' The function standardizes input values, aligns missing data,
#' and assesses mismatches for two sets of alleles per donor and recipient.
#'
#'
#' @param data A data frame containing donor and recipient HLA alleles.
#' @param don_1 Column name for the first donor allele.
#' @param don_2 Column name for the second donor allele.
#' @param recip_1 Column name for the first recipient allele.
#' @param recip_2 Column name for the second recipient allele.
#' @param hmz_cnt Default value is 1. Used to set a homozygosity counter (not used in calculations in this version).
#'
#' @return A modified data frame with mismatches calculated between donor and recipient alleles.
#' @export
#'
#'
#' @import
#' tidyverse
#' utils
#' readr
#' dplyr
#' stringr
#'
#' @examples
#' dat <-  read.csv(system.file("extdata/example", "HLA_Clean_test.csv", package = "tidy_hla"))
#' re <- eval_mism(dat, don_a_1, don_a_2, recip_a_1, recip_a_2)

eval_mism <-
function(data, don_1, don_2, recip_1, recip_2, hmz_cnt=1)
{
    data <- data |>

        # Step 1: Rename columns to simplify processing
        rename(don_1 = {{don_1}}, don_2 = {{don_2}},
               recip_1 = {{recip_1}}, recip_2 = {{recip_2}}) |>

        # Step 2: Convert all columns to character to avoid type mismatch issues
        mutate(across(everything(), as.character),
               hmz_cnt = 1) |>

        # Step 3: Set up for row-wise operations
        rowwise() |>

        # Step 4: Handle missing or empty values for allele columns
        mutate(
            across(c(don_1, don_2), ~ if_else(is.na(.)|. == "", NA, .)),
            across(c(recip_1, recip_2), ~ if_else(is.na(.)|. == "", NA, .)),
            across(c(recip_1, recip_2), ~ if_else(is.na(don_1) & is.na(don_2), NA_character_, .x)),
            across(c(don_1, don_2), ~ if_else(is.na(recip_1) & is.na(recip_2), NA_character_, .x)),

            # Step 5: Ensure donor and recipient alleles are filled where needed
            don_2 = if_else(!is.na(don_1) & is.na(don_2), don_1, don_2),
            don_1 = if_else(!is.na(don_2) & is.na(don_1),don_2,don_1),
            recip_2 = if_else(!is.na(recip_1) & is.na(recip_2), recip_1, recip_2),
            recip_1 = if_else(!is.na(recip_2) & is.na(recip_1), recip_2, recip_1),

            # Step 6: Extract the first two characters from alleles for comparison
            don_1_1 = str_sub(don_1,1,2),
            don_1_1 = if_else(
              substr(don_1_1, 1, 1) == "0",
              substr(don_1_1, 2, nchar(don_1_1)),
              don_1_1),
            recip_1_1 = str_sub(recip_1,1,2),
            recip_1_1 = if_else(
              substr(recip_1_1, 1, 1) == "0",
              substr(recip_1_1, 2, nchar(recip_1_1)),
              recip_1_1),
            don_2_1 = str_sub(don_2,1,2),
            don_2_1 = if_else(
              substr(don_2_1, 1, 1) == "0",
              substr(don_2_1, 2, nchar(don_2_1)),
              don_2_1),
            recip_2_1 = str_sub(recip_2,1,2),
            recip_2_1 = if_else(
              substr(recip_2_1, 1, 1) == "0",
              substr(recip_2_1, 2, nchar(recip_2_1)),
              recip_2_1),

            # Step 7: Calculate mismatch indicators

            # m1_1 =  if_else(don_1_1 %in% c(recip_1_1,recip_2_1), 0,1),
            # m2_1 =  if_else(don_2_1 %in% c(recip_1_1,recip_2_1), 0,1))
            m1_1 = case_when(
              don_1_1 != don_2_1 & don_1_1 == recip_1_1 | don_1_1 != don_2_1 & don_1_1 == recip_2_1  ~ 0,
              don_1_1 != don_2_1 & don_1_1 != recip_1_1 | don_1_1 != don_2_1 & don_1_1 != recip_2_1  ~ 1,
              don_1_1 == don_2_1 & don_1_1 != recip_1_1 & don_1_1 == don_2_1 & don_1_1 != recip_2_1  ~ 1,
              don_1_1 == don_2_1 & (don_1_1 %in% c(recip_1_1,recip_2_1))  ~ 0,
              TRUE ~ 1
            ),
            m2_1 = case_when(
              don_2_1 != don_1_1 & don_2_1 == recip_1_1 | don_2_1 != don_1_1 & don_2_1 == recip_2_1  ~ 0,
              don_2_1 != don_1_1 & don_2_1 != recip_1_1 | don_2_1 != don_1_1 & don_2_1 != recip_2_1  ~ 1,
              don_2_1 == don_1_1 & don_2_1 %in% c(recip_1_1,recip_2_1)  ~0,
              TRUE ~ 1
            ),


            # Step 8: Summing mismatch indicators to get mismatch count
            mism_cnt = if_else(
                is.na(m1_1) & is.na(m2_1),
                NA_real_, # If both are NA, return NA
                sum(c(m1_1, m2_1), na.rm = TRUE)
            )
        ) |>

        # Step 9: Final adjustments to ensure accurate mismatch counts
        mutate(
                mism_cnt =  if_else(
                is.na(don_1_1) & is.na(don_2_1) & is.na(recip_1_1) & is.na(recip_2_1),
                NA_real_, mism_cnt)) |>

        # Step 10: Remove temporary columns and restore original column names
        select(-c(m1_1, m2_1, recip_1_1, recip_2_1, don_1_1, don_2_1, hmz_cnt)) |>
        rename({{don_1}} := don_1 , {{don_2}} := don_2,
               {{recip_1}} := recip_1, {{recip_2}} := recip_2) |>
        ungroup()

}
