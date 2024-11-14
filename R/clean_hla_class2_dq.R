#' @name clean_hla_class2_dq
#' @title Clean and standardize messy HLA Class II - DQ typing data
#' @description This function processes raw HLA Class II typing data,
#' removing inconsistent formatting and unnecessary symbols to ensure a standardized allele format.
#' It also imputes homozygosity at loci where one allele is missing.
#'
#'
#' @param data
#' Data frame containing HLA typing information.
#' @param var_1
#' HLA on allele 1.
#' @param var_2
#' HLA on allele 2.
#' @return Cleaned data frame with standardized HLA Class II data in original columns
#' @export
#'
#'
#' @import
#' tidyverse
#' utils
#' readr
#' tidyr
#' janitor
#'
#' @examples
#' dat <-  read.csv(system.file("extdata/example", "HLA_Clean_test.csv", package = "tidy_hla"))
#' re <- clean_hla_class1(dat, recip_dqa1_1, recip_dqa1_2)
#'
clean_hla_class2_dq <-
function(data, var_1, var_2) {
    data |>
        #* step 0: Rename Columns
        rename(var_1 = {{var_1}}, var_2 = {{var_2}}) |>
        #* step 1: Remove Invalid Suffixes & Replace Missing Values
        mutate(
            across(c(var_1, var_2),
                   ~ ifelse(substr(., nchar(.), nchar(.)) %in% c("N", "n", "p", "g","q"), "", .)),
            across(c(var_1, var_2), ~ if_else(is.na(.), "XX|xx", .))) |>
        #* step 2: Separate Values at Specific Delimiters
        separate(var_1, sep = "[({[]", into = c("var_1", "misc1"), fill = "right", extra = "drop") |>
        separate(var_2, sep = "[({[]", into = c("var_2", "misc2"), fill = "right", extra = "drop") |>
        #* step 3: Remove Unnecessary Columns
        select(-c(misc1, misc2)) |>
        mutate(
            #* step 4: Clean Specific Substrings
            across(c(var_1, var_2), ~ str_remove(., ":bgcwt|:bhjv|:ecaub")),
            #* step 5: Remove Leading Zeros and Characters
            across(c(var_1, var_2), ~ str_replace(., "^[0|*|.]", "")),
            across(c(var_1, var_2), ~ str_replace(., "^[.]", "")),
            #* step 6: Compute String Lengths
            across(c(var_1, var_2),
                   list(c = ~str_length(.)),
                   .names = "{.col}_c"),
            #* step 7: Shorten Long Alleles
            across(c(var_1, var_2),
                   ~ if_else(str_length(.) > 5, str_remove(., ":[0-9]{2}$"), .)),
            #* step 8: Add Missing Colon Separators
            across(c(var_1, var_2), ~ if_else(str_length(.) == 3,
                                              paste0(substr(., 1, nchar(.) - 2), ":", substr(., nchar(.) - 1, nchar(.))),
                                              .)),
            across(c(var_1, var_2), str_length, .names = "{.col}_c"),
            #* step 9: Standardize Length of Single-Digit Alleles
            var_1 = if_else(var_1_c == 1, paste0("0", var_1), var_1),
            var_2 = if_else(var_2_c == 1, paste0("0", var_2), var_2),
            across(c(var_1, var_2), str_length, .names = "{.col}_c"),
            #* step 10: Format Four-Digit Alleles
            across(c(var_1, var_2),
                   ~ if_else(str_length(.) == 4 & !str_detect(., "X|x"),
                             paste0("0", .),
                             .)),
            #* step 11: Fill Missing Values from Paired Column
            var_1 = if_else(str_detect(var_1, "X|x"), var_2, var_1),
            var_2 = if_else(str_detect(var_2, "X|x"), var_1, var_2),
            var_1 = if_else(str_detect(var_2, "X|x|NA_character_| ") & str_detect(var_1, "X|x|NA_character_| "), NA_character_, var_1),
            var_2 = if_else(str_detect(var_2, "X|x|NA_character_| ") & str_detect(var_1, "X|x|NA_character_| "), NA_character_, var_2)
        ) |>
        #* step 12: Remove Redundant Columns
        select(-c(var_1_c, var_2_c))  |>
        #* step 13: Fill Empty Cells Based on Paired Values
        mutate(
            var_1 = if_else(var_1 == "" & var_2 != "", var_2, var_1),
            var_2 = if_else(var_2 == "" & var_1 != "", var_1, var_2)
        ) |>
        #* step 14: Rename Columns to Original Names
        rename({{var_1}} := var_1, {{var_2}} := var_2)
}
