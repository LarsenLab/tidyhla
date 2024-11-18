#' @name clean_hla_class1
#' @title Clean and standardize messy HLA Class I typing data
#' @description This function processes raw HLA Class I typing data,
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
#' @return Cleaned data frame with standardized HLA Class I data in original columns
#' @export
#'
#'
#' @import
#' tidyverse
#' utils
#' readr
#' tidyr
#' janitor
#' dplyr
#'
#'
#'
clean_hla_class1 <-
  function(data, var_1, var_2) {
    data |>
      #* step 0: Rename Columns
      rename(var_1 = {{var_1}}, var_2 = {{var_2}}) |>
      mutate(#* step 1: Remove Invalid Characters, Leading Zeros, Asterisks, and Dots
        across(c(var_1, var_2), ~ if_else(
          substr(., nchar(.), nchar(.)) %in% c("N", "n", "p", "g", "q"), "", .
        )), across(c(var_1, var_2), ~ str_replace_all(., "^[0|*|.]+", ""))) |>
      #* step 2: Separate Extra Information
      tidyr::separate(
        var_1,
        sep = "[({[]",
        into = c("var_1", "misc1"),
        fill = "right",
        extra = "drop"
      ) |>
      tidyr::separate(
        var_2,
        sep = "[({[]",
        into = c("var_2", "misc2"),
        fill = "right",
        extra = "drop"
      ) |>
      mutate(
        #* step 3: Determine Length of Cleaned Alleles
        across(c(var_1, var_2), list(c = ~ str_length(.)), .names = "{.col}_c"),
        #* step 4: Adjust Long Allele Formats
        var_1 = case_when(
          var_1_c > 5 ~ str_remove(var_1, ":[0-9]{2}$"),
          var_1_c == 1 ~ paste0("0", var_1),
          TRUE ~ var_1
        ),
        var_2 = case_when(
          var_2_c > 5 ~ str_remove(var_2, ":[0-9]{2}$"),
          var_2_c == 1 ~ paste0("0", var_2),
          TRUE ~ var_2
        ),
        #* step 5: Replace Incomplete Alleles
        across(c(var_1, var_2), ~ if_else(
          str_detect(.x, "x"), coalesce(var_2, var_1), .x
        )),
        across(c(var_1, var_2), list(c = ~ str_length(.)), .names = "{.col}_c"),
        #* step 6: Add Missing Colon for Four-Character Alleles
        var_1 = if_else(var_1_c == 4 &
                          grepl("^[1-9]:", var_1), paste0("0", var_1), var_1),
        var_2 = if_else(var_2_c == 4 &
                          grepl("^[1-9]:", var_2), paste0("0", var_2), var_2),
        var_2 = if_else(var_2_c == 4, paste0(
          substr(var_2, 1, nchar(var_2) - 2), ":", substr(var_2, nchar(var_2) - 1, nchar(var_2))
        ), var_2),
        var_1 = if_else(var_1_c == 4, paste0(
          substr(var_1, 1, nchar(var_1) - 2), ":", substr(var_1, nchar(var_1) - 1, nchar(var_1))
        ), var_1),
        #* step 7: Remove Unwanted Suffixes
        across(c(var_1, var_2), ~ str_replace_all(., "[:A-Z|*]+$", "")),
        #* step 8: Synchronize Missing Alleles
        var_1 = if_else(var_1 == "" &
                          var_2 != "", var_2, var_1),
        var_2 = if_else(var_2 == "" &
                          var_1 != "", var_1, var_2),
        #* step 9: Replace Double Colons
        across(c(var_1, var_2), ~ str_replace(.x, "::", ":")),
        #* step 10: Handle NA Values and Inconsistencies
        var_1 = if_else(str_detect(var_1, "X|x|NA_character_| "), var_2, var_1),
        var_2 = if_else(str_detect(var_2, "X|x|NA_character_| "), var_1, var_2),
        var_1 = if_else(
          str_detect(var_2, "X|x|NA_character_| ") &
            str_detect(var_1, "X|x|NA_character_| "),
          NA_character_,
          var_1
        ),
        var_2 = if_else(
          str_detect(var_2, "X|x|NA_character_| ") &
            str_detect(var_1, "X|x|NA_character_| "),
          NA_character_,
          var_2
        )
      ) |>
      #* step 11: Finalize Column Names
      rename({{var_1}} := var_1, {{var_2}} := var_2) |>
      #* step 12: Remove Extra Columns
      select(-c("misc1", "misc2"))
  }
