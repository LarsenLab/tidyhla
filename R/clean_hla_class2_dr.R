#' @name clean_hla_class2_dr
#' @title Clean and standardize messy HLA Class II - DR typing data
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
#'
#' @examples
#' dat <-  read.csv(system.file("extdata/example", "HLA_Clean_test.csv", package = "tidy_hla"))
#' re <- clean_hla_class1(dat, recip_dra1_1, recip_dra1_2)

clean_hla_class2_dr <- function(data, var_1, var_2) { 
    data |>
        #* step 0: Rename Columns
        rename(var_1 = {{var_1}}, var_2 = {{var_2}}) |>
        #* step 1: Remove Certain Suffixes & Replace Missing Values
        mutate(
            across(c(var_1, var_2), 
                   ~ ifelse(substr(., nchar(.), nchar(.)) %in% c("N", "n", "p", "g","q"), "", .)),
            across(c(var_1, var_2), ~if_else(is.na(.), "XX|xx", .))
        ) |>
        #* step 2: Separate Extra Information
        separate(var_1, sep = "[({[]", into = c("var_1", "misc1"), fill = "right", extra = "drop") |>
        separate(var_2, sep = "[({[]", into = c("var_2", "misc2"), fill = "right", extra = "drop") |>
        select(-c(misc1, misc2)) |>
        #* step 3: Remove Non-Alphanumeric Characters
        mutate(
            across(c(var_1, var_2), ~ str_remove(., ":bgcwt|:bhjv|:ecaub")),
            across(c(var_1, var_2), ~ str_replace(., "^[0|*|.]", "")),
            across(c(var_1, var_2), ~ str_replace(., "^[.]", "")),
            #* step 4: Calculate String Length
            across(c(var_1, var_2), 
                   list(c = ~str_length(.)), 
                   .names = "{.col}_c"),
            #* step 5: Remove Additional Numeric Suffixes
            across(c(var_1, var_2), ~ case_when(
                str_length(.) > 5 ~ str_remove(., ":[0-9]{2}$"),
                str_length(.) == 1 ~ paste0("0", .),
                TRUE ~ .
            ))) |>
        #* step 6: Format Allele Codes
        mutate(across(c(var_1, var_2), ~ case_when(
            str_length(.) == 3 & substr(., 1, 1) == "1" & substr(., 2, 2) %in% c("0", "1", "2", "3", "4", "5", "6") ~
                paste0(substr(., 1, nchar(.) - 1), ":", substr(., nchar(.), nchar(.))),
            
            str_length(.) == 3 & substr(., 1, 1) == "1" & !substr(., 2, 2) %in% c("0", "1", "2", "3", "4", "5", "6") ~
                paste0(substr(., 1, nchar(.) - 2), ":", substr(., nchar(.) - 1, nchar(.))),
            
            str_length(.) == 3 & !substr(., 1, 1) %in% c("0", "1") ~
                paste0(substr(., 1, nchar(.) - 2), ":", substr(., nchar(.) - 1, nchar(.))),
            
            TRUE ~ .
        ))
        ) |>
        #* step 7: Fix Four-Character Alleles
        mutate(
            across(
                c(var_1, var_2),
                ~ if_else(
                    str_length(.) == 4 & !str_detect(., ":"),
                    str_replace(., "(..)$", ":\\1"),
                    .)),
            #* step 8: Resolve Missing or Placeholder Values
            var_1 = if_else(str_detect(var_1, "X|x"), var_2, var_1),
            var_2 = if_else(str_detect(var_2, "X|x"), var_1, var_2),
            across(c(var_1, var_2), 
                   ~ if_else(str_length(.) == 4 & !str_detect(., "X|x") & !str_detect(., "^[0-9]{2}:"),
                             paste0("0", .),
                             .)),
            var_1 = if_else(str_detect(var_1, "X|x|NA_character_| ") & str_detect(var_1, "X|x|NA_character_| "), NA_character_, var_1),
            var_2 = if_else(str_detect(var_2, "X|x|NA_character_| ") & str_detect(var_1, "X|x|NA_character_| "), NA_character_, var_2)
        ) |>
        select(-c(var_1_c, var_2_c)) |>
        mutate(
            var_1 = if_else(var_1 == "" & var_2 != "", var_2, var_1),
            var_2 = if_else(var_2 == "" & var_1 != "", var_1, var_2)
        ) |>
        #* step 9: Rename Columns Back to Original Names
        rename({{var_1}} := var_1, {{var_2}} := var_2)
}
