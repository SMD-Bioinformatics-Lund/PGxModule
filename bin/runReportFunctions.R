# Functions
risk_duplication <- function(x){
  min_dup_var <- min(abs(c(1/4, 1/3, 2/3, 3/4) - x))
  min_diploid <- min(abs(c(0, 1/2, 1) - x))
  if (min_dup_var > min_diploid){
    return(FALSE)
  }
  return(TRUE)
}

get_haplotypes <- function(ID, haplotype_definitions=haplotype_definitions) {
  haplo <- haplotype_definitions[haplotype_definitions$ID == ID, "HAPLOTYPE"]
  return(paste(haplo, collapse = "/"))
}

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))

  df[,nums] <- round(df[,nums], digits = digits)

  (df)
}

highlight_text <- function(x, color,weight="normal",type="normal") {
  sprintf("<span style='color: %s;font-weight: %s;font-style: %s;'>%s</span>", color,weight,type,x)
}