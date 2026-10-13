library(tools)

package_grep <- function(regex){ grep("revdep", dir("..", pattern = regex, recursive = TRUE, full.names = TRUE), invert = TRUE, value = TRUE) }
.aspell_control <- function(wordlist) c("--master=en_GB", if(missing(wordlist)) NULL else sprintf("--add-wordlists=%s", wordlist))

aspell_R   <- aspell(package_grep("\\.R$"), filter = "R", control = .aspell_control("WORDLIST_R.txt"))
aspell_Rd  <- aspell(package_grep("\\.Rd$"), filter = list("Rd", drop = "\\references"), control = .aspell_control("WORDLIST_Rd.txt"))
aspell_dcf <- aspell("../DESCRIPTION", filter = "dcf", control = .aspell_control())

print(eval(str2expression(sprintf("rbind(%s)", paste(grep("^aspell\\_(R|Rd|pot|dcf|md)$", ls(), value = TRUE), collapse = ", ")))), verbose = TRUE)



# NEWS.Rd --> ../NEWS
Rd2txt_options(
  width            = 80L,
  minIndent        = 2L,
  extraIndent      = 2L,
  sectionIndent    = 2L,
  sectionExtra     = 2L,
  underline_titles = FALSE,
  unicode_symbols  = TRUE,
  itemBullet       = "• "
)
Rd2txt("NEWS.Rd", out="../NEWS", package = "DirichletReg", outputEncoding = "UTF-8")



# NEWS.Rd --> NEWS.pdf
system("R CMD Rd2pdf --no-preview --encoding=UTF-8 --force NEWS.Rd")
compactPDF("NEWS.pdf", gs_quality="printer", verbose = TRUE)
