pkgname <- "leafblower"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('leafblower')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("harvest")
### * harvest

flush(stderr()); flush(stdout())

### Name: harvest
### Title: Generate calibrated weights (drop-in for autumn::harvest)
### Aliases: harvest

### ** Examples

## Not run: 
##D df  <- data.frame(sex = factor(sample(c("M","F"), 500, TRUE)))
##D tgt <- list(sex = c(M = 0.5, F = 0.5))
##D 
##D # Greenkhorn (greedy coordinate-descent IPF)
##D r_grk <- harvest(df, tgt, method = "greenkhorn")
##D 
##D # Logit-distance Newton calibration (Deville-Sarndal 1992)
##D r_logit <- harvest(df, tgt, method = "logit")
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
