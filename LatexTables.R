### Tables
library(foreign)
library(xtable)
library(stargazer)

# Group 1
RES1[[1]] <- round(RES1[[1]], digits = 3)
RES1[[2]] <- round(RES1[[2]], digits = 3)
RES1[[3]] <- round(RES1[[3]], digits = 3)

xtable(RES1[[1]], type = "latex", file = "RES1a.tex")
xtable(RES1[[2]], type = "latex", file = "RES1b.tex")
xtable(RES1[[3]], type = "latex", file = "RES1c.tex")

#Group 2
RES2[[1]] <- round(RES2[[1]], digits = 3)
RES2[[2]] <- round(RES2[[2]], digits = 3)
RES2[[3]] <- round(RES2[[3]], digits = 3)

xtable(RES2[[1]], type = "latex", file = "RES2a.tex")
xtable(RES2[[2]], type = "latex", file = "RES2b.tex")
xtable(RES2[[3]], type = "latex", file = "RES2c.tex")

#Group 3
RES3[[1]] <- round(RES3[[1]], digits = 3)
RES3[[2]] <- round(RES3[[2]], digits = 3)
RES3[[3]] <- round(RES3[[3]], digits = 3)

xtable(RES3[[1]], type = "latex", file = "RES3a.tex")
xtable(RES3[[2]], type = "latex", file = "RES3b.tex")
xtable(RES3[[3]], type = "latex", file = "RES3c.tex")
