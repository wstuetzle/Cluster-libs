input.pathname <- paste(data.dir, "Duong-et-al-2008-CSDA-Fig1-data.txt",
                        sep = dir.sep)
X <- matrix(scan(input.pathname), ncol = 3, byrow = T)
n <- nrow(X)
group.id <- rep(1, n)
assign("X", X, envir = .GlobalEnv)
assign("group.id", group.id, envir = .GlobalEnv)
assign("creation.date", date(), envir = .GlobalEnv)
