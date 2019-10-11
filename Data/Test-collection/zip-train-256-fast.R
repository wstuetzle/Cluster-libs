input.pathname <- paste(data.dir, "zip-train.dat", sep = dir.sep)
zip.train <- matrix(scan(input.pathname), ncol = 257, byrow = T)
## group.id <- zip.train[,1]
group.id <- zip.train[,1] + 1  ## Some code my assume that group.id = 1...k
X <- zip.train[,2:257]
assign("X", X, envir = .GlobalEnv)
assign("group.id", group.id, envir = .GlobalEnv)
assign("creation.date", date(), envir = .GlobalEnv)
