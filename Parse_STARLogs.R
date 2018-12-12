
args = commandArgs( trailingOnly=TRUE )

if (length(args)==0) {
  stop("At least one argument (prefix_name)  must be supplied (input file).n", call.=FALSE)
}

library(stringr)
## Should be in Dir with *Log.final files

logs = list.files()
logs = logs[ grepl("Log.final", logs ) ]

df = c()

for ( i in logs) {
  ID = str_split(i, "Log")[[1]][1]
  val = read.table(i, sep ="|", fill = TRUE)
  val[ ,2] = stringr::str_replace(val[ ,2], "\t", "")
  total = val[5, 2]
  uni = val[9, 2]
  multi = val[24, 2]
  many = val[26, 2]
  tooshort = val[29, 2]
  unmap = val[30, 2]
  
  df = rbind(df, c(ID, total, uni, multi, many, tooshort, unmap) )
}
colnames(df) = c("ID", "Total_Reads", "Unique", "Multimap", "Manymap", "Unmap_short", "Unmap_other")

df = as.data.frame(df)
df[ ,3:ncol(df) ] = mapply( df[ ,3:ncol(df)], FUN = function(x) str_replace(x, "%", "") )
df[ ,3:ncol(df) ] = mapply( df[ ,3:ncol(df)], FUN = function(x) as.numeric(x) )

df$Total_Mapped = df$Unique +df$Multimap +df$Manymap

write.csv(df, paste0("./", args[1], "_STAR_Logs.csv") )
