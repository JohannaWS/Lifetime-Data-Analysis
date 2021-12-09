setwd("C:/Users/johan/Google Drive/UNI/UPC/Lifetime Data Analysis/Project")

###### IMPORT ######
names=list("id", "survival_times", "death", "CD4_cells", "treatment","gender", "AIDS", "AZT")
aids=read.table("AIDS.txt", skip=28, col.names=names)
