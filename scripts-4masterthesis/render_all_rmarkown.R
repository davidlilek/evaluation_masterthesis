# define file path
path = "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/"
# get all Rmd files
files <- list.files(path = path, pattern = "*.Rmd")
# render every Rmd file
for (f in files) rmarkdown::render(f)
