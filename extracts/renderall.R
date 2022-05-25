path = "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/extracts"
files <- list.files(path = path, pattern = "*.Rmd")

for (f in files) rmarkdown::render(f)
