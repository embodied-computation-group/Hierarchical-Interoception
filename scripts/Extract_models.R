get_models <- function(osftoken){
  # get the preanalyzed data from osf:
  if (dir.exists(here::here("results"))) {
    osfr::osf_auth(token = osftoken)
    
    HI <- osfr::osf_retrieve_node("https://osf.io/8yhtp/")
    
    downloaded = HI %>%
      osfr::osf_ls_files(pattern = "fit population") %>%
      osfr::osf_download(path = here::here("results"), recurse = TRUE, conflicts = "overwrite", progress = TRUE)
    
    
    
  }
  
}

get_models(osf_token)

osf_file_path <- here::here("osf","osf.txt")
osf_token <- read_lines(osf_file_path)[1]
osfr::osf_auth(token = osf_token)

HI <- osfr::osf_retrieve_node("https://osf.io/8yhtp/")

HI %>% osfr::osf_upload(
  path = here::here("results","Fit population"),
  recurse = TRUE,
  conflicts = "overwrite" # or "skip" or "ask"
)
