get_models <- function(osftoken){
  # get the preanalyzed data from osf:
  if (!dir.exists(here::here("results"))) {
    osfr::osf_auth(token = osf_token)
    
    HI <- osfr::osf_retrieve_node("https://osf.io/8yhtp/")
    
    downloaded = HI %>%
      osfr::osf_ls_files(pattern = "results") %>%
      osfr::osf_download(path = here::here(""), recurse = TRUE, conflicts = "overwrite", progress = TRUE)
    
    
    
  }
  
}


osf_file_path <- here::here("osf","osf.txt")
osf_token <- read_lines(osf_file_path)[1]
osfr::osf_auth(token = osf_token)

HI <- osfr::osf_retrieve_node("https://osf.io/8yhtp/")

HI %>% osfr::osf_upload(
  path = here::here("results","Fit population"),
  recurse = TRUE,
  conflicts = "overwrite" # or "skip" or "ask"
)
