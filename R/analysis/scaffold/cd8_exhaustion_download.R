# Example fcs file url: https://premium.cytobank.org/cytobank/experiments/154556/fcs_files/2487981/events'

data_dir <- '/Users/eczech/repos/cytometry-experiments/R/analysis/scaffold/data/bengsch_2018'
cookie <- '_ga=GA1.2.1013683102.1539117436; __utmc=51498683; __utmz=51498683.1546976315.1.1.utmcsr=ncbi.nlm.nih.gov|utmccn=(referral)|utmcmd=referral|utmcct=/; __utma=51498683.1013683102.1539117436.1547195094.1547238377.4; __stripe_mid=1256981e-0fd9-4a99-8e16-b91ee76705ff; _gid=GA1.2.2104498739.1547317662; __stripe_sid=2069bd85-ff9b-433a-8bd1-c3ba72815086; user_credentials=3fcbd07a12725f7aed1bdcad5ba6122d9d6b8e806897c572f96f4b9fea98f33f3fc40639acd1ede01f60a8d67dd23a104dc8d53f5df4ce97af5a2a9d3f71c3b3%3A%3A7991; _session_id=899c31d1e39ffa6414d8a9755b14569a'
url_format <- 'https://premium.cytobank.org/cytobank/experiments/154556/fcs_files/%s/events'
fcs_file_ids <- 2487980 + 0:46
fcs_file_names <- c('CD8_HC_001.fcs', 'CD8_HC_002.fcs', 'CD8_HC_003.fcs', 'CD8_HC_004.fcs', 'CD8_HC_005_batch1.fcs', 'CD8_HC_005_batch2.fcs', 'CD8_HC_005_batch3.fcs', 'CD8_HC_006.fcs', 'CD8_HC_007.fcs', 'CD8_HIV_001.fcs', 'CD8_HIV_002.fcs', 'CD8_HIV_003.fcs', 'CD8_HIV_005.fcs', 'CD8_HIV_006.fcs', 'CD8_HIV_009.fcs', 'CD8_HIV_010.fcs', 'CD8_HIV_011.fcs', 'CD8_HIV_015.fcs', 'CD8_HIV_016.fcs', 'CD8_HIV_017.fcs', 'CD8_HIV_021.fcs', 'CD8_HIV_022.fcs', 'CD8_HIV_024.fcs', 'CD8_HIV_026.fcs', 'CD8_HIV_027.fcs', 'CD8_HIV_028.fcs', 'CD8_HIV_029.fcs', 'CD8_LUCA_001_LUNG.fcs', 'CD8_LUCA_001_PBMC.fcs', 'CD8_LUCA_001_TIL.fcs', 'CD8_LUCA_002_Lung.fcs', 'CD8_LUCA_002_PBMC.fcs', 'CD8_LUCA_002_TIL.fcs', 'CD8_LUCA_003_Lung.fcs', 'CD8_LUCA_003_PBMC.fcs', 'CD8_LUCA_003_TIL.fcs', 'CD8_LUCA_004_Lung.fcs', 'CD8_LUCA_004_PBMC.fcs', 'CD8_LUCA_004_TIL.fcs', 'CD8_LUCA_005_Lung.fcs', 'CD8_LUCA_005_PBMC.fcs', 'CD8_LUCA_005_TIL.fcs', 'CD8_LUCA_006_Lung.fcs', 'CD8_LUCA_006_PBMC.fcs', 'CD8_LUCA_006_TIL.fcs', 'CD8_LUCA_010_PBMC.fcs', 'CD8_LUCA_010_TIL.fcs')
urls <- sprintf(url_format, file_ids)

#for (i in 4:length(urls)){
for (i in (length(urls) - c(4, 1, 0))){
  cmd <- sprintf('wget -nv -O %s/%s.txt --no-cookies --header "Cookie: %s" %s', data_dir, fcs_file_names[i], cookie, urls[i])  
  system(cmd)
}

# wget --no-cookies --header "Cookie: _ga=GA1.2.1013683102.1539117436; __utmc=51498683; __utmz=51498683.1546976315.1.1.utmcsr=ncbi.nlm.nih.gov|utmccn=(referral)|utmcmd=referral|utmcct=/; __utma=51498683.1013683102.1539117436.1547195094.1547238377.4; _session_id=5f5eb1cfd93df84ba1f9819e26c4b423; __stripe_mid=1256981e-0fd9-4a99-8e16-b91ee76705ff; __stripe_sid=081358c2-a418-470b-9212-f11ef28a8fe8; user_credentials=c66bc277514a9ec2197a81af255835539ac4e885261843d748ecd80c883035cff644292e0e0e818c24ad36198d032e1ab3d6eec342b749e0eeb822ce51cd99fc%3A%3A7991" https://premium.cytobank.org/cytobank/experiments/154556/fcs_files/2487981/events

# This is what CytobankAPI:authenticate does but it doesn't work for either community or premium
# getURL("http://example.com", httpheader = c(Accept = "application/json", X-My-Token = 666))
# resp <- POST(paste("https://premium.cytobank.org/cytobank/api/v1/authenticate", sep=""),
#              #body=list(username='eczech', password=''),
#              body=list(username='eczech', password=''),
#              encode="json",
#              timeout(5000000)
# )
# str(content(resp))

# library(CytobankAPI)
# cyto_session <- authenticate(site="community", username="eczech", password="")
# downloaded_fcs_file <- fcs_files.download(cyto_session, experiment_id=22, fcs_file_id=10, directory="/my/specified/directory")
