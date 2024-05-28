#######################################################################################################################
## R Workspace
#######################################################################################################################

## Get_Todays_Date
if(exists(todays_date)){
    date_check <- todays_date    
    todays_date <- str_replace_all(string = as.Date(Sys.Date(),
                                                    format = "%m / %d / %Y"), pattern = "-", replacement = "_")
    if(date_check!=todays_date){
        save_counter <- 0
        }
} else{
    todays_date <- str_replace_all(string = as.Date(Sys.Date(),
                                                    format = "%m / %d / %Y"), pattern = "-", replacement = "_")
    save_counter <- 0
}

save_counter <- save_counter + 1

## Name Workspace file name
Workspace <-  paste0("/Users/Aadi/Desktop/PhD_Project/Data/R_Workspaces/Simplified_Model_v1_MH_Shared_All_Adaptive_",
                     todays_date,"_",save_counter,".RData")


## Save Work
save.image(Workspace)
