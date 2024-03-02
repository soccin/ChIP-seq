get_project_number<-function(){
    projNo=grep("^Proj_",strsplit(getwd(),"/")[[1]],value=T)[1]
    if(is.na(projNo)) {
        projNo=Sys.getenv("PROJECT_NUMBER")
        cat("\n")
        cat("Can not find project number from path\n")
        cat("using env var {PROJECT_NUMBER}\n")
        cat("\n\tprojNo =",paste0("[",projNo,"]"),"\n\n")
    }
    projNo
}