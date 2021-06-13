clean_names <- function(.data, unique = FALSE) {
  n <- if (is.data.frame(.data)) colnames(.data) else .data
  # n <- gsub("%+", "_pct_", n)
  # n <- gsub("\\$+", "_dollars_", n)
  # n <- gsub("\\++", "_plus_", n)
  # n <- gsub("-+", "_minus_", n)
  # n <- gsub("\\*+", "_star_", n)
  # n <- gsub("#+", "_cnt_", n)
  # n <- gsub("&+", "_and_", n)
  # n <- gsub("@+", "_at_", n)
  # n <- gsub("[^a-zA-Z0-9_]+", "_", n)
  # n <- gsub("([A-Z][a-z])", "_\\1", n)
  n <- gsub(" ","_", n)
  n <- gsub("\\(","_", n)
  n <- gsub("\\)","_", n)
  n <- tolower(trimws(n))
  
  n <- gsub("(^_+|_+$)", "", n)
  
  n <- gsub("_+", "_", n)
  
  if (unique) n <- make.unique(n, sep = "_")
  
  if (is.data.frame(.data)) {
    colnames(.data) <- n
    .data
  } else {
    n
  }
}


wilcox_lookup <- function(a=NULL, b=NULL, var=NULL, df=NULL){
  
  w <- wilcox.test(value ~ sample_id
                   , data = df %>% 
                     dplyr::filter(sample_id %in% c(a, b) & variable == var) 
                   , exact = FALSE)
  return(w$p.value)
}


metab_change <- function(metab=NULL, df=NULL, a=NULL, b=NULL){
    
    chg_mean <- 
      df %>% filter(variable == metab & sample_id == b) %>% select(mean) / 
      df %>% filter(variable == metab & sample_id == a) %>% select(mean) 
    
    chg_median <- 
      df %>% filter(variable == metab & sample_id == b) %>% select(median) / 
      df %>% filter(variable == metab & sample_id == a) %>% select(median)
    
    return(c(chg_mean %>% as.numeric(), chg_median %>% as.numeric()))
}






# ..device.set.up <- FALSE
# ..current.page <<- 0
# 
# save.bookmark <- function(text,bookmarks=list(),level=1,page=NULL) {
#   if (!..device.set.up) {
#     Cairo.onSave(device = dev.cur(),
#                  onSave=function(device,page){
#                    ..current.page <<- page
#                  })
#     ..device.set.up <<- TRUE
#   }
#   if (missing(page)|| is.null(page)) {
#     page <- ..current.page+1
#   }
#   bookmarks[[length(bookmarks)+1]] <-
#     list(text=text,
#          level=level,
#          page=page)
#   return(bookmarks)
# }
# 
# write.bookmarks <- function(pdf.file,bookmarks=list()) {
#   pdf.bookmarks <- ""
#   for (bookmark in 1:length(bookmarks)) {
#     pdf.bookmarks <-
#       paste0(pdf.bookmarks,
#              "BookmarkBegin\n",
#              "BookmarkTitle: ",bookmarks[[bookmark]]$text,"\n",
#              "BookmarkLevel: ",bookmarks[[bookmark]]$level,"\n",
#              "BookmarkPageNumber: ",bookmarks[[bookmark]]$page,"\n")
#   }
#   temp.pdf <- tempfile(pattern=basename(pdf.file))
#   temp.pdf.info <- tempfile(pattern=paste0(basename(pdf.file),"info_utf8"))
#   cat(file=temp.pdf.info,pdf.bookmarks)
#   system2("pdftk",c(pdf.file,'update_info_utf8',temp.pdf.info,'output',temp.pdf))
#   if (file.exists(temp.pdf)) {
#     file.rename(temp.pdf,pdf.file)
#   } else {
#     warning("unable to properly create bookmarks")
#   }
#   
#   

## check for availability of 'pdftk' and 'convert'
# pdftk <- try(shsystem("pdftk --version", intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE), silent = TRUE)
# magic <- try(shsystem("convert --version", intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE), silent = TRUE)
# if(inherits(pdftk, "try-error")) stop("system requirement 'pdftk' is not available for merging/rotating/splitting PDFs")
# if(inherits(magic, "try-error")) stop("system requirement 'convert' is not available for converting PDF to PNG")
# 
# CairoPDF(file="testing.pdf")
# bookmarks <- list()
# bookmarks <- save.bookmark("First plot",bookmarks)
# plot(1:5,6:10)
# bookmarks <- save.bookmark("Second plot",bookmarks)
# plot(6:10,1:5)
# dev.off()
# write.bookmarks("testing.pdf",bookmarks)

# }

