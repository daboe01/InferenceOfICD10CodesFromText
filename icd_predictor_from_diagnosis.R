# library("devtools")
# install_github("bmschmidt/wordVectors")
library(wordVectors)
library(stringr)
library(plyr)

model = read.vectors('/app/MODEL.bin', binary=1)
words = rownames(model)
words = str_replace_all(words,'&#32;',' ')
rownames(model)=words

icd_clean = read.delim("/app/icd10gm2020syst_kodes.txt", sep=";", stringsAsFactors=F, header=F)
icd = read.delim("/app/diagnosen4.txt", sep="\t", stringsAsFactors=F, header=F)

replace_umlauts <- function(x) {
  umlauts <- "äöü"
  UMLAUTS <- "ÄÖÜ"
  UMLAUT2 <- "ß"
  UMLAUT3 <- "-"
  x <- gsub(pattern = paste0("([", UMLAUTS, "])"), replacement = "\\1E", x)
  x <- gsub(pattern = paste0("([", umlauts, "])"), replacement = "\\1e", x)
  x <- gsub(pattern = paste0("([", UMLAUT2, "])"), replacement = "ss", x)
  x <- gsub(pattern = paste0("([", UMLAUT3, "])"), replacement = " ", x)
  x <- chartr(old = paste0(UMLAUTS, umlauts), new = "AOUaou", x)
  return(x)
}
alpha = read.delim("/app/icd10gm2021_alphaid_edvtxt_20201002.txt", sep="|", stringsAsFactors=F, header=F)
alpha$V5 = ifelse(alpha$V3 == '', alpha$V5, alpha$V3)
alpha$V7 = replace_umlauts (tolower(alpha$V7))
alpha= subset(alpha, select=c(V5, V7))
names(alpha)=c("V1","V2")

icd = unique(rbind(icd, alpha))
icd$V2 = gsub('[^a-z 0-9]+', ' ', icd$V2)
icd$V2 = gsub(' +', ' ', icd$V2)

icd_embeddings_df1 = ddply(icd, .(V1), function(x){
    keywords = unique(x$V2)

    return(model[[keywords, average=T]])
})


icd_embeddings_matrix = as.matrix(subset(rbind(icd_embeddings_df1), select=-c(V1)))
row.names(icd_embeddings_matrix) = c(icd_embeddings_df1$V1)
icd_embeddings = as.VectorSpaceModel(icd_embeddings_matrix)

# inference...
#* @get /icd_from_diagnosis
icd_from_diagnosis= function (input="katarakt")
{
    icd_found = ''
    diagnosis = unlist(strsplit(input, ' '))
   
    num_words = length(diagnosis)

    if (num_words > 1 & input %in% rownames(model))
        diagnosis = c(input, diagnosis)

    diagnosis = unique(diagnosis)
    
    # print(diagnosis)
   
    # direkte treffer bevorzugen
    if (nrow(subset(icd, V2 == input)) > 0)
    {
        icd_found = subset(icd, V2==input)[1,]$V1
    }
    else
    {

        query = model[[diagnosis, average=T]]
        d = cosineDist(query, icd_embeddings)
        icd_found = unlist(dimnames(d))[[which.min(d)]]

        threshold = 0.5
        
        if (d[which.min(d)] > threshold)
            icd_found = ''
    }
    
    if (nchar(icd_found))
    {
        icd_found = gsub('[+]+', '', icd_found)
        icd_found2 = gsub('[A-Z]+$', '', icd_found)
        icd_found2 = gsub('(.+)[A-Z]+[0-9]+$', '\\1', icd_found2)
        return (c(icd_found, subset(icd_clean, V7 == icd_found2)[1,]$V9))
    }
    
    return (NA)
}
