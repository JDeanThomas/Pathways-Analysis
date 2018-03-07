library(xlsx)
library(dplyr)
library(biomaRt)
#require(graphite)
require(KEGGREST)
require(reactome.db)

wb1 <- read.xlsx("/Users/jasonthomas/Downloads/nature13595-s4.xls", 1, as.data.frame = TRUE)
wb2 <- read.xlsx("/Users/jasonthomas/Downloads/nature13595-s4.xls", 2, as.data.frame = TRUE)
colnames(wb1) <- c("PathwayID", "GenensInPathway", "Sig", "ExpSig", "PValue", "CorrP", "PathwayName")
colnames(wb2) <- c("PathwayID", "GenensInPathway", "AssociatedInts", "PValue", "CorrP", "PathwayName")
#wb1$PathwayName <- gsub("_", " ", wb1$PathwayName)
#wb2$PathwayName <- gsub("_", " ", wb2$PathwayName)
df <- bind_rows(wb1, wb2)
df <- df[1:7]
df$GeneID <- ""
df <- df[,c(1, 2, 8, 3, 4, 5, 6, 7)]
df2 <- df
df3 <- df2[1]
colnames(df3) <- "Database"
df3 <- bind_cols(df3, df2)

# regext for DB col

df3$Database <- gsub("(?::).+", "", df3$Database, perl = TRUE)
df3$Database <- gsub("^\\s+(?=\\w)", "", df3$Database, perl = TRUE)
df3$Database <- gsub("(?<=\\s).+", "", df3$Database, perl = TRUE)
# Repeat to remove :6 from Panther:6 in df$Database for some reason
df3$Database <- gsub("(?::).+", "", df3$Database, perl = TRUE)
# Rename PAN-PW (seriously?) to PANTHER
df3$Database <- gsub("PAN-PW", "PANTHER", df3$Database, perl = TRUE)
df3$Database <- gsub("BioCarta", "BIOCARTA", df3$Database, perl = TRUE)

dfPathways <- df3[with(df3, order(Database, PathwayID)),]

rm(df, df2, df3, wb1, wb2)


# regex for Pathway IDs

# wb1 $PathwayID regex fixes 

# GO schema should be GO:*******. Entered replacing leading zeros with whitespace.
# KEGG entries in wrong schema should be hsa*****. Entered without hsa or zeros
#   and whitespace added before ID and missleading additional whitespace
#   before ID number, not just replacing zeros with whitespace
# Panther schema should be p*****. Entered as PAN-PW with leading zeros replaced with whitespace

# NCI: Can't even make sense of NCI ID schema from Google search.
# BIOCARTA
# REACTOME: Can not desipher schema or find pathways after exhaustive search
  # eg: http://www.reactome.org/content/query?q=716+mitotic+g2-g2%2Fm&species=Homo+sapiens&species=Entries+without+species&cluster=true&types=Pathway
  # for entry REACTOME 716 with PathwayName "Mitotic G2-G2"
# MGI schema appears to be MGI:******. Entered as MGI: with leading zeros replaced by whitespace
 

# Stripout leading whitespace from PathwayIDs
dfPathways$PathwayID <- gsub("^\\s+(?=\\w)", "", dfPathways$PathwayID, perl = TRUE)

# KEGG

# Entries from sheet 1
# Replace KEGG with hsa and strip 2 /s incerted before ID number (out of schema)
dfPathways$PathwayID <- gsub("KEGG\\s{2}", "hsa", dfPathways$PathwayID) 
# Replace remaining /s with zeros or correct DB schema
dfPathways$PathwayID <- gsub("\\s", "0", dfPathways$PathwayID, perl = TRUE) # Changes all. Might need below.

# dfPathways$PathwayID <- gsub("(?<=hsa)\\s{3}", "000", dfPathways$PathwayID, perl = TRUE) 
# dfPathways$PathwayID <- gsub("(?<=hsa)\\s{2}", "00", dfPathways$PathwayID, perl = TRUE) 
# dfPathways$PathwayID <- gsub("(?<=hsa)\\s{1}", "0", dfPathways$PathwayID, perl = TRUE) 

# Entries from sheet 2
dfPathways$PathwayID <- gsub("KEGG:HSA", "hsa", dfPathways$PathwayID)

# amiGO

# Entries from sheet 1

# REGG regex 2 from sheet one handles, might change

# Entries from sheet 2
dfPathways$PathwayID <- gsub("(?<=GO:)(?!0)(?!\\d{7})(?=\\d{6})", "0", dfPathways$PathwayID, perl = TRUE) 
dfPathways$PathwayID <- gsub("(?<=GO:)(?!0)(?!\\d{7})(?=\\d{5})", "00", dfPathways$PathwayID, perl = TRUE) 
dfPathways$PathwayID <- gsub("(?<=GO:)(?!0)(?!\\d{7})(?=\\d{4})", "000", dfPathways$PathwayID, perl = TRUE) 
dfPathways$PathwayID <- gsub("(?<=GO:)(?!0)(?!\\d{7})(?=\\d{3})", "0000", dfPathways$PathwayID, perl = TRUE)  
dfPathways$PathwayID <- gsub("(?<=GO:)(?!0)(?!\\d{7})(?=\\d{2})", "00000", dfPathways$PathwayID, perl = TRUE) 

# MGI

# Sheet 1

# Whitespace handled by KEGG regext line 2, sheet 1

dfPathways$PathwayID <- gsub("MGI", "MP", dfPathways$PathwayID) 

# Sheet 2

# Fixed IDs in DB column creation
# Change MP DB column entries from sheet 2 to MGI
dfPathways$Database <- gsub("MP", "MGI", dfPathways$Database) 
# Reorder
dfPathways <- dfPathways[with(dfPathways, order(Database, PathwayID)),]



### IGNORE FOR NOW  ###

# Replace \s not preceeded by :, preceeded by word char with :\s 
text2 <- gsub("\\w(?<!:)\\s", ": ", text, perl = TRUE) 
#alt \s preceeded by word char with :\s
#text3 <- gsub("(?<=\\w)\\s", ": ", text, perl = TRUE) 
#Replace 4 \s preceeded by colon with 0000
text2 <- gsub("(?<=:)\\s{4}", "0000", text2, perl = TRUE) 
#Replace 3 \s preceeded by colon with 000
text2 <- gsub("(?<=:)\\s{3}", "000", text2, perl = TRUE) 
#Replace 4 \s preceeded by colon with 00
text2 <- gsub("(?<=:)\\s{2}", "00", text2, perl = TRUE) 
#Replace 1 \s preceeded by colon with 0
text2 <- gsub("(?<=:)\\s{1}", "0", text2, perl = TRUE) 

# wb2 $PathwayID regex fixes 

# GO, PANTHER, NCI:REACTOME and BioCarta have no preceeding zeros. Check schema 
# MP and NCI appear fine. Check. 
# Check KEGG:HSA scheme (consistant with correct digits, may need lower case hsa)

# going to need to count digits and insert zeros to fit schema

# GO schema should be GO:*******. Entered without leading zeros. Check for overlap with wb1
# KEGG entries in wrong schema should be hsa*****. Entered as KEGG:HSA*****. Check for overlap in wb1
# PANTHER schema should be p*****. Entered as PANTHER:* without leading zeros.
# NCI: Can't even make sense of NCI ID schema from Google search.
# Can't even make sense of NCI:REACTOME ID schema from Google search

### End ignore ####


# Pulling Gene IDs from Pathway IDs

# amiGO

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
goTest <- getBM(attributes=c('go_id', 'hgnc_symbol', 'ensembl_transcript_id'), 
                filters = 'go_id', values = dfPathways$PathwayID[5:350], 
                mart = ensembl, uniqueRows = TRUE)
colnames(goTest) <- c("PathwayID", "HGNC_Symbol", "GeneID (ensembl_transcript_id)")
goTest <- goTest[,c(1, 3, 2)]

write.csv(goTest, file = "./amiGO_IDs.csv")

# KEGG

dfKEGG <- as.data.frame(as.matrix(keggLink("hsa", dfPathways$PathwayID[351:394])))
dfKEGG  <- cbind(PathwayID = rownames(dfKEGG), dfKEGG)
rownames(dfKEGG) <- NULL
colnames(dfKEGG) <- c("PathwayID", "GeneID")
dfKEGG$Database <- "KEGG" 
dfKEGG <- dfKEGG[,c(3, 1, 2)]
dfKEGG$PathwayID <- gsub("path:", "", dfKEGG$PathwayID)
dfKEGG <- dfKEGG[with(dfKEGG, order(PathwayID, GeneID)),]
# Do so keeps GenesInPathway vals
dfPathways2 <- bind_rows(dfPathways, dfKEGG)

# MGI
# (MGI to be decpreciated 10/26/15. Phenotype IDs querried from Mousemine DB)

MGI = useMart("biomart", dataset="markers")
dfMGI <- getBM(attributes=c('mgi_marker_id_att', 'mouse_entrez_gene_id_108',
                            'marker_symbol_107', 'marker_name_107 '), 
               filters = 'mgi_marker_id_att', values = dfPathways$PathwayID[395:500], mart = MGI, uniqueRows = TRUE)



dfMGI <- getBM(attributes=c('mgi_marker_id_att', 'mouse_entrez_gene_id_108', 
                            'marker_symbol_107', 'marker_name_107'), 
               filters = 'entrezGeneid', values = "MP:0000084", mart = MGI, uniqueRows = TRUE)


filters = 'mmarker_id'


MGI = useMart("biomart", dataset="result")
dfMGI <- getBM(attributes=c('markermgiid', 'markersymbol'), 
               filters = 'markermgiid', values = "MP:0000443", mart = MGI, uniqueRows = TRUE)

mart = useMart("biomart", dataset="result")

listDatasets(mart)
listAttributes(MGI)
listFilters(MGI)
listAttributes(mart)

MGI <- getBM(attributes=c('go_id', 'hgnc_symbol', 'marker_symbol_107', 'marker_name_107 '), 
             filters = 'go_id', values = dfPathways$PathwayID[5:350], mart = ensembl, uniqueRows = TRUE)







dfKEGG2 <- as.data.frame(as.matrix(keggLink("hsa", dfPathways$PathwayID[351])))
dfKEGG  <- cbind(PathwayID = rownames(dfKEGG), dfKEGG)
rownames(dfKEGG) <- NULL
colnames(dfKEGG) <- c("PathwayID", "GeneID")
dfKEGG$Database <- "KEGG" 
dfKEGG <- dfKEGG[,c(3, 1, 2)]
dfKEGG$PathwayID <- gsub("path:", "", dfKEGG$PathwayID)
dfKEGG <- dfKEGG[with(dfKEGG, order(PathwayID, GeneID)),]
# Do so keeps GenesInPathway vals
dfPathways2 <- bind_rows(dfPathways, dfKEGG)


dfKEGG2 <- data.frame(, c("PathwayID", "HGNC_Symbol"))
colnames(dfKEGG2) <- c("PathwayID", "HGNC_Symbol")
dfKEGG3 <- keggGet(dfPathways$PathwayID[351:356])
dfKEGG3$PathwayID <- dfKEGG3[[1]]$PATHWAY_MAP
dfKEGG2$HGNC_Symbol <- dfKEGG3[[1]]$GENE

dfKEGG2 <- data.frame("PathwayID", "HGNC_Symbol")
colnames(dfKEGG2) <- c("PathwayID", "HGNC_Symbol")
dfKEGG3 <- keggGet(dfPathways$PathwayID[351])
dfKEGG3$PathwayID <- dfKEGG3[[1]]$PATHWAY_MAP
dfKEGG2$HGNC_Symbol <- dfKEGG3[[1]]$GENE

PathKEGG <- keggGet(dfPathways$PathwayID[351:356])


dfKEGG2 <- data.frame(PathwayID = character(length(dfKEGG$GeneID) * 2), 
                      HGNC_Symbol = character(length(dfKEGG$GeneID) * 2), 
                      stringsAsFactors = FALSE)
idKEGG <- dfKEGG$PathwayID

for i in seq_along(dfKEGG$PathwayID)
  tempKEGG <- data.frame(PathwayID = character(length(unique(dfKEGG$PathwayID) * 2), 
                      HGNC_Symbol = character(length(dfKEGG$GeneID) * 2), 
                      stringsAsFactors = FALSE)
  tempKEGG$PathwayID <- rep(dfPathways$PathwayID[i], length(dfPathways$PathwayID[i]))
  tempKEGG$HGNC_Symbol <- dfKEGG3[[i]]$GENE
  dfKEGG2 <- bind_rows(wb1, wb2)
r <-nrow(dfKEGG2)
dfKEGG2[, !(r%%2==1)]



ConvKEGG <- function(PathKEGG) {
  df <- data.frame()
  for (i in seq_along(PathKEGG)) { 
    Klist <- list()
    #rawData <- fromJSON(file=apiCall[i])
    Klist <- lapply(rawData$response$docs, function(x) {
      x[sapply(x, is.null)] <- NA
      unlist(x)})
    for (i in seq_along(rdParsed)) {
      dfContainer <- data.frame(do.call(unlist("rbind"), rdParsed[i]), stringsAsFactors = FALSE)
      df <- bind_rows(df, dfContainer)
    }
  }
  return(df)
} 


keggList("hsa00601")
keggFind("genes", "hsa00601")

head(keggConv("ncbi-geneid", dfPathways$PathwayID[351]))

head(keggConv("ncbi-geneid", "hsa00601"))

keggLink("hsa00601")



