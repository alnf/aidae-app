# reverse_deps <- function(packagename) {
#   tools::package_dependencies(packagename,
#                               db = installed.packages(),
#                               which = c("Depends", "Imports", "LinkingTo"),
#                               reverse = TRUE)
# }
# 
# reverse_deps("shiny")
# 
# get_ora <- function(gmt_path, t2gh, ccc) {
#   gmt <- CEMiTool::read_gmt(gmt_path)
#   gmt <- merge(gmt, t2gh, by.x="gene", by.y="hsymbol")
#   gmt <- gmt[,c(2,3)]
#   colnames(gmt)[2] <- "gene"
#   cmm <- CEMiTool::mod_ora(ccc, gmt)
#   ora <- cmm@ora
#   return(ora)
# }
# 
# save_ora <- function(ontologies, t2gh, ccc, view) {
#   for (ontology in names(ontologies)) {
#     myora <- get_ora(ontologies[[ontology]], t2gh, ccc)
#     save(myora, file = paste0("data/pathways/", view, "/", ontology, ".RData", sep=""))
#   }
# }
# 
# save_ora(ontologies, t2gh, cem_norm, "norma")
# save_ora(ontologies, t2gh, cem_PE, "PE")


# fetchGenesInfo <- function(t2g, species) {
# 
#   mapping <- mygene::queryMany(
#     unique(t2g$ens_gene[1:20]),
#     scopes = "ensembl.gene",
#     fields = c("entrezgene", "type_of_gene"),
#     species = "rat"
#   )
# 
#   print(mapping)
#   mapping$entrezgene[!is.na(mapping$entrezgene)]
# 
#   summary <- tryCatch({
#     s <- rentrez::entrez_summary(db = "gene", id =   mapping$entrezgene[!is.na(mapping$entrezgene)])
#   }, error = function(e) NA_character_)
# 
#   print(summary)
# 
# 
#   entrez_id = entrez_id,
#   name = summary[["name"]],
#   chromosome = summary[["chromosome"]],
#   map_location = summary[["maplocation"]],
#   summary = summary[["summary"]]  
# 
# }
# 
# fetchGenesInfo(t2g, "rat")
