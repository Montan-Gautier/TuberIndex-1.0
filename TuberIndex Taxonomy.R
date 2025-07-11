# TuberIndex Taxonomic Script
# Author: Montan Gautier
# License: MIT License
# Description: A reproducible pipeline to enrich a taxonomic list of common names with NCBI, GBIF and TAXREF identifiers and metadata for use in ecological analysis.

# =====================
# 1. PACKAGE LOADING
# =====================
library(dplyr)
library(tidyr)
library(taxize)
library(rgbif)
library(readr)
library(purrr)
library(tibble)

# =====================
# 2. PARAMETERS
# =====================
#Required key for interrogate NCBI database
#Create an account on https://www.ncbi.nlm.nih.gov
#Go to 'Settings' > 'API Key Management'
#Generate your key and copy it below
# Replace "YOUR_NCBI_API_KEY" with your actual key
ENTREZ_KEY <- "YOUR_NCBI_API_KEY" 
options(ENTREZ_KEY = ENTREZ_KEY)

# =====================
# 3. FUNCTIONS
# =====================
get_ncbi_metadata <- function(taxa) {
  ncbi_ids <- get_uid(taxa)
  uid_vector <- as.numeric(ncbi_ids)
  valid_idx <- which(!is.na(uid_vector))
  valid_ids <- uid_vector[valid_idx]
  valid_taxa <- taxa[valid_idx]
  tax_data <- ncbi_get_taxon_summary(valid_ids)
  metadata <- data.frame(
    input_name = valid_taxa,
    uid = valid_ids,
    name = tax_data$name,
    rank = tax_data$rank,
    stringsAsFactors = FALSE
  )
  
  return(metadata)
}


get_gbif_metadata <- function(taxa_names) {
  gbif_raw <- get_gbifid_(taxa_names, messages = FALSE)
  gbif_keys <- sapply(gbif_raw, function(df) {
    if (is.data.frame(df) && "status" %in% names(df)) {
      accepted <- df[df$status == "ACCEPTED", ]
      if (nrow(accepted) >= 1) return(accepted$usagekey[1]) else return(NA)
    } else {
      return(NA)
    }
  })
  names_from_gbif <- sapply(gbif_keys, function(tk) {
    res <- tryCatch(name_usage(key = as.numeric(tk)), error = function(e) return(NULL))
    if (!is.null(res) && !is.null(res$data)) return(res$data$canonicalName) else return(NA)
  })
  return(data.frame(
    ncbi_scientific_name = taxa_names,
    gbif_tax_id = gbif_keys,
    gbif_scientific_name = names_from_gbif,
    stringsAsFactors = FALSE
  ))
}

get_phylum_family <- function(tid) {
  if (is.na(tid)) return(c(phylum = NA, family = NA))
  out <- tryCatch(
    classification(tid, db = "ncbi"),
    error = function(e) NULL
  )
  if (is.null(out) || is.null(out[[1]])) return(c(phylum = NA, family = NA))
  classif <- out[[1]]
  phylum <- classif$name[classif$rank == "phylum"]
  family <- classif$name[classif$rank == "family"]
  return(c(phylum = phylum, family = family))
}

# =====================
# 4. DATA IMPORT
# =====================
tax_raw <- read_delim("TuberIndex_Taxonomy_raw.csv", delim = ";") %>%
  distinct(common_name, .keep_all = TRUE) %>%
  filter(!is.na(scientific_name) & trimws(scientific_name) != "") %>%
  mutate(across(contains("scientific_name"), trimws))

# =====================
# 5. ENRICH WITH NCBI
# =====================
# a. ncbi id for scientific_name
# ===============================
#If more than one UID is found for a taxon, press "1" or "2" on the console to make a selection

ncbi_data <- get_ncbi_metadata(tax_raw$scientific_name)
ncbi_data <- ncbi_data %>%
  distinct() %>%
  setNames(c("scientific_name", "ncbi_tax_id", "ncbi_scientific_name", "rank")) %>%
  select(-rank)
tax_raw <- tax_raw %>% left_join(ncbi_data, by = "scientific_name")

#If some taxon identifier codes are not compiled, please enter them manually
tax_raw %>%
  filter(is.na(ncbi_tax_id)) %>%
  print(n = Inf)
tax_raw$ncbi_tax_id[tax_raw$scientific_name == "Petrosedum sediforme"] <- "2726413"
tax_raw$ncbi_tax_id[tax_raw$scientific_name == "Petrosedum ochroleucum"] <- "1239910"
tax_raw$ncbi_tax_id[tax_raw$scientific_name == "Sedum altissimum"] <- "2726413"
tax_raw$ncbi_tax_id[tax_raw$scientific_name == "Potentilla viridis"] <- "1016828"

tax_raw$ncbi_scientific_name[tax_raw$scientific_name == "Petrosedum sediforme"] <- "Sedum crassularia"
tax_raw$ncbi_scientific_name[tax_raw$scientific_name == "Petrosedum ochroleucum"] <- "Sedum ochroleucum"
tax_raw$ncbi_scientific_name[tax_raw$scientific_name == "Sedum altissimum"] <- "Sedum crassularia"
tax_raw$ncbi_scientific_name[tax_raw$scientific_name == "Potentilla viridis"] <- "Potentilla verna"

#If you enter some taxon identifier codes manually, please complete the 'manual_database' field.
manual_ncbi <- c("Petrosedum sediforme", "Petrosedum ochroleucum", "Sedum altissimum", "Potentilla viridis")
tax_raw <- tax_raw %>%
  mutate(
    source_ncbi = if_else(scientific_name %in% manual_ncbi, "manual", "auto")
  )

# b. ncbi id for probable_species
# ================================
#If more than one UID is found for a taxon, press "1" or "2" on the console to make a selection
ncbi_data_probable_species <- get_ncbi_metadata(tax_raw$probable_species)
ncbi_data_probable_species <- ncbi_data_probable_species %>%
  distinct() %>%
  setNames(c("probable_species", "ncbi_tax_id_corrected", "ncbi_scientific_name_corrected", "rank")) %>%
  select(-rank)
tax_raw <- tax_raw %>% left_join(ncbi_data_probable_species, by = "probable_species")

#If some taxon identifier codes are not compiled, please enter them manually
tax_raw %>%
  filter(!is.na(probable_species) & is.na(ncbi_scientific_name_corrected)) %>%
  print(n = Inf)

# c. add phylum and family according to the NCBI taxonomy
# ========================================================
tax_ids_unique <- unique(na.omit(tax_raw$ncbi_tax_id))
tax_class_list <- list()
for (tid in tax_ids_unique) {
  Sys.sleep(0.2)
  res <- get_phylum_family(tid)
  tax_class_list[[as.character(tid)]] <- res
}

tax_class_df <- do.call(rbind, tax_class_list) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ncbi_tax_id") %>%
  mutate(ncbi_tax_id = as.numeric(ncbi_tax_id))

tax_raw <- tax_raw %>% mutate(ncbi_tax_id = as.numeric(ncbi_tax_id)) %>%
  left_join(tax_class_df, by = "ncbi_tax_id")

# =====================
# 6. ENRICH WITH GBIF
# =====================
# a. gbif id for scientific_name
# ===============================
gbif_data <- get_gbif_metadata(ncbi_data$ncbi_scientific_name)

#If some taxon identifier codes are not compiled, please enter them manually
print(subset(gbif_data, is.na(gbif_tax_id)))

corrections <- tibble(
  ncbi_scientific_name = c("Betula", "Carduus", "Lonicera", "Pinus", "Coronilla", "Ulmus", "Lathyrus",
                           "Globularia", "Hymenogaster", "Inula conyza", "Iris", "Lens culinaris", "Melanogaster",
                           "Anagallis arvensis", "Mutarda arvensis", "Morus", "Poterium", "Sanguisorba minor", "Pisum sativum",
                           "Robinia", "Salix", "Salix alba", "Scleroderma", "Setaria", "Silene",
                           "Microthlaspi perfoliatum", "Thymus", "Trifolium", "Trichoderma", "Viola", "Galeopsis",
                           "Fumaria"),
  gbif_tax_id = c(2875008, 3114511, 2888645, 2684241, 2943640, 2984510, 10356062, 
                  3233698, 2521519, 10998833, 2748017, 3973001, 2524217,
                  9631539, 3047598, 2984545, 3001016, 5364245, 5347846,
                  2952066, 3039576, 5372513, 2524936, 2702955, 3085897,
                  7639623, 2926960, 2973363, 7703803, 2874237, 8423109,
                  2888567),
  gbif_scientific_name = c("Betula", "Carduus", "Lonicera", "Pinus", "Coronilla", "Ulmus", "Lathyrus",
                           "Globularia", "Hymenogaster", "Pentanema squarrosum", "Iris", "Vicia lens", "Melanogaster",
                           "Lysimachia arvensis", "Sinapis arvensis", "Morus", "Poterium", "Poterium sanguisorba", "Lathyrus oleraceus",
                           "Robinia", "Salix", "Salix alba", "Scleroderma", "Setaria", "Silene",
                           "Noccaea perfoliata", "Thymus", "Trifolium", "Trichoderma", "Viola", "Galeopsis",
                           "Fumaria")
)


gbif_data <- gbif_data %>%
  left_join(corrections, by = "ncbi_scientific_name") %>%
  mutate(
    gbif_tax_id = if_else(is.na(gbif_tax_id.x), gbif_tax_id.y, gbif_tax_id.x),
    gbif_scientific_name = if_else(is.na(gbif_scientific_name.x), gbif_scientific_name.y, gbif_scientific_name.x)
  ) %>%
  select(-gbif_tax_id.x, -gbif_tax_id.y, -gbif_scientific_name.x, -gbif_scientific_name.y)

gbif_data <- gbif_data %>%
  group_by(ncbi_scientific_name) %>%
  summarise(
    gbif_tax_id = first(na.omit(gbif_tax_id)),
    gbif_scientific_name = first(na.omit(gbif_scientific_name)),
    .groups = "drop"
  )

#If some taxon are poorly compiled, please enter them manually
gbif_data_diff <- gbif_data %>%
  filter(ncbi_scientific_name != gbif_scientific_name)

print(gbif_data_diff, n = Inf)

corrections <- tibble(
  ncbi_scientific_name = c("Avena", "Anethum foeniculum", "Artemisia", "Boletus", "Bromopsis", "Cedrus", "Euphorbia",
                           "Festuca pratensis", "Helleborus", "Helvella crispa", "Prunus dulcis", "Russula lepida",
                           "Sedum rupestre"),
  gbif_tax_id = c(2706017, 7224221, 3120641, 8287374, 2703642, 2685742, 11397237,
                  2706237, 2741682, 2554614, 3020791, 2551240,
                  7334609),
  gbif_scientific_name = c("Lagurus","Foeniculum vulgare", "Artemisia", "Boletus", "Bromus", "Cedrus", "Euphorbia", 
                           "Lolium pratense", "Helleborus", "Helvella crispa", "Prunus avium", "Russula rosea",
                           "Petrosedum rupestre")
)

gbif_data <- gbif_data %>%
  filter(!ncbi_scientific_name %in% corrections$ncbi_scientific_name)
gbif_data <- bind_rows(gbif_data, corrections) %>%
  arrange(ncbi_scientific_name)

tax_raw <- tax_raw %>% left_join(gbif_data, by = "ncbi_scientific_name")

#If you enter some taxon identifier codes manually, please complete the 'manual_database' field.
manual_gbif <- c("Betula", "Carduus", "Lonicera", "Pinus", "Coronilla", "Ulmus", "Lathyrus",
                 "Globularia", "Hymenogaster", "Inula conyza", "Iris", "Lens culinaris", "Melanogaster",
                 "Anagallis arvensis", "Mutarda arvensis", "Morus", "Poterium", "Sanguisorba minor", "Pisum sativum",
                 "Robinia", "Salix", "Salix alba", "Scleroderma", "Setaria", "Silene",
                 "Microthlaspi perfoliatum", "Thymus", "Trifolium", "Trichoderma", "Viola", "Galeopsis",
                 "Fumaria", "Lagurus","Foeniculum vulgare", "Artemisia", "Boletus", "Bromus", "Cedrus", "Euphorbia", 
                 "Lolium pratense", "Helleborus", "Helvella crispa", "Prunus avium", "Russula rosea",
                 "Petrosedum rupestre")

tax_raw <- tax_raw %>%
  mutate(
    source_gbif = if_else(scientific_name %in% manual_gbif, "manual", "auto")
  )

# b. gbif id for probable_species
# ================================
gbif_data_probable_species <- get_gbif_metadata(ncbi_data_probable_species$ncbi_scientific_name_corrected)

#If some taxons identifier codes are not compiled, please enter them manually
gbif_data_probable_species <- gbif_data_probable_species %>%
  distinct() %>%
  setNames(c("ncbi_scientific_name_corrected", "gbif_tax_id_corrected", "gbif_scientific_name_corrected"))

gbif_data_probable_species <- gbif_data_probable_species %>%
  group_by(ncbi_scientific_name_corrected) %>%
  summarise(
    gbif_tax_id_corrected = first(na.omit(gbif_tax_id_corrected)),
    gbif_scientific_name_corrected = first(na.omit(gbif_scientific_name_corrected)),
    .groups = "drop"
  )

tax_raw <- tax_raw %>% left_join(gbif_data_probable_species, by = "ncbi_scientific_name_corrected")

# =====================
# 7. JOIN WITH TAXREF
# =====================
#Required TAXREF V18.0 or more recent version
#Download available at : https://inpn.mnhn.fr/telechargement/referentielEspece/referentielTaxo

taxref <- read.delim("TAXREFv18.txt") %>% filter(REGNE %in% c("Fungi", "Plantae"))

# a. taxref id for scientific_name
# ===============================
gbif_data <- left_join(gbif_data,taxref[, c("LB_NOM", "CD_REF")] %>%
    rename(taxref_tax_scientific_name = LB_NOM, taxref_tax_id = CD_REF),by = c("gbif_scientific_name" = "taxref_tax_scientific_name"))

#If some taxons have multiple CD_REF, please enter them manually according to : https://taxref.mnhn.fr
duplicated_rows <- gbif_data[
  duplicated(gbif_data[, c("gbif_scientific_name", "taxref_tax_id")]) == FALSE &
    duplicated(gbif_data$gbif_scientific_name, fromLast = TRUE), 
]

gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Arabis hirsuta"] <- "83332"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Bothriochloa ischaemum"] <- "86169"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Cardamine hirsuta"] <- "87930"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Carduus crispus"] <- "88104"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Cedrus"] <- "190425"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Crepis vesicaria"] <- "93157"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Daucus carota"] <- "133744"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Festuca glauca"] <- "98258"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Lolium pratense"] <- "121479"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Juncus acutus"] <- "104104"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Limodorum"] <- "194157"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Morchella esculenta"] <- "48691"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Myosotis"] <- "195001"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Orchis purpurea"] <- "110966"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Orchis simia"] <- "110987"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Pinus nigra"] <- "113683"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Prunus armeniaca"] <- "116041"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Pyrus communis"] <- "116574"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Quercus pubescens"] <- "116751"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Russula maculata"] <- "39897"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Russula rosea"] <- "463281"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Salix caprea"] <- "119977"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Scleroderma verrucosum"] <- "40723"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Syringa"] <- "198159"
gbif_data$taxref_tax_id[gbif_data$gbif_scientific_name == "Tuber aestivum"] <- "861413"

gbif_data<-unique(gbif_data)

tax_raw <- tax_raw %>% 
  left_join(gbif_data %>% 
              select(ncbi_scientific_name, taxref_tax_id), by = "ncbi_scientific_name")

taxref_data <- taxref %>% 
  select(NOM_VALIDE, CD_REF) %>% 
  rename(taxref_scientific_name = NOM_VALIDE) %>% 
  unique()

tax_raw <- tax_raw %>% left_join(taxref_data %>% mutate(CD_REF = as.character(CD_REF)), by = c("taxref_tax_id" = "CD_REF"))

#If you enter some taxon identifier codes manually, please complete the 'manual_database' field.
manual_taxref <- c("Arabis hirsuta", "Bothriochloa ischaemum", "Cardamine hirsuta", "Carduus crispus", "Cedrus", "Crepis vesicaria",
                   "Daucus carota", "Festuca glauca", "Lolium pratense", "Juncus acutus", "Limodorum", "Morchella esculenta",
                   "Myosotis", "Orchis purpurea", "Orchis simia", "Pinus nigra", "Prunus armeniaca", "Pyrus communis",
                   "Quercus pubescens", "Russula maculata", "Russula rosea", "Salix caprea", "Scleroderma verrucosum",
                   "Syringa", "Tuber aestivum")

tax_raw <- tax_raw %>%
  mutate(
    source_taxref = if_else(scientific_name %in% manual_taxref, "manual", "auto")
  )


# b. taxref id for probable_species
# ================================
gbif_data_probable_species <- left_join(gbif_data_probable_species,taxref[, c("LB_NOM", "CD_REF")] %>%
                         rename(taxref_tax_scientific_name_corrected = LB_NOM, taxref_tax_id_corrected = CD_REF),by = c("ncbi_scientific_name_corrected" = "taxref_tax_scientific_name_corrected"))

#If some taxons have multiple CD_REF, please enter them manually according to : https://taxref.mnhn.fr
duplicated_rows <- gbif_data_probable_species[
  duplicated(gbif_data_probable_species[, c("gbif_scientific_name_corrected", "taxref_tax_id_corrected")]) == FALSE &
    duplicated(gbif_data_probable_species$gbif_scientific_name_corrected, fromLast = TRUE), 
]

gbif_data_probable_species$taxref_tax_id_corrected[gbif_data_probable_species$gbif_scientific_name_corrected == "Pyrus communis"] <- "116574"
gbif_data_probable_species$taxref_tax_id_corrected[gbif_data_probable_species$gbif_scientific_name_corrected == "Rubus fruticosus"] <- "119097"

gbif_data_probable_species<-unique(gbif_data_probable_species)
tax_raw <- tax_raw %>% left_join(gbif_data_probable_species %>% select(ncbi_scientific_name_corrected, taxref_tax_id_corrected), by = "ncbi_scientific_name_corrected")

taxref_data <- taxref %>% 
  select(NOM_VALIDE, CD_REF) %>% 
  rename(taxref_scientific_name_corrected = NOM_VALIDE) %>% 
  unique()
tax_raw <- tax_raw %>% left_join(taxref_data %>% mutate(CD_REF = as.character(CD_REF)), by = c("taxref_tax_id_corrected" = "CD_REF"))

 # =====================
# 8. EXPORT
# =====================
write_csv(tax_raw, "TuberIndex_Taxonomy.csv")
