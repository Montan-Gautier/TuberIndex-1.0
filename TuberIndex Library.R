# TuberIndex Library Script
# Author: Montan Gautier
# License: MIT License
# Description: This script processes the raw Zotero export of the TuberIndex library by cleaning metadata and generating a cleaned dataset
# =====================
# 1. PACKAGE LOADING
# =====================
library(dplyr)
library(stringr)
library(readr)

# =====================
# 2. DATA IMPORT
# =====================
lib_raw <- read_delim("TuberIndex_Library_raw.csv", delim = ",")

# =====================
# 3. DATA PROCESSING
# =====================

colnames(lib_raw) <- make.names(colnames(lib_raw))

lib_raw <- lib_raw %>% subset(select = c(Key, Item.Type, Publication.Year, Author, Title, Publication.Title,
  Pages, Issue, Volume, Series, Publisher, Language, Type, Extra, Manual.Tags)) %>%
  rename(Doc_ID = Key, Item_Type = Item.Type, Publication_Year = Publication.Year, Publication_Title = Publication.Title,
  Manual_Tags = Manual.Tags, Thesis_Type = Type, Available = Extra) %>%
  subset(Language == "fr") 

colnames(lib_raw) <- tolower(colnames(lib_raw))

lib_raw <- lib_raw %>%
  mutate(source = str_extract_all(manual_tags, "source: *[^;]+") %>% 
      lapply(str_remove_all, "source: *") %>% 
      sapply(paste, collapse = ";"),
  topic = str_extract_all(manual_tags, "topic: *[^;]+") %>%
      lapply(str_remove_all, "topic: *") %>% 
      sapply(paste, collapse = ";")) %>%
  subset(select = c(doc_id, item_type, publication_year, author, title, publication_title,
                                                          pages, issue, volume, series, publisher, language, thesis_type, available, source, topic))

lib_raw$item_type<- str_to_sentence(lib_raw$item_type)

lib_raw <- subset(lib_raw, item_type %in% c("Book", "Journalarticle", "Magazinearticle", "Newspaperarticle", "Report", "Thesis", "Videorecording"))

# =====================
# 4. DATA PROCESSING
# =====================
write_csv(lib_raw, "TuberIndex_Library.csv")
