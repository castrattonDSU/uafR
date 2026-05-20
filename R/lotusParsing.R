.lotus_extract_organism = function(value, context = "") {
  .uaf_first_non_empty_text(.lotus_extract_organisms(value, context))
}

.lotus_extract_organisms = function(value, context = "", markup_text = NA_character_) {
  value = .uaf_squish_text(value)
  context = .uaf_squish_text(context)
  if (is.na(value) || value == "") return(character())

  labeled = .lotus_label_values(value, c("organism", "source organism",
                                         "biological source", "species",
                                         "taxon name"))
  markup_candidates = .lotus_clean_taxon_terms(markup_text)
  reported = .lotus_extract_reported_organisms(value)
  candidates = c(markup_candidates, reported, labeled,
                 .lotus_extract_binomial(value))
  if (.lotus_context_is(context, c("organism", "species", "biological source",
                                   "source organism", "taxon name"))) {
    candidates = c(candidates, value)
  }
  candidates = .lotus_clean_taxon_terms(candidates)
  candidates = candidates[vapply(candidates, .lotus_is_organism_name,
                                 logical(1))]
  unique(candidates)
}

.lotus_extract_taxonomy = function(value, context = "") {
  value = .uaf_squish_text(value)
  context = .uaf_squish_text(context)
  if (is.na(value) || value == "") return(NA_character_)

  terms = character()
  if (.lotus_context_is(context, c("taxonomy", "taxonomic", "lineage",
                                   "classification", "kingdom", "phylum",
                                   "class", "order", "family", "genus"))) {
    terms = c(terms, .lotus_clean_taxon_terms(value))
  }

  for (label in c("taxonomy", "taxonomic lineage", "lineage",
                  "kingdom", "phylum", "class", "order", "family",
                  "genus", "species")) {
    terms = c(terms, .lotus_label_values(value, label))
  }
  terms = .lotus_clean_taxon_terms(terms)
  .pubchem_collapse(terms)
}

.lotus_extract_family = function(text, terms = character()) {
  labeled = .lotus_label_values(text, "family")
  candidates = .lotus_clean_taxon_terms(c(labeled, terms))
  candidates = candidates[grepl("^[A-Z][A-Za-z-]*(aceae|idae|ceae)$",
                                candidates)]
  .uaf_first_non_empty_text(candidates)
}

.lotus_extract_rank = function(text, terms, rank) {
  labeled = .lotus_label_values(text, rank)
  candidate_pool = if (rank %in% c("genus", "species")) {
    labeled
  } else {
    c(labeled, terms)
  }
  if (rank == "species") {
    candidate_pool = c(candidate_pool, .lotus_extract_binomial(text))
  }
  candidates = .lotus_clean_taxon_terms(candidate_pool)
  kingdom_pattern = "^(Plantae|Fungi|Bacteria|Animalia|Archaea|Chromista|Protozoa|Viridiplantae)$"
  pattern = switch(rank,
                   kingdom = kingdom_pattern,
                   phylum = "(phyta|mycota|chordata|arthropoda|proteobacteria|firmicutes|actinobacteria)$",
                   class = "(opsida|phyceae|mycetes|bacilli|clostridia|actinomycetia|alphaproteobacteria|betaproteobacteria|gammaproteobacteria|deltaproteobacteria)$",
                   order = "(ales|formes)$",
                   genus = "^[A-Z][a-z-]{2,}$",
                   species = "^[A-Z][a-z-]{2,}\\s+[a-z][a-z-]{2,}",
                   ".*")
  candidates = candidates[grepl(pattern, candidates, ignore.case = FALSE)]
  if (rank != "kingdom") {
    candidates = candidates[!grepl(kingdom_pattern, candidates,
                                   ignore.case = FALSE)]
  }
  if (rank == "genus") {
    candidates = candidates[!grepl("(aceae|idae|ceae|ales|formes|phyta|mycota)$",
                                   candidates)]
  }
  if (rank == "kingdom") {
    candidates = sub("^Viridiplantae$", "Plantae", candidates)
  }
  .uaf_first_non_empty_text(candidates)
}

.lotus_extract_taxonomy_terms = function(...) {
  values = .uaf_non_empty(unlist(list(...), use.names = FALSE))
  if (length(values) < 1) return(character())
  terms = unlist(lapply(values, .lotus_clean_taxon_terms), use.names = FALSE)
  unique(terms)
}

.lotus_lineage_ranks = function(lineage, organism = NA_character_) {
  terms = .lotus_clean_taxon_terms(lineage)
  genus_species = .normalized_genus_species(organism)
  genus = genus_species$Genus
  species = genus_species$Species
  if (is.na(genus) && length(terms) > 0) {
    family = .lotus_extract_family(lineage, terms)
    family_index = if (!is.na(family)) match(family, terms) else NA_integer_
    genus_pool = if (!is.na(family_index) && family_index < length(terms)) {
      terms[(family_index + 1):length(terms)]
    } else {
      character()
    }
    genus_candidates = genus_pool[grepl("^[A-Z][a-z-]{2,}$", genus_pool)]
    genus_candidates = genus_candidates[
      !grepl("(aceae|idae|ceae|ales|formes|phyta|mycota|opsida|phyceae|mycetes)$",
             genus_candidates)
    ]
    genus = .uaf_first_non_empty_text(utils::tail(genus_candidates, 1))
  }
  domain = .uaf_first_non_empty_text(
    terms[grepl("^(Eukaryota|Bacteria|Archaea)$", terms)]
  )
  list(
    Domain = domain,
    Kingdom = .lotus_extract_rank(lineage, terms, "kingdom"),
    Phylum = .lotus_extract_rank(lineage, terms, "phylum"),
    Class = .lotus_extract_rank(lineage, terms, "class"),
    Order = .lotus_extract_rank(lineage, terms, "order"),
    Family = .lotus_extract_family(lineage, terms),
    Genus = genus,
    Species = species,
    Terms = terms
  )
}

.lotus_clean_taxon_terms = function(x) {
  x = .uaf_non_empty(x)
  if (length(x) < 1) return(character())
  x = gsub("\\b(?:taxonomy|taxonomic lineage|lineage|kingdom|phylum|class|order|family|genus|species|organism|source organism|biological source|taxon name)\\s*[:=]?",
           "", x, ignore.case = TRUE, perl = TRUE)
  pieces = unlist(strsplit(paste(x, collapse = "; "),
                           "\\s*(?:;|\\||>|/|,\\s*(?=[A-Z][a-z]))\\s*",
                           perl = TRUE), use.names = FALSE)
  pieces = .uaf_non_empty(pieces)
  pieces = gsub("^[-:]+\\s*", "", pieces)
  pieces = gsub("\\s*[-:]+$", "", pieces)
  pieces = .uaf_squish_text(pieces)
  pieces = pieces[!vapply(pieces, .lotus_is_junk_taxon, logical(1))]
  unique(pieces)
}

.lotus_label_values = function(text, labels) {
  text = .uaf_squish_text(text)
  if (is.na(text) || text == "") return(character())
  out = character()
  for (label in labels) {
    escaped = .lotus_regex_escape(label)
    pattern = paste0("\\b", escaped, "\\s*[:=]?\\s*([^;|,.]+(?:\\s+[A-Za-z-]+)?)")
    matches = gregexpr(pattern, text, perl = TRUE, ignore.case = TRUE)
    values = regmatches(text, matches)[[1]]
    if (length(values) < 1 || identical(values, character(0))) next
    for (value in values) {
      parts = regmatches(value, regexec(pattern, value, perl = TRUE,
                                        ignore.case = TRUE))[[1]]
      if (length(parts) >= 2) out = c(out, parts[[2]])
    }
  }
  .lotus_clean_taxon_terms(out)
}

.lotus_extract_binomial = function(text) {
  text = .uaf_squish_text(text)
  if (is.na(text) || text == "") return(character())
  pattern = "\\b[A-Z][a-z]{2,}\\s+(?:x\\s+)?[a-z][a-z-]{2,}(?:\\s+(?:subsp\\.|ssp\\.|var\\.|f\\.)\\s+[a-z-]+)?\\b"
  matches = .uaf_extract_pattern(text, pattern)
  matches = matches[!vapply(matches, .lotus_is_junk_taxon, logical(1))]
  matches = matches[vapply(matches, .lotus_is_organism_name, logical(1))]
  unique(matches)
}

.lotus_is_organism_name = function(x) {
  x = .uaf_squish_text(x)
  if (is.na(x) || x == "") return(FALSE)
  if (!grepl("^[A-Z][a-z-]{2,}\\s+(?:x\\s+)?[a-z][a-z-]{2,}", x)) {
    return(FALSE)
  }
  parts = strsplit(x, "\\s+", perl = TRUE)[[1]]
  second = tolower(parts[[min(length(parts), 2)]])
  !second %in% c("has", "was", "were", "is", "are", "and", "or", "with",
                 "without")
}

.lotus_is_junk_taxon = function(x) {
  x = .uaf_squish_text(x)
  if (is.na(x) || x == "") return(TRUE)
  text = tolower(x)
  if (nchar(x) > 80) return(TRUE)
  if (length(strsplit(x, "\\s+", perl = TRUE)[[1]]) > 6) return(TRUE)
  if (grepl("^(none|unknown|not available|null|na|all taxonomy db|cellular organisms)$",
            text)) return(TRUE)
  if (grepl("consolidatedcompoundtaxonomy|pcget_taxonomy|^hid$",
            text)) return(TRUE)
  if (grepl("https?://|www\\.|doi\\b|pubmed|pmid|reference|citation|journal|article",
            text)) return(TRUE)
  if (grepl("\\b(compound|chemical|metabolite|isolated|reported|found|occurrence|database|lotus|natural product class|structure|smiles|inchi|sourceid|record)\\b",
            text)) return(TRUE)
  if (grepl("\\b(?:LTS|LOTUS)[A-Z0-9:-]*[0-9][A-Z0-9:-]*\\b", x)) return(TRUE)
  if (grepl("^[0-9 .:-]+$", x)) return(TRUE)
  FALSE
}

.lotus_extract_reported_organisms = function(text) {
  text = .uaf_squish_text(text)
  if (is.na(text) || text == "") return(character())
  pattern = "\\b(?:reported|found|detected|isolated)\\s+in\\s+(.+?)(?:\\s+with\\s+data\\s+available|\\s+with\\s+|\\.$|$)"
  match = regmatches(text, regexec(pattern, text, perl = TRUE,
                                   ignore.case = TRUE))[[1]]
  if (length(match) < 2) return(character())
  organisms = match[[2]]
  organisms = gsub("\\band\\s+other\\s+organisms\\b.*$", "",
                   organisms, ignore.case = TRUE, perl = TRUE)
  organisms = gsub("\\bother\\s+organisms\\b.*$", "",
                   organisms, ignore.case = TRUE, perl = TRUE)
  pieces = unlist(strsplit(organisms, "\\s*(?:,|;|\\band\\b)\\s*",
                           perl = TRUE), use.names = FALSE)
  pieces = .lotus_clean_taxon_terms(pieces)
  pieces[vapply(pieces, .lotus_is_organism_name, logical(1))]
}

.lotus_context_is = function(context, labels) {
  context = tolower(.uaf_squish_text(context))
  if (is.na(context) || context == "") return(FALSE)
  any(vapply(labels, function(label) {
    grepl(.lotus_regex_escape(label), context, fixed = FALSE)
  }, logical(1)))
}

.lotus_regex_escape = function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}
