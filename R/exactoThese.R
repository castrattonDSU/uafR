#'@title exactoThese
#'
#'@description Takes categorated output as input and makes chemical subsets by the information in it.
#'
#'@details Provides a set of search chemicals for `mzExacto()`. User gets to select
#'chemicals based on information generated from `categorate()`. Database
#'subsetting accepts both historical `Databases` output and current split
#'database tables from `categorate()`.
#'
#'@param categoratedInput the lists output from `categorate()`
#'@param subsetBy specifies the list to subset by (Database, FMCS, or Library)
#'@param subsetArgs additional arguments to specify which values to subset by
#'@param subsetArgs2 additional arguments to specify which values to subset by
#'@param subset_input used when subsetting by FMCS information (e.g. molecular weight)
#'
#'@returns A vector of chemical names that meet each user-specification. Used as input
#'for `mzExacto()`.
#'
#'@examples
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "All")
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "reactives")
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "LOTUS")
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "KEGG")
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "FEMA")
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "FDA_SPL")
#'exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = c("reactives", "FEMA"))
#'exactoThese(standard_categorated, subsetBy = "FMCS", subsetArgs = "MW",
#'subsetArgs2 = "Between", subset_input = c(125, 200))
#'@export

exactoThese = function(categoratedInput, subsetBy = "Database", subsetArgs = "All", subsetArgs2 = NA, subset_input = NA){
  subset_opts = c("Database", "FMCS", "Library")
  if(!(subsetBy %in% subset_opts)){stop("First Subset Argument Unrecognized!\n Your options are: 1. Database, 2. FMCS, 3. Library")}
  if(subsetBy == "Database"){
    db_sets = .exacto_database_sets(categoratedInput)
    db_names = names(db_sets)
    if (identical(subsetArgs, "All")) {
      exactoChems = unique(unlist(db_sets, use.names = FALSE))
    } else {
      if (any(!(subsetArgs %in% db_names))) {
        stop("Database Subset Arguments Unrecognized, Please Try Again!\n Your options are: All, reactives, LOTUS, KEGG, FEMA, FDA_SPL")
      }
      exactoChems = Reduce(intersect, db_sets[subsetArgs])
    }
  }
  if(subsetBy == "FMCS"){
    FMCS_subsetArgs = list("MW", "Rings", "Groups",
                           "Atoms", "NCharges")
    FMCS_subsetArgs2 = list("Equals", "Greater Than", "Less Than", "Between")
    FMCS_subset_val = subset_input


    if(subsetArgs == FMCS_subsetArgs[[1]] & subsetArgs2 == FMCS_subsetArgs2[[1]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$MW != "None" & as.numeric(paste0(categoratedInput$FMCS$MW)) == FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[1]] & subsetArgs2 == FMCS_subsetArgs2[[2]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$MW != "None" & as.numeric(paste0(categoratedInput$FMCS$MW)) > FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[1]] & subsetArgs2 == FMCS_subsetArgs2[[3]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$MW != "None" & as.numeric(paste0(categoratedInput$FMCS$MW)) < FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[1]] & subsetArgs2 == FMCS_subsetArgs2[[4]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$MW != "None" & as.numeric(paste0(categoratedInput$FMCS$MW)) > FMCS_subset_val[1] & as.numeric(paste0(categoratedInput$FMCS$MW)) < FMCS_subset_val[2]]}

    if(subsetArgs == FMCS_subsetArgs[[2]] & subsetArgs2 == FMCS_subsetArgs2[[1]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Rings != "None" & as.numeric(paste0(categoratedInput$FMCS$Rings)) == FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[2]] & subsetArgs2 == FMCS_subsetArgs2[[2]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Rings != "None" & as.numeric(paste0(categoratedInput$FMCS$Rings)) > FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[2]] & subsetArgs2 == FMCS_subsetArgs2[[3]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Rings != "None" & as.numeric(paste0(categoratedInput$FMCS$Rings)) < FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[2]] & subsetArgs2 == FMCS_subsetArgs2[[4]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Rings != "None" & as.numeric(paste0(categoratedInput$FMCS$Rings)) > FMCS_subset_val[1] & as.numeric(paste0(categoratedInput$FMCS$Rings)) < FMCS_subset_val[2]]}

    if(subsetArgs == FMCS_subsetArgs[[3]] & subsetArgs2 == FMCS_subsetArgs2[[1]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Groups != "None" & as.numeric(paste0(categoratedInput$FMCS$GroupCounts)) == FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[3]] & subsetArgs2 == FMCS_subsetArgs2[[2]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Groups != "None" & as.numeric(paste0(categoratedInput$FMCS$GroupCounts)) > FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[3]] & subsetArgs2 == FMCS_subsetArgs2[[3]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Groups != "None" & as.numeric(paste0(categoratedInput$FMCS$GroupCounts)) < FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[3]] & subsetArgs2 == FMCS_subsetArgs2[[4]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Groups != "None" & as.numeric(paste0(categoratedInput$FMCS$GroupCounts)) > FMCS_subset_val[1] & as.numeric(paste0(categoratedInput$FMCS$GroupCounts)) < FMCS_subset_val[2]]}

    if(subsetArgs == FMCS_subsetArgs[[4]] & subsetArgs2 == FMCS_subsetArgs2[[1]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Atom != "None" & as.numeric(paste0(categoratedInput$FMCS$AtomCounts)) == FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[4]] & subsetArgs2 == FMCS_subsetArgs2[[2]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Atom != "None" & as.numeric(paste0(categoratedInput$FMCS$AtomCounts)) > FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[4]] & subsetArgs2 == FMCS_subsetArgs2[[3]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Atom != "None" & as.numeric(paste0(categoratedInput$FMCS$AtomCounts)) < FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[4]] & subsetArgs2 == FMCS_subsetArgs2[[4]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$Atom != "None" & as.numeric(paste0(categoratedInput$FMCS$AtomCounts)) > FMCS_subset_val[1] & as.numeric(paste0(categoratedInput$FMCS$AtomCounts)) < FMCS_subset_val[2]]}

    if(subsetArgs == FMCS_subsetArgs[[5]] & subsetArgs2 == FMCS_subsetArgs2[[1]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$NCharges != "None" & as.numeric(paste0(categoratedInput$FMCS$NCharges)) == FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[5]] & subsetArgs2 == FMCS_subsetArgs2[[2]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$NCharges != "None" & as.numeric(paste0(categoratedInput$FMCS$NCharges)) > FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[5]] & subsetArgs2 == FMCS_subsetArgs2[[3]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$NCharges != "None" & as.numeric(paste0(categoratedInput$FMCS$NCharges)) < FMCS_subset_val]}
    if(subsetArgs == FMCS_subsetArgs[[5]] & subsetArgs2 == FMCS_subsetArgs2[[4]]){exactoChems = categoratedInput$FMCS$Chemical[categoratedInput$FMCS$NCharges != "None" & as.numeric(paste0(categoratedInput$FMCS$NCharges)) > FMCS_subset_val[1] & as.numeric(paste0(categoratedInput$FMCS$NCharges)) < FMCS_subset_val[2]]}
  }
  if(subsetBy == "Library"){
    library_groups = colnames(categoratedInput$FunctionalGroups)
    library_logical = library_groups %in% subsetArgs
    exactoChems = categoratedInput$FunctionalGroups$Chemical[categoratedInput$FunctionalGroups[,library_logical] != "No"]
  }
  tryCatch(return(unique(exactoChems[!is.na(exactoChems)])), error = function(error) {stop("Database Subset Arguments Unrecognized, Please Try Again!\n If subsetting by multiple, their order matters: \n1. reactives, 2. LOTUS, 3. KEGG, 4. FEMA, 5. FDA_SPL")})
}

.exacto_database_sets = function(categoratedInput) {
  database_names = c("reactives", "LOTUS", "KEGG", "FEMA", "FDA_SPL")
  out = stats::setNames(vector("list", length(database_names)), database_names)

  if ("Databases" %in% names(categoratedInput)) {
    databases = categoratedInput$Databases
    column_map = c(reactives = "reactives_df",
                   LOTUS = "LOTUS_df",
                   KEGG = "KEGG_df",
                   FEMA = "FEMA_df",
                   FDA_SPL = "FDA_SPL_df")
    for (database in database_names) {
      out[[database]] = .exacto_non_none_chemicals(databases,
                                                   column_map[[database]])
    }
    return(out)
  }

  for (database in database_names) {
    table = categoratedInput[[database]]
    out[[database]] = .exacto_non_none_chemicals(table, database)
  }
  out
}

.exacto_non_none_chemicals = function(table, value_column) {
  if (!is.data.frame(table) || !("Chemical" %in% colnames(table)) ||
      !(value_column %in% colnames(table))) {
    return(character())
  }
  values = trimws(paste0(table[[value_column]]))
  keep = !is.na(table$Chemical) & !is.na(values) &
    values != "" & values != "None" & values != "NA"
  unique(paste0(table$Chemical[keep]))
}
