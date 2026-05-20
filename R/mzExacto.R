#'@title mzExacto
#'
#'@description Uses the output from `spreadOut()` and a list of chemicals
#'to extract data for. Samples that contain a chemical
#'have all identified area(s) aggregated. The most likely identifications
#'are determined by prioritizing matches of exact chemical names, followed
#'by m/z overlaps within precise retention time windows that are determined
#'by molecular masses published on PubChem.
#'
#'@details Communicates with PubChem to collect information on every search
#'chemical. Uses this information to search the
#'raw `spreadOut()` data for samples where a chemical exists.
#'
#'@param data_in the list output from `spreadOut()`
#'
#'@param chemicals A vector containing chemical names in IUPAC notation
#'
#'@param decontaminate A logical set to remove uncertain chemicals or junk from output
#'
#'@returns A data frame containing chemical names(`chemicals`), their
#'optimal retention time, exact mass, best identified match factor, and
#'aggregated component area across every sample it was identified in.
#'
#'@examples
#'\dontrun{
#'query_chemicals = c("Ethyl hexanoate","Methyl salicylate","Octanal","Undecane")
#'mzExacto(standard_spread, query_chemicals)
#'}
#'
#'@importFrom ChemmineR read.SDFset
#'@importFrom webchem get_cid
#'@importFrom jsonlite fromJSON
#'@importFrom stats lm
#'@importFrom stats median
#'@importFrom stats na.omit
#'@importFrom utils adist
#'@export

mzExacto <- function(data_in, chemicals, decontaminate = T){
  options(digits=22)
  getNCI = function(url_path){
    con = url(url_path)
    on.exit(close(con))
    string_holder = tryCatch(readLines(con, warn = F),
                             error = function(error) {return("None")})
    string_holder
  }

  getPubChem = function(url_path){
    con = url(url_path)
    on.exit(close(con))
    string_holder = tryCatch(ChemmineR::read.SDFset(con, warn = F),
                             error = function(error) {return(T)})
    string_holder
  }

  model_fun_use = function(x, a1, b){
    y = a1*x + b
    return(y)
  }

  custom_min = function(x) {if (length(x)>0) min(x) else 0}
  custom_max = function(x) {if (length(x)>0) max(x) else 0}

  chemicals = unique(chemicals)
  area_matrix = data_in$Area
  probs_matrix = data_in$MatchFactor
  cmp_matrix = data_in$Compounds
  RT_matrix = data_in$RT
  mz_matrix = data_in$MZ
  mass_matrix = data_in$Mass
  rtBYmass_matrix = data_in$rtBYmass

  area_df = as.data.frame(area_matrix)
  probs_df = as.data.frame(probs_matrix)
  cmp_df = as.data.frame(cmp_matrix)
  RT_df = as.data.frame(RT_matrix)
  mz_df = as.data.frame(mz_matrix)
  mass_df = as.data.frame(mass_matrix)
  rtBYmass_df = as.data.frame(rtBYmass_matrix)

  cat("Preparing the data. Stay tuned. \n")
  input_area_long = as.numeric(paste0(area_df[!is.na(area_df)]))
  input_probs_long = as.numeric(paste0(probs_df[!is.na(probs_df)]))
  input_CMPs_long = paste0(cmp_df[!is.na(cmp_df)])
  input_RT_long = as.numeric(paste0(RT_df[!is.na(RT_df)]))
  input_MZs_long = as.numeric(paste0(mz_df[!is.na(mz_df)]))
  input_mass_long = suppressWarnings(as.numeric(paste0(mass_df[!is.na(mass_df)])))
  input_rtBYmass_long = paste0(rtBYmass_df[!is.na(rtBYmass_df)])

  bad_masses = c(58.041864811, 999.99)

  RT_mass_lm = lm(input_RT_long[!(input_mass_long %in% bad_masses) & input_probs_long > max(input_probs_long)-25]~input_mass_long[!(input_mass_long %in% bad_masses) & input_probs_long > max(input_probs_long)-25],
                  na.action = na.omit)

  eqn_coefficients = as.numeric(paste0(RT_mass_lm[[1]]))
  RT_step_set = (max(RT_mass_lm[[5]])-min(RT_mass_lm[[5]]))/length(RT_mass_lm[[5]])

  for (w in 1:nrow(cmp_matrix)){
    input_CMPs_tmp = cmp_matrix[w,][!is.na(cmp_matrix[w,])]
    input_CMPs_long = c(input_CMPs_long, input_CMPs_tmp)

    input_mass_tmp = mass_matrix[w,][!is.na(mass_matrix[w,])]
    input_mass_long = c(input_mass_long, input_mass_tmp)

    input_probs_tmp = probs_matrix[w,][!is.na(probs_matrix[w,])]
    input_probs_long = c(input_probs_long, input_probs_tmp)

    input_RT_tmp = RT_matrix[w,][!is.na(RT_matrix[w,])]
    input_RT_long = c(input_RT_long, input_RT_tmp)
  }

  mz_inputs_tmp = as.data.frame(mz_matrix)
  CMP_inputs_tmp = as.data.frame(cmp_matrix)
  RT_inputs_tmp = as.data.frame(RT_matrix)
  mass_inputs_tmp = as.data.frame(mass_matrix)

  mz_inputs = mz_inputs_tmp[input_CMPs_long %in% chemicals,]
  CMP_inputs = CMP_inputs_tmp[input_CMPs_long %in% chemicals,]
  RT_inputs = RT_inputs_tmp[input_CMPs_long %in% chemicals,]
  mass_inputs = mass_inputs_tmp[input_CMPs_long %in% chemicals,]

  input_MZs_long = c()

  for (x in 1:nrow(mz_inputs)){
    input_MZs_tmp = mz_inputs[x,][!is.na(mz_inputs[x,])]
    input_MZs_long = c(input_MZs_long,
                       input_MZs_tmp)
  }

  df_columns = ncol(area_matrix)+3
  df_rows = length(chemicals)
  all_df_last_merge = data.frame(matrix(nrow = df_rows,
                                        ncol = df_columns))
  num_unique_CMPs = length(chemicals)

  final_list = list()
  searchNames_published = list()

  colnames(all_df_last_merge) = c("RT", "Mass",
                                  "Chemical", colnames(area_matrix))

  chem_mass_set = c()
  mz_mass_list = list()

  all_tentative_list = list()
  all_search_list = list()

  cat("Acquiring relevant info from Pubchem [pubchem.ncbi.nlm.nih.gov]. Please be patient!!")
  for (chem in 1:length(chemicals)){

    alt_trigger = F
    current_CMP = chemicals[chem]
    Sys.sleep(1)
    chem_cid = .uaf_get_cid(chemicals[chem])

    if(is.na(chem_cid)){
      smiles_url = paste0("https://cactus.nci.nih.gov/chemical/structure/",
                          current_CMP,
                          "/smiles")
      inchi_url = paste0("https://cactus.nci.nih.gov/chemical/structure/",
                         current_CMP,
                         "/stdinchikey")
      smiles_url = gsub("\\ ",
                        "%20",
                        smiles_url)
      inchi_url = gsub("\\ ",
                       "%20",
                       inchi_url)
      smile_string = getNCI(smiles_url)
      inchi_string = getNCI(inchi_url)

      if(smile_string != "None"){
        InChiKey = substr(inchi_string,
                          10,
                          nchar(inchi_string))
        smile_cid = .uaf_get_cid(paste0(smile_string),
                                  from = "smiles")

        chem_cid = smile_cid
      }else{}


    }
    chem_cid = paste0(chem_cid)

    if(chem_cid == "0"){chem_cid = "180"}

    chem_cid = gsub("NA",
                    "180",
                    chem_cid)

    chem_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/',
                      chem_cid,
                      '/JSON?heading=GC-MS')

    chem_names_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/',
                            chem_cid,
                            '/JSON?heading=Depositor-Supplied+Synonyms')

    chem_names_info = tryCatch(jsonlite::fromJSON(chem_names_url),
                               error = function(error) {return("None")})
    if(paste0(chem_names_info) == "None"){chem_names_json_dat = data.frame(current_CMP)}else{chem_names_json_dat = data.frame(unlist(chem_names_info$Record$Section))}

    colnames(chem_names_json_dat) = "info"

    if(paste0(chem_names_info) == "None"){
      name_finding = ""
      all_current_chem_names = data.frame(paste0(chem_names_json_dat$info))
    }else{
      name_finding = "Section.Section.Information.Value.StringWithMarkup.String"
      all_current_chem_names = data.frame(paste0(chem_names_json_dat$info[adist(rownames(chem_names_json_dat), name_finding) <= 4]))
    }

    colnames(all_current_chem_names) = "names"

    searchNames_published[[chem]] = unique(all_current_chem_names$names)

    chem_info = tryCatch(jsonlite::fromJSON(chem_url),
                         error = function(error) {return("None")})
    if (chem_info == "None" | chem_cid == 180){
      chem_cid = 180

      chem_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/',
                        chem_cid,
                        '/JSON?heading=GC-MS')

      chem_info = tryCatch(jsonlite::fromJSON(chem_url),
                           error = function(error) {return("None")})
      alt_trigger = T
      } else{}

    mass_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',
                      chem_cid,
                      '/SDF')

    mass_info = tryCatch(ChemmineR::read.SDFset(mass_url),
                         error = function(error) {return("None")})
    mass_tmp = mass_info@SDF[[1]]@datablock[names(mass_info@SDF[[1]]@datablock) == "PUBCHEM_EXACT_MASS"]
    CMP_name_tmp = mass_info@SDF[[1]]@datablock[names(mass_info@SDF[[1]]@datablock) == "PUBCHEM_IUPAC_CAS_NAME"]
    if (length(CMP_name_tmp) <1){
      CMP_name_tmp = current_CMP
    } else {}

    chem_json_dat = data.frame(unlist(chem_info$Record$Section$Section))
    colnames(chem_json_dat) = "info"
    mz_rows_1 = grep("[[:digit:]]+\\.[[:digit:]]\\ [[:digit:]]+\\.[[:digit:]]",
                     chem_json_dat$info)
    mz_rows_2 = grep("[[:digit:]]{2,3}\\ [[:digit:]]+\\.[[:digit:]]",
                     chem_json_dat$info)
    mz_rows_3 = grep("[[:digit:]]{2,3}\\ [[:digit:]]{2,3}",
                     chem_json_dat$info)
    mz_rows_4 = grep("^[[:digit:]]{2,3}$",
                     chem_json_dat$info)

    if (length(mz_rows_1) > 0) {
      mz_rows = mz_rows_1
    } else {
      if (length(mz_rows_2) > 0) {
        mz_rows = mz_rows_2
      } else {
        if(length(mz_rows_3) > 0) {
          mz_rows = mz_rows_3
        } else {
          if(length(mz_rows_4) > 0) {
          mz_rows = mz_rows_4
        } else {
          mz_rows = 0
          alt_trigger = T
        }
        }
      }
    }

    chem_mz_matches = unique(paste0(unlist(strsplit(paste0(chem_json_dat[mz_rows,])," "))))
    chem_mz_matches2 = gsub("\\.0\\>",
                            "",
                            chem_mz_matches)
    chem_mz_matches2 = gsub("\\.[[:digit:]]0\\>",
                            "",
                            chem_mz_matches2)

    mz_count = 1:length(chem_mz_matches2)
    odds = mz_count[mz_count%%2 != 0]
    evens = mz_count[mz_count%%2 == 0]

    mz_primary = chem_mz_matches2[odds]
    mz_secondary = chem_mz_matches2[evens]
    chem_mz_matches2 = mz_primary[order(mz_secondary,
                                        decreasing = T)]

    all_current_names = paste0(unique(all_current_chem_names$names))

    if(alt_trigger == T){
      mz_primary = input_MZs_long[input_CMPs_long == input_CMPs_tmp]
      all_current_names = c(current_CMP,
                            CMP_name_tmp)
      mass_tmp = "NA"}

    mz_print = paste(chem_mz_matches2,
                     collapse = " | ")
    total_chems = length(chemicals)

    step_printer = c(1:num_unique_CMPs*100)
    step_cmper = chem/15

    if(step_cmper %in% step_printer | chem == 1){
      cat(paste0('\n',
                 '[Current/Total]',
                 ' |--Exact Mass--|--CMP Name 1--|--CMP Name 2--|--Top MZ Peaks--|',
                 '\n'))
    } else {}

    names(mass_tmp) = NULL
    names(CMP_name_tmp) = NULL
    all_tmp = paste0(mass_tmp,
                     ';',
                     CMP_name_tmp,
                     ';',
                     current_CMP)
    all_print = paste0(mass_tmp,
                       ' | ',
                       CMP_name_tmp,
                       ' | ',
                       current_CMP,
                       ' | ',
                       mz_print)
    cat(paste0('[',
               chem,
               ' / ',
               num_unique_CMPs,
               ']',
               ': ',
               all_print,
               '\n'))

    mz_mass_list[[chem]] = list(mass_tmp,
                                mz_primary)

    all_search_list[[chem]] = list(all_current_names,
                                   mz_primary,
                                   mass_tmp)
  }

  if(length(mz_mass_list) < length(chemicals)){mz_mass_list[[length(chemicals)]] = "NA"}
  names(mz_mass_list) = chemicals

  if(length(searchNames_published) < length(chemicals)){searchNames_published[[length(chemicals)]] = "NA"}
  names(searchNames_published) = chemicals

  if(length(all_search_list) < length(chemicals)){all_search_list[[length(chemicals)]] = "NA"}
  names(all_search_list) = chemicals

  mass_each_chem = sapply(mz_mass_list,
                          '[[',
                          1)
  mz_each_chem = sapply(mz_mass_list,
                        '[[',
                        -1)
  mass_each_chem[mass_each_chem == "NA"] = 999.999

  mass_orderer = as.numeric(paste0(mass_each_chem))

  masses_ordered = mass_each_chem[order(mass_orderer)]
  mz_mass_ordered = mz_each_chem[order(mass_orderer)]
  chems_ordered = names(masses_ordered)

  CMPs_long = c()
  probs_long = c()
  RTs_long = c()
  masses_long = c()
  areas_long = c()
  MZs_long = c()
  rtBYmass_long = c()

  MZs_published = list()
  names_published = list()

  for (z in 1:nrow(cmp_matrix)){
    alt_trigger = F
    CMPs_tmp = cmp_matrix[z,][!is.na(cmp_matrix[z,])]
    CMPs_long = c(CMPs_long, CMPs_tmp)

    probs_tmp = probs_matrix[z,][!is.na(probs_matrix[z,])]
    probs_long = c(probs_long, probs_tmp)

    RTs_tmp = RT_matrix[z,][!is.na(RT_matrix[z,])]
    RTs_long = c(RTs_long, RTs_tmp)

    masses_tmp = mass_matrix[z,][!is.na(mass_matrix[z,])]
    masses_long = c(masses_long, masses_tmp)

    areas_tmp = area_matrix[z,][!is.na(area_matrix[z,])]
    areas_long = c(areas_long, areas_tmp)

    MZs_tmp = mz_matrix[z,][!is.na(mz_matrix[z,])]
    MZs_long = c(MZs_long, MZs_tmp)

    rtBYmass_tmp = rtBYmass_matrix[z,][!is.na(rtBYmass_matrix[z,])]
    rtBYmass_long = c(rtBYmass_long, rtBYmass_tmp)
  }

  all_tentative_list = data_in$webInfo
  names(all_tentative_list) = paste0(RTs_long,
                                     ' | ',
                                     masses_long)

  all_best_RTs = c()
  previous_areas = c()

  row_find_df = data.frame(cbind(CMPs_long,
                                 as.numeric(paste0(RTs_long)),
                                 as.numeric(paste0(probs_long)),
                                 suppressWarnings(as.numeric(paste0(masses_long))),
                                 MZs_long,
                                 as.numeric(paste0(areas_long))))
  rownames(row_find_df) = NULL
  colnames(row_find_df) = c("Chemical", "RT",
                            "Probs", "Mass",
                            "MZ", "Area")

  probs_found = c()
  RT_found = c()
  area_found = c()
  previous_RTs = c()

  masses_ordered = mass_each_chem[order(mass_orderer)]
  mz_mass_ordered = mz_each_chem[order(mass_orderer)]
  chems_ordered = names(masses_ordered)

  RT_range = c()

  for(RT_rt in 1:length(masses_ordered)){
    RT_range_tmp = model_fun_use(as.numeric(masses_ordered[[RT_rt]]),
                                 a1 = eqn_coefficients[2],
                                 b = eqn_coefficients[1])
    RT_range = c(RT_range,
                 RT_range_tmp)
  }
  names(RT_range) = chems_ordered


  rows_list = list()
  previous_rows = c()

  exact_rows_list = list()
  tentative_rows_list = list()


  for(k in 1:length(all_search_list)){
    print_statement = paste0("Working on ",
                             chems_ordered[k],
                             "\n")
    cat(print_statement)
    current_chem = chems_ordered[k]
    all_rows_current = c()
    exact_probs = c()
    exact_RTs = c()
    exact_rows = c()

    tentative_probs = c()
    tentative_RTs = c()
    tentative_rows = c()
    tentative_mz_counts = c()
    tentative_masses = c()

    previous_remover = names(all_tentative_list) %in% previous_rows
    all_tentative_updated = all_tentative_list[!previous_remover]
    list_spot = names(all_search_list) == current_chem

    for(m in 1:length(all_tentative_list)){

      if(any(paste0(all_search_list[list_spot][[1]][[1]]) %in% paste0(all_tentative_list[[m]][[1]])) & !is.null(all_tentative_list[[m]][[1]])){
        exact_rows_tmp = names(all_tentative_list)[m]
        exact_probs_tmp = probs_long[rtBYmass_long == names(all_tentative_list)[m]]
        exact_RTs_tmp = RTs_long[rtBYmass_long == names(all_tentative_list)[m]]


      }else{exact_rows_tmp = exact_probs_tmp = exact_RTs_tmp = NA}
      exact_probs = c(exact_probs, exact_probs_tmp)
      exact_RTs = c(exact_RTs, exact_RTs_tmp)
      exact_rows = c(exact_rows, exact_rows_tmp)
      if(!is.null(all_tentative_list[[m]][[2]])){
        if(is.na(all_tentative_list[[m]][[3]])){next}
        if(k == length(all_search_list)){RT_spot = k}else{RT_spot = k + 1}
        if(as.numeric(paste0(all_tentative_list[[m]][[4]])) <= RT_range[RT_spot] & any(paste0(all_search_list[[k]][[2]]) %in% paste0(all_tentative_list[[m]][[2]]))){
          tentative_rows_tmp = names(all_tentative_list)[m]
          tentative_probs_tmp = probs_long[rtBYmass_long == names(all_tentative_list)[m] & RTs_long <= RT_range[[k]]]
          tentative_RTs_tmp = RTs_long[rtBYmass_long == names(all_tentative_list)[m] & RTs_long <= RT_range[[k]]]
          tentative_mz_counts_tmp1 = paste0(all_search_list[[k]][[2]]) %in% paste0(all_tentative_list[[m]][[2]])
          tentative_mz_counts_tmp = length(tentative_mz_counts_tmp1[tentative_mz_counts_tmp1 == TRUE])
          tentative_masses_tmp = as.numeric(all_tentative_list[[m]][[3]])
        }else{tentative_masses_tmp = tentative_mz_counts_tmp = tentative_rows_tmp = tentative_probs_tmp = tentative_RTs_tmp = NA}
        tentative_probs = c(tentative_probs, tentative_probs_tmp)
        tentative_RTs = c(tentative_RTs, tentative_RTs_tmp)
        tentative_rows = c(tentative_rows, tentative_rows_tmp)
        tentative_mz_counts = c(tentative_mz_counts, tentative_mz_counts_tmp)
        tentative_masses = c(tentative_masses, tentative_masses_tmp)
      }else{next}
    }

    exact_rows_final = unique(exact_rows)
    tentative_rows_final = unique(tentative_rows)

    exact_RTs_final = unique(exact_RTs)
    tentative_RTs_final = unique(tentative_RTs)

    exact_probs_final = unique(exact_probs)
    tentative_probs_final = unique(tentative_probs)

    tentative_masses_final = unique(tentative_masses)

    all_rows_current = c(unique(exact_rows), unique(tentative_rows))
    previous_rows = c(previous_rows, all_rows_current)

    exact_rows_list[[k]] = list(exact_rows_final,
                                exact_RTs_final,
                                exact_probs_final)
    tentative_rows_list[[k]] = list(tentative_rows_final,
                                    tentative_RTs_final,
                                    tentative_probs_final,
                                    tentative_mz_counts,
                                    tentative_masses_final)
  }

  names(exact_rows_list) = chems_ordered
  names(tentative_rows_list) = chems_ordered

  rtBYmass_previous = c()
  previous_RTs = c()

  final_row_all = c()

  final_df_all = data.frame(matrix(nrow = length(chems_ordered),
                                   ncol = ncol(area_matrix)+4))
  colnames(final_df_all) = c("Compound", "Mass",
                             "RT", "Best Match",
                             colnames(area_matrix))

  final_list = list()

  all_areas_list = list()
  all_RTs_list = list()
  all_probs_list = list()

  exclude_exacts_tmp = as.vector(unlist(sapply(exact_rows_list, '[[', 1)))
  exclude_exacts = exclude_exacts_tmp[!is.na(exclude_exacts_tmp)]

  for(p in seq_along(chems_ordered)){
    exact_trigger = F
    area_rows_first = RT_rows_first = probs_rows_first = NULL
    if(any(!is.na(exact_rows_list[names(exact_rows_list) == chems_ordered[p]][[1]][[1]]))){
      area_rows_first = area_matrix[rtBYmass_long %in% exact_rows_list[p][[1]][[1]],]
      area_rows_first[is.na(area_rows_first)] = 0

      RT_rows_first = RT_matrix[rtBYmass_long %in% exact_rows_list[p][[1]][[1]],]
      RT_rows_first[is.na(RT_rows_first)] = 0

      probs_rows_first = probs_matrix[rtBYmass_long %in% exact_rows_list[p][[1]][[1]],]
      probs_rows_first[is.na(probs_rows_first)] = 0

      exact_trigger = T

    }else{
      area_rows_first = rbind(rep(0,ncol(area_matrix)))
      RT_rows_first = rbind(rep(0,ncol(area_matrix)))
      probs_rows_first = rbind(rep(0,ncol(area_matrix)))
      rtBYmass_first = NA
    }
    if(p > 1){rtBYmass_previous = c(rtBYmass_previous, rtBYmass_tmp)}

    area_rows_second = rbind(rep(0,ncol(area_matrix)))
    RT_rows_second = rbind(rep(0,ncol(area_matrix)))
    probs_rows_second = rbind(rep(0,ncol(area_matrix)))

    area_rows_second_tmp = rbind(rep(0,ncol(area_matrix)))
    RT_rows_second_tmp = rbind(rep(0,ncol(area_matrix)))
    probs_rows_second_tmp = rbind(rep(0,ncol(area_matrix)))

    for(r in 1:ncol(area_rows_first)){

      if(sum(as.numeric(area_rows_first[,r])) == 0){
        best_mz_matches_tmp = rtBYmass_long[RTs_long %in% tentative_rows_list[p][[1]][[2]]][!is.na(tentative_rows_list[p][[1]][[2]])]
        best_mz_matches_tmp = best_mz_matches_tmp[!is.na(best_mz_matches_tmp)]
        tentative_counts = tentative_rows_list[p][[1]][[4]]
        tentative_counts[is.na(tentative_counts)] = 0
        count_test = tentative_counts[RTs_long %in% tentative_rows_list[p][[1]][[2]]][!is.na(tentative_rows_list[p][[1]][[2]])]

        if(suppressWarnings(max(count_test, na.rm = T) <= 1)){next}
        best_mz_matches_code_tmp = best_mz_matches_tmp[count_test >= 2]
        best_mz_matches_code = best_mz_matches_code_tmp[!(best_mz_matches_code_tmp %in% exclude_exacts)]
        if(length(best_mz_matches_code) < 1){next}

        rtBYmass_tmp = sapply(tentative_rows_list[p], '[[', 1)
        RT_tmp = sapply(tentative_rows_list[p],
                        '[[',
                        2)
        if(suppressWarnings(is.na(RT_tmp[[2]]))){RT_tmp = as.numeric(RTs_long[rtBYmass_long %in% rtBYmass_tmp])}
        probs_tmp = sapply(tentative_rows_list[p],
                           '[[',
                           3)
        if(suppressWarnings(is.na(probs_tmp[[2]]))){probs_tmp = as.numeric(probs_long[rtBYmass_long %in% rtBYmass_tmp])}
        mass_tmp = sapply(tentative_rows_list[p],
                          '[[',
                          5)
        max_prob = custom_max(as.numeric(probs_tmp[!is.na(probs_tmp)]))
        median_mass = median(as.numeric(mass_tmp[!is.na(mass_tmp)]))
        best_mass = mean(as.numeric(mass_tmp[rtBYmass_tmp %in% best_mz_matches_code]), na.rm = T)

        best_RT = custom_max(as.numeric(RTs_long[rtBYmass_long %in% best_mz_matches_code][!is.na(RTs_long[rtBYmass_long %in% best_mz_matches_code])])) #RT_tmp[probs_tmp == max_prob][!is.na(RT_tmp)]

        previous_RTs = c(previous_RTs,
                         best_RT)

        mass_RT_test = model_fun_use(median_mass,
                                     eqn_coefficients[2],
                                     eqn_coefficients[1])

        if(p == 1){RT_range = as.numeric(RT_tmp) < best_RT + eqn_coefficients[2]}else{RT_range = as.numeric(RT_tmp) > best_RT - 0.15 & as.numeric(RT_tmp) <= best_RT + 0.15} # RT_step_set

        rtBYmass_tmp = rtBYmass_tmp[RT_range]

        if(p == 1){rtBYmass_tmp = rtBYmass_tmp}else{rtBYmass_tmp = rtBYmass_tmp[!(rtBYmass_tmp %in% rtBYmass_previous)]}

        area_rows_filler_tmp = area_matrix[rtBYmass_long %in% rtBYmass_tmp, r]
        area_rows_filler_tmp[is.na(area_rows_filler_tmp)] = 0
        area_rows_filler = sum(as.numeric(area_rows_filler_tmp))

        RT_rows_filler_tmp = RT_matrix[rtBYmass_long %in% rtBYmass_tmp, r]
        RT_rows_filler_tmp[is.na(RT_rows_filler_tmp)] = 0
        RT_rows_filler = custom_min(as.numeric(RT_rows_filler_tmp))

        probs_rows_filler_tmp = probs_matrix[rtBYmass_long %in% rtBYmass_tmp, r]
        probs_rows_filler_tmp[is.na(probs_rows_filler_tmp)] = 0
        probs_rows_filler = custom_max(as.numeric(probs_rows_filler_tmp))

      }else{
        area_rows_filler = RT_rows_filler = probs_rows_filler = 0
      }
      area_rows_second_tmp[,r] = area_rows_filler
      RT_rows_second_tmp[,r] = RT_rows_filler
      probs_rows_second_tmp[,r] = probs_rows_filler
    }
    area_rows_second = rbind(area_rows_second,
                             area_rows_second_tmp)
    RT_rows_second = rbind(RT_rows_second,
                           RT_rows_second_tmp)
    probs_rows_second = rbind(probs_rows_second,
                              probs_rows_second_tmp)

    colnames(area_rows_first) = NULL
    rownames(area_rows_first) = NULL
    area_rows_first = data.frame(area_rows_first)

    colnames(area_rows_second) = NULL
    rownames(area_rows_second) = NULL
    area_rows_second = data.frame(area_rows_second)

    area_rows_all = rbind(area_rows_first,
                          area_rows_second)

    colnames(RT_rows_first) = NULL
    rownames(RT_rows_first) = NULL
    RT_rows_first = data.frame(RT_rows_first)

    colnames(RT_rows_second) = NULL
    rownames(RT_rows_second) = NULL
    RT_rows_second = data.frame(RT_rows_second)

    RT_rows_all = rbind(RT_rows_first,
                        RT_rows_second)

    colnames(probs_rows_first) = NULL
    rownames(probs_rows_first) = NULL
    probs_rows_first = data.frame(probs_rows_first)

    colnames(probs_rows_second) = NULL
    rownames(probs_rows_second) = NULL
    probs_rows_second = data.frame(probs_rows_second)

    probs_rows_all = rbind(probs_rows_first,
                           probs_rows_second)

    all_areas_list[[p]] = area_rows_all
    all_RTs_list[[p]] = RT_rows_all
    all_probs_list[[p]] = probs_rows_all

    final_list[[p]] = list(area_rows_all,
                           RT_rows_all,
                           probs_rows_all)

    final_row = c()

    for(s in 1:ncol(area_rows_all)){
      final_row_tmp = sum(as.numeric(area_rows_all[,s]))
      final_row = c(final_row,
                    final_row_tmp)
    }
    if(any(is.infinite(unlist(RT_rows_all)))){RT_rows_all[sapply(RT_rows_all, simplify = 'matrix', is.infinite)] = 0}

    final_row_all_tmp = cbind(chems_ordered[p],
                              masses_ordered[p],
                              mean(as.numeric(RT_rows_all[RT_rows_all > 0 & !is.na(RT_rows_all)])),
                              custom_max(as.numeric(probs_rows_all[probs_rows_all > 0 & !is.na(probs_rows_all)])),
                              rbind(final_row))
    colnames(final_row_all_tmp) = c("Compound", "Mass",
                                    "RT", "Best Match",
                                    colnames(area_matrix))

    final_df_all[p,] = final_row_all_tmp
    options(digits=7)
  }
  if(decontaminate == T){final_df_all = final_df_all[final_df_all$RT != "NaN",]}
  if(decontaminate == T){final_df_all = final_df_all[final_df_all$Mass != 999.999,]}
  return(final_df_all)
}
