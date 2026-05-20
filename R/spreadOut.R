#'@title spreadOut
#'
#'@description Prepares input mass spectrometry data for subsequent functions.
#'Relevant information for every instance of every chemical across
#'all samples is placed on the same matrices with no overlaps. Chemicals
#'are also searched on PubChem and all published names, published m/z
#'histograms, and exact mass data are downloaded and stored as a large
#'list.
#'
#'@details Automatically prepares .CSV input containing retention time,
#'detected m/z, component area, compound name, match factor, and sample
#'or file name for downstream uafR functions.
#'
#'@param input A .CSV file containing raw mass spectrometry output with 'Component.RT',
#''Component.Area', 'Base.Peak.MZ', 'File.Name', 'Compound.Name', and 'Match.Factor' as
#'column names/required input
#'
#'@returns Returns a list with matrices storing raw area values, chemical identities,
#'m/z values, match factors, retention times, exact mass data (if published),
#'and a unique code for each data point (retention time pasted to exact mass).
#'The list also contains an additional list with all published chemical names,
#'top m/z peaks, exact mass, and structural fingerprints in SDF format.
#'
#'@examples
#'\dontrun{
#'spreadOut(standard_data)
#'}
#'
#'@importFrom ChemmineR read.SDFset
#'@importFrom webchem get_cid
#'@importFrom jsonlite fromJSON
#'@importFrom utils adist
#'@export

spreadOut = function(input){
   getNCI = function(url_path){
      con = url(url_path)
      on.exit(close(con))
      string_holder = tryCatch(readLines(con, warn = F),
                               error = function(error) {return("None")})
      suppressWarnings(string_holder)
   }

   getPubChem = function(url_path){
      con = url(url_path)
      on.exit(close(con))
      string_holder = tryCatch(ChemmineR::read.SDFset(con, warn = F),
                               error = function(error) {return(T)})
      suppressWarnings(string_holder)
   }
   required_clmns = c("Component.RT", "Component.Area", "Base.Peak.MZ",
                      "File.Name", "Compound.Name", "Match.Factor")
   if(!is.data.frame(input)){
      stop("`input` must be a data frame containing GC/MS output columns.")
   }
   missing_clmns = setdiff(required_clmns, colnames(input))
   if(length(missing_clmns) > 0){
      stop(paste0("Missing required input columns: ",
                  paste(missing_clmns, collapse = ", ")))
   }
   peak_keys = paste(input$Component.RT,
                     input$Component.Area,
                     input$Base.Peak.MZ,
                     sep = " | ")
   if(any(duplicated(peak_keys))){
      stop("Input contains duplicate peak keys based on Component.RT, Component.Area, and Base.Peak.MZ.")
   }
   cat("Welcome to uafR! Preparing your data, please be patient.\n")

   gcms_spread_area = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
   gcms_spread_RT = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
   gcms_spread_prob = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
   gcms_spread_mz = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
   gcms_spread_compound = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
   gcms_spread_sample = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
   gcms_spread_rtBYmass = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))
   gcms_spread_mass = data.frame(matrix(ncol = length(unique(input$File.Name)), nrow = length(input$Component.RT)))

   row.names(gcms_spread_area) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
   colnames(gcms_spread_area) = unique(input$File.Name)

   row.names(gcms_spread_RT) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
   colnames(gcms_spread_RT) = unique(input$File.Name)

   row.names(gcms_spread_prob) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
   colnames(gcms_spread_prob) = unique(input$File.Name)

   row.names(gcms_spread_mz) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
   colnames(gcms_spread_mz) = unique(input$File.Name)

   row.names(gcms_spread_compound) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
   colnames(gcms_spread_compound) = unique(input$File.Name)

   row.names(gcms_spread_sample) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
   colnames(gcms_spread_sample) = unique(input$File.Name)

   row.names(gcms_spread_rtBYmass) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
   colnames(gcms_spread_rtBYmass) = unique(input$File.Name)

   row.names(gcms_spread_mass) = paste0(input$Component.RT, input$Component.Area, input$Base.Peak.MZ)
   colnames(gcms_spread_mass) = unique(input$File.Name)

   tentative_identities = as.character(paste0(input$Compound.Name))
   tentative_RTs = as.numeric(paste0(input$Component.RT))
   num_tentative = length(tentative_identities)

   mz_mass_list = list()

   all_tentative_list = list()
   all_search_list = list()

   tentative_masses = c()
   searchNames_published = list()

   for (chem in 1:num_tentative){
      Sys.sleep(1)
      alt_trigger = F
      current_CMP = tentative_identities[chem]
      current_RT = tentative_RTs[chem]
      chem_cid = .uaf_get_cid(tentative_identities[chem])

      if(is.na(chem_cid)){
         smiles_url = paste0("https://cactus.nci.nih.gov/chemical/structure/",current_CMP,"/smiles")
         inchi_url = paste0("https://cactus.nci.nih.gov/chemical/structure/",current_CMP,"/stdinchikey")
         smiles_url = gsub("\\ ", "%20", smiles_url)
         inchi_url = gsub("\\ ", "%20", inchi_url)
         smile_string = getNCI(smiles_url)
         inchi_string = getNCI(inchi_url)

         if(smile_string != "None"){
            InChiKey = inchi_string
            smile_cid = .uaf_get_cid(paste0(smile_string), from = "smiles")

            chem_cid = smile_cid
         }else{}
      }
      chem_cid = paste0(chem_cid)

      if(chem_cid == "0"){chem_cid = "180"}

      chem_cid = gsub("NA",
                      "180",
                      chem_cid)
      if(chem_cid == "180"){alt_trigger = T}

      chem_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/', chem_cid,'/JSON?heading=GC-MS')

      chem_names_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/', chem_cid,'/JSON?heading=Depositor-Supplied+Synonyms')

      chem_names_info = tryCatch(jsonlite::fromJSON(chem_names_url),
                                 error = function(error) {return("None")})
      if(chem_names_info == "None"){next}

      chem_names_json_dat = data.frame(unlist(chem_names_info$Record$Section))
      colnames(chem_names_json_dat) = "info"

      name_finding = "Section.Section.Information.Value.StringWithMarkup.String"
      all_current_chem_names = data.frame(paste0(chem_names_json_dat$info[adist(rownames(chem_names_json_dat), name_finding) <= 4]))
      colnames(all_current_chem_names) = "names"

      searchNames_published[[chem]] = unique(all_current_chem_names$names)

      chem_info = tryCatch(jsonlite::fromJSON(chem_url),
                           error = function(error) {return("None")})
      if (chem_info == "None" | chem_cid == 180){
         chem_cid = 180

         chem_url = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/', chem_cid, '/JSON?heading=GC-MS')

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
      mz_rows_1 = grep("[[:digit:]]+\\.[[:digit:]]\\ [[:digit:]]+\\.[[:digit:]]", chem_json_dat$info)
      mz_rows_2 = grep("[[:digit:]]{2,3}\\ [[:digit:]]+\\.[[:digit:]]", chem_json_dat$info)
      mz_rows_3 = grep("[[:digit:]]{2,3}\\ [[:digit:]]{2,3}", chem_json_dat$info)
      mz_rows_4 = grep("^[[:digit:]]{2,3}$", chem_json_dat$info)

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
      chem_mz_matches2 = gsub("\\.0\\>","",chem_mz_matches)
      chem_mz_matches2 = gsub("\\.[[:digit:]]0\\>","",chem_mz_matches2)

      mz_count = 1:length(chem_mz_matches2)
      odds = mz_count[mz_count%%2 != 0]
      evens = mz_count[mz_count%%2 == 0]

      mz_primary = chem_mz_matches2[odds]
      mz_secondary = chem_mz_matches2[evens]
      chem_mz_matches2 = mz_primary[order(mz_secondary, decreasing = T)]

      all_current_names = paste0(unique(all_current_chem_names$names))

      if(alt_trigger == T){
         mz_primary = input$Base.Peak.MZ[input$Compound.Name == current_CMP]
         all_current_names = c(current_CMP, CMP_name_tmp)
         mass_tmp = "NA"}

      if(length(chem_mz_matches2) > 5){chem_mz_matches2 = chem_mz_matches2[1:5]}
      mz_print = paste(chem_mz_matches2, collapse = " | ")

      total_chems = length(tentative_identities)

      step_printer = c(1:num_tentative*100)
      step_cmper = chem/25

      if(step_cmper %in% step_printer | chem == 1){
         cat(paste0('\n', '[Current/Total]', ' |--Exact Mass--|--CMP Name 1--|--CMP Name 2--|--Top MZ Peaks--|', '\n'))
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
      cat(paste0('[', chem, ' / ', num_tentative, ']', ': ', all_print, '\n'))

      mz_mass_list[[chem]] = list(mass_tmp, mz_primary)

      tentative_masses = c(tentative_masses, as.numeric(paste0(mass_tmp)))
      all_search_list[[chem]] = list(all_current_names, mz_primary, as.numeric(paste0(mass_tmp)), current_RT)
   }
   if(length(mz_mass_list) < length(tentative_identities)){mz_mass_list[[length(tentative_identities)]] = "NA"}
   names(mz_mass_list) = tentative_identities

   if(length(searchNames_published) < length(tentative_identities)){searchNames_published[[length(tentative_identities)]] = "NA"}
   names(searchNames_published) = tentative_identities

   if(length(all_search_list) < length(tentative_identities)){all_search_list[[length(tentative_identities)]] = "NA"}
   names(all_search_list) = tentative_identities

   sample_names = unique(input$File.Name)
   compound_names = unique(input$Compound.Name)
   compound_area = unique(input$Component.Area)
   compound_RT = unique(input$Component.RT)
   compound_match = unique(input$Match.Factor)
   print("Spreading out your data, one datum at a time.")

   for (j in 1:nrow(gcms_spread_area)){
      area_tmp = paste0(input$Component.Area[j])
      compound_tmp = paste0(input$Compound.Name[j])
      sample_tmp = paste0(input$File.Name[j])
      mz_tmp = paste0(input$Base.Peak.MZ[j])
      prob_tmp = paste0(input$Match.Factor[j])
      RT_tmp = paste0(input$Component.RT[j])
      mass_tmp = paste0(tentative_masses[j])
      rtBYmass_tmp = paste0(RT_tmp,' | ', mass_tmp)

      gcms_spread_area[rownames(gcms_spread_area) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_area) == paste0(sample_tmp)] = area_tmp
      gcms_spread_compound[rownames(gcms_spread_compound) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_compound) == paste0(sample_tmp)] = compound_tmp
      gcms_spread_sample[rownames(gcms_spread_sample) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_sample) == paste0(sample_tmp)] = sample_tmp
      gcms_spread_mz[rownames(gcms_spread_mz) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_mz) == paste0(sample_tmp)] = mz_tmp
      gcms_spread_prob[rownames(gcms_spread_prob) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_prob) == paste0(sample_tmp)] = prob_tmp
      gcms_spread_RT[rownames(gcms_spread_RT) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_RT) == paste0(sample_tmp)] = RT_tmp
      gcms_spread_rtBYmass[rownames(gcms_spread_rtBYmass) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_rtBYmass) == paste0(sample_tmp)] = rtBYmass_tmp
      gcms_spread_mass[rownames(gcms_spread_mass) == paste0(RT_tmp, area_tmp, mz_tmp), colnames(gcms_spread_mass) == paste0(sample_tmp)] = mass_tmp
   }

   print("Losing rownames (data ID codes).")
   row.names(gcms_spread_area) = NULL
   row.names(gcms_spread_RT) = NULL
   row.names(gcms_spread_prob) = NULL
   row.names(gcms_spread_mz) = NULL
   row.names(gcms_spread_compound) = NULL
   row.names(gcms_spread_sample) = NULL
   row.names(gcms_spread_mass) = NULL
   row.names(gcms_spread_rtBYmass) = NULL

   gcms_list = list(gcms_spread_area, gcms_spread_compound,
                    gcms_spread_mz, gcms_spread_prob,
                    gcms_spread_RT, gcms_spread_mass,
                    gcms_spread_rtBYmass, all_search_list)
   names(gcms_list) = c("Area", "Compounds",
                        "MZ", "MatchFactor",
                        "RT", "Mass",
                        "rtBYmass", "webInfo")
   return(gcms_list)
}
