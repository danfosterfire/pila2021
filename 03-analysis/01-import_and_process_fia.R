
library(here)
library(tidyverse)
library(DBI)
library(RSQLite)


#### read in data from SQLite databases ########################################

# sqlite dbs downloaded from 
# https://apps.fs.usda.gov/fia/datamart/datamart_sqlite.html
# on 12/14/2021

fiadb.ca = 
  dbConnect(RSQLite::SQLite(),
            here::here('02-data',
                       '00-source',
                       'fia',
                       'FIADB_CA.db'))

fiadb.nv = 
  dbConnect(RSQLite::SQLite(),
            here::here('02-data',
                       '00-source',
                       'fia',
                       'FIADB_NV.db'))
  
fiadb.or = 
  dbConnect(RSQLite::SQLite(),
            here::here('02-data',
                       '00-source',
                       'fia',
                       'FIADB_OR.db'))

fiadb.wa = 
  dbConnect(RSQLite::SQLite(),
            here::here('02-data',
                       '00-source',
                       'fia',
                       'FIADB_WA.db'))

# they're not all the same which is super annoying
dbListTables(fiadb.ca) == dbListTables(fiadb.nv) &
  dbListTables(fiadb.ca) == dbListTables(fiadb.or) &
  dbListTables(fiadb.ca) == dbListTables(fiadb.wa)

dbListTables(fiadb.ca)

dbReadTable(fiadb.ca, 'REF_FIADB_VERSION')

fia = 
  lapply(X =  c('COND', 'PLOT', 'REF_FOREST_TYPE',
                'REF_FOREST_TYPE_GROUP', 'REF_SPECIES', 'SEEDLING', 'SUBPLOT', 
                'SUBP_COND', 'TREE'),
         FUN = function(tname){
           dbReadTable(fiadb.ca, tname) %>%
             bind_rows(dbReadTable(fiadb.nv, tname)) %>%
             bind_rows(dbReadTable(fiadb.or, tname)) %>%
             bind_rows(dbReadTable(fiadb.wa, tname)) %>%
             as_tibble()
         })

names(fia) = c('COND', 'PLOT', 'REF_FOREST_TYPE',
               'REF_FOREST_TYPE_GROUP', 'REF_SPECIES', 'SEEDLING', 'SUBPLOT', 
               'SUBP_COND', 'TREE')
         
dbDisconnect(fiadb.ca)
dbDisconnect(fiadb.or)
dbDisconnect(fiadb.nv)
dbDisconnect(fiadb.wa)

# on this, ans presumably many other plt_cns, the plot code is 
# inconsistent between the PLOT table and the TREE table in the FIAPNW
# database. Is it also broken here? No.
bad_cn = "290008535489998"
fia$PLOT %>% filter(CN==bad_cn)
fia$TREE %>% filter(PLT_CN==bad_cn)

#### create crosswalks for codes ###############################################


# copying the relevant crosswalks from the FIA documentation
cond_dstrbcds = 
  data.frame(
    code = 
      c(0, 10, 11, 12, 20, 21, 22, 30, 31, 32, 40, 41, 42, 43, 44,
        45, 46, 50, 51, 52, 53, 54, 60, 70, 80, 90, 91, 92, 93, 94,
        95),
    cond_dstrbdesc = 
      c('no visible disturbance',
        'insect damage', 'insect damage to understory vegetation', 
        'insect damage to rees, including seedlings and saplings',
        'disease damage', 'disease damage to understory vegetation',
        'disease damage to trees, including seedlings and saplings',
        'fire damage (from crown and ground fire, either rx or natural',
        'ground fire damage', 'crown fire damage',
        'animal damage', 'beaver (includes flooding caused by beaver',
        'porcupine', 'deer-ungulate', 'bear', 'rabbit', 
        'domestic animal / livestock',
        'weather damage', 'ice', 'wind (includes hurricane, tornado)', 
        'flooding (weather induced)', 'drought', 
        'vegetation (suppression, competition, vines)',
        'unknown / not sure / other',
        'human-induced damage - any significant level of human-caused damage not described in the disturbance codes or the treatment codes',
        'geologic disturbances', 'landslide', 'avalanche track', 'volcanic blast zone', 
        'other geologic event', 'earth movement / avalanches'),
    cond_dstrbtype = 
      c('none', 'insects', 'insects (understory)', 'insects (trees)',
        'disease', 'disease (understory)', 'disease (trees)', 
        'fire', 'fire (ground)', 'fire (crown)', 
        'animal', 'animal', 'animal', 'animal', 'animal', 'animal', 'animal',
        'weather (other)', 'weather (other)', 'weather (other)', 'weather (other)',
        'drought', 'vegetation', 'unknown / other', 'human (other)', 
        'geologic', 'geologic', 'geologic', 'geologic', 'geologic', 'geologic'))  %>%
  mutate(cond_dstrbdesc = as.character(cond_dstrbdesc),
         cond_dstrbtype = as.character(cond_dstrbtype)) %>%
  as_tibble()

cond_trtcds = 
  data.frame(
    code = 
      c('0', '10', '20', '30', '40', '50'),
    cond_trtdesc = 
      c('no observable treatment',
        'cutting - removal of one or more trees',
        'site prep - clearing, slash burning, chopping, disking, bedding, or other practices clearly intended to prepare a site for either natural or artificial regeneration',
        'artificial regeneration - following a disturbance or treatment (usually cutting), a new stand where at least 50 percent of the live trees present resulted from planting or direct seeding',
        'natural regeneration - following a disturbance or treatment (usually cutting), a new stand where at least 50 percent of the live trees present (of any size) were established through the growth of existing trees and/or natural regen or sprouting',
        'other silvicultural treatment - the use of fertilizers, herbicides, girdling, purning, or other activities (not covered by codes 10-40) designed to improve the commercial value of the residual stand; or chaining, which is a practice used on woodlands to encourage wildlife forage'),
    cond_trttype = 
      c('none',
        'cutting',
        'siteprep',
        'regen (artificial)',
        'regen (natural)',
        'other')
  ) %>%
  mutate(code = as.integer(as.character(code)),
         cond_trtdesc = as.character(cond_trtdesc),
         cond_trttype = as.character(cond_trttype)) %>%
  as_tibble()

cond_statuscds = 
  data.frame(
    code = c(1, 2, 3, 4, 5),
    cond_statustype = 
      c('accessible_forest',
        'nonforest',
        'noncensus_water',
        'census_water',
        'nonsampled_potentialforest')
  ) %>%
  mutate(cond_statustype = as.character(cond_statustype)) %>%
  as_tibble()

plot_kindcds = 
  data.frame(
    code = c(0, 1, 2, 3, 4),
    plot_kindtype = 
      c('periodic', 'national_initial', 'national_remeasure', 'national_replacement',
        'modeled_periodic')
  ) %>%
  mutate(across(plot_kindtype, as.character)) %>%
  as_tibble()

plot_statuscds = 
  data.frame(
    code = c(1, 2, 3),
    plot_statustype = 
      c('sampled_accessible_forest',
        'sampled_no_forest',
        'nonsampled')
  ) %>%
  mutate(across(plot_statustype, as.character)) %>%
  as_tibble()

plot_designcds = 
  data.frame(
    code = c(1, 501, 502), 
    plot_designtype = 
      c('standard', 'macro_24', 'macro_30')
  ) %>%
  mutate(plot_designtype = as.character(plot_designtype))

subp_statuscds = plot_statuscds

tree_statuscds = 
  data.frame(
    code = c('0', '1', '2', '3'),
    tree_statustype = 
      c('outofsample', 'live', 'dead', 'harvested')
  ) %>%
  mutate(code = as.integer(as.character(code)),
         tree_statustype = as.character(tree_statustype)) %>%
  as_tibble()

tree_cclasscds = 
  data.frame(
    code = c('1', '2', '3', '4', '5'),
    cclctype = 
      c('opengrown', 'dominant', 'codominant', 'intermediate', 'overtopped')
  ) %>%
  mutate(code = as.integer(as.character(code)),
         cclctype = as.character(cclctype)) %>%
  as_tibble()


tree_reconcilecds = 
  data.frame(
    code = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    tree_reconciletype = 
      c('ingrowth (other)',
        'ingrowth (growth)',
        'missed live',
        'missed dead',
        'shrank',
        'missing (moved)',
        'previous error',
        'procedural change',
        'area nonsampled')
  ) %>%
  mutate(across(tree_reconciletype, as.character)) %>%
  as_tibble()


tree_snagdiscds = 
  data.frame(
    code = 
      c(2, 3, 4, 5, 6),
    tree_snagdistype = 
      c('fell_present',
        'fell_missing',
        'cut_present',
        'cut_missing',
        'shrank')
  ) %>%
  mutate(across(tree_snagdistype, as.character)) %>%
  as_tibble()

cond_physclcds = 
  data.frame(
    code = 
      c(11, 12, 13, 19, 
        21, 22, 23, 24, 25, 29,
        31, 32, 33, 34, 35, 39),
    physcldesc = 
      c('dry_tops', 'dry_slopes', 'deep_sands', 'other_xeric',
        'flatwoods', 'rolling_uplands', 'moist_slopes', 'narrow_bottomlands',
        'broad_bottomlands', 'other_mesic',
        'swampsbogs', 'small_drains', 'bays_pocosins', 'beaver_ponds',
        'cypress_ponds', 'other_hydric'),
    physcltype = 
      c(rep('xeric', 4), rep('mesic', 6), rep('hydric', 6))
  ) %>%
  mutate(across(c(physcldesc, physcltype), as.character)) %>%
  as_tibble()

# ecoregion subsection codes; disabled for now
# join the descriptive names with the codes
#plot_ecosubcds = 
#  st_read(here::here('02-data','00-source',
#                     'usda', 'ecoregions',
#                     'S_USA.EcoMapProvinces.shp')) %>%
#  as.data.frame() %>%
#  dplyr::select(province_code = MAP_UNIT_S,
#                province_name = MAP_UNIT_N) %>%
#  left_join(st_read(here::here('02-data','00-source',
#                               'usda', 'ecoregions',
#                               'S_USA.EcomapSections.shp')) %>%
#              as.data.frame() %>%
#              select(section_code = MAP_UNIT_S,
#                     section_name = MAP_UNIT_N) %>%
#              mutate(province_code = stringr::str_sub(section_code, end = -2)),
#            by = 'province_code') %>%
#  left_join(st_read(here::here('02-data','00-source',
#                                 'usda', 'ecoregions',
#                                 'S_USA.EcomapSubsections.shp')) %>%
#              as.data.frame() %>%
#              select(subsection_code = MAP_UNIT_S,
#                     subsection_name = MAP_UNIT_N) %>%
#              mutate(section_code = stringr::str_sub(subsection_code, end = -2))) %>%
#  filter(!is.na(section_code)) %>%
#  mutate_all(as.character)

# join the descriptive names to the FIA plot data
#fia_plots = 
#  fia_plots %>%
# get some of the more-coarse ecoregions
#  mutate(ecosub = stringr::str_remove(ecosub, ' '),
#         province = stringr::str_sub(ecosub, end = -3),
#         section = stringr::str_sub(ecosub, end = -2),
#         subsection = ecosub) %>%
#  select(-ecosub) %>%
#  left_join(ecoregion_defs,
#            by = c('province' = 'province_code',
#                   'section' = 'section_code',
#                   'subsection' = 'subsection_code')) %>%
#  as_tibble()

plot_topopositioncds = 
  data.frame(
    code = 
      c(1:9),
    plot_topopositiontype = 
      c('ridgetop', 'ridgetop_narrow', 'sidehill_upper', 'sidehill_middle',
        'sidehill_lower', 'bottom_narrow', 'bench', 'bottom_alluvialflat', 
        'swamp_wetflat')
  ) %>%
  mutate(across(everything(), as.character)) %>%
  as_tibble()

plot_nonsamplereasncds = 
  data.frame(
    code = c(2, 3, 8),
    plot_nonsamplereasntype = 
      c('denied_access', 'hazard', 'skipped')
  ) %>%
  mutate(plot_nonsamplereasntype = as.character(plot_nonsamplereasntype)) %>%
  as_tibble()

#### prepare subplots data #####################################################

# one row per observation of a subplot, with columns indicating the disturbances 
# present at that observation, the spatial location, etc.
# the goal is to have a table with one row per observation of a condition (stand 
# on a subplot at a specific inventory date)
# context
subplots = 
  
  # start with the plot table, with one row per plot observation
  fia$PLOT %>%
  as_tibble() %>%
  rename(PLT_CN = CN) %>%
  select(PLT_CN, PREV_PLT_CN, INVYR, STATECD, UNITCD, COUNTYCD, PLOT, 
         PLOT_STATUS_CD, PLOT_NONSAMPLE_REASN_CD, MEASYEAR, 
         MEASMON, MEASDAY, KINDCD, DESIGNCD, MANUAL, MACRO_BREAKPOINT_DIA,
         LAT, LON, ELEV, ECOSUBCD) %>%
  
  # join in the subplot info, with one row per subplot observation
  left_join(fia$SUBPLOT %>%
              rename(SUBP_CN = CN) %>%
              select(SUBP_CN, PLT_CN, 
                     STATECD, UNITCD, COUNTYCD, PLOT, SUBP),
            by = c('PLT_CN' = 'PLT_CN',
                   'STATECD' = 'STATECD',
                   'UNITCD' = 'UNITCD',
                   'COUNTYCD' = 'COUNTYCD',
                   'PLOT' = 'PLOT')) %>%
  
  # join the subplot conditions; this creates one row per condition:subplot:plot:time
  left_join(fia$SUBP_COND %>%
              rename(SUBPCOND_CN = CN) %>%
              select(PLT_CN, SUBPCOND_CN, STATECD, UNITCD,
                     COUNTYCD, PLOT, SUBP, CONDID),
            by = c('PLT_CN' = 'PLT_CN', 
                   'STATECD' = 'STATECD',
                   'UNITCD' = 'UNITCD',
                   'COUNTYCD' = 'COUNTYCD',
                   'PLOT' = 'PLOT',
                   'SUBP' = 'SUBP')) %>%
  
  # join the condition details
  left_join(fia$COND %>%
              rename(CND_CN = CN) %>%
              select(PLT_CN, CND_CN, STATECD, UNITCD, 
                     COUNTYCD, PLOT, CONDID, 
                     #FORTYPCD,
                     DSTRBCD1, DSTRBCD2, DSTRBCD3,
                     TRTCD1, TRTCD2, TRTCD3),
            by = c('PLT_CN' = 'PLT_CN',
                   'STATECD' = 'STATECD',
                   'UNITCD' = 'UNITCD',
                   'COUNTYCD' = 'COUNTYCD',
                   'PLOT' = 'PLOT',
                   'CONDID' = 'CONDID')) %>%
  
  # join the code meanings
  left_join(cond_dstrbcds %>%
              select(DSTRBCD1 = code,
                     DSTRBCD1.DESC = cond_dstrbdesc,
                     DSTRBCD1.TYPE = cond_dstrbtype)) %>%
  
  left_join(cond_dstrbcds %>%
              select(DSTRBCD2 = code,
                     DSTRBCD2.DESC = cond_dstrbdesc,
                     DSTRBCD2.TYPE = cond_dstrbtype)) %>%
  
  left_join(cond_dstrbcds %>%
              select(DSTRBCD3 = code,
                     DSTRBCD3.DESC = cond_dstrbdesc,
                     DSTRBCD3.TYPE = cond_dstrbtype)) %>%
  
  left_join(cond_trtcds %>%
              select(TRTCD1 = code,
                     TRTCD1.DESC = cond_trtdesc,
                     TRTCD1.TYPE = cond_trttype)) %>%
  
  left_join(cond_trtcds %>%
              select(TRTCD2 = code,
                     TRTCD2.DESC = cond_trtdesc,
                     TRTCD2.TYPE = cond_trttype)) %>%
  
  left_join(cond_trtcds %>%
              select(TRTCD3 = code,
                     TRTCD3.DESC = cond_trtdesc,
                     TRTCD3.TYPE = cond_trttype)) %>%
  
  left_join(plot_designcds %>%
              select(DESIGNCD = code,
                     DESIGNCD.TYPE = plot_designtype)) %>%
  
  left_join(plot_kindcds %>%
              select(KINDCD = code,
                     KINDCD.TYPE = plot_kindtype)) %>%
  
  left_join(plot_statuscds %>%
              select(PLOT_STATUS_CD = code,
                     PLOT_STATUS_CD.TYPE = plot_statustype)) %>%
  
  left_join(plot_nonsamplereasncds %>%
              select(PLOT_NONSAMPLE_REASN_CD = code,
                     PLOT_NONSAMPLE_REASN_CD.TYPE = plot_nonsamplereasntype))# %>%
  
  #left_join(fia$REF_FOREST_TYPE %>%
  #            select(FORTYPCD = VALUE,
  #                   FORTYPCD.TYPE = MEANING))



plots =
  subplots %>%
  
  # convert the disturbance and treatment codes to presence/absence
  # treat NA disturbance and treatment codes as 'none'
  mutate(
    disturb1 = ifelse(is.na(DSTRBCD1.TYPE), 'none', DSTRBCD1.TYPE),
    disturb2 = ifelse(is.na(DSTRBCD2.TYPE), 'none', DSTRBCD2.TYPE),
    disturb3 = ifelse(is.na(DSTRBCD3), 'none', DSTRBCD3.TYPE),
    treat1 = ifelse(is.na(TRTCD1), 'none', TRTCD1),
    treat2 = ifelse(is.na(TRTCD2), 'none', TRTCD2),
    treat3 = ifelse(is.na(TRTCD3), 'none', TRTCD3)
  ) %>%
  
  # scan the disturbance columns and turn them into relevant logicals
  mutate(
    
    # first, scan disturb1, 2, and 3 separately for fire, and fill in fire year
    # disturb_yr is well covered for fire, not so much for others
    fire_1 = is.element(disturb1, c('fire', 'fire (crown)', 'fire (ground)')),
    crownfire_1 = disturb1=='fire (crown)',
    fire_2 = is.element(disturb2, c('fire', 'fire (crown)', 'fire (ground)')),
    crownfire_2 = disturb2=='fire (crown)',
    fire_3 = is.element(disturb3, c('fire', 'fire (crown)', 'fire (ground)')),
    crownfire_3 = disturb3=='fire (crown)',
    
    fire = fire_1|fire_2|fire_3,
    
    insects = 
      is.element(disturb1, c('insects (trees)', 'insects')) | 
      is.element(disturb2, c('insects (trees)', 'insects')) | 
      is.element(disturb3, c('insects (trees)', 'insects')),
    
    disease = 
      is.element(disturb1, c('disease (trees)', 'disease')) | 
      is.element(disturb2, c('disease (trees)', 'disease')) | 
      is.element(disturb3, c('disease (trees)', 'disease')),
    
    drought = 
      disturb1=='drought'|disturb2=='drought'|disturb3=='drought',
    
    suppression = 
      disturb1=='vegetation'|disturb2=='vegetation'|disturb3=='vegetation',
    
    other = 
      (is.element(disturb1, 
                  c('weather (other)', 'human (other)', 'unknown / other',
                    'animal', 'geologic', 'insects (understory)', 
                    'disease (understory)')))| 
      (is.element(disturb2, 
                  c('weather (other)', 'human (other)', 'unknown / other',
                    'animal', 'geologic', 'insects (understory)', 
                    'disease (understory)')))|
      (is.element(disturb3, 
                  c('weather (other)', 'human (other)', 'unknown / other',
                    'animal', 'geologic', 'insects (understory)', 
                    'disease (understory)'))),
    
    cutting = 
      treat1=='cutting'|treat2=='cutting'|treat3=='cutting',
    
    siteprep = 
      treat1=='siteprep'|treat2=='siteprep'|treat3=='siteprep',
    
    artificialregen = 
      treat1=='regen (artificial)'|treat2=='regen (artificial)'|treat3=='regen (artificial)',
    
    naturalregen = 
      treat1=='regen (natural)'|treat2=='regen (natural)'|treat3=='regen (natural)'
    ) %>%
  
  # aggregate the conds together to subplots
  group_by(PLT_CN, PREV_PLT_CN, INVYR,
           STATECD, UNITCD, COUNTYCD, PLOT,
           MEASYEAR, MEASMON, MEASDAY,
           KINDCD.TYPE, DESIGNCD.TYPE, MANUAL,
           PLOT_STATUS_CD.TYPE, PLOT_NONSAMPLE_REASN_CD.TYPE,
           ELEV, LAT, LON, ECOSUBCD, MACRO_BREAKPOINT_DIA) %>%
  summarise(
    fire = any(fire),
    insects = any(insects),
    disease = any(disease),
    drought = any(drought),
    suppression = any(suppression),
    cutting = any(cutting),
    siteprep = any(siteprep),
    artificialregen = any(artificialregen),
    naturalregen = any(naturalregen)
    
  ) %>%
  ungroup() %>%
  
  mutate(
    plot_id = paste(STATECD, UNITCD, COUNTYCD, PLOT, 
                    sep = '-'),
    invdate = as.Date(paste(MEASYEAR, MEASMON, MEASDAY, sep = '-'))
  ) %>%
  
  select(
    plt_cn = PLT_CN,
    prev_plt_cn = PREV_PLT_CN,
    invyr = INVYR,
    plot_id, 
    invdate,
    inv_kind = KINDCD.TYPE,
    inv_design = DESIGNCD.TYPE,
    inv_manual = MANUAL,
    plot_status = PLOT_STATUS_CD.TYPE,
    plot_nonsamp = PLOT_NONSAMPLE_REASN_CD.TYPE,
    macro_break = MACRO_BREAKPOINT_DIA,
    elev_ft = ELEV, 
    lat = LAT,
    lon = LON,
    ecosubcd = ECOSUBCD,
    fire, insects, disease, drought, suppression, cutting, siteprep, 
    artificialregen, naturalregen
  )
  

#### prepare treelist data #####################################################

treelist = 
  
  fia$TREE %>%
  as_tibble() %>%
  rename(TRE_CN = CN) %>%
  left_join(tree_cclasscds %>%
              select(CCLCD = code,
                     CCLCD.TYPE = cclctype)) %>%
  
  left_join(tree_reconcilecds %>%
              select(RECONCILECD = code,
                     RECONCILECD.TYPE = tree_reconciletype)) %>%
  left_join(tree_snagdiscds %>%
              select(SNAG_DIS_CD_PNWRS = code,
                     SNAG_DIS_CD_PNWRS.TYPE = tree_snagdistype)) %>%
  left_join(tree_statuscds %>%
              select(STATUSCD = code,
                     STATUSCD.TYPE = tree_statustype)) %>%
  # join in the species labels
  left_join(fia$REF_SPECIES %>%
              select(SPCD,
                     SPECIES = SPECIES_SYMBOL) %>%
              group_by(SPCD, SPECIES) %>%
              summarise() %>%
              ungroup()) %>%
  
  mutate(plot_id = 
           paste(STATECD, UNITCD, COUNTYCD, PLOT,
                 sep = '-'),
         subp_id = 
           paste(STATECD, UNITCD, COUNTYCD, PLOT, SUBP,
                 sep = '-'),
         tree_id = 
           paste(STATECD, UNITCD, COUNTYCD, PLOT, SUBP, TREE,
                 sep = '-'),
         wpbr = 
           (!is.na(DMG_AGENT1_CD_PNWRS) & DMG_AGENT1_CD_PNWRS==36)|
           (!is.na(DMG_AGENT2_CD_PNWRS) & DMG_AGENT2_CD_PNWRS==36)|
           (!is.na(DMG_AGENT3_CD_PNWRS) & DMG_AGENT3_CD_PNWRS==36)|
           (!is.na(DAMAGE_AGENT_CD1) & DAMAGE_AGENT_CD1 == 26001)|
           (!is.na(DAMAGE_AGENT_CD2) & DAMAGE_AGENT_CD2 == 26001)|
           (!is.na(DAMAGE_AGENT_CD3) & DAMAGE_AGENT_CD3 == 26001)) %>%
  
  select(tre_cn = TRE_CN,
         plt_cn = PLT_CN,
         prev_tre_cn = PREV_TRE_CN,
         plot_id, subp_id, tree_id,
         #cond_id = CONDID,
         tree_status = STATUSCD.TYPE,
         species = SPECIES,
         dbh_in = DIA,
         tpa_unadj = TPA_UNADJ,
         cclass = CCLCD.TYPE,
         reconcile = RECONCILECD.TYPE,
         wpbr) %>%
  
  mutate(species = ifelse(is.element(species,
                                     c('ABCO', 'CADE27', 'PILA', 'PIPO',
                                       'PSME', 'QUKE')),
                          species,
                          'OTHER'))


         


#### filter plots ##############################################################


# aspatial filtering first
plots = 
  plots %>%
  
  # within CA or OR
  filter(is.element(gsub(x = plot_id, pattern = '-.*$', replacement = ''), 
                    c('6', '41'))) %>%
  
  # invyr in correct range; drops the phase 3 "off subpanel"
  filter(invyr >= 2000 & invyr <= 2021) %>%
  
  # only sampled plots
  filter(plot_status != 'nonsampled') %>%
  
  # plots where PILA is present
  left_join(
    fia$PLOT %>%
      left_join(fia$TREE %>%
                  select(PLT_CN, STATUSCD, SPCD),
                by = c('CN' = 'PLT_CN')) %>%
      mutate(live_pila = STATUSCD==1 & SPCD==117,
             plot_id = paste(STATECD, UNITCD, COUNTYCD, PLOT,
                             sep = '-')) %>%
      select(plot_id, live_pila) %>%
      group_by(plot_id) %>%
      summarise(live_pila = any(live_pila, na.rm = TRUE)) %>%
      ungroup()
  ) %>%
  filter(live_pila)


#### filter treelist ###########################################################

# want only trees from included plots, and no trees which were 
# excluded on remeasurement, and no trees with NA remeasure DBH and remasure
# status=='live' (theres 17 of the latter)
treelist = 
  treelist %>%
  filter(is.element(plt_cn, plots$plt_cn)) %>%
  filter(tree_status != 'outofsample' & 
           !is.element(tre_cn, 
                       treelist %>%
                         filter(tree_status == 'outofsample') %>%
                         pull(prev_tre_cn))) %>%
  filter(!(tree_status=='live'&
             !is.na(prev_tre_cn)&
             is.na(dbh_in)))


#### build subplot covariates data frame #######################################


# all three models - growth, mortality, and recruitment - use the 
# same plot-level explanatory variables:
# presence of fire/insects/disease/cutting flag at remeasurement
# basal area at initial measurement
# max CWD between initial and remeasurement
# for a plot-level analysis, also want TPH and BA of live PILA at initial
# and remeasurement

# need a dataframe with one row per plot
plot_data = 
  
  plots %>%
  
  # just the remeasurement plots which have an initial measure associated
  filter(!is.na(prev_plt_cn) & inv_kind == 'national_remeasure') %>%
  select(plt_cn, prev_plt_cn, plot_id, 
         elev_ft, lat, lon, ecosubcd,
         invdate, inv_manual, macro_break,
         fire, insects, disease, cutting) %>%
  
  # get the invdate and inv manual for the initial measurement
  left_join(plots %>%
              select(plt_cn, plot_id, invdate, inv_manual),
            suffix = c('.re', '.init'),
            by = c('prev_plt_cn' = 'plt_cn',
                   'plot_id' = 'plot_id')) %>%
  
  # get the live BA for the initial measurement
  left_join(treelist %>%
              filter(tree_status == 'live') %>%
              mutate(ba_ft2ac = (pi*((0.5*(dbh_in/12))^2))*tpa_unadj) %>%
              group_by(plt_cn, plot_id) %>%
              summarise(ba_ft2ac = sum(ba_ft2ac, na.rm = TRUE)) %>%
              ungroup(),
            by = c('prev_plt_cn' = 'plt_cn',
                   'plot_id' = 'plot_id')) %>%
  mutate(ba_ft2ac = 
           ifelse(is.na(ba_ft2ac), 0, ba_ft2ac)) %>%
  
  # get the presence or absence of WPBR for the initial measurement
  left_join(treelist %>%
              group_by(plt_cn, plot_id) %>%
              summarise(wpbr = any(wpbr)) %>%
              ungroup(),
            by = c('prev_plt_cn' = 'plt_cn',
                   'plot_id' = 'plot_id')) %>%
  mutate(wpbr = 
           ifelse(is.na(wpbr), FALSE, wpbr)) %>%
  
  # get the live TPH and BA of pila at initial measurement
  left_join(treelist %>%
              filter(tree_status=='live'& species=='PILA') %>%
              mutate(pila_ba_m2ha = (pi*((0.5*(dbh_in*2.54)/100)**2)*(tpa_unadj/0.404686)),
                     pila_tph = tpa_unadj/0.404686) %>%
              group_by(plt_cn, plot_id) %>%
              summarise(pila_ba_m2ha = sum(pila_ba_m2ha, na.rm = TRUE),
                        pila_tph = sum(pila_tph, na.rm = TRUE)) %>%
              ungroup(),
            by = c('prev_plt_cn' = 'plt_cn',
                   'plot_id' = 'plot_id')) %>%
  # get the live TPH and BA of pila at remeasurement
  left_join(treelist %>%
              filter(tree_status=='live'& species=='PILA') %>%
              mutate(pila_ba_m2ha = (pi*((0.5*(dbh_in*2.54)/100)**2)*(tpa_unadj/0.404686)),
                     pila_tph = tpa_unadj/0.404686) %>%
              group_by(plt_cn, plot_id) %>%
              summarise(pila_ba_m2ha = sum(pila_ba_m2ha, na.rm = TRUE),
                        pila_tph = sum(pila_tph, na.rm = TRUE)) %>%
              ungroup(),
            by = c('plt_cn' = 'plt_cn',
                   'plot_id' = 'plot_id'),
            suffix = c('.init', '.re'))  %>%
  # fill in missing TPH and BA with 0
  mutate(pila_ba_m2ha.init = ifelse(is.na(pila_ba_m2ha.init), 0, pila_ba_m2ha.init),
         pila_ba_m2ha.re = ifelse(is.na(pila_ba_m2ha.re), 0, pila_ba_m2ha.re),
         pila_tph.init = ifelse(is.na(pila_tph.init), 0, pila_tph.init),
         pila_tph.re = ifelse(is.na(pila_tph.re), 0, pila_tph.re))
  
  

#### get timestep distribution #################################################

head(plot_data)

invdate_diffs = 
  as.numeric(plot_data$invdate.re-plot_data$invdate.init) / 365

summary(invdate_diffs)

quantile(invdate_diffs, c(0, .05, 0.5, 0.95, 1))

#### pull in CWD data ##########################################################

library(terra)

head(plot_data)

plot_data.sf = 
  plot_data %>%
  sf::st_as_sf(coords = c('lon', 'lat'),
           crs = 'EPSG:4269')


plots_bbox = 
  list('lat_min' = min(plot_data$lat)-1,
       'lat_max' = max(plot_data$lat)+1,
       'lon_min' = min(plot_data$lon)-1,
       'lon_max' = max(plot_data$lon)+1)

plots_bbox

cwd_growseason_means = 
  lapply(X = 2000:2020,
         FUN = function(y){
           
           cwd_year = 
             rast(here::here('02-data',
                             '00-source',
                             'terraclimate',
                             paste0('TerraClimate_def_',y,'.nc')))
           
           cwd_year = 
             crop(cwd_year,
                  c(plots_bbox$lon_min, plots_bbox$lon_max,
                    plots_bbox$lat_min, plots_bbox$lat_max))
           
           cwd_year = 
             mean(cwd_year[[5:10]])
           
           return(cwd_year)
           
         }) %>%
  rast()

names(cwd_growseason_means) = 
  paste0('growseasonmean_',as.character(2000:2020))

plot(cwd_growseason_means)

# get departure from "normal" (20 year mean) CWD for each 
# year on each location
cwd_departure = 
  cwd_growseason_means - mean(cwd_growseason_means)

names(cwd_departure) = paste0('cwddeparture_', as.character(2000:2020))

plot(cwd_departure)

head(plot_data)

cwd_departures = 
  extract(cwd_departure, plot_data[,c('lon', 'lat')]) %>%
  bind_cols('plot_id' = plot_data$plot_id,
            'year_begin' = lubridate::year(plot_data$invdate.init),
            'year_end' = lubridate::year(plot_data$invdate.re))

cwd_departure_span = 
  sapply(X = 1:nrow(cwd_departures),
         FUN = function(i){
           
           year_range = seq(from = cwd_departures$year_begin[i],
                            to = cwd_departures$year_end[i],
                            by = 1)
           
           columns_to_select = 
             paste0('cwddeparture_', year_range)
           
           values_in_span = 
             as.numeric(cwd_departures[i,columns_to_select])
           
           # return the 90th percentile CWD departure
           return(as.numeric(quantile(values_in_span, probs = 0.9)))
           
         })

head(cwd_departure_span)

plot_data$cwd_departure90 = cwd_departure_span

plot_data$cwd_mean = 
  extract(mean(cwd_growseason_means),
          plot_data[,c('lon', 'lat')])$mean


#### make individual mortality data frame ######################################

mort_data = 
  
  # start with the treelist
  treelist %>%
  
  select(tre_cn, prev_tre_cn, plt_cn, plot_id, tree_status, species, dbh_in) %>%
  
  # filter to only trees which were recorded at remeasurement
  filter(!is.na(prev_tre_cn)) %>%
  
  # join in the initial size and status
  left_join(.,
            treelist %>%
              select(tre_cn, plt_cn, tree_status, dbh_in),
            by = c('prev_tre_cn' = 'tre_cn'),
            suffix = c('.re', '.init')) %>%
  
  # filter to only trees which were alive at initial measurement
  filter(tree_status.init == 'live') %>%
  
  # create a column for survival 
  mutate(survived = ifelse(tree_status.init=='live'&tree_status.re=='live',
                           TRUE,
                           FALSE))


mort_data %>%
  filter(species=='PILA'&dbh_in.init>70) %>%
  print(width = Inf)

#### make individual growth data frame #########################################

growth_data = 
  
  # start with the treelist
  treelist %>%
  
  select(tre_cn, prev_tre_cn, plt_cn, plot_id, tree_status, species, dbh_in) %>%
  
  # filter to only trees which were recorded at remeasurement
  filter(!is.na(prev_tre_cn)) %>%
  
  # join in the initial size and status
  left_join(.,
            treelist %>%
              select(tre_cn, plt_cn, tree_status, dbh_in),
            by = c('prev_tre_cn' = 'tre_cn'),
            suffix = c('.re', '.init')) %>%
  
  # filter to only trees which were alive at initial measurement and the 
  # remeasurement
  filter(tree_status.init == 'live' & tree_status.re=='live') 


#### recruits as >= 1" #########################################################

# want a size distribution for new recruits, and count of new recruits on 
# each plot

new_recruits = 
  # start with the treelist
  treelist %>%
  
  # filter to only included remeasurement surveys, and only trees which werent 
  # present at initial measurement
  filter(is.element(plt_cn, plot_data$plt_cn) & 
           is.na(prev_tre_cn)) %>%
  
  # pull in the invdates to check
  left_join(plot_data %>%
              select(plt_cn, prev_plt_cn, plot_id, invdate.re, invdate.init)) %>%
  
  # only trees < 5"; assume others are not true new recruits
  filter(dbh_in < 5)

summary(new_recruits) # inv dates look good

hist(new_recruits$dbh_in)

new_recruits %>% filter(species=='PILA') %>% nrow()

# size distribution of new recruits
r = 
  new_recruits %>%
  filter(species=='PILA') %>%
  mutate(dbh_class = cut(dbh_in, 
                         breaks = c(1, 2, 3, 4, 5), 
                         include.lowest = TRUE,
                         labels = FALSE)) %>%
  group_by(dbh_class) %>%
  summarise(tpa_unadj = sum(tpa_unadj)) %>%
  mutate(tpa_total = sum(pull(., tpa_unadj)),
         prop = tpa_unadj / tpa_total) %>%
  pull(prop)

# count of new recruits on each plot
# start with a frame of all the plots and size classes
recruits_data = 
  plot_data %>%
  
  select(plt_cn, prev_plt_cn, plot_id) %>%
  
  # expand it to get 1 row per species 
  tidyr::expand(nesting(plt_cn, prev_plt_cn, plot_id),
         species = c('ABCO', 'CADE27', 'PILA', 'PIPO', 'PSME', 'QUKE', 'OTHER')) %>%

  # add in the counts from the bigger size classes from the trees data 
  # for the remeasurement, **BUT ONLY INCLUDING UNTAGGED TREES**
  left_join(new_recruits %>%
              filter(is.na(prev_tre_cn) & tree_status=='live') %>%
              mutate(count = 1) %>%
              select(plt_cn, plot_id, species, count) %>%
              group_by(plt_cn, plot_id, species) %>%
              summarise(count = sum(count, na.rm = TRUE)) %>%
              ungroup(),
            by = c('plt_cn' = 'plt_cn',
                   'plot_id' = 'plot_id',
                   'species' = 'species')) %>%
  
  # fill in the missings with 0s
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  
  select(plot_id, species, count)

recruits_data

#### make size distribution data frame #########################################

# start with the subplots data
sizedist_data = 
  
  plot_data %>%
  
  # filter to only subplots where both initial and remeasurement had manual 
  # greater than or equal to 2; no longer necessary if we're only looking at 
  # new recruits >= 1" DBH
  #filter(inv_manual.init >= 2.0 & inv_manual.re >= 2.0) %>%
  
  select(plt_cn, prev_plt_cn, plot_id) %>%
  
  # expand it to get 1 row per size bin and species (should be 
  # 17232 * 7 spp * 100 bins = 12062400 rows)
  tidyr::expand(nesting(plt_cn, prev_plt_cn, plot_id),
         species = c('ABCO', 'CADE27', 'PILA', 'PIPO', 'PSME', 'QUKE', 'OTHER'),
         dbh_class = 
             cut(seq(from = 1.5, to = 99.5, by = 1),
                 breaks = seq(from = 1, to = 100, by = 1),
                 labels = FALSE,
                 right = FALSE)) %>%
  
  # add in the counts for each size classe from an aggregated 
  # treelist 
  left_join(treelist %>%
              filter(tree_status == 'live')  %>%
              mutate(dbh_class = cut(dbh_in,
                                    breaks = 
                                      seq(from = 1, to = 100, by = 1),
                                    labels = FALSE,
                                    right = FALSE)) %>%
              select(plt_cn, plot_id, species, dbh_class, tpa_unadj) %>%
              group_by(plt_cn, plot_id, species, dbh_class) %>%
              # note that TPA-unadj is for summed subplots; to get tpa for 
              # an individual subplot need to multiply by four
              summarise(tpa_unadj = sum(tpa_unadj, na.rm = TRUE)) %>%
              ungroup(),
            by = c('plt_cn' = 'plt_cn',
                   'plot_id' = 'plot_id',
                   'species' = 'species',
                   'dbh_class' = 'dbh_class')) %>%
  
  left_join(treelist %>%
              filter(tree_status == 'live') %>%
              mutate(dbh_class = cut(dbh_in,
                                     breaks = 
                                       seq(from = 1, to = 100, by = 1),
                                     labels = FALSE,
                                     right = FALSE)) %>%
              select(plt_cn, plot_id, species, dbh_class, tpa_unadj) %>%
              group_by(plt_cn, plot_id, species, dbh_class) %>%
              summarise(tpa_unadj= sum(tpa_unadj, na.rm = TRUE)) %>%
              ungroup(),
            by = c('prev_plt_cn' = 'plt_cn',
                   'plot_id' = 'plot_id',
                   'species' = 'species',
                   'dbh_class' = 'dbh_class'),
            suffix = c('.re', '.init')) %>%
  
  
  
  # fill in the NAs with 0s
  mutate(tpa_unadj.init = ifelse(is.na(tpa_unadj.init),
                                   0,
                                   tpa_unadj.init),
         tpa_unadj.re = ifelse(is.na(tpa_unadj.re),
                                0,
                                tpa_unadj.re)) %>%
  select(plt_cn, prev_plt_cn, plot_id, species, dbh_class,
         tpa_unadj.init, tpa_unadj.re) %>%
  arrange(plot_id, species, dbh_class)


sizedist_data$plot_id[1]

treelist %>% filter(plot_id == '41-0-47-92433') %>% arrange(species, dbh_in, plt_cn)


test = 
  sizedist_data %>%
  filter(species == 'PILA' &
           is.element(plot_id,
                      mort_data %>%
                        filter(species=='PILA') %>%
                        pull(plot_id))&
           is.element(plot_id,
                      growth_data %>%
                        filter(species=='PILA') %>%
                        pull(plot_id))) %>%
  group_by(plt_cn, prev_plt_cn, plot_id) %>%
  summarise(tpa_unadj.init = sum(tpa_unadj.init),
            tpa_unadj.re = sum(tpa_unadj.re)) %>%
  ungroup() %>%
  filter(tpa_unadj.init==0)

test

treelist %>% filter(plt_cn =='24988722010900' & species=='PILA')  %>% print(width = Inf)

sizedist_data %>% filter(plot_id=='6-2-105-93474' & species=='PILA')%>% print(n = Inf)


 
#### size classes metadata #####################################################

# the model also needs metadata about the size bins for the bins included 
# in the recruitment model; specifically the upper 
# and lower bounds for each bin, the midpoint of each bin, and the plot size 
# of each bin
size_metadata = 
  data.frame(bin_midpoint = 
               seq(from = 1.5, to = 99.5, by = 1)) %>%
  mutate(bin_id = cut(bin_midpoint,
                      breaks = seq(from = 1, to = 100, by = 1),
                      labels = FALSE,
                      right = FALSE),
         bin_lower = seq(from = 1, to = 99, by = 1),
         bin_upper = c(seq(from = 2, to = 99, by = 1), 400),
         
         # <5" dbh are measured on a 6.8' radius (.00333ac) microcplot
         # >= 5" dbh measured on a 24' radius (0.0415ac) subplot
         # including the macroplots here because we want to use the whole a
         # vector and stan won't allow NAs, but we only use the areas for the 
         # smalles classes so it doesn't matter that the macro breakpoint 
         # diameter is inconsistent, because we only need this info for the 
         # size classes included as responses in the recruitment submodel; 
         # min macroplot dbh is 24"
         plot_area_ac = 
           c(rep(pi*(6.8**2)*2/43560, times = 4), # saplings on microplot
             rep(pi*(24**2)*4/43560, times = 19), # small trees on subplot
             rep(pi*(58.9**2)*4/43560, times = 76)) # big trees macroplot; 
         # incoorectly assuming that the macroplot dbh is 24" everywhere (it 
         # varies) but it doesn't matter because this data doesn't get used 
         # anywhere; just cant have it be NA because stan wants the vector for 
         # the small classes
         )


# to avoid evicting big trees, the upper bound for the largest size class 
# needs to be really high
size_metadata[size_metadata$bin_id==100,'bin_upper'] = 400

# join in the recruit size distribution
size_metadata$r = 
  c(r, rep(0, times = 99-length(r)))

size_metadata

#### write results #############################################################


head(mort_data)
head(growth_data)
head(plot_data)
head(sizedist_data)

write.csv(mort_data,
          here::here('02-data',
                     '01-preprocessed',
                     'mort_data.csv'),
          row.names = FALSE)
saveRDS(mort_data,
        here::here('02-data',
                   '01-preprocessed',
                   'mort_data.rds'))

write.csv(growth_data,
          here::here('02-data',
                     '01-preprocessed',
                     'growth_data.csv'),
          row.names = FALSE)
saveRDS(growth_data,
        here::here('02-data',
                   '01-preprocessed',
                   'growth_data.rds'))

write.csv(sizedist_data,
          here::here('02-data',
                     '01-preprocessed',
                     'sizedist_data.csv'),
          row.names = FALSE)
saveRDS(sizedist_data,
        here::here('02-data',
                   '01-preprocessed',
                   'sizedist_data.rds'))

write.csv(plot_data,
          here::here('02-data',
                     '01-preprocessed',
                     'plot_data.csv'),
          row.names = FALSE)
saveRDS(plot_data,
        here::here('02-data',
                   '01-preprocessed',
                   'plot_data.rds'))

write.csv(recruits_data,
          here::here('02-data', 
                     '01-preprocessed', 
                     'recruits_data.csv'),
          row.names = FALSE)
saveRDS(recruits_data,
        here::here('02-data',
                   '01-preprocessed',
                   'recruits_data.rds'))

write.csv(size_metadata,
          here::here('02-data',
                     '01-preprocessed',
                     'size_metadata.csv'),
          row.names = FALSE)
saveRDS(size_metadata,
        here::here('02-data',
                   '01-preprocessed',
                   'size_metadata.rds'))

