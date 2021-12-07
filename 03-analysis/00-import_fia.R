
#### setup #####################################################################

library(here)
library(tidyverse)
library(DBI)
library(RSQLite)

# read in the relevant tables;
# https://www.fs.usda.gov/pnw/tools/pnw-fiadb-forest-inventory-and-analysis-database
# as of 12/6/21 (they change the database ffs)
fiadb = dbConnect(RSQLite::SQLite(),
                  here::here('02-data',
                             '00-source',
                             'FIAPNW',
                             'FIAPNW.db'))
tablenames = dbListTables(fiadb)

tablenames

fiatables = 
  lapply(X = 
           c('COND', 'COND_PNW', 'PLOT', 'PLOT_PNW', 'REF_FOREST_TYPE',
             'REF_FOREST_TYPE_GROUP', 'REF_SPECIES', 'SEEDLING', 'SUBPLOT', 
             'SUBPLOT_PNW', 'SUBP_COND', 'SUBP_COND_CHNG_MTRX', 'TREE', 
             'TREE_PNW'),
         FUN = function(t){dbReadTable(fiadb,t)})

names(fiatables) = 
  c('COND', 'COND_PNW', 'PLOT', 'PLOT_PNW', 'REF_FOREST_TYPE',
    'REF_FOREST_TYPE_GROUP', 'REF_SPECIES', 'SEEDLING', 'SUBPLOT', 
    'SUBPLOT_PNW', 'SUBP_COND', 'SUBP_COND_CHNG_MTRX', 'TREE', 
    'TREE_PNW')
  
dbDisconnect(fiadb)

#### initial read of FIA tables ################################################

db.PLOT = 
  fiatables$PLOT %>%
  as_tibble() %>%
  
  # create inv_date
  mutate(inv_date = 
           as.Date(paste(MEASYEAR, MEASMON, MEASDAY,
                         sep = '/'))) %>%
  # keep only the still-useful columns
  select(CN, STATECD, INVYR, UNITCD, COUNTYCD, PLOT, 
       inv_date, KINDCD, DESIGNCD, PLOT_STATUS_CD,
       PLOT_NONSAMPLE_REASN_CD, PLOT_NONSAMPLE_REASN_CD_PNWRS,
       PREV_PLT_CN, ECOSUBCD, ELEV, LAT, LON,
       MACRO_BREAKPOINT_DIA, TOPO_POSITION_PNW)

plot_plot_cats = 
  function(response){
    ggplot(data = db.PLOT,
           aes(x = as.character(INVYR), fill = as.character(.data[[response]])))+
      geom_bar()+
      theme_minimal()+
      labs(fill = response,
           title = response)
  }

purrr::map(names(db.PLOT)[c(8,9,10)],
           ~plot_plot_cats(.x))


db.COND = 
  fiatables$COND %>%
  select(CN, PLT_CN, STATECD, PLOT, CONDID,
         COND_STATUS_CD, 
         CONDPROP_UNADJ,
         ASPECT, 
         
         # even coverage
         DSTRBCD1, DSTRBCD2, DSTRBCD3,
         
         # even coverage
         DSTRBYR1, DSTRBYR2, DSTRBYR3,
         
         FORTYPCDCALC, MIXEDCONFCD, PHYSCLCD,
         PRESNFCD, SICOND, SLOPE, STDORGCD,
         
         #2001-2010
         #STND_COND_CD_PNWRS, 
         
         # 2001-2010 
         #STUMP_CD_PNWRS,
         
         # even coverage
         TRTCD1, TRTCD2, TRTCD3,
         
         # even coverage
         TRTYR1, TRTYR2, TRTYR3,
         
         LIVE_CANOPY_CVR_PCT,
         LIVE_MISSING_CANOPY_CVR_PCT) %>%
  
  # keep only the inventory events (PLT_CNs) we kept 
  # above
  #filter(is.element(PLT_CN, db.PLOT$CN)) %>%
  
  # join with the cond PNW table
  left_join(.,
            fiatables$COND_PNW %>%
              select(CND_CN,
                     CONDPROP_ADJ_2010,
                     FORESTLAND_YN,
                     
                     # only good in 2010-2011
                     #FIRE_CD_PNWRS,
                     
                     # only after 2010, then good
                     LAND_STATUS_CD,
                     
                     # only after 2008, much better coverage after 2010, 
                     # lots of '0' to NA after 2017 in cd2 and cd3
                     TRTCD1_PNWRS,TRTCD2_PNWRS,TRTCD3_PNWRS),
            
            # all NA
            #TRTYR1_PNWRS,TRTYR2_PNWRS,TRTYR3_PNWRS,
            
            # sparse except for 2010-2011
            #HIST_DSTRBCD1_PNWRS,HIST_DSTRBCD2_PNWRS,HIST_DSTRBCD3_PNWRS,
            
            # very sparse except 2010-2011
            #HIST_DSTRBYR1_PNWRS,HIST_DSTRBYR2_PNWRS,HIST_DSTRBYR3_PNWRS,
            
            # present only after 2009, sparse except for 2010-2011
            #HIST_TRTCD1_PNWRS,HIST_TRTCD2_PNWRS,HIST_TRTCD3_PNWRS,
            
            # present only after 2009, very sparse except for 2010-2011
            #HIST_TRTYR1_PNWRS,HIST_TRTYR2_PNWRS,HIST_TRTYR3_PNWRS),
            by = c('CN' = 'CND_CN')) %>%
  
  # join in the forest type code and forest type group code
  left_join(.,
            fiatables$REF_FOREST_TYPE %>%
              select(REF_FORTYPCD, forest_type = MEANING, TYPGRPCD),
            by = c('FORTYPCDCALC' = 'REF_FORTYPCD')) %>%
  left_join(.,
            fiatables$REF_FOREST_TYPE_GROUP %>%
              select(TYPGRPCD, forest_type_group = MEANING, for_type_grp = ABBR),
            by = c('TYPGRPCD' = 'TYPGRPCD'))


ggplot(data = db.COND %>% left_join(db.PLOT, by = c('PLT_CN'='CN')) %>% filter(INVYR<=2021),
       aes(x = INVYR, fill = is.na(SLOPE)))+
  geom_bar()+
  theme_minimal()

# lets see how good the coverage for condition code data is over time...
testcond = 
  db.COND %>% 
  left_join(db.PLOT,
            by = c('PLT_CN' = 'CN')) %>%
  filter(INVYR<=2021)
plot_condition_cats = 
  function(response){
    ggplot(data = testcond,
           aes(x = INVYR, 
               fill = ifelse(is.na(.data[[response]]),
                             'NA',
                             ifelse(.data[[response]]=='0',
                                    '0',
                                    'other'))))+
      geom_bar()+
      theme_minimal()+
      labs(fill = response)
  }

purrr::map(names(testcond)[c(9:14, 22:29, 31:35)],
           ~plot_condition_cats(.x))

db.SUBPLOT = 
  fiatables$SUBPLOT %>%
  select(CN,
         PLT_CN,
         SUBP,
         SUBP_STATUS_CD,
         CONDLIST,
         MACRCOND,
         SUBPCOND,
         MICRCOND,
         ASPECT,
         SLOPE) %>%
  left_join(.,
            fiatables$SUBPLOT_PNW %>%
              select(SBP_CN,
                     BURN_ASSESS_CD_PNWRS,
                     MECH_ASSESS_CD_PNWRS,),
            by = c('CN' = 'SBP_CN'))


db.SUBP_COND = 
  fiatables$SUBP_COND %>%
  select(CN, PLT_CN, SUBP, CONDID,
         MACRCOND_PROP, SUBPCOND_PROP, MICRCOND_PROP)


db.TREE = 
  fiatables$TREE %>%
  select(CN,
         PLT_CN,
         CONDID,
         INVYR,
         STATECD,
         UNITCD,
         COUNTYCD,
         PLOT,
         SUBP,
         TREE,
         SPCD, # always there
         SPGRPCD, # always there
         STATUSCD,
         TPA_UNADJ,
         DIA, # a bunch of NA diameters in the remeasurements post 2010
         HT,
         CCLCD,
         CR,
         AGENTCD, # a few starting 2004, ramping up to many after 2010
         AGENTCD_PNWRS, # a few starting 2004, ramping up to many after 2010
         
         # always present, but always '0'
         #DAMAGE_AGENT_CD1_PNWRS, DAMAGE_AGENT_CD2_PNWRS, DAMAGE_AGENT_CD3_PNWRS, 
         
         # 2013 onwards
         DAMAGE_AGENT_CD1, DAMAGE_AGENT_CD2, DAMAGE_AGENT_CD3,
         
         # many 2000-2004, some 2005-2012; looks like the difference 
         # was that they stopped recording '0' and started leaving it 
         # as NA; the amount of non-zero non-NA values stays consistent
         # through 2012; this is what i want
         DMG_AGENT1_CD_PNWRS, DMG_AGENT2_CD_PNWRS, DMG_AGENT3_CD_PNWRS,
         
         #DAMTYP1, DAMTYP2, # 2000-2003 
         #MORTCD, # never populated
         SNAG_DIS_CD_PNWRS, # c. 2010 onwards
         MORTYR, MORTYR_PNWRS, # base 2004-2019, better coverage after 2010 
         #OLD_TREE_NO_PNWRS, # 2000-2009, no idea what this is
         PREV_CONDID_PNWRS, PREV_HT_PNWRS, PREV_TRE_CN, 
         PREV_STATUS_CD, PREV_STATUS_CD_PNWRS, # 2004-2009 pnwrs; 2006-2019 base
         PREVCOND, PREVSUBC, 
         PREVDIA, 
         
         RECONCILECD, # RECONCILECD_P2A, RECONCILECD_PNWRS, 
         #2004-2018 in pnwrs, 2005-2019 in base
         
         SEVERITY1_CD_PNWRS, SEVERITY2_CD_PNWRS, SEVERITY3_CD_PNWRS, 
         # damage severity in 2001-2004
         
         SEVERITY1A_CD_PNWRS, SEVERITY1B_CD_PNWRS, 
         SEVERITY2A_CD_PNWRS, SEVERITY2B_CD_PNWRS, # damage severity in 2005-2012
         #MIST_CL_CD,
         #MIST_CL_CD_PNWRS 
  ) %>%
  
  left_join(.,
            fiatables$TREE_PNW %>%
              select(TRE_CN,
                     CND_CN,
                     TPA_ADJ,
                     TPA_ADJ_2010,
                     INC10YR_PNWRS,
                     INC5YR_PNWRS,
                     INC5YRHT_PNWRS),
            by = c('CN' = 'TRE_CN'))


# doing some basic data exploration: what's the distribution of these 
# different categories across years?
plot_across_years = 
  function(response){
    
    ggplot(data = db.TREE,
           aes(x = as.character(INVYR), 
               fill = ifelse(is.na(.data[[response]]),
                             'NA',
                             ifelse(.data[[response]]=='0',
                                    '0',
                                    'other'))))+
      geom_bar()+
      theme_minimal()+
      labs(fill = response)
    
  }

checking_categorical = 
  purrr::map(names(db.TREE),
             ~plot_across_years(.x))

checking_categorical



db.SEEDLING = 
  fiatables$SEEDLING %>%
  select(CN,
         PLT_CN,
         SPCD,
         TPA_UNADJ,
         TREECOUNT,
         TREECOUNT_CALC) %>%
  as_tibble()

db.REF_SPECIES = 
  fiatables$REF_SPECIES %>%
  select(SPCD, SPECIES, SPECIES_SYMBOL)

db.REF_FOREST_TYPE = 
  fiatables$REF_FOREST_TYPE %>%
  select(REF_FORTYPCD, MEANING, TYPGRPCD)

db.REF_FOREST_TYPE_GROUP = 
  fiatables$REF_FOREST_TYPE_GROUP %>%
  select(TYPGRPCD, MEANING, ABBR)

#### create crosswalks for codes ###############################################

# copying the relevant crosswalks from the FIA documentation
cond_dstrbcds = 
  data.frame(
    code = 
      c(0, 10, 11, 12, 20, 21, 22, 30, 31, 32, 40, 41, 42, 43, 44,
        45, 56, 50, 51, 52, 53, 54, 60, 70, 80, 90, 91, 92, 93, 94,
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


cond_trtpnwrscds = 
  data.frame(
    code = c(0, 
             10:17,
             20,
             30:33,
             40,
             50:52,
             60),
    cond_trtpnwrsdesc = 
      c('none',
        'cutting',
        'clearcut',
        'partialcut_heavy',
        'partialcut_light',
        'firewoodcut',
        'incidentalcut',
        'pct',
        'improvementcut',
        'siteprep',
        'artificial_regen',
        'planting_throughout',
        'planting_gaps',
        'underplanting',
        'natural_regen',
        'other_silv',
        'conversion',
        'clean_release',
        'chaining'),
    cond_trtpnwrstype = 
      c('none',
        rep('cutting', 8),
        'siteprep',
        rep('regen (artificial)', 4),
        'regen (natural)',
        'other',
        'other',
        'cutting',
        'other')
  ) %>%
  mutate(across(c(cond_trtpnwrsdesc, cond_trtpnwrstype), as.character)) %>%
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


cond_storgcds = 
  data.frame(
    code = c(0, 1),
    cond_storgtype = c('natural', 'artificial')
  ) %>%
  mutate(across(cond_storgtype, as.character)) %>%
  as_tibble()

tree_agentcd = 
  data.frame(
    code = c('0', '10', '20', '30', '40', '50', '60', '70', '80'),
    tree_agenttype = 
      c('none', 'insect', 'disease', 'fire', 'animal', 'weather', 'vegetation',
        'unknown/other', 'silviculture/landclearing')) %>%
  mutate(code = as.integer(as.character(code)),
         tree_agenttype = as.character(tree_agenttype)) %>%
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



# damage codes to cover:
damageagentcds = 
  unique(c(db.TREE$DAMAGE_AGENT_CD1, 
           db.TREE$DAMAGE_AGENT_CD2, 
           db.TREE$DAMAGE_AGENT_CD3))[
             order(unique(c(
               db.TREE$DAMAGE_AGENT_CD1,
               db.TREE$DAMAGE_AGENT_CD2,
               db.TREE$DAMAGE_AGENT_CD3
             )))
           ]

damageagentcds[!is.na(damageagentcds)]

tree_damageagentcds = 
  data.frame(
    code = damageagentcds[!is.na(damageagentcds)],
    tree_damageagentdesc = 
      c('none',
        'insects', 'barkbeetles', 'spruce_beetle',
        'defoliators', 'western_spruce_budworm', 'budworm',
        'chrewing_insects', 'sucking_insects', 'baslsam_wooly_adelgid', 'boring_insects',
        'general_disease', rep('root_rot', 11), rep('cankers', 5),
        'stem_decay', rep('parasitic_plants', 4), 'decline_complexes', 
        rep('foliage_disease', 2), 'stem_rust', 'wpbr', 'stem_rust',
        'broom_rust', 'fire', rep('wild_animals', 11), 'domestic_animals',
        rep('abiotic', 4), 'suppression', rep('human', 3), 'harvest',
        rep('other', 13)),
    tree_damageagenttype = 
      c('none',
        rep('insects', 10),
        rep('disease', 1+11+5+1+4+1+2+3+1),
        'fire',
        rep('animal', 12),
        rep('abiotic', 4),
        'suppression',
        rep('animal', 3),
        'harvest', 
        rep('other', 13))
  ) %>%
  mutate(across(c(tree_damageagentdesc, tree_damageagenttype), as.character)) %>%
  as_tibble()

tree_dmgagentpnwrscds = 
  data.frame(
    code = 
      c('1', '2', '3', '4', '5', '6', '7', '8', '9', '26',
        '10', '11', '12' ,'13', '14', '15', '16', '17', '18', '19',
        '60', '61', '62', '63', '65', '66',
        '36',
        '31',
        '20', '21', '22', '23', '24', '25',
        '33', '40', '41', '42', '43', '44', '45',
        '32',
        '46', '47', '48', '49',
        '50', '51',
        '55', '56', '57', '58', '59',
        '70', '71', '72', '73', '74', '75', '76', '77', '78',
        '80', '81', '82', '83', '84', '85', '86', '87',
        '90', '91', '92', '93', '94',
        '95', '96', '97', '98', '99'),
    tree_dmgagentpnwrsdesc = 
      c(rep('barkbeetles', 10), 
        rep('defoliators', 10),
        rep('rootdis', 6),
        rep('wpbr', 1),
        rep('suddenoakdeath', 1),
        rep('otherinsects', 6),
        rep('cankers', 7),
        rep('pitchcanker', 1),
        rep('stemdecay', 4),
        rep('suppression', 1),
        rep('deformed', 1),
        rep('foliardis', 5),
        rep('animal', 9),
        rep('weather', 5), 'drought', rep('weather', 2),
        rep('physicalinjury', 2), 'fire', rep('physicalinjury', 2),
        rep('defect', 5)),
    tree_dmgagentpnwrstype = 
      c(rep('insects', 10),
        rep('insects', 10),
        rep('disease', 6),
        rep('disease', 1),
        rep('disease', 1),
        rep('insects', 6),
        rep('disease', 7+1+4),
        rep('suppression', 1),
        rep('other', 1),
        rep('disease', 5),
        rep('animal', 9),
        rep('abiotic', 5), 'drought', rep('abiotic', 2),
        rep('other', 2), 'fire', rep('other', 2),
        rep('other', 5))
  ) %>%
  mutate(across(c(tree_dmgagentpnwrsdesc, tree_dmgagentpnwrstype), as.character),
         code = as.integer(as.character(code))) %>%
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

#### prepare context data ######################################################

# the goal is to have a table with one row per observation of a condition (stand 
# on a subplot at a specific inventory date)
# context
context_full = 
  
  # start with the plot table, with one row per plot observation
  db.PLOT %>%
  as_tibble() %>%
  rename(PLT_CN = CN) %>%
  
  # join in the subplot info, with one row per subplot observation
  left_join(db.SUBPLOT %>%
              rename(SUBP_CN = CN),
            by = c('PLT_CN' = 'PLT_CN'),
            suffix = c('.PLOT', '.SUBPLOT')) %>%
  
  # join the subplot conditions; this creates one row per condition:subplot:plot:time
  left_join(db.SUBP_COND %>%
              rename(SUBPCOND_CN = CN),
            by = c('PLT_CN' = 'PLT_CN', 
                   'SUBP' = 'SUBP'),
            suffix = c('.context', '.SUBP_COND')) %>%
  
  # join the condition details
  left_join(db.COND %>%
              rename(CND_CN = CN),
            by = c('PLT_CN' = 'PLT_CN',
                   'CONDID' = 'CONDID',
                   'STATECD' = 'STATECD',
                   'PLOT' = 'PLOT'),
            suffix = c('.context', '.COND'))

#### join code descriptions ####################################################

context_full = 
  context_full %>%
  
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
  
  left_join(cond_physclcds %>%
              select(PHYSCLCD = code,
                     PHYSCLCD.DESC = physcldesc,
                     PHYSCLCD.TYPE = physcltype)) %>%
  
  left_join(cond_statuscds %>%
              select(COND_STATUS_CD = code,
                     COND_STATUS_CD.TYPE = cond_statustype)) %>%
  
  left_join(cond_storgcds %>%
              select(STDORGCD = code,
                     STDORGCD.TYPE = cond_storgtype)) %>%
  
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
  
  left_join(cond_trtpnwrscds %>%
              select(TRTCD1_PNWRS = code,
                     TRTCD1_PNWRS.DESC = cond_trtpnwrsdesc,
                     TRTCD1_PNWRS.TYPE = cond_trtpnwrstype)) %>%
  
  left_join(cond_trtpnwrscds %>%
              select(TRTCD2_PNWRS = code,
                     TRTCD2_PNWRS.DESC = cond_trtpnwrsdesc,
                     TRTCD2_PNWRS.TYPE = cond_trtpnwrstype)) %>%
  
  left_join(cond_trtpnwrscds %>%
              select(TRTCD3_PNWRS = code,
                     TRTCD3_PNWRS.DESC = cond_trtpnwrsdesc,
                     TRTCD3_PNWRS.TYPE = cond_trtpnwrstype)) %>%
  
  left_join(plot_designcds %>%
              select(DESIGNCD = code,
                     DESIGNCD.type = plot_designtype)) %>%
  
  left_join(plot_kindcds %>%
              select(KINDCD = code,
                     KINDCD.type = plot_kindtype)) %>%
  
  left_join(plot_statuscds %>%
              select(PLOT_STATUS_CD = code,
                     PLOT_STATUS_CD.TYPE = plot_statustype)) %>%
  
  left_join(plot_topopositioncds %>%
              select(TOPO_POSITION_PNW = code,
                     TOPO_POSITION_PNW.TYPE = plot_topopositiontype)) %>%
  
  left_join(subp_statuscds %>%
              select(SUBP_STATUS_CD = code,
                     SUBP_STATUS_CD.TYPE = plot_statustype)) %>%
  
  left_join(plot_nonsamplereasncds %>%
              select(PLOT_NONSAMPLE_REASN_CD = code,
                     PLOT_NONSAMPLE_REASN_CD.TYPE = plot_nonsamplereasntype))
  

treelist_full = 
  
  db.TREE %>%
  as_tibble() %>%
  rename(TRE_CN = CN) %>%
  left_join(tree_agentcd %>%
              select(AGENTCD = code,
                     AGENTCD.TYPE = tree_agenttype)) %>%
  
  left_join(tree_cclasscds %>%
              select(CCLCD = code,
                     CCLCD.TYPE = cclctype)) %>%
  
  left_join(tree_damageagentcds %>%
              select(DAMAGE_AGENT_CD1 = code,
                     DAMAGE_AGENT_CD1.DESC = tree_damageagentdesc,
                     DAMAGE_AGENT_CD1.TYPE = tree_damageagenttype)) %>%
  
  left_join(tree_damageagentcds %>%
              select(DAMAGE_AGENT_CD2 = code,
                     DAMAGE_AGENT_CD2.DESC = tree_damageagentdesc,
                     DAMAGE_AGENT_CD2.TYPE = tree_damageagenttype)) %>%
  
  left_join(tree_damageagentcds %>%
              select(DAMAGE_AGENT_CD3 = code,
                     DAMAGE_AGENT_CD3.DESC = tree_damageagentdesc,
                     DAMAGE_AGENT_CD3.TYPE = tree_damageagenttype)) %>%
  
  left_join(tree_dmgagentpnwrscds %>%
              select(DMG_AGENT1_CD_PNWRS = code,
                     DMG_AGENT1_CD_PNWRS.DESC = tree_dmgagentpnwrsdesc,
                     DMG_AGENT1_CD_PNWRS.TYPE = tree_dmgagentpnwrstype)) %>%
  
  left_join(tree_dmgagentpnwrscds %>%
              select(DMG_AGENT2_CD_PNWRS = code,
                     DMG_AGENT2_CD_PNWRS.DESC = tree_dmgagentpnwrsdesc,
                     DMG_AGENT2_CD_PNWRS.TYPE = tree_dmgagentpnwrstype)) %>%
  
  left_join(tree_dmgagentpnwrscds %>%
              select(DMG_AGENT3_CD_PNWRS = code,
                     DMG_AGENT3_CD_PNWRS.DESC = tree_dmgagentpnwrsdesc,
                     DMG_AGENT3_CD_PNWRS.TYPE = tree_dmgagentpnwrstype)) %>%
  
  
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
  left_join(db.REF_SPECIES %>%
              select(SPCD,
                     spp = SPECIES_SYMBOL))

#### select useful columns from context and treelist ###########################

names(context_full)[order(names(context_full))]

context_full = 
  context_full %>%
  select(
    
    ## observation IDs ##
    PLT_CN, # UID FOR THE PLOT-LEVEL SAMPLING EVENT
    SUBP_CN,
    SUBPCOND_CN, # UID FOR THE SUBPLOT CONDITION TABLE
    CND_CN,  # exists in the tree table, not consistently recorded there
    STATECD, # CODE FOR THE STATE
    UNITCD,
    COUNTYCD,
    PLOT, # ID FOR THE PLOT; CONFUSINGLY NON-UNIQUE (needs state and countycd), 
    # and also not always consistent with the PLOT in db.TREE for a given PLT_CN, 
    # which is weird...
    CONDID, # CONDITION CODE FOR THE PLT_CN:SUBP CONDLIST THAT THE TREE IS IN;
    # ALWAYS PRESENT
    SUBP, # SUBPLOT id
    CONDLIST, # LIST OF CONDITIONS (UP TO 4) ON THE SUBPLOT
    
    ## observation metadata ##
    INVYR, # NOMINAL YEAR OF THE INVENTORY
    inv_date, # actual measurement date
    KINDCD, KINDCD.type, # TYPE OF SAMPLING: 
    DESIGNCD, DESIGNCD.type, # PLOT SAMPLING DESIGN; 
    PLOT_STATUS_CD, PLOT_STATUS_CD.TYPE, # STATUS CODE
    PLOT_NONSAMPLE_REASN_CD, PLOT_NONSAMPLE_REASN_CD.TYPE,
    MACRO_BREAKPOINT_DIA, # MINIMUM SIZE FOR TREES IN THE MACROPLOT; VARIES
    SUBP_STATUS_CD, SUBP_STATUS_CD.TYPE,
    CONDPROP_UNADJ, # the unadjust proportion of the plot that is 
    # in the condition
    
    ## matching across time ##
    PREV_PLT_CN, # PREVIOUS CONTROL NUMBER; ONLY POPULATED FOR KINDCD = 2
    
    ## environmental data (plot or cond-level)
    PHYSCLCD, PHYSCLCD.DESC, PHYSCLCD.TYPE, # PHYSIOGRAPHIC CLASS CODE
    ECOSUBCD, # ECOREGION SUBSECTION
    ELEV, LAT, LON, # ELEV, LAT, AND LON (NAD83 DATUM) OF THE FUZZED/SWAPPED
    # COORDS; NOT SURE IF ELEV IS TRUE OR FUZZ/SWAP
    TOPO_POSITION_PNW, TOPO_POSITION_PNW.TYPE, # topo position, usually but not always filled
    ASPECT = ASPECT.context, #ASPECT.COND, # aspect, almost always filled for both
    DSTRBCD1, DSTRBCD2, DSTRBCD3, 
    DSTRBCD1.DESC, DSTRBCD2.DESC, DSTRBCD3.DESC,
    DSTRBCD1.TYPE, DSTRBCD2.TYPE, DSTRBCD3.TYPE, 
    # disturbance since the last measurement, 
    # or in the last 5 years for new plots; mort or damage to >25% of the trees
    # in the condition is required; 
    DSTRBYR1, DSTRBYR2, DSTRBYR3, # field-est. year of disturbance; 9999 means 
    # continuous or over time; NA means dstrbcd = 0
    TRTCD1, TRTCD2, TRTCD3, # TYPE OF STAND TREATMENT AFFECTING CONDITION 
    # SINCE LAST REAMEASURE (OR LAST 5 YRS FOR NEW PLOTS); MUST BE 1AC TOTAL 
    # IN SIZE
    TRTCD1.DESC, TRTCD2.DESC, TRTCD3.DESC,
    TRTCD1.TYPE, TRTCD2.TYPE, TRTCD3.TYPE,
    TRTYR1, TRTYR2, TRTYR3,
    LIVE_CANOPY_CVR_PCT, # CANOPY COVER AT TIME OF SAMPLING
    LIVE_MISSING_CANOPY_CVR_PCT, # ESTIMATED CANOPY COVER PRE-DISTURBANCE
    TRTCD1_PNWRS, TRTCD2_PNWRS, TRTCD3_PNWRS, # detailed treatment codes
    TRTCD1_PNWRS.DESC, TRTCD2_PNWRS.DESC, TRTCD3_PNWRS.DESC,
    TRTCD1_PNWRS.TYPE, TRTCD2_PNWRS.TYPE, TRTCD3_PNWRS.TYPE,
    
    ## useful for filtering ##
    COND_STATUS_CD, COND_STATUS_CD.TYPE, # whether the condition is forest, 
    # nonforest, nonsampled, etc 
    STDORGCD, STDORGCD.TYPE, # STAND ORIGIN CODE; 0 IS NATURAL 1 IS ARTIFICIAL
    FORTYPCDCALC, forest_type, TYPGRPCD, forest_type_group, for_type_grp 
    # forest type codes
    
    ## codes to ditch ##
    #TPA_ADJ, TPA_ADJ_2010, # TPA ADJUSTED FOR NONSAMPLING BIAS IN THE STRATA,
    # ONLY COMMON AFTER 2010
    # INC10YR_PNWRS, INC5YR_PNWRS, INC5YRHT_PNWRS, # growth increment, mostly NA
    #BURN_ASSESS_CD_PNWRS, # almost always NA
    #MECH_ASSESS_CD_PNWRS, # almost always NA
    # CONDPROP_ADJ_2010
    # FORESTLAND_YN, # only after 2010
    # MIXEDCONFCD, some pnw forest type flag, usually NA before 2012 and always 
    # NA after
    # PRESNFCD ALMOST ALWAYS NA
    #SICOND, # SITE CONDITOIN CLASS
    # MACRCOND, SUBPCOND, MICRCOND, # CONDITION NUMBER AT THE CENTER OF THE 
    # MACRO, SUB, AND MICROPLOT, not that useful
    #MACRCOND_PROP, SUBPCOND_PROP, MICRCOND_PROP, # THE PROPORTION OF THE 
    # MACRO, SUB, OR MICROPLOT OCCUPIED BY THE CONDITION PRESENT AT EACH'S CENTER
    #LAND_STATUS_CD, # ONLY AFTER 2010; reserved, timberland, etc. 
  )


names(treelist_full)[order(names(treelist_full))]
treelist_full = 
  treelist_full %>%
  select(
    
    ## observation IDs ##
    TRE_CN, # UID for the tree-level sampling event
    PLT_CN, # UID FOR THE PLOT-LEVEL SAMPLING EVENT
    CONDID, # CONDITION CODE FOR THE PLT_CN:SUBP CONDLIST THAT THE TREE IS IN;
    # ALWAYS PRESENT
    STATECD, # CODE FOR THE STATE
    UNITCD,
    COUNTYCD,
    PLOT, # ID FOR THE PLOT; CONFUSINGLY NON-UNIQUE and not always matching the context
    SUBP, # SUBPLOT id FOR THE TREE
    TREE, # TREE IDENTIFIER; UNIQUE IN A SUBPLOT, CONSITENT ACROSS TIME WHEN 
    # DESIGNCD IS CONSISTENT
    # CND_CN, # was not consistently recorded in the tree table
    
    ## observation metadata ##
    INVYR, # NOMINAL YEAR OF THE INVENTORY
    
    ## tree-level data ##
    SPCD, SPGRPCD, spp, # species and species-group codes for the tree, use lookup
    STATUSCD, STATUSCD.TYPE, # status code for the tree;
    TPA_UNADJ, # UNADJUSTED SCALING FACTOR; na for ~10% of reamasure trees
    DIA, # diameter (in)
    HT,
    CCLCD, CCLCD.TYPE, # crown class position code
    CR, # crown ratio
    AGENTCD, AGENTCD.TYPE, # cause of death code
    # (Core: all remeasured plots when the tree was alive at the previous 
    # visit and at revisit is dead or removed OR the tree is standing dead in 
    # the current inventory and the tree is ingrowth, through growth, or a 
    # missed live tree; Core optional: all initial plot visits when tree 
    # qualifies as a mortality tree)
    DAMAGE_AGENT_CD1, DAMAGE_AGENT_CD2, DAMAGE_AGENT_CD3,
    DAMAGE_AGENT_CD1.DESC, DAMAGE_AGENT_CD2.DESC, DAMAGE_AGENT_CD3.DESC,
    DAMAGE_AGENT_CD1.TYPE, DAMAGE_AGENT_CD2.TYPE, DAMAGE_AGENT_CD3.TYPE,
    # DAMAGE AGENT CODE; ONLY PRESENT 2013 ONWARDS
    DMG_AGENT1_CD_PNWRS, DMG_AGENT2_CD_PNWRS, DMG_AGENT3_CD_PNWRS,
    DMG_AGENT1_CD_PNWRS.DESC, DMG_AGENT2_CD_PNWRS.DESC, DMG_AGENT3_CD_PNWRS.DESC,
    DMG_AGENT1_CD_PNWRS.TYPE, DMG_AGENT2_CD_PNWRS.TYPE, DMG_AGENT3_CD_PNWRS.TYPE,
    # DAMAGE AGENT CODE; ONYL PRESENT UP TO 2012; MANY 0S BECOME NAS AFTER 2005
    MORTYR, # MORTYR_PNWRS, estimated year of mortality
    
    ## matching across time ##
    PREV_TRE_CN, # (remeasure only) UID for this tree at the previous sample
    PREV_STATUS_CD, # 1 = LIVE 2 = STANDING DEAD
    PREVCOND,  # PREVIOUS CONDITION CLASS ID
    PREVSUBC,
    RECONCILECD, RECONCILECD.TYPE, # reason why the tree appears / disappears 
    # from the inventory at a remeasurement
    SNAG_DIS_CD_PNWRS, SNAG_DIS_CD_PNWRS.TYPE # remeasures only, 
    # snag reason for disappearance 
    
    ## codes to ditch ##
    #TPA_ADJ, TPA_ADJ_2010, # TPA ADJUSTED FOR NONSAMPLING BIAS IN THE STRATA,
    # ONLY COMMON AFTER 2010
    # INC10YR_PNWRS, INC5YR_PNWRS, INC5YRHT_PNWRS, # growth increment, mostly NA
    #AGENTCD_PNWRS, # PNWRS specific agent code, look sabout as common; dropped 
    # because i can't find a crosswalk for it in either the FIA or FIA-PNW docs
    #SEVERITY1_CD_PNWRS, SEVERITY2_CD_PNWRS, SEVERITY3_CD_PNWRS,
    #SEVERITY1A_CD_PNWRS, SEVERITY1B_CD_PNWRS, SEVERITY2A_CD_PNWRS, SEVERITY2B_CD_PNWRS,
    # damage severity codes; 2001-2004 and then 2005-2012; too much info to 
    # deal with atm
    #PREV_CONDID_PNWRS, MOSTLY BLANK
    # PREV_HT_PNWRS, PREVDIA
  )

#### make streamlined versions of each table ###################################

# don't usually want or need all of the ID codes, full descriptions, etc, and 
# i hate the capitalized column names. make streamlined versions of each table 
# for actually using them:
context = 
  context_full %>%
  mutate(plot_id = 
           paste(STATECD,
                 UNITCD,
                 COUNTYCD,
                 PLOT,
                 sep = '-'),
         subp_id = 
           paste(STATECD,
                 UNITCD,
                 COUNTYCD,
                 PLOT,
                 SUBP,
                 sep = '-')) %>%
  select(
    plt_cn = PLT_CN,
    subp_cn = SUBP_CN,
    subpcond_cn = SUBPCOND_CN,
    cnd_cn = CND_CN,
    state_id = STATECD,
    plot_id,
    subp_id,
    cond_id = CONDID,
    condlist = CONDLIST,
    inv_year_nominal = INVYR,
    inv_date = inv_date,
    inv_kind = KINDCD.type,
    inv_design = DESIGNCD.type,
    inv_macrobrk = MACRO_BREAKPOINT_DIA,
    inv_plotstatus = PLOT_STATUS_CD.TYPE,
    inv_nonsamp = PLOT_NONSAMPLE_REASN_CD.TYPE,
    inv_subpstatus = SUBP_STATUS_CD.TYPE,
    condprop_unadj = CONDPROP_UNADJ,
    inv_condstatus = COND_STATUS_CD.TYPE,
    
    prev_plt_cn = PREV_PLT_CN,
    physiographic = PHYSCLCD.TYPE,
    ecosub_cd = ECOSUBCD,
    elev_ft = ELEV,
    lat = LAT, lon = LON,
    aspect = ASPECT,
    topographic = TOPO_POSITION_PNW.TYPE,
    disturb1 = DSTRBCD1.TYPE, disturb2 = DSTRBCD2.TYPE, disturb3 = DSTRBCD3.TYPE,
    disturb1_yr = DSTRBYR1, disturb2_yr = DSTRBYR2, disturb3_yr = DSTRBYR3,
    treat1 = TRTCD1.TYPE, treat2 = TRTCD2.TYPE, treat3 = TRTCD3.TYPE,
    treatpnw1 = TRTCD1_PNWRS.TYPE, treatpnw2 = TRTCD2_PNWRS.TYPE, 
    treatpnw3 = TRTCD3_PNWRS.TYPE,
    treat1_yr = TRTYR1, treat2_yr = TRTYR2, treat3_yr = TRTYR3,
    stand_orig = STDORGCD.TYPE,
    forest_type_group = for_type_grp,
    forest_type
    
  )

names(treelist_full)

treelist = 
  treelist_full %>%
  
  mutate(plot_id = 
           paste(STATECD,
                 UNITCD,
                 COUNTYCD,
                 PLOT,
                 sep = '-'),
         subp_id = 
           paste(STATECD,
                 UNITCD,
                 COUNTYCD,
                 PLOT,
                 SUBP,
                 sep = '-'),
         tree_id = 
           paste(subp_id, TREE)) %>%
  select(
    tre_cn = TRE_CN,
    plt_cn = PLT_CN,
    cond_id = CONDID,
    state_id = STATECD,
    plot_id,
    subp_id,
    tree_id,
    inv_year_nominal = INVYR,
    
    spp,
    tree_status = STATUSCD.TYPE,
    dbh_in = DIA,
    height_ft = HT,
    cclass = CCLCD.TYPE,
    crownrat = CR,
    damage1 = DAMAGE_AGENT_CD1.TYPE, damage2 = DAMAGE_AGENT_CD2.TYPE,
    damage3 = DAMAGE_AGENT_CD3.TYPE,
    damage1_pnw = DMG_AGENT1_CD_PNWRS.TYPE, damage2_pnw = DMG_AGENT2_CD_PNWRS.TYPE,
    damage3_pnw = DMG_AGENT3_CD_PNWRS.TYPE,
    death_cause = AGENTCD.TYPE,
    death_year = MORTYR,
    tpa_unadj = TPA_UNADJ,
    
    prev_tre_cn = PREV_TRE_CN,
    prev_status = PREV_STATUS_CD,
    prev_cond = PREVCOND,
    prev_subc = PREVSUBC,
    reconcile = RECONCILECD.TYPE,
    snag_missing = SNAG_DIS_CD_PNWRS.TYPE
  )

#### write results #############################################################

# has all the columns and all rows
write.csv(treelist_full,
          here::here('02-data', '01-preprocessed', 'fia_treelist_full.csv'),
          row.names = FALSE)
saveRDS(treelist_full,
        here::here('02-data', '01-preprocessed', 'fia_treelist_full.rds'))

# only the more useful columns
write.csv(treelist,
          here::here('02-data', '01-preprocessed', 'fia_treelist.csv'),
          row.names = FALSE)
saveRDS(treelist,
        here::here('02-data', '01-preprocessed', 'fia_treelist.rds'))


# has all the columns and all rows
write.csv(context_full,
          here::here('02-data', '01-preprocessed', 'fia_context_full.csv'),
          row.names = FALSE)
saveRDS(context_full,
        here::here('02-data', '01-preprocessed', 'fia_context_full.rds'))

# only the more useful columns
write.csv(context,
          here::here('02-data', '01-preprocessed', 'fia_context.csv'),
          row.names = FALSE)
saveRDS(context,
        here::here('02-data', '01-preprocessed', 'fia_context.rds'))

