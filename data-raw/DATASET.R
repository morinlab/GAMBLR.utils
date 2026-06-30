#!/usr/bin/env Rscript

library(usethis)
library(readr)
library(tidyr)
library(dplyr)


aliases = read_tsv("inst/extdata/aliases.tsv")

usethis::use_data(aliases,overwrite = T)

hotspot_regions_grch37 =  read_tsv("inst/extdata/hotspot_regions.grch37.tsv") %>%
  dplyr::filter(include == TRUE)
usethis::use_data(hotspot_regions_grch37,overwrite = T)

dcol_ord = colnames(hotspot_regions_grch37)

hotspot_regions_grch37_bed = hotspot_regions_grch37 %>%
  select(chrom,start,end,strand,any_of(colnames(hotspot_regions_grch37))) %>%
  mutate(fake_e = ifelse(end==300000000,TRUE,FALSE)) %>%
  mutate(end=ifelse(end==300000000,start+2,end)) %>%
  mutate(fake_s = ifelse(start==1,TRUE,FALSE)) %>%
  mutate(start=ifelse(start==1,end-2,start))

hotspot_regions_hg38 = liftover(hotspot_regions_grch37_bed,
                                target_build = "hg38",
                                mode="bed")

hotspot_regions_hg38 = mutate(hotspot_regions_hg38,
                              end=ifelse(fake_e,300000000,end),
                              start=ifelse(fake_s,1,start))

hotspot_regions_hg38 = select(hotspot_regions_hg38,all_of(col_ord))

#lost during liftover due to fake start-end
missing_hg38 = filter(hotspot_regions_grch37,
                      !alias %in% hotspot_regions_hg38$alias)


#add prefix
missing_hg38 = mutate(missing_hg38,chrom = paste0("chr",chrom))


hotspot_regions_hg38 = bind_rows(hotspot_regions_hg38,
                                 missing_hg38)


usethis::use_data(hotspot_regions_hg38,overwrite = T)
