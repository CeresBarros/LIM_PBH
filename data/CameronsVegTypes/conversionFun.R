dendro.cov.BA.synthesis.f <- function(x, dry.sp, moist.sp, pure.cutoff) {
  #Calculate the total BA and percent of plot BA by species
  x <- mutate(x, BA.total=sum(BA.weight))
  cov.BA <- ddply(x, .(Species), summarize, Prop.BA.sp=sum(BA.weight)/unique(BA.total)) #cov.BA is used for defining rulesets since there is no repitition of values, making coding simpler
  x <- ddply(x, .(Species), mutate, Prop.BA.sp=sum(BA.weight)/unique(BA.total)) #repeating the Prop.BA.sp calculation through a mutate function is simpler than merging x and cov.BA and ensures that Prop.BA.sp is included in outputs

  #Cover classification ruleset
  if (sum(cov.BA[cov.BA$Species %in% c('QUGA'),]$Prop.BA.sp) >= pure.cutoff &
      sum(cov.BA[cov.BA$Species %in% c('PIPO', 'PSME', 'PIED', 'PIMO2','JUSC', 'JUOC', 'JUOS'),]$Prop.BA.sp) < .05) {
    #Oak woodlands are dominated by oaks with no more dominant tree stature species
    x$Cover.dendro <- 'Oak'
  } else if (sum(cov.BA[cov.BA$Species %in% c('PIED', 'PIMO2','JUSC', 'JUOC', 'JUOS', 'QUGA'),]$Prop.BA.sp) >= pure.cutoff &
             sum(cov.BA[cov.BA$Species %in% c('PIPO', 'PSME'),]$Prop.BA.sp) < .05) {
    #P-J woodlands are dominated by pinyon juniper trees with no more dominant tree stature species
    x$Cover.dendro <- 'PJ'
  }  else if (sum(cov.BA[cov.BA$Species %in% c('PIPO'),]$Prop.BA.sp) >= pure.cutoff &
              sum(cov.BA[cov.BA$Species %in% c('PSME', 'PIFL', 'PIED', 'PIMO2', 'JUSC', 'JUOC', 'JUOS', 'QUGA'),]$Prop.BA.sp) < .30) {
    #Pure PIPO if PIPO is heavily dominant and accompanied by small amount of other species
    x$Cover.dendro <- 'Pure-PIPO'
  } else if(sum(cov.BA[cov.BA$Species %in% c('PIPO'),]$Prop.BA.sp) >= .10 &
            sum(cov.BA[cov.BA$Species %in% c('PIPO'),]$Prop.BA.sp) < pure.cutoff &
            sum(cov.BA[cov.BA$Species %in% c('PSME', 'PIPO', 'PIFL', 'PIED', 'PIMO2', 'JUSC', 'JUOC', 'JUOS', 'QUGA'),]$Prop.BA.sp) >= .50) {
    #DMC if PIPO present at >= 10% but less than 70% and other species are all dry site species
    x$Cover.dendro <- 'DMC-PIPO'
  } else if(sum(cov.BA[cov.BA$Species %in% c('PSME'),]$Prop.BA.sp) >= pure.cutoff &
            sum(cov.BA[cov.BA$Species %in% moist.sp,]$Prop.BA.sp) < .10) {
    if (sum(cov.BA[cov.BA$Species %in% c('JUSC', 'JUOC', 'JUOS', 'PIFL', 'PIED', 'PIMO2', 'QuGA'),]$Prop.BA.sp) > 0.05) {
      #If PSME is dominant and dry site species are present
      x$Cover.dendro <- 'dry-PSME'
    } else {
      #If PSME is dominant and dry site species are absent
      x$Cover.dendro <- 'PSME'
    }
  } else if(sum(cov.BA[cov.BA$Species %in% dry.sp,]$Prop.BA.sp) > .50 &
            sum(cov.BA[cov.BA$Species %in% moist.sp,]$Prop.BA.sp) < .10) {
    #DMC if ponderosa pine not present, dry site species are dominant but may be micov.BAed with some other species (e.g. POTR, LAOC, PICO), and few moist site species are present in significant numbers
    x$Cover.dendro <- 'DMC-PSME'
  } else if(sum(cov.BA[cov.BA$Species == 'PICO',]$Prop.BA.sp) >= .50) {
    #PICO if PICO dominates stand
    x$Cover.dendro <- 'PICO'
  } else if (sum(cov.BA[cov.BA$Species %in% c('PIEN', 'PIEN/PIGL', 'PIGL', 'ABLA'),]$Prop.BA.sp) > .50) {
    #lowland PIEN if dominated by spruce
    x$Cover.dendro <- 'PIEN'
  } else if (sum(cov.BA[cov.BA$Species %in% c('POTR', 'POTR5', 'POBA', 'BEPA'),]$Prop.BA.sp) > .25) {
    if (sum(cov.BA[cov.BA$Species %in% c('POTR', 'POTR5', 'POBA', 'BEPA'),]$Prop.BA.sp) >= pure.cutoff) {
      x$Cover.dendro <- 'Broadleaf'
    } else {
      x$Cover.dendro <- 'Mixedwood'
    }
  } else {
    x$Cover.dendro <- 'MMC'
  }

  #Convert cover types to factor
  # x$Cover.dendro <- as.factor(x$Cover.dendro)
  return(x)
}
