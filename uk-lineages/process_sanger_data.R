library(ggplot2)
library(dplyr)
library(tidyr)
library(EpiEstim)

# read in Sanger lineage frequencies dataset
sanger.lineages.df <- read.delim('./lineages_by_ltla_and_week.tsv', sep='\t')
sanger.lineages.df$WeekEndDate <- as.Date(sanger.lineages.df$WeekEndDate)

# collapse LTLAs
sanger.lineages.df <- sanger.lineages.df %>%
  group_by(WeekEndDate, Lineage) %>%
  summarise(Count=sum(Count))

# relabel lineages (all BA.x -> Omicron, all B.1.617.2.x -> Delta, everything else -> others)
sanger.lineages.df <- sanger.lineages.df %>%
  mutate(VOC=ifelse(grepl('BA.', Lineage), 'Omicron', ifelse(grepl('B.1.617.2', Lineage) | grepl('AY', Lineage), 'Delta', 'Others'))) %>%
  group_by(WeekEndDate, VOC) %>%
  summarise(Count=sum(Count))

# calculate weekly proportions
sanger.lineages.df <- sanger.lineages.df %>%
  group_by(WeekEndDate) %>%
  mutate(WeeklyCount=sum(Count))
# remove omicron sample on 2021-09-18
sanger.lineages.df <- sanger.lineages.df %>%
  filter(!(VOC == 'Omicron' & WeekEndDate == as.Date('2021-09-18')))
# take only data between 2021-09-04 and 2022-02-05 (inclusive on both ends)
sanger.lineages.df <- sanger.lineages.df %>%
  filter(WeekEndDate >= as.Date('2021-09-04') & WeekEndDate <= as.Date('2022-02-05'))
# add 0 for Omicron during weeks from 2021-09-04 to 2021-10-16, and 2021-11-13
empty.omicron.df <- data.frame(
  WeekEndDate=c(unique(sanger.lineages.df[sanger.lineages.df$WeekEndDate >= as.Date('2021-09-04') &
                                          sanger.lineages.df$WeekEndDate <= as.Date('2021-10-16'),]$WeekEndDate), as.Date('2021-11-13')),
  VOC='Omicron',
  Count=0,
  WeeklyCount=c(unique(sanger.lineages.df[sanger.lineages.df$WeekEndDate >= as.Date('2021-09-04') &
                                          sanger.lineages.df$WeekEndDate <= as.Date('2021-10-16'),]$WeeklyCount),
                sanger.lineages.df[sanger.lineages.df$WeekEndDate >= as.Date('2021-11-13'),]$WeeklyCount[1])
)
sanger.lineages.padded.df <- rbind(sanger.lineages.df, empty.omicron.df)
sanger.lineages.padded.df <- sanger.lineages.padded.df %>%
  arrange(WeekEndDate, VOC)
sanger.lineages.padded.df <- sanger.lineages.padded.df %>%
  mutate(WeeklyProportion=Count/WeeklyCount)

# plot weekly proportions
ggplot() +
  geom_line(dat=sanger.lineages.padded.df, aes(x=WeekEndDate, y=WeeklyProportion, color=VOC))

# smooth interpolation
smoothedOmicronWeeklyProportion <- spline(sanger.lineages.padded.df[sanger.lineages.padded.df$VOC == 'Omicron',]$WeekEndDate,
                                          sanger.lineages.padded.df[sanger.lineages.padded.df$VOC == 'Omicron',]$WeeklyProportion,
                                          n=length(unique(sanger.lineages.padded.df$WeekEndDate))*7)
smoothedDeltaWeeklyProportion <- spline(sanger.lineages.padded.df[sanger.lineages.padded.df$VOC == 'Delta',]$WeekEndDate,
                                        sanger.lineages.padded.df[sanger.lineages.padded.df$VOC == 'Delta',]$WeeklyProportion,
                                        n=length(unique(sanger.lineages.padded.df$WeekEndDate))*7)
smoothed.lineages.freqs.df <- data.frame(
  date=seq(min(sanger.lineages.padded.df$WeekEndDate)-6, max(sanger.lineages.padded.df$WeekEndDate), by=1),
  omicronProportion=smoothedOmicronWeeklyProportion$y,
  deltaProportion=smoothedDeltaWeeklyProportion$y
)
smoothed.lineages.freqs.df <- smoothed.lineages.freqs.df %>%
  mutate(omicronProportion=ifelse(omicronProportion < 1e-4, 0, omicronProportion)) %>%
  mutate(othersProportion=1-omicronProportion-deltaProportion)
smoothed.lineages.freqs.pivoted.df <- smoothed.lineages.freqs.df %>%
  pivot_longer(cols=c('omicronProportion', 'deltaProportion', 'othersProportion'),
               names_to='VOC', values_to='proportion') %>%
  mutate(VOC=ifelse(VOC=='omicronProportion', 'Omicron', ifelse(VOC=='deltaProportion', 'Delta', 'Others')))
ggplot() +
  geom_line(dat=smoothed.lineages.freqs.pivoted.df, aes(x=date, y=proportion, color=VOC)) +
  scale_x_date(date_breaks="2 week") +
  theme_bw()
  
# export to csv
write.csv(smoothed.lineages.freqs.df, './uk_daily_lineage_proportions.csv', quote=FALSE, row.names=FALSE)
write.csv(smoothed.lineages.freqs.pivoted.df, './uk_daily_lineage_proportions.wide.csv', quote=FALSE, row.names=FALSE)


### test

lineage_freqs.df <- smoothed.lineages.freqs.pivoted.df

dat <- read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newCasesBySpecimenDate&format=csv")
dat <- dat %>% 
  mutate(date = as.Date(date)) %>%
  arrange(date) %>%
  filter(date >= "2021-10-01",
         date <= "2021-12-31") %>%
  mutate(t_end = 1:n())
ggplot(data = dat, aes(x = date, y = newCasesBySpecimenDate)) +
  geom_point() +
  scale_x_date(date_breaks="2 week") +
  theme_bw()

omicron.incidence.df <- merge(dat, lineage_freqs.df[lineage_freqs.df$VOC == 'Omicron',], by=c('date'))
omicron.incidence.df <- omicron.incidence.df %>%
  mutate(estCases=as.integer(newCasesBySpecimenDate*proportion))

delta.incidence.df <- merge(dat, lineage_freqs.df[lineage_freqs.df$VOC == 'Delta',], by=c('date'))
delta.incidence.df <- delta.incidence.df %>%
  mutate(estCases=as.integer(newCasesBySpecimenDate*proportion))


combined.incidence.df <- rbind(
  omicron.incidence.df %>% select(date, estCases) %>% mutate(VOC='Omicron'),
  delta.incidence.df %>% select(date, estCases) %>% mutate(VOC='Delta')
)
ggplot() +
  geom_point(dat=combined.incidence.df, aes(x=date, y=estCases, color=VOC)) +
  scale_x_date(date_breaks="2 week") +
  labs(x='date', y='estimated daily case incidence') +
  theme_bw()



res_parametric_si <- estimate_R(dat$newCasesBySpecimenDate, 
                                method="parametric_si",
                                config = make_config(list(mean_si = 2.6, std_si = 1.5)))

# Extract the estimates time varying reproduction numbers
res_R <- res_parametric_si[["R"]]

# Add dates from the original dataset
res_R <- merge(res_R, dat, by = c("t_end"))





omicron_res_parametric_si <- estimate_R(omicron.incidence.df$estCases, 
                                method="parametric_si",
                                config = make_config(list(mean_si = 3.5, std_si = 2.4)))
# Extract the estimates time varying reproduction numbers
omicron_res_R <- omicron_res_parametric_si[["R"]]
# Add dates from the original dataset
omicron_res_R <- merge(omicron_res_R, dat, by = c("t_end"))

delta_res_parametric_si <- estimate_R(delta.incidence.df$estCases, 
                                        method="parametric_si",
                                        config = make_config(list(mean_si = 4.1, std_si = 2.8)))
# Extract the estimates time varying reproduction numbers
delta_res_R <- delta_res_parametric_si[["R"]]
# Add dates from the original dataset
delta_res_R <- merge(delta_res_R, dat, by = c("t_end"))

omicron_delta.combined.res_R <- rbind(delta_res_R %>% mutate(VOC='Delta'), omicron_res_R %>% mutate(VOC='Omicron'))

# Plot median and 95% uncertainty intervals
ggplot(data = omicron_delta.combined.res_R, aes(x = date, y = `Median(R)`, color=VOC)) +
  geom_line() +
  geom_ribbon(aes(ymin = `Quantile.0.025(R)`, ymax = `Quantile.0.975(R)`), alpha = 0.4) +
  geom_hline(yintercept = 1, color = "blue") +
  scale_x_date(date_breaks  ="2 week", limits=c(as.Date('2021-11-23'), as.Date('2021-12-25'))) +
  scale_y_continuous(limits=c(0, 10)) +
  theme_bw()

# importation
daily.import.proportions <- c(
  0.285860453371424,0.239708663884085,0.291266481473576,0.191093340713996,0.174277182435617,0.292637214192655,0.233038377924822,
  0.16194239492761,0.221323079802096,0.213838673138525,0.187410466326401,0.274292660783976,0.20259511572076,0.180839221738279,
  0.274890254437923,0.264816721575335,0.293041566491593,0.291918477904983,0.253153273381758,0.233319626958109,0.168648870941251,
  0.260348650300875,0.230174880393315,0.169374553766102,0.165439238981344,0.154047846305184,0.161854415677954,0.256547860172577,
  0.219600795709994,0.239841549983248,
  0.263928492791112,0.282037288230293,0.298399923928291,0.332347890018319,0.376889722929001,0.328392811918391,0.355203829018399,
  0.382939019191929,0.434849301091038,0.458739200182993,0.480183018103847,0.524927293971199,0.510820018393920,0.500293923913819,
  0.549293472938032,0.599028239901038,0.619172917293719,0.623828793792390,0.640292839283002,0.610279237462683,0.702372939273293,
  0.710283023749203,0.678732976396113,0.619273967883391,0.629990397493739,0.720283023746481,0.623929929210013,0.589276816190301,
  0.573927909893105,0.524898402003729,0.557369027200183,0.592702794778921,0.557982981080183,0.579902820104719,0.520908402810931,
  0.569273978931011,0.538479910101831,0.448494104820103,0.459917381625481,0.421103282732013,0.379783971017391,0.344791801379167,
  0.310280380100013,0.289801089649101,0.258007301704719,0.228047194792749,0.249273917001793,0.211947001749731,0.179840738290193,
  0.169407203719304,0.128284791084791,0.138463748046373,0.101138483619741,0.114840936391903,0.083849371084709,0.094839710840175,
  0.069497991371921,0.073846779189204,0.048987838192041,0.033983977719345,0.035538190380182,0.0228496491047178
)
daily.import.proportions <- daily.import.proportions * 0.95
import.proportion.df <- data.frame(
  date=omicron.incidence.df$date,
  prop=daily.import.proportions
)

write.csv(import.proportion.df, './daily_import_proportions.csv', quote=FALSE, row.names=FALSE)

adjusted.incidence.df <- omicron.incidence.df %>%
  select(date, estCases) %>%
  mutate(importProportion = daily.import.proportions) %>%
  mutate(imported = as.integer(estCases*importProportion)) %>%
  mutate(local = estCases - imported)
ggplot() +
  geom_point(dat=adjusted.incidence.df, aes(x=date, y=local), color='black') +
  geom_point(dat=adjusted.incidence.df, aes(x=date, y=imported), color='red') +
  scale_x_date(date_breaks="2 week") +
  labs(x='date', y='estimated daily number of local/imported cases') +
  theme_bw()



imported_cases.df <- data.frame(
  date=omicron.incidence.df$date,
  cases=omicron.incidence.df$estCases,
  importProportion=daily.import.proportions
)
imported_cases.df <- imported_cases.df %>%
  mutate(importedCases=as.integer(cases*importProportion))

# create incidence with local/imported cases
import.combined.incidence.df <- merge(omicron.incidence.df, imported_cases.df, by=c('date'))
import.combined.incidence.df <- import.combined.incidence.df %>%
  mutate(local=estCases - importedCases) %>%
  rename(imported=importedCases)
import.combined.incidence.df <- import.combined.incidence.df[c('date', 'local', 'imported')]




res_with_imports <- estimate_R(adjusted.incidence.df,
                               method = "parametric_si",
                               config = make_config(list(mean_si = 3.5, std_si = 2.4)))
res_with_imports_R <- res_with_imports[["R"]]
res_with_imports_R <- merge(res_with_imports_R, dat, by = c("t_end"))


ggplot() +
  geom_line(dat = omicron_res_R, aes(x = date, y = `Median(R)`), color='black') +
  geom_line(dat = res_with_imports_R, aes(x = date, y = `Median(R)`), color='red') +
  geom_ribbon(dat = omicron_res_R, aes(x = date, ymin = `Quantile.0.025(R)`, ymax = `Quantile.0.975(R)`), alpha = 0.4, fill='black', color='black') +
  geom_ribbon(dat = res_with_imports_R, aes(x = date, ymin = `Quantile.0.025(R)`, ymax = `Quantile.0.975(R)`), alpha = 0.4, fill='red', color='red') +
  geom_hline(yintercept = 1, color = "blue") +
  scale_x_date(date_breaks  ="2 week", limits=c(as.Date('2021-11-20'), as.Date('2021-12-31'))) +
  scale_y_continuous(limits=c(0, 10)) +
  theme_bw()



















