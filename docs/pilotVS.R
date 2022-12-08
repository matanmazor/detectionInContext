library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)


N_perm <- 1000;
bootstrap_error <- function(x, N_perm) {
  N <- length(x)
  medians = c();
  for (i in 1:N_perm) {
    medians = c(medians,sample(x,replace=TRUE,size=N)%>%median())
  };
  return(sd(medians))
}

FlippedZsPilot.df <- read.csv('..\\experiments\\pilots\\flippedZsSameDistractors\\data\\jatos_results_batch1.txt',na.strings=c(""," ","NA")) %>%
  filter(trial_type=='p5vs_yn' & test_part=='main') %>%
  mutate(subj_id=subject_identifier,
         RT=as.numeric(RT),
         set_size=as.numeric(set_size),
         correct = correct=='true',
         target_present = target_present=='true',
         block=as.numeric(block))%>%
  dplyr::select(subj_id,flipped_first,RT,correct,set_size,target_present,search_type, block)

FlippedZsPilot.exclude <- FlippedZsPilot.df %>%
  group_by(subj_id,search_type)%>%
  summarise(acc=mean(correct))%>%
  group_by(subj_id)%>%
  summarise(acc=min(acc))%>%
  filter(acc<0.8)%>%
  pull(subj_id)%>%
  unique()

FlippedZsPilot.RTplot <- FlippedZsPilot.df %>%
  filter(correct & !subj_id %in% FlippedZsPilot.exclude) %>%
  group_by(target_present, 
           search_type,set_size, 
           subj_id) %>%
  summarise(RT=median(RT)) %>%
  group_by(target_present,
           search_type,
           set_size)%>%
  summarise(sem_RT=bootstrap_error(RT,N_perm),
            median_RT=median(RT))%>%
  mutate(target_present=factor(ifelse(target_present,'present','absent'),
                               levels=c('present','absent')))%>%
  ggplot(aes(x=set_size,y=median_RT,color=search_type,fill=search_type)) +
  geom_line(size=1) +
  # geom_point(shape=21,size=4)+
  facet_grid(cols=vars(target_present))+
  scale_fill_manual(values = c("black","#377eb8"))+
  scale_color_manual(values = c("black","#377eb8"))+
  scale_linetype_manual(values=c("solid","21"))+
  geom_errorbar(aes(ymin=median_RT-sem_RT,ymax=median_RT+sem_RT),linetype="solid", width=1.2, color="black") +
  # facet_grid(cols = vars(test_part),
  #            labeller = labeller(test_part = block_names))+
  labs(x='set size',y='median RT (ms)', title='Pilot: all trials') +
  theme_bw()+
  scale_x_continuous(limits=c(2,13),breaks = c(3,12))+
  theme(legend.position='none',
        legend.background = element_rect(fill=NA))+
  guides(color = FALSE, linetype=FALSE)

ggsave('figures/FlippedZsPilot_all_trials.png',FlippedZsPilot.RTplot,width=5,height=3)

FlippedZsPilot.slopes <- FlippedZsPilot.df %>%
  group_by(subj_id,search_type, target_present) %>%
  do(model=lm(RT~set_size,data=.)) %>%
  mutate(tidys=list(broom::tidy(model))) %>%
  unnest(tidys) %>%
  filter(term=='set_size') %>%
  dplyr::select(subj_id,search_type,target_present,estimate) %>%
  spread(search_type,estimate) %>%
  mutate(targetEffect = searchZ-searchFlippedZ)


FlippedZsPilot2.df <- read.csv('..\\experiments\\pilots\\flippedZsSameDistractors\\data\\jatos_results_batch2.txt',na.strings=c(""," ","NA")) %>%
  filter(trial_type=='p5vs_yn' & test_part=='main') %>%
  mutate(subj_id=subject_identifier,
         RT=as.numeric(RT),
         set_size=as.numeric(set_size),
         correct = correct=='true',
         target_present = target_present=='true',
         block=as.numeric(block))%>%
  dplyr::select(subj_id,flipped_first,RT,correct,set_size,target_present,search_type, block)

FlippedZsPilot2.exclude <- FlippedZsPilot2.df %>%
  group_by(subj_id,search_type)%>%
  summarise(acc=mean(correct))%>%
  group_by(subj_id)%>%
  summarise(acc=min(acc))%>%
  filter(acc<0.8)%>%
  pull(subj_id)%>%
  unique()

FlippedZsPilot2.RTplot <- FlippedZsPilot2.df %>%
  filter(correct & !subj_id %in% FlippedZsPilot.exclude) %>%
  group_by(target_present, 
           search_type,set_size, 
           subj_id) %>%
  summarise(RT=median(RT)) %>%
  group_by(target_present,
           search_type,
           set_size)%>%
  summarise(sem_RT=bootstrap_error(RT,N_perm),
            median_RT=median(RT))%>%
  mutate(target_present=factor(ifelse(target_present,'present','absent'),
                               levels=c('present','absent')))%>%
  ggplot(aes(x=set_size,y=median_RT,color=search_type,fill=search_type)) +
  geom_line(size=1) +
  # geom_point(shape=21,size=4)+
  facet_grid(cols=vars(target_present))+
  scale_fill_manual(values = c("black","#377eb8"))+
  scale_color_manual(values = c("black","#377eb8"))+
  scale_linetype_manual(values=c("solid","21"))+
  geom_errorbar(aes(ymin=median_RT-sem_RT,ymax=median_RT+sem_RT),linetype="solid", width=1.2, color="black") +
  # facet_grid(cols = vars(test_part),
  #            labeller = labeller(test_part = block_names))+
  labs(x='set size',y='median RT (ms)', title='Pilot: all trials') +
  theme_bw()+
  scale_x_continuous(limits=c(2,13),breaks = c(3,12))+
  theme(legend.position='none',
        legend.background = element_rect(fill=NA))+
  guides(color = FALSE, linetype=FALSE)

ggsave('figures/FlippedZsPilot_all_trials2.png',FlippedZsPilot2.RTplot,width=5,height=3)

FlippedZsPilot.slopes <- FlippedZsPilot.df %>%
  group_by(subj_id,search_type, target_present) %>%
  do(model=lm(RT~set_size,data=.)) %>%
  mutate(tidys=list(broom::tidy(model))) %>%
  unnest(tidys) %>%
  filter(term=='set_size') %>%
  dplyr::select(subj_id,search_type,target_present,estimate) %>%
  spread(search_type,estimate) %>%
  mutate(targetEffect = searchZ-searchFlippedZ)


# Evs3Pilot.df <- read.csv('..\\experiments\\pilots\\EsAndThrees\\data\\jatos_results_batch1.txt',na.strings=c(""," ","NA")) %>%
# Evs3Pilot.df <- read.csv('..\\experiments\\pilots\\EsAndThrees\\data\\jatos_results_longPause.txt',na.strings=c(""," ","NA")) %>%
Evs3Pilot.df <- read.csv('..\\experiments\\pilots\\EsAndThrees\\data\\jatos_results_symmetric.txt',na.strings=c(""," ","NA")) %>%
  filter(trial_type=='p5vs_yn' & test_part=='main') %>%
  mutate(subj_id=subject_identifier,
         RT=as.numeric(RT),
         set_size=as.numeric(set_size),
         correct = correct=='true',
         target_present = target_present=='true',
         block=as.numeric(block),
         search_type=ifelse(search_type=='threeInNumm','threeInNum',search_type))%>%
  dplyr::select(subj_id,three_first,RT,correct,set_size,target_present,search_type, block)


Evs3Pilot.exclude <- Evs3Pilot.df %>%
  group_by(subj_id,search_type)%>%
  summarise(acc=mean(correct))%>%
  group_by(subj_id)%>%
  summarise(acc=min(acc))%>%
  filter(acc<0.7)%>%
  pull(subj_id)%>%
  unique()

Evs3Pilot.RTplot <- Evs3Pilot.df %>%
  filter(correct & !subj_id %in% Evs3Pilot.exclude) %>%
  group_by(target_present, 
           search_type,set_size, 
           subj_id) %>%
  summarise(RT=mean(RT)) %>%
  group_by(target_present,
           search_type,
           set_size)%>%
  summarise(sem_RT=bootstrap_error(RT,N_perm),
            median_RT=mean(RT))%>%
  mutate(target=factor(substr(search_type,0,1),
                       levels=c('E','t'), labels=c('E','3')),
         distractors=factor(ifelse(search_type %in% c('EinLet', 'threeInLet'),
                                   'letters', 'numbers'),
                            levels=c('letters','numbers')),
         target_present=factor(ifelse(target_present,'present','absent'),
                               levels=c('present','absent'))
         )%>%
  ggplot(aes(x=set_size,y=median_RT,color=target,fill=target)) +
  geom_line(size=1) +
  geom_point(shape=21,size=4)+
  facet_grid(cols=vars(target_present), rows=vars(distractors))+
  scale_fill_manual(values = c("black","#377eb8"))+
  scale_color_manual(values = c("black","#377eb8"))+
  scale_linetype_manual(values=c("solid","21"))+
  geom_errorbar(aes(ymin=median_RT-sem_RT,ymax=median_RT+sem_RT),linetype="solid", width=1.2, color="black") +
  # facet_grid(cols = vars(test_part),
  #            labeller = labeller(test_part = block_names))+
  labs(x='set size',y='median RT (ms)', title='Pilot: all trials') +
  theme_bw()+
  scale_x_continuous(limits=c(2,13),breaks = c(3,12))+
  theme(legend.position='none',
        legend.background = element_rect(fill=NA))+
  guides(color = FALSE, linetype=FALSE)

ggsave('figures/Esand3sPilotSymmetric_all_trials.png',Evs3Pilot.RTplot,width=5,height=4.5)

Evs3Pilot.accPlot <- Evs3Pilot.df %>%
  filter(!subj_id %in% Evs3Pilot.exclude) %>%
  group_by(target_present, 
           search_type,set_size, 
           subj_id) %>%
  summarise(acc=mean(correct)) %>%
  group_by(target_present,
           search_type,
           set_size)%>%
  summarise(sd_acc=sd(acc),
            mean_acc=mean(acc))%>%
  mutate(target=factor(substr(search_type,0,1),
                       levels=c('E','t'), labels=c('E','3')),
         distractors=factor(ifelse(search_type %in% c('EinLet', 'threeInLet'),
                                   'letters', 'numbers'),
                            levels=c('letters','numbers')),
         target_present=factor(ifelse(target_present,'present','absent'),
                               levels=c('present','absent'))
  )%>%
  ggplot(aes(x=set_size,y=mean_acc,color=target,fill=target)) +
  geom_line(size=1) +
  geom_point(shape=21,size=4)+
  facet_grid(cols=vars(target_present), rows=vars(distractors))+
  scale_fill_manual(values = c("black","#377eb8"))+
  scale_color_manual(values = c("black","#377eb8"))+
  scale_linetype_manual(values=c("solid","21"))+
  # geom_errorbar(aes(ymin=median_RT-sem_RT,ymax=median_RT+sem_RT),linetype="solid", width=1.2, color="black") +
  # facet_grid(cols = vars(test_part),
  #            labeller = labeller(test_part = block_names))+
  labs(x='set size',y='mean accuracy', title='Pilot: all trials') +
  theme_bw()+
  scale_x_continuous(limits=c(2,13),breaks = c(3,12))+
  theme(legend.position='none',
        legend.background = element_rect(fill=NA))+
  guides(color = FALSE, linetype=FALSE)

ggsave('figures/Esand3sAccPilotSymmetric_all_trials.png',Evs3Pilot.accPlot,width=5,height=4.5)

Evs3Pilot.slopes <- Evs3Pilot.df %>%
  filter(!subj_id %in% Evs3Pilot.exclude) %>%
  group_by(subj_id,search_type, target_present) %>%
  do(model=lm(RT~set_size,data=.)) %>%
  mutate(tidys=list(broom::tidy(model))) %>%
  unnest(tidys) %>%
  filter(term=='set_size') %>%
  dplyr::select(subj_id,search_type,target_present,estimate) %>%
  spread(search_type,estimate) %>%
  mutate(targetEffectInLetters = EinLet-threeInLet,
         targetEffectInNumbers = EinNum-threeInNum)

Evs3Pilot.congSlopes <- Evs3Pilot.df %>%
  filter(!subj_id %in% Evs3Pilot.exclude) %>%
  mutate(cong=ifelse(search_type %in% c('EinLet','threeInNum'),'cong','incong')) %>%
  group_by(subj_id,cong, target_present) %>%
  do(model=lm(RT~set_size,data=.)) %>%
  mutate(tidys=list(broom::tidy(model))) %>%
  unnest(tidys) %>%
  filter(term=='set_size') %>%
  dplyr::select(subj_id,cong,target_present,estimate) %>%
  spread(cong,estimate) %>%
  mutate(diff=cong-incong)

exclude <- pilot.df %>%
  group_by(subj_id,search_type)%>%
  summarise(acc=mean(correct))%>%
  filter(acc<0.95)%>%
  pull(subj_id)%>%
  unique()

N_perm <- 1000;
bootstrap_error <- function(x, N_perm) {
  N <- length(x)
  medians = c();
  for (i in 1:N_perm) {
    medians = c(medians,sample(x,replace=TRUE,size=N)%>%median())
  };
  return(sd(medians))
}

pilot.median_search_times <- pilot.df %>%
  filter(correct & !(subj_id %in% exclude)) %>%
  group_by(subj_id,target_present,set_size,search_type) %>%
  summarise(RT=median(RT)) %>%
  group_by(target_present,search_type,set_size)%>%
  summarise(median_RT=median(RT),
            sem_RT = bootstrap_error(RT,N_perm))%>%
  mutate(target_present = factor(ifelse(target_present,'present','absent'), levels=c('present','absent')))

pilot.median_search_times_first_block <- pilot.df %>%
  filter(correct & (block==1 | block==4) & !(subj_id %in% exclude)) %>%
  group_by(subj_id,target_present,set_size,search_type) %>%
  summarise(RT=median(RT)) %>%
  group_by(target_present,search_type,set_size)%>%
  summarise(median_RT=median(RT),
            sem_RT = bootstrap_error(RT,N_perm))%>%
  mutate(target_present = factor(ifelse(target_present,'present','absent'), levels=c('present','absent')))

pilot.first_search_times <- pilot.df %>%
  filter(correct & !subj_id %in% exclude) %>%
  group_by(subj_id,target_present,set_size,search_type) %>%
  summarise(RT=RT[1]) %>%
  group_by(target_present,search_type,set_size)%>%
  summarise(median_RT=median(RT),
            sem_RT = bootstrap_error(RT,N_perm))%>%
  mutate(target_present = factor(ifelse(target_present,'present','absent'), levels=c('present','absent')))

pilot.mean_accuracy <- pilot.df %>%
  group_by(subj_id,target_present,set_size,search_type) %>%
  summarise(acc=mean(correct)) %>%
  group_by(target_present,search_type,set_size)%>%
  summarise(acc=mean(acc)) %>%
  mutate(target_present = factor(ifelse(target_present,'present','absent'), levels=c('present','absent')))

RTplot <- ggplot(data=pilot.median_search_times,
                 aes(x=set_size, y=median_RT, linetype=target_present, color=search_type)) +
  geom_line(size=1) +
  # geom_point(aes(shape = search_type), size=4, color="black",stroke=1.5, alpha=0.8) +
  # scale_shape_manual(values=c(4,21))+
  scale_fill_manual(values = c("black","#377eb8"))+
  scale_color_manual(values = c("black","#377eb8"))+
  scale_linetype_manual(values=c("solid","21"))+
  facet_grid(cols = vars(target_present))+
  # geom_errorbar(aes(ymin=median_RT-sem_RT,ymax=median_RT+sem_RT),linetype="solid", width=1.2, color="black") +
  # facet_grid(cols = vars(test_part),
  #            labeller = labeller(test_part = block_names))+
  labs(x='set size',y='median RT (ms)', title='Pilot: all trials') +
  theme_bw()+
  scale_x_continuous(limits=c(2,13),breaks = c(3,12))+
  theme(legend.position='none',
        legend.background = element_rect(fill=NA))+
  guides(color = FALSE, linetype=FALSE)

ggsave('figures/pilotVS_all_trials.png',RTplot,width=5,height=4.5)
ggsave('figures/pilotVS_all_trials.svg',RTplot,width=5,height=4.5)

FirsBlockRTplot <- ggplot(data=pilot.median_search_times_first_block,
                 aes(x=set_size, y=median_RT, linetype=target_present, color=search_type)) +
  geom_line(size=1) +
  # geom_point(aes(shape = search_type), size=4, color="black",stroke=1.5, alpha=0.8) +
  # scale_shape_manual(values=c(4,21))+
  scale_fill_manual(values = c("black","#377eb8"))+
  scale_color_manual(values = c("black","#377eb8"))+
  scale_linetype_manual(values=c("solid","21"))+
  facet_grid(cols = vars(target_present))+
  # geom_errorbar(aes(ymin=median_RT-sem_RT,ymax=median_RT+sem_RT),linetype="solid", width=1.2, color="black") +
  # facet_grid(cols = vars(test_part),
  #            labeller = labeller(test_part = block_names))+
  labs(x='set size',y='median RT (ms)', title='Pilot: first block only') +
  theme_bw()+
  scale_x_continuous(limits=c(2,7),breaks = c(3,6))+
  theme(legend.position='none',
        legend.background = element_rect(fill=NA))+
  guides(color = FALSE, linetype=FALSE)

ggsave('figures/pilotVS_first_block.png',FirsBlockRTplot,width=5,height=4.5)
ggsave('figures/pilotVS_first_block.svg',FirsBlockRTplot,width=5,height=4.5)

FirstRTplot <- ggplot(data=pilot.first_search_times,
                          aes(x=set_size, y=median_RT, linetype=target_present, color=search_type)) +
  geom_line(size=1) +
  # geom_point(aes(shape = search_type), size=4, color="black",stroke=1.5, alpha=0.8) +
  # scale_shape_manual(values=c(4,21))+
  scale_fill_manual(values = c("black","#377eb8"))+
  scale_color_manual(values = c("black","#377eb8"))+
  scale_linetype_manual(values=c("solid","21"))+
  facet_grid(cols = vars(target_present))+
  # geom_errorbar(aes(ymin=median_RT-sem_RT,ymax=median_RT+sem_RT),linetype="solid", width=1.2, color="black") +
  # facet_grid(cols = vars(test_part),
  #            labeller = labeller(test_part = block_names))+
  labs(x='set size',y='median RT (ms)', title='Pilot: first trials') +
  theme_bw()+
  scale_x_continuous(limits=c(2,7),breaks = c(3,6))+
  theme(legend.position='none',
        legend.background = element_rect(fill=NA))+
  guides(color = FALSE, linetype=FALSE)

ggsave('figures/pilotVS_first_trials.png',FirstRTplot,width=5,height=4.5)
ggsave('figures/pilotVS_first_trials.svg',FirstRTplot,width=5,height=4.5)



pilot2.df <- read.csv('..\\experiments\\pilots\\flippedZsFirstTrials\\data\\jatos_results_batch2.txt',na.strings=c(""," ","NA")) %>%
  filter(trial_type=='p5vs_yn' & test_part %in% c('absence1','mixed2')) %>%
  mutate(subj_id=subject_identifier,
         RT=as.numeric(RT),
         set_size=as.numeric(set_size),
         correct = correct=='true',
         target_present = target_present=='true',
         block=as.numeric(block))%>%
  dplyr::select(subj_id,flipped_first,RT,correct,set_size,target_present,search_type, block,test_part)

pilot2.summary_stats = pilot2.df %>%
  group_by(subj_id) %>%
  summarise(nerrors=sum(!correct),
            nerrors_first = sum(!correct & test_part=='absence1'),
            RT25=quantile(RT,0.25),
            RT75=quantile(RT,0.75))

pilot2.included = pilot2.summary_stats %>%
  filter(nerrors<5 &
           RT25>100 &
           RT75<5000) %>%
  pull(subj_id)

N_perm <- 1000;
bootstrap_error <- function(x, N_perm) {
  N <- length(x)
  medians = c();
  for (i in 1:N_perm) {
    medians = c(medians,sample(x,replace=TRUE,size=N)%>%median())
  };
  return(sd(medians))
}

pilot2.RT_by_position <- pilot2.df %>%
  group_by(subj_id)%>%
  mutate(i=seq_along(RT))%>%
  group_by(i) %>%
  summarise(RT=mean(RT))

pilot2.mean_RT <- pilot2.df$RT %>% mean()

pilot2.pp_df <- pilot2.df %>%
  group_by(subj_id) %>%
  mutate(i=seq_along(RT)) %>%
  rowwise() %>%
  mutate(meanRT_per_i = pilot2.RT_by_position$RT[pilot2.RT_by_position$i==i],
         meanRT_per_s = pilot2.RT_by_subj_id$RT[pilot2.RT_by_subj_id$subj_id==subj_id],
         RTcorrected = RT-meanRT_per_i+pilot2.mean_RT)

pilot2.median_search_times <- pilot2.pp_df %>%
  filter(correct & RT>100 & RT<5000 & subj_id %in% pilot2.included) %>%
  group_by(subj_id,target_present,set_size,search_type,test_part) %>%
  summarise(RT=median(RTcorrected)) %>%
  group_by(target_present,search_type,set_size,test_part)%>%
  summarise(median_RT=median(RT),
            sem_RT = bootstrap_error(RT,N_perm))%>%
  mutate(target_present = factor(ifelse(target_present,'present','absent'), levels=c('present','absent')))

pilot2.slopes_by_part <- pilot2.pp_df %>%
  filter(correct & RT>100 & RT<5000 & subj_id %in% pilot2.included) %>%
  group_by(subj_id,search_type,test_part, target_present) %>%
  do(model=lm(RTcorrected~set_size,data=.)) %>%
  mutate(tidys=list(broom::tidy(model))) %>%
  unnest(tidys) %>%
  # we are interested in the slope, i.e., the effect of set size.
  filter(term=='set_size') %>%
  dplyr::select(subj_id,search_type,test_part,target_present,estimate)%>%
  rename(slope=estimate)%>%
  mutate(target_present = factor(ifelse(target_present,'present','absent'), levels=c('present','absent')),
         search_type=ifelse(search_type=='searchZ','Z','flipZ'),
         test_part=ifelse(test_part=='absence1','a1','m2'),
         condition=paste(test_part,search_type,target_present,sep='_')) %>%
  dplyr::select(subj_id,condition,slope) %>%
  spread(condition,slope)


pilot2.pooled_slopes <- pilot2.pp_df %>%
  filter(correct & RT>100 & RT<5000 & subj_id %in% pilot2.included) %>%
  group_by(subj_id,search_type,target_present) %>%
  do(model=lm(RTcorrected~set_size,data=.)) %>%
  mutate(tidys=list(broom::tidy(model))) %>%
  unnest(tidys) %>%
  # we are interested in the slope, i.e., the effect of set size.
  filter(term=='set_size') %>%
  dplyr::select(subj_id,search_type,target_present,estimate)%>%
  rename(slope=estimate)%>%
  mutate(target_present = factor(ifelse(target_present,'present','absent'), levels=c('present','absent')),
         search_type=ifelse(search_type=='searchZ','Z','flipZ'),
         condition=paste('pooled',search_type,target_present,sep='_')) %>%
  dplyr::select(subj_id,condition,slope) %>%
  spread(condition,slope)

pilot2.last_slopes <- pilot2.pp_df %>%
  filter(correct & RT>100 & RT<5000 & subj_id %in% pilot2.included) %>%
  group_by(subj_id,search_type,target_present,set_size) %>%
  arrange(-i) %>%
  summarise(RTcorrected=RTcorrected[1]) %>%
  group_by(subj_id,search_type,target_present) %>%
  do(model=lm(RTcorrected~set_size,data=.)) %>%
  mutate(tidys=list(broom::tidy(model))) %>%
  unnest(tidys) %>%
  # we are interested in the slope, i.e., the effect of set size.
  filter(term=='set_size') %>%
  dplyr::select(subj_id,search_type,target_present,estimate)%>%
  rename(slope=estimate)%>%
  mutate(target_present = factor(ifelse(target_present,'present','absent'), levels=c('present','absent')),
         search_type=ifelse(search_type=='searchZ','Z','flipZ'),
         condition=paste('last',search_type,target_present,sep='_')) %>%
  dplyr::select(subj_id,condition,slope) %>%
  spread(condition,slope)

pilot2.first_slopes <- pilot2.pp_df %>%
  filter(correct & RT>100 & RT<5000 & subj_id %in% pilot2.included) %>%
  group_by(subj_id,search_type,target_present,set_size) %>%
  arrange(i) %>%
  summarise(RTcorrected=RTcorrected[1]) %>%
  group_by(subj_id,search_type,target_present) %>%
  do(model=lm(RTcorrected~set_size,data=.)) %>%
  mutate(tidys=list(broom::tidy(model))) %>%
  unnest(tidys) %>%
  # we are interested in the slope, i.e., the effect of set size.
  filter(term=='set_size') %>%
  dplyr::select(subj_id,search_type,target_present,estimate)%>%
  rename(slope=estimate)%>%
  mutate(target_present = factor(ifelse(target_present,'present','absent'), levels=c('present','absent')),
         search_type=ifelse(search_type=='searchZ','Z','flipZ'),
         condition=paste('first',search_type,target_present,sep='_')) %>%
  dplyr::select(subj_id,condition,slope) %>%
  spread(condition,slope)

pilot2.first_search_times <- pilot2.df %>%
  filter(correct & RT>100 & RT<5000 & subj_id %in% pilot2.included) %>%
  group_by(subj_id,target_present,set_size,search_type) %>%
  summarise(RT=RT[1]) %>%
  group_by(target_present,search_type,set_size)%>%
  summarise(median_RT=median(RT),
            sem_RT = bootstrap_error(RT,N_perm))%>%
  mutate(target_present = factor(ifelse(target_present,'present','absent'), levels=c('present','absent')))

pilot2.slopes <- pilot2.slopes_by_part %>%
  merge(pilot2.pooled_slopes) %>%
  merge(pilot2.last_slopes) %>%
  merge(pilot2.first_slopes) %>%
  mutate(H1 = pooled_Z_present-pooled_flipZ_present,
         H2 = pooled_Z_absent-pooled_flipZ_absent,
         H3 = a1_Z_absent-a1_flipZ_absent,
         H4 = last_Z_absent-last_flipZ_absent,
         H5 = H3-H4,
         H6 = H1-H2,
         H7a = first_Z_present-first_flipZ_present,
         H7b = first_Z_absent-first_flipZ_absent,
         H7 = H7a-H7b)


pilot2.mean_position <- pilot2.df %>%
  group_by(subj_id)%>%
  mutate(i=seq_along(RT))%>%
  group_by(target_present,set_size,search_type,test_part) %>%
  summarise(meani=mean(i))

pilot.first_search_times <- pilot.df %>%
  filter(correct) %>%
  group_by(subj_id,target_present,set_size,search_type) %>%
  summarise(RT=RT[1]) %>%
  group_by(target_present,search_type,set_size)%>%
  summarise(median_RT=median(RT),
            sem_RT = bootstrap_error(RT,N_perm))%>%
  mutate(target_present = factor(ifelse(target_present,'present','absent'), levels=c('present','absent')))

pilot.mean_accuracy <- pilot.df %>%
  group_by(subj_id,target_present,set_size,search_type) %>%
  summarise(acc=mean(correct)) %>%
  group_by(target_present,search_type,set_size)%>%
  summarise(acc=mean(acc)) %>%
  mutate(target_present = factor(ifelse(target_present,'present','absent'), levels=c('present','absent')))

RTplot <- ggplot(data=pilot.median_search_times,
                 aes(x=set_size, y=median_RT, linetype=target_present, color=search_type)) +
  geom_line(size=1) +
  # geom_point(aes(shape = search_type), size=4, color="black",stroke=1.5, alpha=0.8) +
  # scale_shape_manual(values=c(4,21))+
  scale_fill_manual(values = c("black","#377eb8"))+
  scale_color_manual(values = c("black","#377eb8"))+
  scale_linetype_manual(values=c("solid","21"))+
  facet_grid(cols = vars(target_present))+
  # geom_errorbar(aes(ymin=median_RT-sem_RT,ymax=median_RT+sem_RT),linetype="solid", width=1.2, color="black") +
  # facet_grid(cols = vars(test_part),
  #            labeller = labeller(test_part = block_names))+
  labs(x='set size',y='median RT (ms)', title='Pilot: all trials') +
  theme_bw()+
  scale_x_continuous(limits=c(2,7),breaks = c(3,6))+
  theme(legend.position='none',
        legend.background = element_rect(fill=NA))+
  guides(color = FALSE, linetype=FALSE)

ggsave('figures/pilotVS_all_trials.png',RTplot,width=5,height=4.5)
ggsave('figures/pilotVS_all_trials.svg',RTplot,width=5,height=4.5)

RTplot2 <- ggplot(data=pilot2.median_search_times %>% filter(test_part=='mixed2'),
                 aes(x=set_size, y=median_RT, linetype=target_present, color=search_type)) +
  geom_line(size=1) +
  # geom_point(aes(shape = search_type), size=4, color="black",stroke=1.5, alpha=0.8) +
  # scale_shape_manual(values=c(4,21))+
  scale_fill_manual(values = c("black","#377eb8"))+
  scale_color_manual(values = c("black","#377eb8"))+
  scale_linetype_manual(values=c("solid","21"))+
  facet_grid(cols = vars(target_present))+
  # geom_errorbar(aes(ymin=median_RT-sem_RT,ymax=median_RT+sem_RT),linetype="solid", width=1.2, color="black") +
  # facet_grid(cols = vars(test_part),
  #            labeller = labeller(test_part = block_names))+
  labs(x='set size',y='median RT (ms)', title='Pilot: all trials') +
  theme_bw()+
  scale_x_continuous(limits=c(2,7),breaks = c(3,6))+
  theme(legend.position='none',
        legend.background = element_rect(fill=NA))+
  guides(color = FALSE, linetype=FALSE)

ggsave('figures/pilotVS2_all_trials.png',RTplot2,width=5,height=4.5)
ggsave('figures/pilotVS2_all_trials.svg',RTplot2,width=5,height=4.5)

FirsBlockRTplot <- ggplot(data=pilot.median_search_times_first_block,
                          aes(x=set_size, y=median_RT, linetype=target_present, color=search_type)) +
  geom_line(size=1) +
  # geom_point(aes(shape = search_type), size=4, color="black",stroke=1.5, alpha=0.8) +
  # scale_shape_manual(values=c(4,21))+
  scale_fill_manual(values = c("black","#377eb8"))+
  scale_color_manual(values = c("black","#377eb8"))+
  scale_linetype_manual(values=c("solid","21"))+
  facet_grid(cols = vars(target_present))+
  # geom_errorbar(aes(ymin=median_RT-sem_RT,ymax=median_RT+sem_RT),linetype="solid", width=1.2, color="black") +
  # facet_grid(cols = vars(test_part),
  #            labeller = labeller(test_part = block_names))+
  labs(x='set size',y='median RT (ms)', title='Pilot: first block only') +
  theme_bw()+
  scale_x_continuous(limits=c(2,7),breaks = c(3,6))+
  theme(legend.position='none',
        legend.background = element_rect(fill=NA))+
  guides(color = FALSE, linetype=FALSE)

ggsave('figures/pilotVS_first_block.png',FirsBlockRTplot,width=5,height=4.5)
ggsave('figures/pilotVS_first_block.svg',FirsBlockRTplot,width=5,height=4.5)

FirstRTplot <- ggplot(data=pilot.first_search_times,
                      aes(x=set_size, y=median_RT, linetype=target_present, color=search_type)) +
  geom_line(size=1) +
  # geom_point(aes(shape = search_type), size=4, color="black",stroke=1.5, alpha=0.8) +
  # scale_shape_manual(values=c(4,21))+
  scale_fill_manual(values = c("black","#377eb8"))+
  scale_color_manual(values = c("black","#377eb8"))+
  scale_linetype_manual(values=c("solid","21"))+
  facet_grid(cols = vars(target_present))+
  # geom_errorbar(aes(ymin=median_RT-sem_RT,ymax=median_RT+sem_RT),linetype="solid", width=1.2, color="black") +
  # facet_grid(cols = vars(test_part),
  #            labeller = labeller(test_part = block_names))+
  labs(x='set size',y='median RT (ms)', title='Pilot: first trials') +
  theme_bw()+
  scale_x_continuous(limits=c(2,7),breaks = c(3,6))+
  theme(legend.position='none',
        legend.background = element_rect(fill=NA))+
  guides(color = FALSE, linetype=FALSE)

ggsave('figures/pilotVS_first_trials.png',FirstRTplot,width=5,height=4.5)
ggsave('figures/pilotVS_first_trials.svg',FirstRTplot,width=5,height=4.5)

FirstRTplot2 <- ggplot(data=pilot2.first_search_times,
                      aes(x=set_size, y=median_RT, linetype=target_present, color=search_type)) +
  geom_line(size=1) +
  # geom_point(aes(shape = search_type), size=4, color="black",stroke=1.5, alpha=0.8) +
  # scale_shape_manual(values=c(4,21))+
  scale_fill_manual(values = c("black","#377eb8"))+
  scale_color_manual(values = c("black","#377eb8"))+
  scale_linetype_manual(values=c("solid","21"))+
  facet_grid(cols = vars(target_present))+
  # geom_errorbar(aes(ymin=median_RT-sem_RT,ymax=median_RT+sem_RT),linetype="solid", width=1.2, color="black") +
  # facet_grid(cols = vars(test_part),
  #            labeller = labeller(test_part = block_names))+
  labs(x='set size',y='median RT (ms)', title='Pilot: first trials') +
  theme_bw()+
  scale_x_continuous(limits=c(2,7),breaks = c(3,6))+
  theme(legend.position='none',
        legend.background = element_rect(fill=NA))+
  guides(color = FALSE, linetype=FALSE)

ggsave('figures/pilotVS2_first_trials.png',FirstRTplot2,width=5,height=4.5)
ggsave('figures/pilotVS2_first_trials.svg',FirstRTplot2,width=5,height=4.5)
