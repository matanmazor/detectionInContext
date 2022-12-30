Ntrials <- 3000;

mu_absence_cong <- 0;
mu_absence_incong <- 0;

mu_presence_cong <- 0.5;
mu_presence_incong <- 0.3;

sd_presence <- 1.2;
sd_absence <- 1;

upper_boundary <- function(t) {
  # return(20-20/(1+exp(-0.001*(t-60))))
  return(10)
}

lower_boundary <- function(t) {
  # return(-13+20/(1+exp(-0.05*(t-60))))
  return(-10+t*0.2)
}

simulation.df <- data.frame(trial_number=integer(),target=character(),
                            context=character(), condition=character(),
                            t=integer(), x=double(), ub = double(), lb = double())


for (i in seq(Ntrials)) {
  t = 0;
  target = sample(c('absence','presence'),1);
  context = sample(c('cong','incong'),1);
  drift_rate = ifelse(target=='absence' & context=='cong', mu_absence_cong,
                      ifelse(target=='absence' & context=='incong', mu_absence_incong,
                        ifelse(target=='presence' & context=='cong', mu_presence_cong,
                             mu_presence_incong)))
  drift_sd = ifelse(target=='absence',sd_absence,sd_presence)
  x = 0
  simulation.df[nrow(simulation.df)+1,] = list(i,target,context,paste(target,context,sep=''),t,x,upper_boundary(t),lower_boundary(t))
  while (x<upper_boundary(t) & x>lower_boundary(t)) {
    t=t+1
    x=x+rnorm(n=1,mean=drift_rate,sd=drift_sd)
    simulation.df[nrow(simulation.df)+1,] = list(i,target,context,paste(target,context,sep=''),t,x,upper_boundary(t),lower_boundary(t))
  }
}

trial.df = simulation.df %>%
    group_by(trial_number) %>%
    filter(t==max(t)) %>%
    mutate(response=factor(ifelse(x>ub,1,0),levels=c(1,0)),
           presence=factor(ifelse(target=='presence',1,0),levels=c(1,0)),
           correct=response==presence)

p <- trial.df %>%
  group_by(response,context) %>%
  summarize(t=mean(t)) %>%
  ggplot(aes(x=context,fill=factor(response),y=t)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#377eb8","#e41a1c")) +
  coord_flip()+
  theme_classic()
ggsave('figures/sim_t_by_cong_resp.png',p,width=6,height=2)


p <- trial.df %>%
  group_by(presence,context) %>%
  summarize(rate=mean(response==1)) %>%
  ggplot(aes(x=context,fill=factor(presence),y=rate)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#377eb8","#e41a1c")) +
  coord_flip()+
  theme_classic()
ggsave('figures/sim_resp_by_cong_presence.png',p,width=6,height=2)
