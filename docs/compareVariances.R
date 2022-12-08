E1.true_sd_diff <- E1.RT_by_resp$respFALSE%>%sd()-E1.RT_by_resp$respTRUE%>%sd()

E1.sd_null = c();

for (i in seq(100)) {

  shuffled_sd_diff <- E1.df %>%
    filter((test_part=='test1' | test_part=='test2') & RT>100) %>%
    transform(subj_id=sample(subj_id))%>%
    group_by(subj_id,resp) %>%
    summarise(RT=median(RT))%>%
    group_by(resp)%>%
    summarise(sd=sd(RT))%>%
    spread(resp,sd,sep='')%>%
    mutate(diff=respFALSE-respTRUE)%>%
    pull(diff)

  E1.sd_null = c(E1.sd_null,shuffled_sd_diff)
}

E2.true_sd_diff <- E2.RT_by_resp$respFALSE%>%sd()-E2.RT_by_resp$respTRUE%>%sd()

E2.sd_null = c();

for (i in seq(100)) {

  shuffled_sd_diff <- E2.df %>%
    filter((test_part=='test1' | test_part=='test2') & RT>100) %>%
    transform(subj_id=sample(subj_id))%>%
    group_by(subj_id,resp) %>%
    summarise(RT=median(RT))%>%
    group_by(resp)%>%
    summarise(sd=sd(RT))%>%
    spread(resp,sd,sep='')%>%
    mutate(diff=respFALSE-respTRUE)%>%
    pull(diff)

  E2.sd_null = c(E2.sd_null,shuffled_sd_diff)
}
