library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)
library(mgcv)
library(gratia)
library(RColorBrewer)
library(rstan)

# fit a generalized additive model (GAM) to the data
fitGam <- function(data, sbst, n)
{
  # initialize
  knts = list(tim=c(0,0,0,1,15,30,45,60,100,101,101,101))
  gamma = 0.5 # 1 - 0.5 
  max.deriv.se.sbst = 3 # something greater than 2
  fit = NULL
  # fitting
  while (max.deriv.se.sbst>2 & gamma<10) {
    gamma = gamma+0.5
    fit.check = gam(reg~s(tim, bs='bs', k=8), gamma=gamma, data=list(reg=data, tim=1:n), family=poisson(), knots=knts)
    max.deriv.se.sbst = max(derivatives(fit.check, n=n)$se[sbst])
    if (!is.null(fit) & max.deriv.se.sbst<0.01) break
    fit = fit.check
  }
  if (gamma==10 & max.deriv.se.sbst>2) fit = NULL
  return(fit)
}

# the primary load data function to run the model. Does all the neuron screening
LoadDataFiltered <- function(n=100, mltplyr=2)
{
  print(paste('RUNNING LoadAndScreenData FOR SESSION ID:', session_id))
  filter.strt = 30/mltplyr
  filter.end = 130/mltplyr
  sbst = filter.strt:filter.end
  data = list()
  configs.path = file.path(data.path, 'input', paste('units_data', session_id, sep='_'))
  configs = list.files(configs.path)
  for (config in configs)
  {
    data[[config]] = list()
    regions.path = file.path(configs.path, config)
    regions = intersect(c('V1','LM','RL','AL','AM','PM','TH'), list.files(regions.path))
    for (region in regions)
    {
      data[[config]][[region]] = list()
      units.path = file.path(regions.path, region)
      units = list.files(units.path)
      region.sum = 0
      region.data.for.fitlering = list()
      for (unit in units)
      {
        unit_key = strsplit(unit, ".csv", fixed = TRUE)[[1]][1]
        print(paste(config, region, unit_key))
        csv.data = read.csv(file = file.path(units.path, unit), nrows=n)
        unit.data = data.frame(matrix(vector('double'), n, 0))
        for (trial in names(csv.data))
        {
          trial.data = csv.data[[trial]]
          if (sum(trial.data[sbst])==0) next
          unit.data[[trial]] = trial.data
        }
        if (ncol(unit.data)==0) next
        region.data.for.fitlering[[unit_key]] = unit.data
      }
      if (length(region.data.for.fitlering)==0) {
        data[[config]][[region]] = NULL
        next
      }
      # Calculate the threshold for the top 60%
      thresholds.set = sapply(region.data.for.fitlering,
                              function(df) {
                                subset_df = df[sbst, ]
                                thissum = sum(as.matrix(subset_df))
                                unit.sum = rowSums(df)
                                fit = fitGam(unit.sum, sbst, n)
                                if (is.null(fit)) return(c(0,0,0,0))
                                else {
                                  rng = derivatives(fit, n=n)$derivative
                                  dff = sign(rng[(filter.strt+1):filter.end])-sign(rng[filter.strt:(filter.end-1)]) 
                                  crit.pt.indx = which(dff==-2)
                                  return(c(thissum, length(crit.pt.indx), max(rng[sbst]), max(fit$fitted.values[sbst])))
                                }
                              })
      cutoffs = apply(thresholds.set, 1, function(row) {
        quantile(row, 0.4)
      })
      cutoffs[2] = 0
      selected.df.cols = apply(thresholds.set, 2, function(col) {
        all(col > cutoffs)
      })
      selected.dfs = region.data.for.fitlering[selected.df.cols]
      for (unit_key in names(selected.dfs))
      {
        unit.data = selected.dfs[[unit_key]]
        if (ncol(unit.data)<10) next
        data[[config]][[region]][[unit_key]] = unit.data
        region.sum = region.sum+rowSums(unit.data)
      }
      if (length(data[[config]][[region]])<2) {
        data[[config]][[region]] = NULL
        next
      }
      fit = fitGam(region.sum, sbst, n)
      rng = derivatives(fit, n=n)$derivative
      dff = sign(rng[(filter.strt+1):filter.end])-sign(rng[filter.strt:(filter.end-1)]) 
      crit.pt.indx = which(dff==-2)
      if (length(crit.pt.indx)==0)  {
        data[[config]][[region]] = NULL
        next
      }
    }
    if (length(data[[config]])==0) data[[config]] = NULL
  }
  return(data)
}

# extract the population peak times for each trial in each region 
selectTrials <- function(data, n=100, mltplyr=2)
{
  print(paste('RUNNING selectTrials FOR SESSION ID:', session_id))
  filter.strt = 20/mltplyr
  filter.end = 130/mltplyr
  sbst = filter.strt:filter.end
  theNames = c('config','region','trial','num.units', 'trial.spikeCt', 'edf', 'deriv', 'num.crit.pts', 'trial.bl', 'fit.max', 'increase',
               'gain.time', 'reg.gain.time', 'relative.time', 'sd', 'num.trials', 'reg.spikeCt','reg.edf',
               'reg.deriv','reg.increase','reg.crit.pts')
  data.readout = data.frame(matrix(vector('double'), 0, length(theNames)))
  colnames(data.readout) = theNames
  ct = 1
  data.structure = list()
  loop.configs=names(data)
  inits = data.frame(matrix(vector('double'), 0, 0))
  membership.count = 1
  for (config in loop.configs) {
    loop.regions=names(data[[config]])
    count.reg = 0
    for (region in loop.regions) {
      loop.units=names(data[[config]][[region]])
      trials = list()
      for (unit in loop.units) {
        for (trl in names(data[[config]][[region]][[unit]])) {
          if (!(trl %in% names(trials))) trials[[trl]] = data.frame(matrix(vector(), n, 0))
          trials[[trl]][[unit]] = data[[config]][[region]][[unit]][[trl]][1:n]
        }
      }
      reg.dat = 0
      gain.time.vec = trial.sum.spikes.vec = trial.increase.vec = vector()
      for (trl in names(trials)) {
        print(paste(config,region,trl,sep='_'))
        num.units = ncol(trials[[trl]])
        if (num.units<5) next
        dat = rowSums(trials[[trl]])
        sum.spikes = sum(dat[sbst])
        if (sum.spikes<4) next
        fit = fitGam(dat, sbst, n)
        if (is.null(fit)) next
        rng = derivatives(fit, n=n)$derivative
        max.fd.rng = max(rng[sbst])
        if (max.fd.rng<=0.05) next
        edf = pen.edf(fit)
        if (edf<=2) next
        dff = sign(rng[(filter.strt+1):filter.end])-sign(rng[filter.strt:(filter.end-1)]) 
        crit.pt.indx = as.integer(names(which(dff==-2)))
        if (length(crit.pt.indx)==0) next
        if (length(crit.pt.indx)==1)
          max.fit.peak.indx = crit.pt.indx
        else {
          pos.loc = which(fit$fitted.values[crit.pt.indx]==max(fit$fitted.values[crit.pt.indx]))[1]
          if (pos.loc>1) {
            if ((fit$fitted.values[crit.pt.indx[pos.loc-1]]/fit$fitted.values[crit.pt.indx[pos.loc]])>=0.7) {
              max.fit.peak.indx = crit.pt.indx[pos.loc-1]
            }
            else
              max.fit.peak.indx = crit.pt.indx[pos.loc]
          }
          else
            max.fit.peak.indx = crit.pt.indx[pos.loc]
        }
        fit.bl = mean(fit$fitted.values[1:filter.strt])
        fit.max = fit$fitted.values[max.fit.peak.indx]
        increase = fit.max-fit.bl
        if (increase<=0.1) next
        gain.time = mltplyr*max.fit.peak.indx
        gain.time.vec[paste(config,trl,sep='_')] = gain.time
        trial.sum.spikes.vec[paste(config,trl,sep='_')] = sum.spikes 
        trial.increase.vec[paste(config,trl,sep='_')] = increase
        reg.dat = reg.dat + dat
        
        data.readout[paste(config,region,trl,sep='_'),] = c(config, region, trl, num.units, sum.spikes, 
                                                            edf, max.fd.rng, length(crit.pt.indx), fit.bl, 
                                                            fit.max, increase, gain.time, rep(0,length(theNames)-12))
      }
      if (length(reg.dat)==1) next
      sub.ct = length(gain.time.vec)
      reg.fit = fitGam(reg.dat, sbst, n)
      reg.sum.spikes = sum(reg.dat[sbst])
      reg.edf = pen.edf(reg.fit)
      rng = derivatives(reg.fit, n=n)$derivative
      reg.deriv = max(rng[sbst])
      dff = sign(rng[(filter.strt+1):filter.end])-sign(rng[filter.strt:(filter.end-1)]) 
      crit.pt.indx = as.integer(names(which(dff==-2)))
      if (length(crit.pt.indx)==0)
        max.fit.peak.indx = filter.end
      else if (length(crit.pt.indx)==1)
        max.fit.peak.indx = crit.pt.indx
      else
        max.fit.peak.indx = which(reg.fit$fitted.values==max(reg.fit$fitted.values[crit.pt.indx]))
      
      reg.fit.bl = mean(reg.fit$fitted.values[1:filter.strt])
      reg.fit.max = max(reg.fit$fitted.values[sbst])
      reg.increase = reg.fit.max-reg.fit.bl
      reg.gain.time = mltplyr*max.fit.peak.indx
      rl.tm = gain.time.vec-reg.gain.time
      sd = sqrt(sum((gain.time.vec-reg.gain.time)^2)/(length(gain.time.vec)-1))
      if (length(gain.time.vec)==1) sd = 1
      data.readout[ct:(ct+sub.ct-1),c('relative.time','reg.gain.time','sd','num.trials','reg.spikeCt',
                                      'reg.edf','reg.deriv','reg.increase','reg.crit.pts')] =
        c(rl.tm, rep(c(reg.gain.time,sd,sub.ct,reg.sum.spikes,reg.edf,reg.deriv,
                       reg.increase,length(crit.pt.indx)),each=sub.ct))
      
      ct = ct+sub.ct
      if (reg.edf<=2) next
      if (reg.deriv<=0.05) next
      if (length(crit.pt.indx)==0) next
      if (reg.increase<3) next
      if (reg.sum.spikes<100) next
      if (!(region %in% names(data.structure))) {
        cnms = c('trial.tm', 'trial.spikeCt', 'trial.increase', 'reg.tm', 'reg.spikeCt', 'reg.increase', 'membership')
        data.structure[[region]] = data.frame(matrix(vector('double'), 0, length(cnms)))
        colnames(data.structure[[region]]) = cnms
      }
      filtr = abs(gain.time.vec-reg.gain.time)>ceiling(1.6*sd)
      trial.increase.vec[filtr] = (trial.increase.vec/(1+0.1*abs(gain.time.vec-reg.gain.time)))[filtr]
      data.structure[[region]] = bind_rows(data.structure[[region]],
                                           data.frame(trial.tm=gain.time.vec, trial.spikeCt=trial.sum.spikes.vec,
                                                      trial.increase=trial.increase.vec, reg.tm=reg.gain.time,
                                                      reg.spikeCt=reg.sum.spikes, reg.increase=reg.increase,
                                                      membership=membership.count))
      
      inits[membership.count,region] = sd
      count.reg = 1
    }
    membership.count = membership.count+count.reg
  }
  
  for (region in names(data.structure)) {
    if (nrow(data.structure[[region]])<100)
      data.structure[[region]] = NULL
  }
  regs.ordrd = intersect(c('V1','LM','RL','AL','AM','PM','TH'), names(data.structure))
  data.structure = data.structure[regs.ordrd]
  inits = inits[regs.ordrd]
  return(list(readout=data.readout,structure=data.structure,inits=inits))
}

# compute the standard error of the peak time estimates for each trial in each region using bootstrap
BoostrapTrials <- function(structure, data, n=100, mltplyr=2, G=15)
{
  print(paste('RUNNING BoostrapTrials FOR SESSION ID:', session_id))
  filter.strt = 20/mltplyr
  filter.end = 130/mltplyr
  sbst = filter.strt:filter.end
  
  regn.names = names(structure)
  trial.bootstrap.samples = vector("list", length(regn.names))
  region.bootstrap.samples = vector("list", length(regn.names))
  names(trial.bootstrap.samples) = names(region.bootstrap.samples) = regn.names
  
  for (regn in regn.names) {
    regn.trials = rownames(structure[[regn]])
    num.regn.trials = length(regn.trials)
    trial.bootstrap.samples[[regn]] = matrix(NA, nrow=num.regn.trials, ncol=G)
    region.bootstrap.samples[[regn]] = matrix(NA, nrow=num.regn.trials, ncol=G)
    colnames(trial.bootstrap.samples[[regn]]) = colnames(region.bootstrap.samples[[regn]]) = paste0('g', 1:G)
    rownames(trial.bootstrap.samples[[regn]]) = rownames(region.bootstrap.samples[[regn]]) = regn.trials
    
    region.bootstrap.data = region.bootstrap.data = lapply(vector("list", G), function(x) data.frame(matrix(NA, n, 15)))
    config.track = NULL
    trials.so.far = character(0)
    ct = 0
    for (trial in regn.trials) {
      ct = ct + 1
      items = strsplit(trial, "_", fixed = TRUE)[[1]]
      config = items[1]
      trl = items[2]
      
      if (ct == 1) config.track = config
      else if (config!=config.track)
      {
        region.bootstrap.est = numeric(G)
        for (g in 1:G)
        {
          region.bootstrap.data[[g]] = region.bootstrap.data[[g]][, !apply(is.na(region.bootstrap.data[[g]]), 2, any)] %>% 
            rowwise %>% 
            mutate(sm=sum(c_across(everything())))
          reg.dat = region.bootstrap.data[[g]]$sm
          reg.fit = fitGam(reg.dat, sbst, n)
          rng = derivatives(reg.fit, n=n)$derivative
          dff = sign(rng[(filter.strt+1):filter.end])-sign(rng[filter.strt:(filter.end-1)]) 
          crit.pt.indx = as.integer(names(which(dff==-2)))
          if (length(crit.pt.indx)==0)
            max.fit.peak.indx = filter.end
          else if (length(crit.pt.indx)==1)
            max.fit.peak.indx = crit.pt.indx
          else
            max.fit.peak.indx = which(reg.fit$fitted.values==max(reg.fit$fitted.values[crit.pt.indx]))
          
          region.bootstrap.est[g] = mltplyr*max.fit.peak.indx
        }
        region.bootstrap.samples[[regn]][trials.so.far,] = matrix(rep(region.bootstrap.est, length(trials.so.far)), ncol=G, byrow=TRUE)
        region.bootstrap.data = region.bootstrap.data = lapply(vector("list", G), function(x) data.frame(matrix(NA, n, 15)))
        trials.so.far = character(0)
        config.track = config
      }
      print(paste(config, regn, trl,sep='_'))
      trials.so.far = c(trials.so.far, trial)
      trial.units = names(data[[config]][[regn]])[!sapply(data[[config]][[regn]], function(x) is.null(x[[trl]]))]
      num.units = length(trial.units)
      trial.bootstrap.est = numeric(G)
      for (g in 1:G) {
        cont = TRUE
        while (cont) {
          trials.all.units = do.call(cbind, lapply(sample(trial.units, size=num.units, replace=TRUE), function(unit) data[[config]][[regn]][[unit]][[trl]][1:n]))
          dat = rowSums(trials.all.units)
          fit = fitGam(dat, sbst, n)
          if (is.null(fit))
          {
            print('FAIL fit!')
            next
          }
          rng = derivatives(fit, n=n)$derivative
          dff = sign(rng[(filter.strt+1):filter.end])-sign(rng[filter.strt:(filter.end-1)])
          crit.pt.indx = as.integer(names(which(dff==-2)))
          if (length(crit.pt.indx)==0)
          {
            print('FAIL crit.pt.indx!')
            next
          }
          cont = FALSE
          if (length(crit.pt.indx)==1)
            max.fit.peak.indx = crit.pt.indx
          else {
            pos.loc = which(fit$fitted.values[crit.pt.indx]==max(fit$fitted.values[crit.pt.indx]))[1]
            if (pos.loc>1) {
              if ((fit$fitted.values[crit.pt.indx[pos.loc-1]]/fit$fitted.values[crit.pt.indx[pos.loc]])>=0.7) {
                max.fit.peak.indx = crit.pt.indx[pos.loc-1]
              }
              else
                max.fit.peak.indx = crit.pt.indx[pos.loc]
            }
            else
              max.fit.peak.indx = crit.pt.indx[pos.loc]
          }
          gain.time = mltplyr*max.fit.peak.indx
          trial.bootstrap.est[g] = gain.time
          region.bootstrap.data[[g]][[length(trials.so.far)]] = dat
        }
      }
      trial.bootstrap.samples[[regn]][trial,] = trial.bootstrap.est
      if (ct ==length(regn.trials))
      {
        region.bootstrap.est = numeric(G)
        for (g in 1:G)
        {
          region.bootstrap.data[[g]] = region.bootstrap.data[[g]][, !apply(is.na(region.bootstrap.data[[g]]), 2, any)] %>% 
            rowwise %>% 
            mutate(sm=sum(c_across(everything())))
          reg.dat = region.bootstrap.data[[g]]$sm
          reg.fit = fitGam(reg.dat, sbst, n)
          rng = derivatives(reg.fit, n=n)$derivative
          dff = sign(rng[(filter.strt+1):filter.end])-sign(rng[filter.strt:(filter.end-1)]) 
          crit.pt.indx = as.integer(names(which(dff==-2)))
          if (length(crit.pt.indx)==0)
            max.fit.peak.indx = filter.end
          else if (length(crit.pt.indx)==1)
            max.fit.peak.indx = crit.pt.indx
          else
            max.fit.peak.indx = which(reg.fit$fitted.values==max(reg.fit$fitted.values[crit.pt.indx]))
          
          region.bootstrap.est[g] = mltplyr*max.fit.peak.indx
        }
        region.bootstrap.samples[[regn]][trials.so.far,] = matrix(rep(region.bootstrap.est, length(trials.so.far)), ncol=G, byrow=TRUE)
      }
    }
    trial.bootstrap.samples[[regn]] = trial.bootstrap.samples[[regn]] %>% as.data.frame %>% rowwise %>% mutate(sd=sd(c_across(everything())))
    region.bootstrap.samples[[regn]] = region.bootstrap.samples[[regn]] %>% as.data.frame %>% rowwise %>% mutate(sd=sd(c_across(everything())))
    rownames(trial.bootstrap.samples[[regn]]) = rownames(region.bootstrap.samples[[regn]]) = regn.trials
  }
  return(list(trial.bootstrap.samples=trial.bootstrap.samples, region.bootstrap.samples=region.bootstrap.samples))
}

# format the trial data for the stan model
formatTrialData <- function(structure, bootstrap.trials, regions=NULL)
{
  print(paste('RUNNING formatTrialData FOR SESSION ID:', session_id))
  trial.bootstrap.samples = bootstrap.trials$trial.bootstrap.samples
  region.bootstrap.samples = bootstrap.trials$region.bootstrap.samples
  if (is.null(regions))
    regions = names(structure)
  for (regn in regions)
  {
    structure[[regn]] = structure[[regn]] %>% 
      rename(trial.sd=trial.spikeCt, reg.sd=reg.spikeCt) %>%
      mutate(trial.sd=trial.bootstrap.samples[[regn]]$sd, reg.sd=region.bootstrap.samples[[regn]]$sd)
  }
  
  data.table = structure[regions] %>%
    map2(regions, ~.x %>%
           rename_with(~paste0(.x,'.',.y),everything(),.y) %>%
           rownames_to_column('rn')) %>%
    reduce(full_join, by='rn') %>%
    rowwise() %>%
    mutate(mem.raw=first(na.omit(c_across(contains('membership'))))) %>%
    ungroup() %>%
    mutate(mem.rank=dense_rank(mem.raw)) %>%
    dplyr::select(!contains(c('membership','mem.raw'))) %>%
    arrange(mem.rank) %>% 
    group_by(mem.rank) %>%
    mutate(mem.row=row_number())
  
  intermediate = data.table %>%
    pivot_longer(contains(c('trial','reg')),
                 names_to=c('.value','region'),
                 names_pattern=paste0('(.+)\\.(',paste(regions,collapse='|'),')')) %>%
    group_by(mem.rank, mem.row) %>%
    mutate(mem.col=row_number()) %>%
    ungroup() %>%
    filter(!is.na(trial.tm))
  
  X_cr = intermediate %>%
    dplyr::select(contains('reg'),mem.rank,mem.col,-region) %>% #col is region, mem.rank is config
    distinct() %>%
    arrange(mem.rank, mem.col) %>%
    mutate(row.indx=row_number())
  
  X_ctr_data = intermediate$trial.tm
  X_cr_data = X_cr$reg.tm 
  X_ctr_sd = replace(intermediate$trial.sd, intermediate$trial.sd==0, 0.01)
  X_cr_sd = replace(X_cr$reg.sd, X_cr$reg.sd==0, 0.01)
  
  intermediate_1 = X_cr %>%
    dplyr::select(mem.rank,mem.col,row.indx)
  
  intermediate_2 = intermediate %>%
    full_join(intermediate_1, by=c('mem.rank','mem.col'))
  
  X_ctr_info = intermediate_2 %>%
    dplyr::select(mem.rank,mem.row,mem.col,row.indx)
  
  input_data = list(Nctr=length(X_ctr_data), Ncr=length(X_cr_data),
                    c=max(X_cr$mem.rank), t=max(intermediate$mem.row), r=length(regions),
                    reg_index=ifelse(('V1' %in% regions),which(regions=='V1'),1),
                    X_ctr_data=X_ctr_data, X_ctr_sd=X_ctr_sd, X_ctr_info=X_ctr_info,
                    X_cr_data=X_cr_data, X_cr_sd=X_cr_sd)
  
  return(input_data)
}

# initialize the stan model parameters
paramsInitials <- function(inits, regions, num.chains=4)
{
  inits = inits[regions] %>%
    rowwise() %>%
    mutate(mask=first(na.omit(c_across(everything())))) %>%
    ungroup() %>%
    filter(!is.na(mask)) %>%
    dplyr::select(!mask)
  parInits = list(gamma_cr_raw = unname(unlist(inits)) %>% replace_na(1))
  parInits$gamma_cr_raw[!parInits$gamma_cr_raw] = 1
  initials = list()
  for (i in 1:num.chains) {
    initials[[i]] = parInits
  }
  return(initials)
}

# from the stan model output, extract the posterior peak time estimates for each region
writeLagsSingle <- function(fit_cp, input_data, regns, alpha=0.05)
{
  unpermuted = rstan::extract(fit_cp, pars='P_ctr_unstructured', permute=FALSE, inc_warmup=FALSE)
  dms = dim(unpermuted)
  num.chains = dms[2]
  trace = data.frame(matrix(NA,0,0))
  for (chain in 1:num.chains)
    trace = bind_rows(trace, as.data.frame(unpermuted[,chain,]))
  
  num.pars = dms[3]
  num.iters = nrow(trace)
  structure = data.frame(matrix(NA, input_data$c*input_data$t*num.iters, (input_data$r+1)))
  for (indx in 1:num.pars) {
    row = input_data$X_ctr_info[indx,]
    c = row$mem.rank
    t = row$mem.row
    r = row$mem.col
    structure[(1:num.iters)+(input_data$t*(c-1)+t-1)*num.iters, r] = trace[,indx]
    structure[(1:num.iters)+(input_data$t*(c-1)+t-1)*num.iters, (input_data$r+1)] = input_data$t*(c-1)+t #trial
  }
  structure$iter = rep(1:num.iters, input_data$t*c)
  qntls = abs(c(0,-1)+0.5*alpha)
  
  #single
  data = data.frame(matrix(NA,0,0))
  for (i in 1:input_data$r) {
    diffs = structure %>% 
      mutate(dif=.[[i]]) %>%
      select(dif,input_data$r+1,iter) %>%
      filter(!is.na(dif)) %>%
      group_by(.[2]) %>% # grouped by trial
      mutate(trl.post.se2=var(dif), trl.post.cntr=mean(dif)) %>%
      ungroup %>%
      mutate(post.cntr.se2=var(unique(trl.post.cntr)), 
             wts=1/(trl.post.se2+post.cntr.se2), 
             sum.wts=1/sum(wts)) %>%
      group_by(iter) %>%
      summarise(sample=sum(dif*wts*sum.wts)) %>% .$sample
    diffs = diffs*length(diffs) #rescaling because they were scaled within and across iterations
    # it is really unweighting but every sample has the same weight
    lag = mean(diffs)
    ci = unname(quantile(diffs,qntls))
    data = bind_rows(data, data.frame(mouse.id=session_id, row.reg=regns[i],
                                      lag=lag, lci=ci[1], uci=ci[2]))
  }
  filePath = file.path(output.path, 'output', 'writeCorrMatrixInfoLagsSingle')
  dir.create(filePath)
  filePathNow = file.path(filePath, paste0(session_id, '_lags.csv'))
  write.csv(data, filePathNow, row.names=F)
}

# from the stan model output, extract the posterior correlation estimates for each region
writeCorrs <- function(fit_cp, input_data, regns, alpha=0.05)
{
  unpermuted = rstan::extract(fit_cp, pars='Omega', permute=FALSE, inc_warmup=FALSE)
  dms = dim(unpermuted)
  num.chains = dms[2]
  trace = data.frame(matrix(NA,0,0))
  for (chain in 1:num.chains) {
    t = as.data.frame(unpermuted[,chain,])
    select.pars = names(t) %>%
      str_match(paste0('^Omega.([0-9]).([0-9]).$')) %>%
      as.data.frame() %>%
      mutate(V2=as.numeric(V2), V3=as.numeric(V3), V4=V3>V2) %>%
      filter(V4) %>% .$V1
    trace = bind_rows(trace, t[select.pars])
  }
  qntls = abs(c(0,-1)+0.5*alpha)
  
  data = data.frame(matrix(NA,0,0))
  corr.names = names(trace)
  for (nm in corr.names) {
    loc.dat = trace[[nm]]
    corr = mean(loc.dat)
    ci = unname(quantile(loc.dat, qntls))
    loc = strsplit(nm,'\\[|,|\\]')[[1]]
    data = bind_rows(data, data.frame(mouse.id=session_id, 
                                      row.reg=regns[as.integer(loc[2])], 
                                      col.reg=regns[as.integer(loc[3])], 
                                      corr=corr, lci=ci[1], uci=ci[2]))
  }
  filePath = file.path(output.path, 'output', 'writeCorrMatrixInfoMarginals')
  dir.create(filePath)
  filePathNow = file.path(filePath, paste0(session_id, '_marginals.csv'))
  write.csv(data, filePathNow, row.names=F)
}

# from the stan model output, extract the posterior partial correlation estimates for each region
writePCorrs <- function(trials, num.chains, iter, alpha=0.05)
{
  regns.full = names(trials$structure)[names(trials$structure)!='V1']
  end = length(regns.full)
  if (end<2) return()
  data = data.frame(matrix(NA,0,0))
  for (i in 1:(end-1)) {
    for (j in (i+1):end) {
      regns = c('V1', regns.full[c(i,j)])
      print(paste(regns,collapse='_'))
      print(paste('RUNNING stan fit FOR SESSION ID:', session_id))
      initials = paramsInitials(trials$inits, regns, num.chains)
      input_data = formatTrialData(trials$structure, bootstrap.trials, regns)
      fit_cp = stan(file='corr_posterior.stan', data=input_data, chains=num.chains, iter=iter,
                    init=initials, pars=c('Z_ctr_raw','Z_cr_raw','L_Omega','gamma_cr_raw'), include=FALSE)
      
      unpermuted = rstan::extract(fit_cp, pars='Pcorr', permute=FALSE, inc_warmup=FALSE)
      dms = dim(unpermuted)
      num.chains = dms[2]
      trace = data.frame(matrix(NA,0,0))
      for (chain in 1:num.chains) {
        trace = trace %>%
          bind_rows(as.data.frame(unpermuted[,chain,]) %>% 
                      select('Pcorr[2,3]'))
      }
      qntls = abs(c(0,-1)+0.5*alpha)
      
      loc.dat = trace[[1]]
      pcorr = mean(loc.dat)
      ci = unname(quantile(loc.dat, qntls))
      data = data %>% 
        bind_rows(data.frame(mouse.id=session_id, given='V1',
                             row.reg=regns[2], col.reg=regns[3],  
                             pcorr=pcorr, lci=ci[1], uci=ci[2]))
    }
  }
  filePath = file.path(output.path, 'output', 'writeCorrMatrixInfoPartials')
  dir.create(filePath)
  filePathNow = file.path(filePath, paste0(session_id, '_partials.csv'))
  write.csv(data, filePathNow, row.names=F)
}

# write the naieve correlation estimates for each region
writeCorrsNoisy <- function(trials, regns.full)
{
  end = length(regns.full)
  data = data.frame(matrix(vector(), 0, 0))
  for (i in 1:(end-1)) {
    for (j in (i+1):end) {
      regns = regns.full[c(i,j)]
      print(paste(regns,collapse='_'))
      rel.reg = trials$structure[c(regns)]
      info = rel.reg[[1]] %>% select(trial.tm) %>% 
        merge(rel.reg[[2]]%>%select(trial.tm), by='row.names')
      corrs = cor(info[[2]],info[[3]])
      row.reg = regns[1]
      col.reg = regns[2]
      data = bind_rows(data, data.frame(mouse.id=session_id, row.reg=row.reg, col.reg=col.reg, 
                                        noisycorr=corrs, lci=corrs, uci=corrs))
    }
  }
  filePath = file.path(output.path, 'output', 'writeCorrMatrixInfoNoisyMarginals')
  dir.create(filePath)
  filePathNow = file.path(filePath, paste0(session_id, '_noisymarginals.csv'))
  write.csv(data, filePathNow, row.names=F)
}
