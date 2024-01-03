
args = commandArgs(trailingOnly=TRUE)
session_id = toString(args[1])

data.path = path.expand('~/data/r/dataset')
output.path = path.expand('~/projects/r/ADAproject')
project.path = path.expand('~/projects/r/ADAproject/src')
setwd(project.path)
source('functions.R')
iter = 20000
num.chains = 6
options(mc.cores=12)

Sys.getpid()
print('')
print(paste('RUNNING analysis FOR SESSION ID:', session_id))
data = LoadDataFiltered()
trials = selectTrials(data)
bootstrap.trials = BoostrapTrials(trials$structure, data, G=15)

regns = names(trials$structure)
if (length(regns)<2) {
  print('Returned less than 2 areas')
} else {
  print(paste('RUNNING stan fit FOR SESSION ID:', session_id))
  initials = paramsInitials(trials$inits, regns, num.chains)
  input_data = formatTrialData(trials$structure, bootstrap.trials)
  fit_cp = stan(file='corr_posterior.stan', data=input_data, chains=num.chains, iter=iter,
                init=initials, pars=c('Z_ctr_raw','Z_cr_raw','L_Omega','gamma_cr_raw'), include=FALSE)
  
  print(paste('RUNNING post processing FOR SESSION ID:', session_id))
  writeLagsSingle(fit_cp, input_data, regns)
  writeCorrs(fit_cp, input_data, regns)
  writeCorrsNoisy(trials, regns)
}

if (length(regns)<3 | !('V1' %in% regns)) {
  print('Returned insufficient areas to compute pcorr given V1.')
} else {
  writePCorrs(trials, num.chains, iter)
}
print('DONE!')
