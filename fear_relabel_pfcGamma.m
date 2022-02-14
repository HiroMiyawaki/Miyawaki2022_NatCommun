function fear_relabel_pfcGamma(basename)
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with date of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basename '.pfcGamma.events.mat'])
load([basename '.pfclowGamma.events.mat'])

pfcRipple=pfcGamma(4);
pfcFastGamma=pfcGamma(1);
pfcSlowGamma=pfcLowGamma;
 
pfcRipple.generator=mfilename;
pfcRipple.generatedate=datestr(now,'yyyy-mm-dd');

pfcFastGamma.generator=mfilename;
pfcFastGamma.generatedate=datestr(now,'yyyy-mm-dd');

pfcSlowGamma.generator=mfilename;
pfcSlowGamma.generatedate=datestr(now,'yyyy-mm-dd');

save([basename '.pfcRipple.events.mat'],'pfcRipple','-v7.3')
save([basename '.pfcFastGamma.events.mat'],'pfcFastGamma','-v7.3')
save([basename '.pfcSlowGamma.events.mat'],'pfcSlowGamma','-v7.3')
