function [particlelist] = phaseI(params)
tic
target_volfrac = params.volfrac;
box_vol=params.boxlength^3;
vol_prefactor=4/3*pi;
mean_vol=vol_prefactor*params.mean_R^3;

%Total number of particles
particlelist.N = round(2*target_volfrac/mean_vol*box_vol); %factor of 2 is for good measure.

%Params for lognormal distribution
R_logmu=log(params.mean_R^2/sqrt(params.mean_R^2+params.sd_R^2));
gamma_logmu=log(params.mean_gamma^2/sqrt(params.mean_gamma^2+params.sd_gamma^2));
R_logsigma=sqrt(log(1+params.sd_R^2/params.mean_R^2));
gamma_logsigma=sqrt(log(1+params.sd_gamma^2/params.mean_gamma^2));

%Generate random samples of R 
Ri = lognormrandvar([particlelist.N 1],R_logmu,R_logsigma);
partialvoli = vol_prefactor.*Ri.^3/box_vol;
actual_volfrac = sum(partialvoli);
while actual_volfrac <= target_volfrac
    Ri_extra = lognormrandvar([particlelist.N 1],R_logmu,R_logsigma);
    Ri = [Ri;Ri_extra]; %#ok<AGROW> 
    voli_extra = vol_prefactor.*Ri_extra.^3/box_vol;
    partialvoli = [partialvoli;voli_extra]; %#ok<AGROW> 
    actual_volfrac = sum(partialvoli);
end
if actual_volfrac > target_volfrac
    particlelist.N=find(cumsum(partialvoli)>target_volfrac,1,"first")-1;
    Ri=Ri(1:particlelist.N);
end
%Generate random samples of gamma 
gammai = lognormrandvar([particlelist.N 1],gamma_logmu,gamma_logsigma);

%Set spheroid axis lengths
ai=Ri./gammai.^(1/3);
ci=Ri.*gammai.^(2/3);
particlelist.ac = [ai ci];
Lambda= generate_quat(params.W,params.omega); 
particlelist.quat = samplequat(particlelist.N,Lambda,params.kappa);
disp(['Phase I completed! N is ' num2str(particlelist.N) '. The time elapsed during Phase I is ' num2str(toc) ' seconds.']);
end

function result = lognormrandvar(size,logmu,logsigma)
temp=-sqrt(2)*erfcinv(2*rand(size)); %std normal random variable
result = exp(logmu+temp*logsigma);
end