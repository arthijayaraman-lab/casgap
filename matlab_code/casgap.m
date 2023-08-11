function casgap(params,restart,output)
%Initialize random number generator
currentseed=rng(params.seed);
mkdir(restart.path);
mkdir(output.path);
%If restarting, skip phase I
if restart.flag    
	load([restart.path restart.file]);
    rng(currentseed);
else
	particlelist=phaseI(params);
    particlelist.polyhedra=struct('vertices',[],'faces',{},'center',[]);
    particlelist.Nprime=0; %Current population
    %Initialize coordinates
    particlelist.xyz=(rand(particlelist.N,3)-0.5).*params.boxlength;
    loopstart = 1;
    loopend = particlelist.N;
    currentseed=rng;
    save([restart.path restart.file]);
end

%Use particle list from restart or phase I as input for phase II
tic
for i = loopstart:loopend
    [successflag,particlelist] = phaseII(particlelist,params);
    if ~successflag
        break;
    end
    %Restart
    if ~mod(i,restart.freq)
        currentseed=rng;
        loopstart=i;
        save([restart.path restart.file]);
        disp(['Writing Restart File! Current population is ' num2str(i) '. The time elapsed is ' num2str(toc) ' seconds.']);
    end
    %Output
    if ~mod(i,output.freq)
        write_data(output,particlelist,params);
    end
end
%Final restart
currentseed=rng;
save([restart.path restart.file]);
%Final output
write_data(output,particlelist,params);

partialvoli = 4/3*pi.*particlelist.ac(1:particlelist.Nprime,1).^2.*particlelist.ac(1:particlelist.Nprime,2);
actual_volfrac = sum(partialvoli)/params.boxlength^3;
if ~successflag
    disp(['Not all ellipsoids could be added. Final population is ' num2str(particlelist.Nprime) '. Final volume fraction is ' num2str(actual_volfrac) '. The time elapsed is ' num2str(toc) ' seconds.']);
else
    disp(['All ellipsoids were successfully added! Final population is ' num2str(particlelist.Nprime) '. Final volume fraction is ' num2str(actual_volfrac) '. The time elapsed is ' num2str(toc) ' seconds.']);
end