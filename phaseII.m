function [successflag,particlelist] = phaseII(particlelist,params)
if particlelist.Nprime
    %Check and adjust coords so that the new ellipsoid does not overlap
    %with any previous ellipsoids.
    [successflag,coords]=check_overlap(particlelist,params);
    if ~successflag
        return;
    end
    %Increment population
    population = particlelist.Nprime+1;
    particlelist.xyz(population,:)=coords;
else
    % It is the first particle, there can't be any overlap!
    successflag=1;
    population = 1;
    coords=particlelist.xyz(population,:);
end

%Append polyhedra
axislengths=particlelist.ac(population,:);
quats=particlelist.quat(population,:);
particlelist.polyhedra(:,:,population)=discretize_ellipsoid(axislengths,quats,coords);
particlelist.Nprime=population;
end