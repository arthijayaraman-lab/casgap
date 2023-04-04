function [successflag,coords] = check_overlap(particlelist,params)
successflag=0;
%
population=particlelist.Nprime;
newind=population+1;
ac=particlelist.ac(newind,:);
q=particlelist.quat(newind,:);
boundingradius=max(ac);
newpolyhedra_unshifted=discretize_ellipsoid(ac,q,[0 0 0]);
%
maxattempts=10000;
perturbfreq=100; %resample new coords at this rate
maxradii=max(particlelist.ac(1:population,:),[],2);
pos=particlelist.xyz(1:population,:);

EPAclearance=0.1;
for numattempts = 1:maxattempts
    if mod(numattempts,perturbfreq)==1
        %this assumes that previous attempts at shifting coords were not too successful
        %sample new coords
        if numattempts<perturbfreq %first sampling is uniformly random
            coords=(rand(1,3)-0.5)*params.boxlength;
        else
            randvars=rand(1,3);
            zbins = (-5:5)/10*params.boxlength;
            znum = histcounts(pos(:,3),zbins);
            maxznum=max(znum);
            minznum=min(znum);
            if maxznum==minznum
                zcoord = (randvars(1)-0.5)*params.boxlength;
            else
                newznum = (maxznum-znum)/(maxznum-minznum);
                cumznum=[0 cumsum(newznum)]+(0:10)*1e-10; %Small constants added to make them unique
                cumznum=cumznum/cumznum(end);
                zcoord = interp1(cumznum,zbins,randvars(3));
            end
            zbinwidth=zbins(2)-zbins(1);
            zfilter=(pos(:,3)>zcoord-zbinwidth/2) & (pos(:,3)<=zcoord+zbinwidth/2);
            ybins = (-5:5)/10*params.boxlength;
            ynum = histcounts(pos(zfilter,2),ybins);
            maxynum=max(ynum);
            minynum=min(ynum);
            if ~sum(zfilter) || maxynum==minynum
                ycoord = (randvars(2)-0.5)*params.boxlength;
            else
                newynum = (maxynum-ynum)/(maxynum-minynum);
                cumynum=[0 cumsum(newynum)]+(0:10)*1e-10; %Small constants added to make them unique
                cumynum=cumynum/cumynum(end);
                ycoord = interp1(cumynum,ybins,randvars(2));
            end
            ybinwidth=ybins(2)-ybins(1);
            yzfilter=(pos(zfilter,2)>ycoord-ybinwidth/2) & (pos(zfilter,2)<=ycoord+ybinwidth/2);
            xbins = (-5:5)/10*params.boxlength;
            xnum = histcounts(pos(yzfilter,1),xbins);
            maxxnum=max(xnum);
            minxnum=min(xnum);
            if ~sum(yzfilter) || maxxnum==minxnum
                xcoord = (randvars(2)-0.5)*params.boxlength;
            else
                newxnum = (maxxnum-xnum)/(maxxnum-minxnum);
                cumxnum=[0 cumsum(newxnum)]+(0:10)*1e-10; %Small constants added to make them unique
                cumxnum=cumxnum/cumxnum(end);
                xcoord = interp1(cumxnum,xbins,randvars(1));
            end
            coords=[xcoord ycoord zcoord];
        end
    end
    newpolyhedra=newpolyhedra_unshifted;
    newpolyhedra.vertices=newpolyhedra.vertices+ones(size(newpolyhedra.vertices,1),1)*coords;
    dist=sqrt(sum((ones(population,1)*coords-pos).^2,2));
    potentialoverlap=(maxradii+boundingradius-dist>0);
    if all(~potentialoverlap)
        successflag=1;
        break;
    end
    indices=find(potentialoverlap);
    numintersections=0;
    cummulativeshift=[0 0 0];
    for i=1:length(indices)
        [intersectionflag,intersectionsimplex] = gjk_simplex(newpolyhedra,particlelist.polyhedra(indices(i)));
        if intersectionflag
            numintersections=numintersections+1;
            [shiftvector,shiftdist] = expandingpolytope_shift(newpolyhedra,particlelist.polyhedra(indices(i)),intersectionsimplex);
            cummulativeshift = -shiftvector*(shiftdist+EPAclearance);
        end
    end
    if any(cummulativeshift)
        if ~mod(numattempts,100)
            disp(['Warning: Previous attempts failed to converge after ' num2str(numattempts) ' attemps. Current ellipsoid has ' num2str(numintersections) ' intersections. The population is ' num2str(population) '.']);
        end
        coords=coords+cummulativeshift;
        coords=mod(coords+params.boxlength/2,params.boxlength)-params.boxlength/2;
    else
        successflag=1;
        break;
    end
end
end