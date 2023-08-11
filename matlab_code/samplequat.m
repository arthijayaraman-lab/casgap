function result=samplequat(num,pref_q,kappa)
v=-sqrt(2)*erfcinv(2*rand(num,2));
vsum=sqrt(sum(v.^2,2));
v(:,1)=v(:,1)./vsum;
v(:,2)=v(:,2)./vsum;
%scatter3(v(:,1),v(:,2),v(:,3));
if kappa
    randnum=rand(num,1);
    w=1+1/kappa*log(randnum+(1-randnum)*exp(-2*kappa));
else
    w=2*rand(num,1)-1;
end

orgvec = [zeros(num,2) ones(num,1)]; 
newvec = [w sqrt(1-w.^2).*v(:,1) sqrt(1-w.^2).*v(:,2)];

axisvec = cross(orgvec,newvec,2);
axisvec = axisvec./sum(axisvec.^2,2);
axistheta = acos(dot(orgvec,newvec,2));

quats=[cos(axistheta/2) sin(axistheta/2).*axisvec];
if kappa
    Mumat= zeros(4,4);
    Mumat(:,1)=pref_q';
    [Q,R]=qr(Mumat);
    if R(1,1) < 0
        Q=-Q;
    end
    result=(Q*quats')';
else
    result=quats;
end