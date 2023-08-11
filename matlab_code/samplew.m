function sampledvalues = samplew(size,kappa)
%size=100;
%kappa=1000;
u=rand(size,1);
w=(-1:0.0001:1)';
%kappa=100;
first_two_terms=1+(w.*sqrt(1-w.^2)-acos(w))/pi;
%scaled version of the bessel function is used to avoid overflow for large
%kappa and since the bessel functions are present in a ratio, the scaling
%factors will exactly cancel out.
%mplus2_term=@(m)-2/pi*(besseli(m-2,kappa,1)-2*besseli(m,kappa,1)+besseli(m+2,kappa,1))./(besseli(-2,kappa,1)-2*besseli(0,kappa,1)+besseli(2,kappa,1)).*sin(m*acos(w))/m;
mplus2_term=@(m)(m+1)/pi*(sin((m+2)*acos(w))/(m+2)-sin(m*acos(w))/m)*besseli(m+1,kappa,1)/besseli(1,kappa,1);
Fw=first_two_terms;
mvalue=1;
if kappa
    while 1
        newterm=mplus2_term(mvalue);
        newterm_magnitude=sqrt(sum(newterm.^2));
        if isnan(newterm_magnitude)
            error("Nan values detected, kappa is too high!");
        elseif abs(newterm_magnitude)<1e-10
            %disp(['mvalue = ' num2str(mvalue)]);
            break;
        end
        Fw=Fw+newterm;
        mvalue=mvalue+1;
    end
end
%check for monotonic behavior
Fw(abs(Fw)<1e-10)=0;
if ~all(diff(Fw)>=0)
    error("Cummulative function is not monotonic!");
end
%at this point the data is at least monotonic, to make it strictly
%monotonic add a small increasing arithmetic progression
smallseries=linspace(0,1e-15,length(w))';
data=[w (Fw+smallseries)/(Fw(end)+smallseries(end))];
sampledvalues=interp1(data(:,2),data(:,1),u);
%plot(data(:,1),data(:,2));