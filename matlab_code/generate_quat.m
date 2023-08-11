function q = generate_quat(vec,alpha) %angle in degree
%alpha=60;
%vec=[1 0 0];
norm=vec/sqrt(sum(vec.^2));
alpha_by_2=alpha/2*pi/180; %Also convert degree to radian
q=[cos(alpha_by_2)...
    norm(1)*sin(alpha_by_2)...
    norm(2)*sin(alpha_by_2)...
    norm(3)*sin(alpha_by_2)];
end