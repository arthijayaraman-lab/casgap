function [intersection_flag,simplex] = gjk_simplex(polyhedron1,polyhedron2)
direction=polyhedron1.center-polyhedron2.center;
support = supportpoint(polyhedron1,polyhedron2,direction);
simplexvertices = support;
direction = -support;
while 1
    support = supportpoint(polyhedron1,polyhedron2,direction);
    if all(~direction) || dot(support,direction)<=0
        intersection_flag=0;
        simplex=struct('vertices',[]);
        break;
    end
    simplexvertices = [support;simplexvertices]; %#ok<AGROW>
    [flag,simplexvertices,direction] = nextsimplex(simplexvertices);

    if flag
        if size(simplexvertices,1)<4
            intersection_flag=0;
            simplex=struct('vertices',[]);
        else
            intersection_flag=1;
            simplexfaces=gettetrahedronfaces(simplexvertices,[0 0 0]);
            simplex=struct('vertices',simplexvertices,'faces',{simplexfaces});
        end
        break;
    end
end

end

function [flag,vertices,direction] = nextsimplex(vertices)
switch size(vertices,1)
    case 2
        [flag,vertices,direction]=linesimplex(vertices);
    case 3
        [flag,vertices,direction]=trianglesimplex(vertices);
    case 4
        [flag,vertices,direction]=tetrahedronsimplex(vertices);
    otherwise
        flag=0;
end
end

function [flag,vertices,direction]=linesimplex(vertices)
direction=[0 0 0];
point_a = vertices(1,:);
point_b = vertices(2,:);
vec_ab = point_b-point_a;
vec_ao = -point_a;
norm_abo = cross(vec_ab,vec_ao);
if norm(norm_abo) %check for collinearity
    direction = cross(norm_abo,vec_ab);
    flag = 0;
else
    flag = 1;
end
end

function [flag,vertices,direction]=trianglesimplex(vertices)
direction=[0 0 0];
point_a = vertices(1,:);
point_b = vertices(2,:);
point_c = vertices(3,:);

vec_ab = point_b-point_a;
vec_ac = point_c-point_a;
vec_ao = -point_a;

facenorm_abc = cross(vec_ab,vec_ac);

if norm(facenorm_abc) %check if origin is coplanar
    if dot(facenorm_abc,vec_ao)>0
        direction = facenorm_abc;
        flag = 0;
    else
        direction = -facenorm_abc;
        flag = 0;

    end
else
    flag = 1;
end
end

function [flag,vertices,direction]=tetrahedronsimplex(vertices)
direction=[0 0 0];
point_a = vertices(1,:);
point_b = vertices(2,:);
point_c = vertices(3,:);
point_d = vertices(4,:);

vec_ab = point_b-point_a;
vec_ac = point_c-point_a;
vec_ad = point_d-point_a;
vec_ao = -point_a;

facenorm_abc = cross(vec_ab,vec_ac);
if dot(facenorm_abc,vec_ad)>0
    facenorm_abc=-facenorm_abc;
end
facenorm_acd = cross(vec_ac,vec_ad);
if dot(facenorm_acd,vec_ab)>0
    facenorm_acd=-facenorm_acd;
end
facenorm_abd = cross(vec_ab,vec_ad);
if dot(facenorm_abd,vec_ac)>0
    facenorm_abd=-facenorm_abd;
end

outsideabc=(dot(facenorm_abc,vec_ao)>0);
outsideacd=(dot(facenorm_acd,vec_ao)>0);
outsideabd=(dot(facenorm_abd,vec_ao)>0);

if (outsideabc && outsideacd) || (outsideacd && outsideabd) ||(outsideabd && outsideabc)
    flag=0;
elseif outsideabc
    vertices = [point_a;point_b;point_c];
    direction = facenorm_abc;
    flag=0;
elseif outsideacd
    vertices = [point_a;point_c;point_d];
    direction = facenorm_acd;
    flag=0;
elseif outsideabd
    vertices = [point_a;point_b;point_d];
    direction = facenorm_abd;
    flag=0;
else
    flag=1;
end
end

function [faces]=gettetrahedronfaces(vertices, interiorpoint)
if size(vertices,1)~=4
    error('Input should be a set of tetrahedral vertices')
end
faces={[1 2 3]; [1 4 2]; [1 3 4]; [2 4 3]};
for i=1:4 %number of faces in the simplex
    face=faces{i};
    point_a = vertices(face(1),:);
    point_b = vertices(face(2),:);
    point_c = vertices(face(3),:);

    vec_ab = point_b-point_a;
    vec_bc = point_c-point_b;
    facenorm_abc = cross(vec_ab,vec_bc);
    if(dot(facenorm_abc,point_b-interiorpoint)<0)
        faces{i}=[face(1) face(3) face(2)];
    end
end
end