function [shiftvector,shiftdist] = expandingpolytope_shift(polyhedron1,polyhedron2,simplex)
    polytope = struct('vertices',simplex.vertices,'faces',{simplex.faces});
    facenorms = getfacenormals(polytope);
    plotminkowskiflag=0;
    tolerance=1e-6;
    while 1
        [minInd,mindist] = findnearestface(polytope,facenorms);
        nearestnormal = facenorms(minInd,:);
        support = supportpoint(polyhedron1,polyhedron2,nearestnormal);
        supportdist = dot(nearestnormal,support);
        if plotminkowskiflag
            plotminkowskiepa;
        end
        if (supportdist - mindist) > tolerance
            numnormals=size(facenorms,1);
            uniqueedges=[];
            newfaces={};
            newnormals=[];
            for i = 1:numnormals
                face=polytope.faces{i};
                facevertex1=polytope.vertices(face(1),:);
                if dot(facenorms(i,:),support-facevertex1)>tolerance
                    faceedges = [face(1:end-1)' face(2:end)'; face(end) face(1)];
                    uniqueedges = finduniqueedges(uniqueedges,faceedges);
                else
                    newfaces=[newfaces;polytope.faces(i,:)]; %#ok<AGROW> 
                    newnormals=[newnormals;facenorms(i,:)]; %#ok<AGROW>
                end
            end
            %add the new support point to the polytope
            polytope.vertices = [polytope.vertices;support];
            supportind = size(polytope.vertices,1);
            numuniqueedges=size(uniqueedges,1);
            for i = 1:numuniqueedges
                [newaddedface,newaddedfacenorm] = createnewpolytopeface(polytope,uniqueedges(i,:),supportind);
                newfaces=[newfaces;newaddedface]; %#ok<AGROW>
                newnormals=[newnormals;newaddedfacenorm]; %#ok<AGROW>
            end
            polytope.faces=newfaces;
            facenorms=newnormals;
        else
            shiftvector = nearestnormal;
            shiftdist = mindist;
            break;
        end
    end
end

function facenorms=getfacenormals(polytope) %Assumes right hand rule to get the face normal
faces=polytope.faces;
numfaces = size(faces,1);
facenorms=zeros(numfaces,3);
for i=1:numfaces %get the first three points of the face
    face=faces{i};
    point_a = polytope.vertices(face(1),:);
    point_b = polytope.vertices(face(2),:);
    point_c = polytope.vertices(face(3),:);
    vec_ab = point_b-point_a;
    vec_bc = point_c-point_b;
    facenorm_abc = cross(vec_ab,vec_bc);
    facenorms(i,:)=facenorm_abc./norm(facenorm_abc); %These are unit vectors
end
end

function uniqueedges = finduniqueedges(uniqueedges,faceedges)
if ~size(uniqueedges,1)
    uniqueedges = faceedges;
else
    for i=1:size(faceedges,1)
        %check if the reverse of a uniqueedge is in the faceedge. If there
        %is, then delete.
        uniqueedge_vert1 = uniqueedges(:,1)-faceedges(i,2);
        uniqueedge_vert2 = uniqueedges(:,2)-faceedges(i,1);
        uniqueness_check = (~uniqueedge_vert1) & (~uniqueedge_vert2);
        if any(uniqueness_check)
            uniqueedges(uniqueness_check,:)=[];
        else
            uniqueedges = [uniqueedges;faceedges(i,:)]; %#ok<AGROW> 
        end
    end
end
end

function [face,facenorm] = createnewpolytopeface(polytope,edge,pointind) %assumes that the normal is away from the origin
    point_a = polytope.vertices(edge(1),:);
    point_b = polytope.vertices(edge(2),:);
    point_c = polytope.vertices(pointind,:);
    tolerance=1e-6;
    vec_ab = point_b-point_a;
    vec_bc = point_c-point_b;
    facenorm_abc = cross(vec_ab,vec_bc);
    facenorm=facenorm_abc./norm(facenorm_abc); %These are unit vectors
    facedist = dot(point_a,facenorm);
    if facedist> tolerance
        face = {[edge(1) edge(2) pointind]};
    else
        face = {[edge(2) edge(1) pointind]};
        facenorm = -facenorm;
    end
end

function [minInd,mindist] = findnearestface(polytope,facenorms)
firstel=@(x)x(1); %define an inline function to get the first element of an array
distances = dot(polytope.vertices(cellfun(firstel,polytope.faces),:),facenorms,2);
[mindist,minInd]=min(distances);
end