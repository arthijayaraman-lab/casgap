function [point,furthestpoint1,furthestpoint2] = supportpoint(polyhedron1,polyhedron2,direction)
furthestpoint1=furthestpoint(polyhedron1.vertices,direction);
furthestpoint2=furthestpoint(polyhedron2.vertices,-direction);
point = furthestpoint1-furthestpoint2;
end

function maxpoint = furthestpoint(vertices,direction)
maxpoint=[0 0 0];
maxprojection=-inf;
for i=1:size(vertices,1)
    projection=dot(vertices(i,:),direction);
    if maxprojection<projection
        maxprojection=projection;
        maxpoint=vertices(i,:);
    end
end
end