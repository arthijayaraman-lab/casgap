function write_data(output,particlelist,params)
population=particlelist.Nprime;
fileID = fopen([output.path output.file],'w');
fprintf(fileID,'ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n%d\n',population);
fprintf(fileID,'ITEM: BOX BOUNDS pp pp pp\n%f %f\n%f %f\n%f %f\n',-params.boxlength/2,params.boxlength/2,-params.boxlength/2,params.boxlength/2,-params.boxlength/2,params.boxlength/2);
fprintf(fileID,'ITEM: ATOMS id type x y z a b c qw qx qy qz\n');
Alldata=[(1:population)' particlelist.xyz(1:population,:) particlelist.ac(1:population,[1 1 2]) particlelist.quat(1:population,:)];
fprintf(fileID,'%d 1 %f %f %f %f %f %f %f %f %f %f\n',Alldata');
fclose(fileID);
end
