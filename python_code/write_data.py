def write_data(output, particlelist, params):
    population = particlelist['Nprime']
    fileID = open(output['path'] + output['file'], 'w')
    fileID.write('ITEM: TIMESTEP\n0\n')
    fileID.write('ITEM: NUMBER OF ATOMS\n%d\n' % population)
    fileID.write('ITEM: BOX BOUNDS pp pp pp\n%f %f\n%f %f\n%f %f\n' % (-params['boxlength']/2, params['boxlength']/2, -params['boxlength']/2, params['boxlength']/2, -params['boxlength']/2, params['boxlength']/2))
    fileID.write('ITEM: ATOMS id type x y z a b c qw qx qy qz\n')
    Alldata = [(i+1, particlelist['xyz'][i][0], particlelist['xyz'][i][1], particlelist['xyz'][i][2], particlelist['ac'][i][0], particlelist['ac'][i][0], particlelist['ac'][i][1], particlelist['quat'][i][0], particlelist['quat'][i][1], particlelist['quat'][i][2], particlelist['quat'][i][3]) for i in range(population)]
    for data in Alldata:
        fileID.write('%d 1 %f %f %f %f %f %f %f %f %f %f\n' % data)
    
    fileID.close()
