from boto.s3.connection import S3Connection
import os
import re
import numpy as np
import iris

def list_files(bucket, prefix='prods', n=100):
    '''List the first n files in an S3 bucket, with the given prefix.'''
    os.environ['S3_USE_SIGV4'] = 'True'
    conn = S3Connection(host='s3.eu-west-2.amazonaws.com')
    bucket = conn.get_bucket(bucket)
    results = []
    keys = iter(bucket.list(prefix=prefix))
    for i in range(n):
        try:
            results.append(next(keys))
        except StopIteration:
            pass
    return ['{}'.format(k.key) for k in results]


def get_info(fname):
    '''
    Extract the sections of a Mogreps filename. 
    Returns [date string, runtime (hour), ensemble member, lead time (hours)]
    
    e.g. get_info('prods_op_mogreps-g_20160107_00_00_003.nc') 
     -> ['20160107', 0, 0, 3]
    '''
    segs = fname.split('_')
    segs[-1] = segs[-1].split('.')[0]
    info = [segs[-4]] + [int(s) for s in segs[-3:]]
    return info


def get_ground_level(cube):
    '''
    Extract the level with either minimum height or maximum pressure.
    If neither is found, return the entire cube.
    Returns (ground-level cube, coord name)
    '''
    if any([re.match('height.*', d.name()) for d in cube.coords()]):
        for d in cube.coords():
            match = re.match('height.*', d.name())
            if match:
                name = match.group(0)
                vert_coord = cube.coord(name)
        ground_level = cube.extract(iris.Constraint(**{name: min(vert_coord.points)}))
        
    elif any([re.match('pressure.*', d.name()) for d in cube.coords()]):
        for d in cube.coords():
            match = re.match('pressure.*', d.name())
            if match:
                name = match.group(0)
                vert_coord = cube.coord(name)
        ground_level = cube.extract(iris.Constraint(**{name: max(vert_coord.points)}))
    
    else:
        ground_level = cube
        name = None
    
    return (ground_level, name)


def unrotate(uk_cube, g_cube):
    '''Translate uk_cube from a rotated pole grid to an unrotated grid.'''
    glat = uk_cube.coord('grid_latitude').points
    glon = uk_cube.coord('grid_longitude').points
    x, y = np.meshgrid(glon, glat)
    cs = uk_cube.coord_system()
    lons, lats = iris.analysis.cartography.unrotate_pole(x, y, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)
    clons = iris.coords.DimCoord(lons[0,:], standard_name='longitude', var_name='longitude', 
                                 units=g_cube.coord('longitude').units, coord_system=iris.coord_systems.GeogCS(6371229.0))
    clats = iris.coords.DimCoord(lats[:,0], standard_name='latitude', var_name='latitude', units=g_cube.coord('latitude').units, 
                                 coord_system=iris.coord_systems.GeogCS(6371229.0))

    uk_unrotate = uk_cube.copy()
    uk_unrotate.remove_coord('grid_latitude')
    uk_unrotate.add_dim_coord(clats, 0)
    uk_unrotate.remove_coord('grid_longitude')
    uk_unrotate.add_dim_coord(clons, 1)
    return uk_unrotate


def get_uk_from_global(g_cube, uk_cube):
    '''Extract the (unrotated) UK area from g_cube.'''
    minlat = uk_cube.coord('latitude').points[0]
    maxlat = uk_cube.coord('latitude').points[-1]
    minlon = uk_cube.coord('longitude').points[0]
    maxlon = uk_cube.coord('longitude').points[-1]
    
    g_cube_uk = g_cube.intersection(latitude=(minlat,maxlat)).intersection(longitude=(minlon,maxlon))
    return g_cube_uk


def get_uk_global_pairs(param, uk_files, g_files, uk_path='../data/mogreps-uk/', g_path='../data/mogreps-gg/'):
    '''Return two CubeLists of matching UK and global cubes, each containing one timestep, for the given param.'''
    uks = iris.cube.CubeList([])
    gs = iris.cube.CubeList([])
    
    for fname in uk_files:
        info = get_info(fname)
        uk_fname = fname
        g_fname = 'prods_op_mogreps-g_{}_{:02d}_{:02d}_{:03d}.nc'.format(info[0], info[1]-3, info[2], info[3]+3)

        if g_fname in g_files:
            try:
                g_cubes = iris.load(g_path + g_fname)
                uk_cubes = iris.load(uk_path + uk_fname)
                uk_cube = uk_cubes[[c.name() for c in uk_cubes].index(param)]
                g_cube = g_cubes[[c.name() for c in g_cubes].index(param)]
                
                uk_cube = uk_cube.extract(iris.Constraint(time=uk_cube.coord('time').points[-1]))
                g_cube = g_cube.extract(iris.Constraint(time=uk_cube.coord('time').points[-1]))
                print(uk_cube.coord('time'))
                
                uk_cube, uk_h = get_ground_level(uk_cube)
                if uk_h: uk_cube.remove_coord(uk_h)
                g_cube, g_h = get_ground_level(g_cube)
                if g_h: g_cube.remove_coord(g_h)
                    
                uk_unrotate = unrotate(uk_cube, g_cube)
                g_cube_uk = get_uk_from_global(g_cube, uk_unrotate)

                uks.append(uk_unrotate)
                gs.append(g_cube_uk)
            except (IndexError, AttributeError):
                pass
    return (uks, gs)