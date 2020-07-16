import numpy as np
import scipy
from scipy.interpolate import interp2d

from tqdm import tqdm 

def trace_image_between_rows(row_start,row_end, image):
    return np.sum(image[row_start:row_end,:], axis=0)

def trace_image_between_cols(col_start,col_end, image):
    return np.sum(image[:,col_start:col_end], axis=1)

def trace_image_between_lines(gradeint, upper_line_intersection, lower_line_intersection,   image, vertical_oversample=10):
    dc = upper_line_intersection - lower_line_intersection
    dist = dc / np.sqrt(gradeint**2 + 1)
    xs = np.sqrt(dc**2 - dist**2)

    # Now need to trace the top line's start and end positions
    top_start = [0,upper_line_intersection]

    xstep = np.sqrt(xs**2 / (gradeint**2 + 1))
    top_end = [image.shape[1]  - xstep, upper_line_intersection + gradeint*(image.shape[1]  - xstep)]
    R = int(np.floor(np.hypot(top_end[0] - top_start[0], top_end[1] - top_start[1])))

    # Define the arrays
    xtrace = np.arange(R)
    flux_trace = np.arange(R)

    x_image, y_image = np.zeros(R), np.zeros(R)
    dx = np.sqrt(1 / (gradeint**2 + 1))
    dy = np.sqrt(1 / ((1/gradeint**2) + 1))

    # Coordinates of the top line
    position_top_line = np.transpose(( np.arange(top_start[0],top_end[0], dx ), np.arange(top_start[1],top_end[1], dy )  ))

    # calculate the corrosponding coordinates of the lower line
    position_lower_line = np.copy(position_top_line)
    position_lower_line[:,0] = position_top_line[:,0] + np.sqrt(dist**2 / (1 / (gradeint**2) + 1))
    position_lower_line[:,1] = position_top_line[:,1] - np.sqrt(dist**2 / (1 + gradeint**2))

    positions_mask = (position_top_line[:,0] > image.shape[1]-1) | (position_lower_line[:,0] > image.shape[1]-1) | (position_top_line[:,1] > image.shape[0]-1) | (position_lower_line[:,1] > image.shape[0]-1)
    position_top_line = position_top_line[~positions_mask]
    position_lower_line = position_lower_line[~positions_mask]

    # Now we have the lines, create th eimage interpolater
    image_interpolater = interp2d(np.arange(image.shape[1]), np.arange(image.shape[0]), image)

    # Interpolate the line joingin the upper line and the lower line preserving the angle
    # Making sure to conserve flux
    flux = np.zeros(position_top_line.shape[0])
    for i in range(flux.shape[0]): 
        x = np.linspace(position_top_line[i,0], position_lower_line[i,0], vertical_oversample )
        y = np.linspace(position_top_line[i,1], position_lower_line[i,1], vertical_oversample )
        flux[i] = np.sum(image_interpolater(x,y))*dist / vertical_oversample

    return flux


















    '''
    f = plt.figure()
    plt.imshow(image,  aspect='auto', norm=norm, origin='lower', cmap='Greys')
    plt.scatter(*position_top_line.T, c='r')
    plt.scatter(*position_lower_line.T, c='b')
    for i in range(position_lower_line.shape[0])[::10]: plt.plot([position_top_line[i,0],position_lower_line[i,0] ],   [position_top_line[i,1],position_lower_line[i,1] ], c='r')
    '''



    '''
    f = plt.figure()
    plt.plot(np.arange(flux.shape[0]), flux)
    plt.show()
    '''
'''
#zi = scipy.ndimage.map_coordinates(z, np.vstack((x,y)))

from astropy.io import fits 
import matplotlib.pyplot as plt 
from astropy.visualization import simple_norm

image = fits.open('a1341001.fits')[0].data
norm = simple_norm(image, 'sqrt', min_percent=5, max_percent=95)

gradeint = 0.001
upper_line_intersection = 72
lower_line_intersection = 66
star_flux = trace_image_between_lines(gradeint, upper_line_intersection, lower_line_intersection,   image, vertical_oversample=100)
upper_sky = trace_image_between_lines(gradeint, upper_line_intersection+6, lower_line_intersection+6,   image, vertical_oversample=100)
lower_sky = trace_image_between_lines(gradeint, upper_line_intersection-6, lower_line_intersection-6,   image, vertical_oversample=100)

mean_sky = (upper_sky + lower_sky) / 2

f = plt.figure()
plt.plot(np.arange(star_flux.shape[0]), star_flux -mean_sky )
plt.show()
'''