import numpy as np
import scipy
from scipy.interpolate import interp2d
from scipy.signal import find_peaks
from tqdm import tqdm 
from lmfit.models import GaussianModel

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

def walk_to_threashold(x, start, threshold, step=1):
    while (start > 0) and (start < x.shape[0]-1):
        start += step
        if x[start] < threshold : break
    return start

def find_science_image_bounds(gradeint, upper_line_intersection, lower_line_intersection,   image, vertical_oversample=10, horizontal_oversample=10):
    gradeint = np.clip(gradeint, 1e-5, np.inf)
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
    print(top_start, top_end,dy, dx)
    position_top_line = np.transpose(( np.arange(top_start[0],top_end[0], dx ), np.arange(top_start[1],top_end[1], dy )  ))
    print(position_top_line)

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
    yoffsets = np.linspace(0, upper_line_intersection-lower_line_intersection, vertical_oversample)
    flux = np.zeros(yoffsets.shape[0])

    import matplotlib.pyplot as plt
    for i in range(flux.shape[0]): 
        #plt.plot(position_top_line.T[0], position_top_line.T[1] - yoffsets[i], c='r')
        flux[i] = np.sum(image_interpolater(position_top_line.T[0], position_top_line.T[1] - yoffsets[i]))
    flux = flux - np.percentile(flux, 20)
    flux = flux / np.max(flux)


    #plt.figure()
    #plt.scatter(yoffsets, flux)


    peaks, meta = find_peaks(flux, distance = 5, height = np.percentile(flux, 70))
    #plt.plot(yoffsets[peaks], flux[peaks], 'r+')

    fitted_offets = []
    for peak in peaks:
        gaus = GaussianModel()
        pars = gaus.guess(flux, x=yoffsets, center=yoffsets[peak])


        out = gaus.fit(flux, pars, x=yoffsets)

        fitted_offets.append([out.params['center'] - 1*out.params['fwhm'], out.params['center'] + 1*out.params['fwhm']])

    #    plt.plot(yoffsets, out.init_fit, 'r')
    #    plt.plot(yoffsets, out.best_fit, 'g')

    #for fitted_offet in fitted_offets:
    #    plt.axvline(fitted_offet[0], c='r')
    #    plt.axvline(fitted_offet[1], c='r')

    #plt.show()

    fitted_c_values = [[ upper_line_intersection - fitted_offet[0], upper_line_intersection - fitted_offet[1]] for fitted_offet in fitted_offets]
    return fitted_c_values














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


def plot_rectangle(x, y, ax, c='r'):
    # Plot a rectangle 
    # assume corrdinates are
    # (x1,y1)     (x2,y2)
    #
    # (x4,y4)     (x3,y3)
    for i in range(x.shape[0]) : 
        i1 = (i+1) % (x.shape[0])
        ax.plot([x[i],x[i1]], [y[i],y[i1]], c=c)

def trace_rectangle(image, x, y, upsamplex=5, upsampley=5, ax = None):
    # Trace a rectangle over an image
    # assume corrdinates are
    # (x1,y1)     (x2,y2)
    #
    # (x4,y4)     (x3,y3)

    # The idea is that we will start at (x1,y1) and work towards (x2,y2)
    # First, we will create the upsampled top lines
    Nsample = int(upsamplex*(x[1] - x[0]))
    Nsampley = int(upsampley*(y[1]-y[3]))
    top_line_x = np.linspace(x[0], x[1], Nsample)
    top_line_y = np.linspace(y[0], y[1], Nsample)

    # Then we will create the upsampled bottom line
    bottom_line_x = np.linspace(x[3], x[2], Nsample)
    bottom_line_y = np.linspace(y[3], y[2], Nsample) 

    if ax is not None:
        for i in range(top_line_x.shape[0]) : ax.plot([top_line_x[i], bottom_line_x[i]], [top_line_y[i], bottom_line_y[i]], 'y')

    # Now create the image interpolater
    image_interpolater = interp2d(np.arange(image.shape[1]), np.arange(image.shape[0]), image)

    # Now interpolate
    traced_flux = np.zeros(top_line_x.shape[0])
    for i in range(top_line_x.shape[0]) : 
        # Noew construct line linking the relavent coordinates
        x_ = np.linspace(top_line_x[i], bottom_line_x[i], Nsampley)
        y_ = np.linspace(top_line_y[i], bottom_line_y[i], Nsampley)
        traced_flux[i] = np.sum(image_interpolater(x_,y_)) / x_.shape
    return top_line_x, traced_flux


'''
import matplotlib.pyplot as plt
image = np.ones((100,100))
image[45:55,45:50] = 2
plt.imshow(image, aspect='auto', origin='lower')
x = np.array([20,60, 65, 25])
y = np.array([60,65,35,30])
plot_rectangle(x, y, plt.gca(), c='r')
trace_rectangle(image, x, y, upsamplex=2, upsampley=5, ax = plt.gca())
plt.show()
'''