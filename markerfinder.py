# Written by: Anton Cecic
# Date: April 2023

from matplotlib import pyplot as plt
import skimage
from skimage.measure import profile_line
from skimage.io import imread
import numpy as np
from scipy.signal import find_peaks
from scipy import ndimage
from skimage.draw import line
import pylab as plb
from scipy.optimize import curve_fit

def gauss(x, a, x0, sigma, c):
    '''
    Returns a gaussian distribution
    
    INPUT:
        x: array of x-values
        a: max height
        x0: mean
        sigma: 
        c: height offset
    
    OUTPUT:
        Gaussian distribution array.
    '''
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + c


def gauss_fit(x, y): 
    '''
    Performs a gaussian fit:
    
    INPUT:
        x: array of x-values
        y: array of height values
    
    OUTPUT:
        popt:  optimized fit parameters
    '''
    # Try to guess approximate values for the fit parameters
    
    a_guess = max(y)
    x0_guess = sum(x * y) / sum(y)
    sigma_guess = np.sqrt(sum(y * (x - x0_guess) ** 2) / sum(y))
    c_guess = min(y)
    
    popt, pcov = curve_fit(gauss, x, y, p0=[a_guess, x0_guess, sigma_guess, c_guess])
    return popt


def GetCross(image, canvas, s, h, d, w, rh, nump, deriv):
    
    '''
    Analyzes the image and determines the centers of the markers.
    
    INPUT:
        image:     image in 2d array form
        canvas:    ...
        s:         array of coordinates for a rectangular window around markers
        h:         min height for find_peaks
        d:         min distance for find_peaks
        w:         min width for find_peaks
        rh:        relative height for find peaks
        nump:      number of peaks
        deriv:     derivative (0 or 1)
    
    OUTPUT:
        xbox:      ...
        ybox:      ...
        canvas:    ...
        
    '''

    blankcanvas = np.zeros(image.shape[0:2])
    
    x1 = s[0]
    y1 = s[1]
    x2 = s[2]
    y2 = s[3]
    
    pt1 = (x1, y1)
    pt2 = (x2, y2)

    xbox = (pt1[0], pt2[0], pt2[0], pt1[0], pt1[0])
    ybox = (pt1[1], pt1[1], pt2[1], pt2[1], pt1[1])

    for j in range(pt2[1] - pt1[1]):

        start = (pt1[1] + j, pt1[0])
        end = (pt1[1] + j, pt2[0])

        profile = skimage.measure.profile_line(image, start, end)

        prof = []
        for i in profile:
            prof.append(i[1])
            
        if deriv==1:

            prof = np.abs(np.gradient(prof))

        peaks, properties = find_peaks(prof, height=h, threshold=None, distance=d, 
                                prominence=None, width=w, wlen=None, rel_height=rh, plateau_size=None)

        if len(peaks)==nump:
           
            try:
            
                if nump==2:

                    peaks_center = int((peaks[1]-peaks[0])/2 + peaks[0])
                    arr0 = prof[0:peaks_center]
                    arr1 = prof[peaks_center:-1]

                    fit_params0 = gauss_fit(np.arange(0, len(arr0), 1), arr0)
                    fit_params1 = gauss_fit(np.arange(0, len(arr1), 1), arr1)

                    peaks = [int(fit_params0[1]), int(len(arr0)) + int(fit_params1[1])]

                    for i in range(len(peaks)):

                        blankcanvas[pt1[1] + j, pt1[0] + peaks[i]] = 100

                if nump==1:

                    fit_params0 = gauss_fit(np.arange(0, len(prof), 1), prof)
                    peak = int(fit_params0[1])
                    blankcanvas[pt1[1] + j, pt1[0] + peak] = 100
                    
            except Exception as e:
                print(e)
                
        
    for j in range(pt2[0] - pt1[0]):
        
        start = (pt1[1], pt1[0] + j)
        end = (pt2[1], pt1[0] + j)

        profile = skimage.measure.profile_line(image, start, end)

        prof = []
        for i in profile:
            prof.append(i[1])
            
        if deriv==1:

            prof = np.abs(np.gradient(prof))

        peaks, properties = find_peaks(prof, height=h, threshold=None, distance=d, 
                                prominence=None, width=w, wlen=None, rel_height=rh, plateau_size=None)


        if len(peaks)==nump:
           
            try:
            
        
                if nump==2:

                    peaks_center = int((peaks[1]-peaks[0])/2 + peaks[0])
                    arr0 = prof[0:peaks_center]
                    arr1 = prof[peaks_center:-1]

                    fit_params0 = gauss_fit(np.arange(0, len(arr0), 1), arr0)
                    fit_params1 = gauss_fit(np.arange(0, len(arr1), 1), arr1)

                    peaks = [int(fit_params0[1]), int(len(arr0)) + int(fit_params1[1])]

                    for i in range(len(peaks)):

                        blankcanvas[pt1[1] + peaks[i], pt1[0] + j] = 100


                if nump==1:

                    fit_params0 = gauss_fit(np.arange(0, len(prof), 1), prof)
                    peak = int(fit_params0[1])
                    blankcanvas[pt1[1] + peak, pt1[0] + j] = 100
                    
            except Exception as e:
                print(e)
    
    CM = ndimage.measurements.center_of_mass(blankcanvas)
    CM = np.round(CM)
    print(CM)
    
    try:
        canvas[int(CM[0]), int(CM[1])] = 100
        rr, cc = line(int(CM[0])-50, int(CM[1])-50, int(CM[0])+50, int(CM[1])+50)
        canvas[rr, cc] = 100
        rr, cc = line(int(CM[0])+50, int(CM[1])-50, int(CM[0])-50, int(CM[1])+50)
        canvas[rr, cc] = 100
    except:
        CM = [0,0]

    return xbox, ybox, canvas
            
def GetImage(path_and_name, s1, s2, s3, s4, h, d, w, rh, nump, deriv):
    
    image = skimage.io.imread(path_and_name)
    canvas = np.zeros(image.shape[0:2])
    
    xbox1, ybox1, canvas = GetCross(image, canvas, s1, h, d, w, rh, nump, deriv)
    xbox2, ybox2, canvas = GetCross(image, canvas, s2, h, d, w, rh, nump, deriv)
    xbox3, ybox3, canvas = GetCross(image, canvas, s3, h, d, w, rh, nump, deriv)
    xbox4, ybox4, canvas = GetCross(image, canvas, s4, h, d, w, rh, nump, deriv)
    
    plt.plot(xbox1, ybox1, 'r')
    plt.plot(xbox2, ybox2, 'r')
    plt.plot(xbox3, ybox3, 'r')
    plt.plot(xbox4, ybox4, 'r')


    return canvas, image
    