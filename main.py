"""
Exercise 2 Full Code

Created on Sun Nov 10 14:25:51 2024

@author: Sarah Straw

"""

#==============================================================================

# Imports

from scipy.integrate import dblquad
import matplotlib.pyplot as plt
import numpy as np

#==============================================================================

# Functions

def monte(Fresnel2dreal,Fresnel2dimag,xp1,xp2,N_samples, _x, _y, R):
   
    """
    Function
    ----------
    Calculates the intensity of diffracted waves using the monte carlo
    integration method to solve the Fresnel diffraction equation.
    
    
    Parameters
    ----------
    Fresnel2dreal : Function
        The real part of function to be integrated.
    Fresnel2dimag : Function
        The imaginary part of function to be integrated.
    xp1 : float
        Minimum value of range for integration.
    xp2 : float
        Maximum value of range for integration.
    N _samples : int
        Number of MC points.
    _x : float
        x value of the point on the screen intensity is being evalucated at
    _y : float
        y value of the point on the screen intensity is being evalucated at
    R : float
        radius of circular aperture
        
    Returns
    -------
    intensity : float
        The value of the intensity
    error : float
        The absolute error in the intensity.

    """
    
    # mask set up
    
    samples_xp = np.random.uniform(low=xp1,high=xp2,size=N_samples) # random x values of samples on the aperture
    samples_yp = np.random.uniform(low=yp1,high=yp2,size=N_samples) # random y values of samples on the aperture
    samples_rp = np.sqrt( samples_xp**2 + samples_yp**2 ) # array of displacements of samples from the origin 
    
    mask = samples_rp < R # mask for all sample points within circular aperture
    
    N_mask = sum(mask) # now there are N_mask number of sample points we will use, not N_samples
    
    good_samples_xp = samples_xp[mask] # all values in samples_rp which align for TRUE in mask boolian will be taken into good_samples_xp 
    good_samples_yp = samples_yp[mask]
    
    # now only points within circle are passed through function so no need for a loop 'if', 'else' later
    
    # sum of the real and imaginary parts of the function within the electric field integral
    real_values_sum = sum(Fresnel2dreal(good_samples_yp, good_samples_xp, _y, _x, k, z))
    imag_values_sum = sum(Fresnel2dimag(good_samples_yp, good_samples_xp, _y, _x, k, z))
    
    real_mean = real_values_sum / N_mask # mean of real parts of function across all sample points
    imag_mean = imag_values_sum / N_mask # mean of imaginary parts of function across all sample points
    
    area_aperture = np.pi*R**2 # m^2, area of the aperture
    
    real_E = coeff*area_aperture*real_mean # real parts of electric field following the MC integral formula
    imag_E = coeff*area_aperture*imag_mean # imaginary parts of electric field following the MC integral formula
    
    intensity = ep0*c*( real_E**2 + imag_E**2 ) # intensity calculation 
   
    # error calculation 
        
    squared_real_values_sum = sum((Fresnel2dreal(good_samples_yp, good_samples_xp, _y, _x, k, z))**2)
    squared_imag_values_sum = sum((Fresnel2dimag(good_samples_yp, good_samples_xp, _y, _x, k, z))**2)
    
    meansq = (squared_real_values_sum + squared_imag_values_sum) / N_mask
    
    # done for intensity, cut out middle man of electric field array, now compute <f**2> from real mean and imag mean

    error = (area_aperture) * np.sqrt(( meansq- (real_mean**2 + imag_mean**2) ) / N_mask )
    
    return intensity, error

def Fresnel2dreal(yp, xp, y, x, k, z): # real part of inner integral

    """
    Function
    ----------
    Real component on function being integrated in the Fresnel diffraction equation
    
    Parameters
    ----------
    yp : float
        y coordinate of point on aperture being evaluated
    xp : float
        x coordinate of point on aperture being evaluated
    y : float
        y coordinate of point on screen being evaluated
    x : float
        x coordinate of point on screen being evaluated
    k : float
        2 pi over wavelength
    z : float
        distance from screen to aperture
        
    Returns
    -------
    function within Fresnel diffraction integral
    
    """
    
    return np.cos((k/(2*z)) * ( (x-xp)**2 + (y-yp)**2 ))


def Fresnel2dimag(yp, xp, y, x, k, z): # imaginary part

    """
    Function
    ----------
    Imaginary component on function being integrated in the Fresnel diffraction equation
    
    Parameters
    ----------
    yp : float
        y coordinate of point on aperture being evaluated
    xp : float
        x coordinate of point on aperture being evaluated
    y : float
        y coordinate of point on screen being evaluated
    x : float
        x coordinate of point on screen being evaluated
    k : float
        2 pi over wavelength
    z : float
        distance from screen to aperture
        
    Returns
    -------
    function within Fresnel diffraction integral
    
    """
    
    return np.sin((k/(2*z)) * ( (x-xp)**2 + (y-yp)**2 ))


def yp1func(xp): 

    """
    Function
    ----------
    Function of inner integral y limits for a circular aperture
    
    Parameters
    ----------
    xp : float
        x coordinate of point along aperture being evaluated
        
    Returns
    -------
    yp lower limit to inner integral
    
    """

    return - np.sqrt(R*R - xp*xp)

def yp2func(xp): 
    
    """
    Function
    ----------
    Function of inner integral y limits for a circular aperture
    
    Parameters
    ----------
    xp : float
        x coordinate of point along aperture being evaluated
        
    Returns
    -------
    yp upper limit to inner integral
    
    """

    return np.sqrt(R*R - xp*xp)

#=================================================================================================================

# constants 

ep0 = 8.854187817e-12 # Fm^-1, permativity of free space
c = 299792458 # ms^-1, speed of light
E0 = 1 # NC^-1, electric field of incident light on the aperture 

k = 2e6 # m, wavenumber = 1 / wavelength
z = 0.175 # m, distance from screen to apature 

xp1 = -1e-5 # m,  limits of the aperture centered 
yp1 = xp1
xp2 = -xp1
yp2 = -xp1


coeff = (k*E0) / (2*(np.pi)*z) # coefficient of real and imaginary parts outside of integral

#=================================================================================================================

# Menu for 4 parts of exersise

MyInput = '0'
while MyInput != 'q':
    
    print('')
    print('"1" for part (1): 1-d diffraction from 2-d aperture')
    print('"2" for 2-d diffraction from 2-d rectangular aperture')
    print('"3" for 2-d diffraction from 2-d circular aperture')
    print('"4" for 2-d circular diffraction using Monte Carlo')
    print('or "q" to quit')
    MyInput = input('Enter a choice here:')
        
    print('You entered the choice: ',MyInput)
    
    if MyInput == '1':
        print('You have chosen part (1): 1-d diffraction from 2-d aperture')
        
        #=================================================================================================================
        
        # constants and variables
        
        k = 2e6 # m, wavenumber = 1 / wavelength

        z = 0.05 # m, distance from screen to apature 
        
        xp1 = -1e-5 # m,  limits of the aperture centered 
        yp1 = xp1
        xp2 = -xp1
        yp2 = -xp1
        
        num = 200 # number of data points for the 1-d diffraction pattern
        
        #=================================================================================================================
        
        # array set up

        realparts = np.zeros(num) # to store the real part of the solved double integral
        imagparts = np.zeros(num)

        I_xarray = np.zeros(num) # to store the array of intensity values along the cross section

        real_error = np.zeros(num) # to store real error in the integral
        imag_error = np.zeros(num) # to store imaginary error in the integral

        err = np.zeros(num) # to store error in intensity across the cross section

        #=================================================================================================================

        # main code

        # Iterate over each value in xarray, calculate the integral

        y = 0 # 1-d diffraction pattern is cross section of intensity at a specific y value, 

        xarray = np.linspace(-4e-2, 4e-2, num) # array of x values to solve double integral to


        for i in range(num): # loop for every data point

            x = xarray[i] # running through xarray setting x to each value in turn
            realpart, realerr = dblquad(Fresnel2dreal, xp1, xp2, yp1, yp2, args=(y, x, k, z)) # solving integral with the specific x value
            realparts[i] = realpart  # Store the result in thearray realparts

            real_error[i] = realerr

            imagpart, imagerr = dblquad(Fresnel2dimag, xp1, xp2, yp1, yp2, args=(y, x, k, z))
            imagparts[i] = imagpart # Store the result in the array imagparts
            imag_error[i] = imagerr

            #use these to calculate intensity at each x value and store in array I_xarray

            I = ep0*c*((coeff*realpart)**2 + (coeff*imagpart)**2) # intensity equation from real and imaginary parts
            I_xarray[i] = I # storing intensity values 


            err_ = ep0*c*((coeff*realerr)**2 + (coeff*imagerr)**2) # intensity equation from real and imaginary parts
            err[i] = err_

        #=================================================================================================================

        # plotting results

        plt.plot(xarray, I_xarray)
        plt.xlabel("Screen Coordinate (m)")
        plt.ylabel("Relative Intensity (Wm$^{-2}$)")
        plt.title('1D Diffraction Pattern')
        plt.show()
        
    elif MyInput == '2':
        
        print('You have chosen part (2): 2-d diffraction from 2-d rectangular aperture')
        print('Plot loading time: 1.5 minutes')
        
        #=================================================================================================================
        
        # constants and variables
                
        k = 2e6 # m, wavenumber = 1 / wavelength
        z = 0.1 # m, distance from screen to apature 
        
        #num = 100 # number of data points(smooth)
        
        x1=-0.015 # m,  x and y limits for the plot
        x2=-x1
        y1=0.005
        y2=-y1
        
        num=100
        
        #==============================================================================
        
        # array set up
        
        # empty arrays to store results of integral for the array of x and y values
        real_parts = np.zeros((num,num)) 
        imag_parts = np.zeros((num,num))
        intensity = np.zeros((num,num))
        
        # array of x and y values to integrate with 
        xarray = np.linspace(-4e-2, 4e-2, num) # array of x values to solve double integral to
        yarray = np.linspace(-4e-2, 4e-2, num) # array of y values to solve double integral to
        
        #==============================================================================
        
        # main code
        
        for i in range(num): # loop for every data point
        
            x = xarray[i] # running through xarray setting x to each value in turn
            
            for j in range (num):
                
                y = yarray[j]
                
                # solving integral for array of y values in yarray[j] at specifc x value xarray[i]
                # storing results in array where real_parts[:,0] is from array of x values and real_parts[:,1] is from array of y values
                
                real_parts[i,j], realerr = dblquad(Fresnel2dreal, xp1, xp2, yp1, yp2, args=(y, x, k, z))
                
                # repeat for imaginary part of the equation
                
                imag_parts[i,j], imagerr = dblquad(Fresnel2dimag, xp1, xp2, yp1, yp2, args=(y, x, k, z))
                
                # solve for intensity and store into an array
                
                intensity[i,j] = ep0*c*((coeff*real_parts[i,j])**2 + (coeff*imag_parts[i,j])**2)
        
        #==============================================================================
        
        # plotting results
        
        extents = (x1,x2,y1,y2) # Sets the limits for the plot
        plt.imshow(intensity,vmin=0.0,vmax=1.0*intensity.max(),extent=extents,\
                   origin="lower",cmap="gist_ncar_r") # Try different
        #plt.xlabel("x (m)")
        #plt.ylabel("y (m)")
        #plt.title('Rectangular diffraction;\nz = {:4.2f}m'.format(z))
        plt.figure(figsize=(6,1))
        plt.colorbar()
        plt.show()
        
        
    elif MyInput == '3':
        
        print('You have chosen part 3): 2-d diffraction from 2-d circular aperture')
        print('Plot loading time: 3.5 minutes')
        
        R = 2e-5 # radius of aperture

        xp1 = -R # limits of aperture in terms of radius
        xp2 = R

        #==============================================================================

        # arguments and variables

        k = 2e6 # m, wavenumber = 1 / wavelength
        z = 0.175 # m, distance from screen to apature 

        coeff = (k*E0) / (2*(np.pi)*z) # coefficient of real and imaginary parts outside of integral

        num = 200  # number of data points(smooth)

        x1=-0.01  # m, x and y limits for the plot
        x2=-x1
        y1=x1
        y2=x2

        #==============================================================================

        # array set up

        real_parts = np.zeros((num,num)) #empty arrays to store results of integral for the array of x and y values
        imag_parts = np.zeros((num,num))
        intensity = np.zeros((num,num))

        # array of x and y values to integrate with 

        xarray = np.linspace(-4e-2, 4e-2, num) # array of x values to solve double integral to
        yarray = np.linspace(-4e-2, 4e-2, num) # array of y values to solve double integral to

        #==============================================================================

        # main code

        # loop of x and y values 

        for i in range(num): # loop for every data point

            x = xarray[i] # running through xarray setting x to each value in turn
            
            for j in range (num):
                
                y = yarray[j]
                
                # solving integral for array of y values in yarray[j] at specifc x value xarray[i]
                
                # storing results in array where real_parts[:,0] is from array of x values and real_parts[:,1] is from array of y values
                
                real_parts[i,j], realerr = dblquad(Fresnel2dreal,xp1,xp2,yp1func,yp2func,args=(y, x, k, z), epsabs=1e-10, epsrel=1e-10)

                # repeat for imaginary part of the equation

                imag_parts[i,j], imagerr = dblquad(Fresnel2dimag,xp1,xp2,yp1func,yp2func,args=(y, x, k, z), epsabs=1e-10, epsrel=1e-10)

                # solve for intensity and store into an array

                intensity[i,j] = ep0*c*((coeff*real_parts[i,j])**2 + (coeff*imag_parts[i,j])**2)

        #==============================================================================

        # plotting results

        extents = (x1,x2,y1,y2) # Sets the limits for the plot
        plt.imshow(intensity,vmin=0.0,vmax=1.0*intensity.max(),extent=extents,\
                   origin="lower",cmap="gist_ncar_r") 
        #plt.xlabel("x (m)")
        #plt.ylabel("y (m)")
        #plt.title('Circular aperture diffraction;\nz = {:4.2f}m'.format(z))
        plt.figure(figsize=(6,1))

        plt.colorbar()
        plt.show()
        
    elif MyInput =='4':
        
        print('You have chosen part (4): 2-d circular diffraction using Monte Carlo')
        print('Plot loading time: 2 minutes')
        
        #=================================================================================================================
        
        # constants and variables
        
        # x and y limits for the plot

        x1=-0.01
        x2=-x1
        y1=x1
        y2=x2

        R = 2e-5 # m, radius of aperture

        # limits of aperture

        xp1 = -R
        xp2 = R
        yp1 = -R
        yp2 = R

        #==============================================================================

        # arguments and variables

        x1=-0.01 # x and y limits for the plot
        x2=-x1
        y1=x1
        y2=x2

        R = 2e-5 # m, radius of aperture

        xp1 = -R # limits of aperture
        xp2 = R
        yp1 = -R
        yp2 = R

        k = (np.pi*2) / 5e-7 # m, wavenumber = 2pi / wavelength
        z = 0.8 # m, distance from screen to apature 

        coeff = (k*E0) / (2*(np.pi)*z) # coefficient of real and imaginary parts outside of integral

        N_samples = 3000 # number of samples per point on screen, the higher this is the more acurate the intensity is 

        N = 200 # number of points along x and y for 40000 points on screen (smooth plot)

        #==============================================================================

        # array set up

        intensity_array = np.zeros((N,N)) # empty array to store results of integral for the array of x and y values
        error_array = np.zeros((N,N))

        xarray = np.linspace(-4e-2, 4e-2, N) # array of x values to evaluate the Fresnel expression to 
        yarray = np.linspace(-4e-2, 4e-2, N) # array of y values to evaluate the Fresnel expression to

        #==============================================================================

        # main code

        for i in range(N): # double loop for every point on the screen from array of x and y values
            
            x = xarray[i]
            
            for j in range(N):
                
                y = yarray[j] 
                
                # for point on screen (xarray[i], yarray[j]):
                
                intensity , error = monte(Fresnel2dreal,Fresnel2dimag,xp1,xp2,N_samples, x, y, R)
                intensity_array[i,j] = intensity
                error_array[i,j] = error
               
        #==============================================================================

        # plotting results
            
        plt.imshow(intensity_array,\
                   origin="lower",cmap="nipy_spectral_r") # Try different
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.title('Diffraction of a Circular Aperture MC integration;\nz = {:4.2f}m'.format(z))
        plt.colorbar()
        plt.show()
        
        
    elif MyInput != 'q':
        print('This is not a valid choice')
print('You have chosen to finish - goodbye.')
