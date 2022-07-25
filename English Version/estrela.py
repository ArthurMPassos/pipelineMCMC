__author__ = "Adriana Valio, Beatriz Duque, Felipe Pereira Pinho"
__copyright__ = "..."
__credits__ = ["Universidade Presbiteriana Mackenzie, CRAAM"]
__license__ = ""
__version__ = ""
__maintainer__ = ""
__email__ = "biaduque7@hotmail.com"
__status__ = "Production"

'''
This program simulates the plot of a star with spots, through parameters such as radius, intensity, dimming
of limbo, etc.
The imported libraries are:
math
matplotlib
numpy
verify: function created to validate inputs, for example non-float/int or negative numbers
'''


import math
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from verify import Validar
import random
from ctypes import *
from numpy.ctypeslib import ndpointer
import time
import sys
import platform


class Estrela:
    '''
    The star class receives as object the radius, maximum intensity, limb darkening coefficients.
    The star is formatted in an array of default size 856.
    The parameters passed to the spot are objects belonging to the class, such as: radius, intensity, longitude and latitude
    in relation to the star.
    ************ STAR PARAMETERS***************
    :raio parameter: The radius of the star in pixels
    :raioSun parameter: The radius of the star in units of Rsun
    :intensidadeMaxima parameter: Intensity of the center of the star
    :coeficienteHum parameter: Limbo darkening coefficient
    :coeficienteDois parameter: Limbo darkening coefficient
    :tamanhoMatriz parameter: Size of the matrix in which the star will be built
    :estrela parameter: Star constructed with the limbo darkening coefficients
    '''
   

    def __init__(self,raio,raioSun,intensidadeMaxima,coeficienteHum,coeficienteDois,tamanhoMatriz):
        
        self.raio=raio #em pixel
        self.raioSun = raioSun
        self.intensidadeMaxima=intensidadeMaxima
        self.coeficienteHum=coeficienteHum
        self.coeficienteDois=coeficienteDois
        self.tamanhoMatriz=tamanhoMatriz
        #self.colors = ["gray","pink","hot"]
        error=0
        
        #start = time.time()
        # Verifica se o Python é 32 ou 64 bit
        if(platform.architecture()[0] == "32bit"):
            my_func = WinDLL('scripts/func32.dll', winmode = 0x8)
        elif(platform.architecture()[0] == "64bit"):
            my_func = WinDLL('scripts/func64.dll', winmode = 0x8)

        my_func.criaEstrela.restype = ndpointer(dtype=c_int, ndim=2, shape=(self.tamanhoMatriz,self.tamanhoMatriz))
        self.estrela = my_func.criaEstrela(self.tamanhoMatriz,self.tamanhoMatriz,self.tamanhoMatriz,c_float(self.raio),c_float(self.intensidadeMaxima),c_float(self.coeficienteHum),c_float(self.coeficienteDois))

        del my_func
        #print(self.estrela)

        self.error=error
        self.Nx = self.tamanhoMatriz
        self.Ny = self.tamanhoMatriz
        #self.color = random.choice(self.colors)
        self.color = "hot"
        #Plotar(self.tamanhoMatriz,self.estrela)
        #end = time.time()
        #print(end - start)
    
    #######  inserir manchas


    def manchas(self,r,intensidadeMancha,lat,longt):
        '''
        Function where the star spot(s) is created. all parameters
         are related to the size of the star, the user being able to choose values
         or select the default option.
         ********START OF STAIN PARAMETERS*******
         :raioMancha parameter: Radius of the spot in relation to the radius of the star
         :intensidadeMancha parameter: Spot intensity as a function of the star's maximum intensity
         :latitudeMancha parameter: Latitude coordinate of the spot in relation to the star
         :longitudeMancha parameter: Longitude coordinate of the spot in relation to the star
        '''
        # Spots parameter
        # #r=0.05 (for tests)
        # intensidadeMancha=0.5 (for tests)
        # spot coordinates in degrees
        #test latitude=-30
        #test longitude=20

        self.raioMancha = self.raio * r  # radius as a function of the star's radius in pixels
        self.intensidadeMancha = intensidadeMancha # spot intensity as a function of the star's maximum intensity

        #spot positioning coordinates in degrees


        degreeToRadian = np.pi/180. #A read-only variable containing the floating-point value used to convert degrees to radians.
        self.latitudeMancha  = lat * degreeToRadian 
        self.longitudeMancha =  longt * degreeToRadian

        #position of the spot in pixels with respect to the center of the star
        ys=self.raio*np.sin(self.latitudeMancha)  
        xs=self.raio*np.cos(self.latitudeMancha)*np.sin(self.longitudeMancha)
        anguloHelio=np.arccos(np.cos(self.latitudeMancha)*np.cos(self.longitudeMancha))

        
        # projection effect by the spot being at a heliocentric angle from the center of the star - ellipticity
        yy = ys + self.Ny/2 #pixel position with respect to the matrix origin
        xx = xs + self.Nx/2 

        kk = np.arange(self.Ny * self.Nx)
        vx = kk-self.Nx*np.int64(1.*kk/self.Nx) - xx
        vy = kk/self.Ny - yy

        # spot rotation angle
        anguloRot=np.abs(np.arctan(ys/xs))    # in radians
        if self.latitudeMancha*self.longitudeMancha > 0: anguloRot=-anguloRot

        ii, = np.where((((vx*np.cos(anguloRot)-vy*np.sin(anguloRot))/np.cos(anguloHelio))**2+(vx*np.sin(anguloRot)+vy*np.cos(anguloRot))**2) < self.raioMancha**2)
        
        spot = np.zeros(self.Ny * self.Nx) + 1
                
        spot[ii]=self.intensidadeMancha
        spot = spot.reshape([self.Ny, self.Nx])
    
        self.estrela= self.estrela * spot
        plt.axis([0,self.Nx,0,self.Ny])
                
        #self.estrelaManchada= estrelaManchada
        
        #Plotar(self.tamanhoMatriz,self.estrela)
        error=0
        self.error=error
        return self.estrela #returns the decision: whether there are stains or not

    #inserção de flares
    def faculas(self,estrela,count): 
        
        #receives as parameter the updated star
        #not implemented 
        '''
        Function where the faculae of the star are created. all parameters
        are related to the size of the star, and the user can choose values
        or select the default option.
        ---Parameters not yet defined
        ********START OF FACULA PARAMETERS*******
        :parameter
        :parameter
        :parameter
        :parameter
        '''
        error=0
        self.error=error
        self.estrela=estrela
        #returns the decision: whether there is facula or not
        return self.estrela
    
    def flares(self,estrela,count): 
        #receives as parameter the updated star
        #not implemented
        '''
        Function where the star's flares are created. all parameters
        are related to the size of the star, and the user can choose values
        or select the default option.
        ---Parameters not yet defined
        ********START OF FLARES PARAMETERS*******
        :parameter
        :parameter
        :parameter
        :parameter
        '''
        error=0
        self.error=error
        self.estrela=estrela
        #returns the decision: whether there is flares or not
        return self.estrela


    def getNx(self):
        '''
        Returns Nx parameter, required by Eclipse.        
        '''
        return self.Nx
    def getNy(self):
        '''
        Returns Ny parameter, required by Eclipse.
        '''
        return self.Ny

    def getRaioStar(self):
        '''
        Returns the radius of the star in pixels, necessary for the Eclipse program, since the radius of the planet is given in
        with respect to the radius of the star.
        '''
        return self.raio
    def getEstrela(self):
        '''
        Returns the star, plotted without the spots, necessary if the user chooses the plot without spots.
        '''
        return self.estrela

    def getu1(self):
        return self.coeficienteHum
    
    def getu2(self):
        return self.coeficienteDois

    def getTamanhoMatriz(self):
        return self.tamanhoMatriz
    
    def getRaioSun(self):
        return self.raioSun

    def getIntensidadeMaxima(self):
        return self.intensidadeMaxima

    def getError(self):
        '''
        Returns error value. If there are no errors, the variable will assume 0. If there are errors, the program will keep
        the source value of the variable (which is -1).
        '''
        return self.error
    def setStarName(self,starName):
        self.starName = starName

    def getStarName(self):
        return self.starName

    def setCadence(self,cadence):
        self.cadence = cadence

    def getCadence(self):
        return self.cadence

    def Plotar(self,tamanhoMatriz,estrela):
        Nx = tamanhoMatriz
        Ny = tamanhoMatriz
        plt.axis([0,Nx,0,Ny])
        plt.imshow(estrela,self.color)
        plt.show()