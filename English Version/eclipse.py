__author__ = "Adriana Valio, Beatriz Duque, Felipe Pereira Pinho"
__copyright__ = "..."
__credits__ = ["Universidade Presbiteriana Mackenzie, CRAAM"]
__license__ = ""
__version__ = ""
__maintainer__ = ""
__email__ = "biaduque7@hotmail.com"
__status__ = "Production"

'''
Program that simulates the eclipse and the light curve of a planet when transiting
your host star.
In this program, the light curve of the star is calculated in relation to the parameters of the planet added
by the user.
***Imported Libraries***
numpy: big calculations
matplotlib: gaphic plots
star: program file where the star parameters are calculated, given user inputs (radius, intensity, etc.)
verify: function created to validate inputs, for example non-float/int or negative numbers
'''

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import pyplot
from estrela import Estrela
from moon import Moon
from verify import Validar
#from keplerAux import keplerfunc  #helper library in case the kepler library doesn't work
import matplotlib.animation as animation
import kepler # for calculating eccentric orbits (pip install kepler)
import os
from ctypes import *
from numpy.ctypeslib import ndpointer
import time
import gc
import sys
import platform

class Eclipse:

   
    def __init__(self,Nx,Ny,raioEstrelaPixel,estrelaManchada):
        

        '''
        : parameter Nx and Ny: size of the star array
        : RayStarPixel parameter: radius of the star in pixels
        : Smearedstar parameter: STAR object passed as Smudgedstar after inserting spots
        '''
        self.Nx = Nx
        self.Ny = Ny
        self.raioEstrelaPixel = raioEstrelaPixel
        self.estrelaManchada = estrelaManchada
        
        # OUTPUT
        curvaLuz =[ 1.0 for i in range(self.Nx)]
        self.curvaLuz = curvaLuz

    def geraTempoHoras(self):
        '''
        Function called in Main to calculate the Transit time in Hours
        '''
        x=int(input("Time range=1. Do you want to change? 1. YES | 2. NO:"))
        if x ==1:
            self.intervaloTempo=float(input('Enter the time range in minutes:'))
        elif x==2:
            self.intervaloTempo = 1.   # in Minuts


        self.tamanhoMatriz= self.Nx #Nx or Ny
        tempoHoras = (np.arange(self.tamanhoMatriz)-self.tamanhoMatriz/2)*self.intervaloTempo/60.   # in Hours
        self.tempoHoras= tempoHoras

    def setTempoHoras(self,intervalo):
        self.intervaloTempo = intervalo   # in Minuts
        self.tamanhoMatriz= self.Nx #Nx or Ny
        tempoHoras = (np.arange(self.tamanhoMatriz)-self.tamanhoMatriz/2)*self.intervaloTempo/60.   # in Hours
        self.tempoHoras= tempoHoras

    #from the moment lua is instantiated in main, these objects become objects of the class with 'self'.
    def criarLua(self, raioM, massM, raioPlanetaPixel, raioStar,tempoHoras,anguloInclinacao,periodo,distancia):
        moon = Moon(raioM, massM, self.raioEstrelaPixel,anguloInclinacao ,periodo, raioPlanetaPixel, self.tempoHoras,distancia)
        moon.moonOrbit(raioStar)
        Rmoon = moon.getRmoon()

        #gathering necessary data to eclipse construction
        self.xxm = moon.getxm()
        self.yym = moon.getym()
        self.Rmoon = Rmoon #em pixel 
        self.massM = massM
        self.tamanhoMatriz= self.Nx
        #gathering moon data
        self.ppMoon = moon.getppMoon(self.tamanhoMatriz)
        self.xl = moon.getxl()
        self.yl = moon.getyl()
        return moon
        


    def criarEclipse(self,semiEixoRaioStar,semiEixoUA, raioPlanetaRstar, raioPlanJup ,periodo,anguloInclinacao,lua,ecc,anom,anim=True,plot=True):

        '''
        Creation of the eclipse class, which will return the light curve of the planet's transit around the star

        ****parameters assigned to the planet****
        :periodo parameter: rotation period of the planet
        :semiEixoRaioStar parameter: semi axis of the planet in relation to the radius of the star
        :semiEixoUA parameter: semiaxis of the planet in UA
        :anguloInclinacao parameter: tilt angle of the planet
        :raioPlanetaRstar parameter: radius of the planet in relation to the radius of the star
        :raioPlanJup parameter: radius of the planet in relation to the radius of Jupiter
        :lua parameter: moon orbiting the planet (enters as True or False)
        :ecc parameter: eccentricity of the planet's orbit
        :anom parameter: planet orbit anomaly
        :anim parameter: checks if the animation will be shown to the user (True by default)
        :plot parameter: checks if the light curve graph will be shown to the user (True by default)
        '''

        intervaloTempo = self.intervaloTempo
        tamanhoMatriz = self.tamanhoMatriz
        self.semiEixoRaioStar = semiEixoRaioStar
        self.semiEixoUA = semiEixoUA
        self.raioPlanetaRstar = raioPlanetaRstar
        self.raioPlanJup = raioPlanJup
        self.periodo = periodo
        self.anguloInclinacao = anguloInclinacao
        self.lua = lua
        self.ecc = ecc
        self.anom = anom 

        dtor = np.pi/180.
        semiEixoPixel = self.semiEixoRaioStar * self.raioEstrelaPixel

        '''Beginning of the calculation of the TOTAL transit time through the parameters passed to the planet.'''

        #ecc = 0. #default
        #anom = 0.  #default

        #calculando obliquidade

        '''
        orbit parameters
        :parameter xplaneta: x in the matrix that will project the planet
        :parameter yplaneta: y in the matrix that will project the planet
        '''

        nk=2*np.pi/(self.periodo*24)    # em horas^(-1)
        Tp=self.periodo*anom/360.*24. # tempo do pericentro (em horas)
        m = nk*(self.tempoHoras-Tp)     # em radianos

        # calculando a anomalia excentrica em radianos
        eccanom = kepler.solve(m,ecc)  # subrotina em anexo
        xs=semiEixoPixel*(np.cos(eccanom)-ecc)
        ys=semiEixoPixel*(math.sqrt(1-(ecc**2))*np.sin(eccanom))

        ang=anom*dtor-(np.pi/2)
        xp=xs*np.cos(ang)-ys*np.sin(ang)
        yp=xs*np.sin(ang)+ys*np.cos(ang)

        ie, = np.where(self.tempoHoras == min(abs(self.tempoHoras)))

        xplaneta=xp-xp[ie[0]]
        yplaneta=yp*np.cos(self.anguloInclinacao*dtor)

        pp, = np.where((abs(xplaneta) < 1.2 * tamanhoMatriz/2) & (abs(yplaneta) < tamanhoMatriz/2)) #rearranges the vector with only the points needed for the light curve analysis
        xplan = xplaneta[pp] + tamanhoMatriz/2
        yplan = yplaneta[pp] + tamanhoMatriz/2

        raioPlanetaPixel = self.raioPlanetaRstar * self.raioEstrelaPixel

       

        '''
         Beginning of the time calculation in Hours and the Light curve in the matrix
         :parameter nn: calculation of the number of points on the light curve
         :tamanhoMatriz parameter: receives the smeared star and then plots the planet
         :tempoHoras parameter: calculates the transit time in hours, transforming it into an object of the Eclipse class
         :curvaLuz parameter: calculates the light curve of the transit of the planet when eclipsing the star, it also becomes
         Eclipse object
        '''
        latitudeTransito = -np.arcsin(self.semiEixoRaioStar*np.cos(self.anguloInclinacao*dtor))/dtor # latitude Sul (arbitraria)
        # duracao do transito em horas
        duracaoTransito=2 * (90.-np.arccos((np.cos(latitudeTransito*dtor))/self.semiEixoRaioStar)/dtor)*self.periodo/360*24. 
        tempoTotal = 3 * duracaoTransito
        self.tempoTotal= tempoTotal

        
        # calculo do numero de pontos na curva de luz
        nn=np.fix(tempoTotal*60./intervaloTempo)

        #seleciona a maior orbita para que a curva de luz seja plotada de maneira correta (observando ela inteira)
        if(lua == True):
            if (len(pp)>len(self.ppMoon)):
                rangeloop = pp
            else: 
                rangeloop = self.ppMoon
                xplan = xplaneta[self.ppMoon] + tamanhoMatriz/2 #x plan e y plan se alteram caso haja o acrescimo de luas 
                yplan = yplaneta[self.ppMoon] + tamanhoMatriz/2
        else:
            rangeloop = pp


        ''''
        Light curve and intensity normalization
        '''
        # maximo da curva de luz, usado na normalizacao da intensidade
        maxCurvaLuz = np.sum(self.estrelaManchada)

        # definição de variaveis para utilizacao da função de calculo da curva de luz em C
        tamanho = self.tamanhoMatriz*self.tamanhoMatriz

        # Matriz em auxiliar para ser passada como parametro para o script em C
        em = (c_double*tamanho)()
        for j in range (self.tamanhoMatriz):
            for i in range(self.tamanhoMatriz):
                index = i*self.tamanhoMatriz + j
                num = self.estrelaManchada[i][j]
                num = (c_double)(num)
                em[index] = num

        # Matriz plan auxiliar para ser passada como parametro para o script em C
        plan = (c_double*tamanho)()
        for i in range(tamanho):
            num = 1.
            num = (c_double)(num)
            plan[i] = num

        kk=np.arange(tamanhoMatriz*tamanhoMatriz)

        # Matriz kk auxiliar para ser passada como parametro para o script em C
        kk2 = (c_double * len(kk))(*kk)

        # Verifica se o Python é 32 ou 64 bit
        if(platform.architecture()[0] == "32bit"):
            my_func = WinDLL('scripts/func32.dll', winmode = 0x8)
        elif(platform.architecture()[0] == "64bit"):
            my_func = WinDLL('scripts/func64.dll', winmode = 0x8)

        # Prepara os tipos de cada variável dos argumentos e do retorno da função do calculo da curva de luz
        my_func.curvaLuz.restype = c_double
        my_func.curvaLuz.argtypes = c_double,c_double,c_int,c_int,POINTER(c_double),POINTER(c_double),c_double
        my_func.curvaLuzLua.restype = c_double
        my_func.curvaLuzLua.argtypes = c_double,c_double,c_double,c_double,c_double,c_int,c_int,POINTER(c_double),POINTER(c_double),c_double

        raioPlanetaPixel = int(raioPlanetaPixel)

        '''
        Matrix generation for ploting
        '''
        if(anim):
            #criacao de variaveis para plotagem da animacao 
            fig, (ax1, ax2) = plt.subplots(2,1)
            ims = []
            plota = True #variavel FLAG que indica quando armazenar a imagem do PLOT 
            numAux = 0 #variavel FLAG que indica quantidade de imagens no vetor de PLOT

            print("\nWait a moment please. Eclipse animation is being generated.\n")
            #Start of loops for plotting and calculating traffic
            #start = time.time()
            intervalo = math.ceil(len(rangeloop)/400)
            if (lua == False):
                for i in range(0,len(rangeloop)):

                                x0 = xplan[i]
                                y0 = yplan[i]

                                self.curvaLuz[rangeloop[i]]=my_func.curvaLuz(x0,y0,self.tamanhoMatriz,raioPlanetaPixel,em,kk2,maxCurvaLuz)

                                if(plota and self.curvaLuz[rangeloop[i]] != 1 and numAux<200):
                                    plan = np.zeros(tamanhoMatriz*tamanhoMatriz)+1.
                                    ii = np.where(((kk/tamanhoMatriz-y0)**2+(kk-tamanhoMatriz*np.fix(kk/tamanhoMatriz)-x0)**2 <= raioPlanetaPixel**2))
                                    plan[ii]=0.
                                    plan = plan.reshape(self.tamanhoMatriz, self.tamanhoMatriz) #position included in the matrix
                                    plt.axis([0,self.Nx,0,self.Ny])
                                    im = ax1.imshow(self.estrelaManchada*plan,cmap="hot", animated = True)
                                    ims.append([im]) #store in the animation the graphic's spots (image)
                                    numAux+=1
                                plota = not(plota) #auxiliary variable that select the correct time interval for plotting
            else:
                for i in range(0,len(rangeloop)):

                                x0 = xplan[i] 
                                y0 = yplan[i]

                                ### adding moons ###
                                xm = x0-self.xxm[i]         
                                ym = y0-self.yym[i]  
                                
                                self.curvaLuz[rangeloop[i]]=my_func.curvaLuzLua(x0,y0,xm,ym,self.Rmoon,self.tamanhoMatriz,raioPlanetaPixel,em,kk2,maxCurvaLuz)

                                if(plota and self.curvaLuz[rangeloop[i]] != 1 and numAux<200):
                                    plan = np.zeros(tamanhoMatriz*tamanhoMatriz)+1.
                                    ii = np.where(((kk/tamanhoMatriz-y0)**2+(kk-tamanhoMatriz*np.fix(kk/tamanhoMatriz)-x0)**2 <= raioPlanetaPixel**2))
                                    ll = np.where((kk/tamanhoMatriz-ym)**2+(kk-tamanhoMatriz*np.fix(kk/tamanhoMatriz)-xm)**2 <= self.Rmoon**2)
                                    plan[ii]=0.
                                    plan[ll]=0.
                                    plan = plan.reshape(self.tamanhoMatriz, self.tamanhoMatriz) #position included in the matrix
                                    plt.axis([0,self.Nx,0,self.Ny])
                                    im = ax1.imshow(self.estrelaManchada*plan,cmap="hot", animated = True)
                                    ims.append([im]) #store in the animation the graphic's spots (image)
                                    numAux+=1
                                plota = not(plota) #auxiliary variable that select the correct time interval for plotting

            #end = time.time()
            #print(end-start)
            ax2.plot(self.tempoHoras,self.curvaLuz)
            ax2.axis([-self.tempoTotal/2,self.tempoTotal/2,min(self.curvaLuz)-0.001,1.001])
            ani =animation.ArtistAnimation(fig, ims, interval=50, blit=True,repeat_delay=0.1)
            
            plt.show()
            #ani.save('animacao_transito.gif',writer="PillowWriter") #save the generated gif in the root of the file, for use by the user
        else:
            #Start of loops for plotting and calculating traffic
            start = time.time()
            if (lua == False):
                for i in range(0,len(rangeloop)):

                                x0 = xplan[i]
                                y0 = yplan[i]

                                self.curvaLuz[rangeloop[i]]=my_func.curvaLuz(x0,y0,self.tamanhoMatriz,raioPlanetaPixel,em,kk2,maxCurvaLuz)
            else:
                for i in range(0,len(rangeloop)):

                                x0 = xplan[i] 
                                y0 = yplan[i]

                                ### adicionando luas ###
                                xm = x0-self.xxm[i]         
                                ym = y0-self.yym[i]  
                                
                                self.curvaLuz[rangeloop[i]]=my_func.curvaLuzLua(x0,y0,xm,ym,self.Rmoon,self.tamanhoMatriz,raioPlanetaPixel,em,kk2,maxCurvaLuz)
            if(plot):
                end = time.time()
                print(end-start)
                plt.plot(self.tempoHoras,self.curvaLuz)
                plt.axis([-self.tempoTotal/2,self.tempoTotal/2,min(self.curvaLuz)-0.001,1.001])
                plt.show()

        locals().clear # Limpa qualquer possível sujeira de memória
        del my_func

        error=0
        self.error=error

    '''Calling the objects assigned to the Eclipse class.'''
    def getTempoTransito(self):
        '''Returns the tempoTotal parameter, representing the planet's transit time on its host star.'''
        return self.tempoTotal
    def getTempoHoras(self):
        '''Returns the tempoHoras parameter, representing the planet's transit time on its host star in Hours.'''
        return self.tempoHoras
    def getCurvaLuz(self):
        '''Returns the curvaLuz parameter, representing the light curve of the star that has a planet orbiting it.'''
        return self.curvaLuz

    def getRaioPlan(self):
        return self.raioPlanetaRstar

    def getRplanJup(self):
        return self.raioPlanJup

    def getSemiEixo(self):
        return self.semiEixoUA
    
    def getsemiEixoRaioStar(self):
        return self.semiEixoRaioStar

    def getPeriodo(self):
        return self.periodo

    def getInc(self):
        return self.anguloInclinacao

    def getEccAnom(self):
        return self.ecc, self.anom 

    def getLua(self):
        return self.lua

    def getError(self):
        '''
        Returns the error value, whether or not one occurs. If there is no error, it receives 0. If there is, the variable will have
        its start value (which is -1)
        '''
        return self.error
    def setEstrela(self,estrela):
        '''
        with this function, it is possible to pass the updated star to the eclipse that is forming, in case more spots are added.
        '''
        self.estrelaManchada = estrela