import numpy as np
from matplotlib import pyplot
from estrela import Estrela
from eclipse import Eclipse
#from model import Modelo
from verify import Validar,calSemiEixo,calculaLat

'''
main programmed for professionals and students familiar with the area
--- star ---
raio parameter: The radius of the star in pixels
raioSun parameter: The radius of the star in units of Rsun
intensidadeMaxima parameter: Intensity of the center of the star
coeficienteHum parameter: Limbo darkening coefficient
coeficienteDois parameter: Limbo darkening coefficient
tamanhoMatriz parameter: Size of the matrix in which the star will be built
estrela parameter: Star constructed with the limbo darkening coefficients
star_ object :: is the star object where the star array is stored according to the parameters. Class function calls
star are made through it.
--- planet ---
:periodo parameter: rotation period of the planet
:semiEixoRaioStar parameter: semi axis of the planet in relation to the radius of the star
:semiEixoUA parameter: semiaxis of the planet in UA
:anguloInclinacao parameter: tilt angle of the planet
:raioPlanetaRstar parameter: radius of the planet in relation to the radius of the star
:raioPlanJup parameter: radius of the planet in relation to the radius of Jupiter
---  spot --- 
latsugerida parameter:: suggested latitude for the spot
fa parameter:: vector with the area of each spot
fi parameter :: vector with the intensity of each spot
li parameter:: vector with the length of each spot
quantidade parameter:: variable that stores the amount of stains
r parameter:: radius of the spot in relation to the radius of the star
intensidadeMancha parameter:: stain intensity in relation to the star's intensity
lat parameter:: latitude of the spot
longt parameter:: longitude of the spot
raioMancha parameter:: actual smear radius
area parameter:: area of the spot
--- eclipse ---
eclipse parameter:: variable that holds the object of the eclipse class that generates the light curve. Class function calls
    Eclipse() are done through it.
tempoTransito parameter:: transit time of the planet
curvaLuz parameter:: light curve matrix that will be plotted as a graph
tempoHoras parameter:: transit time in hours
'''

raio= 373. #default (pixel)
intensidadeMaxima=240 #default
tamanhoMatriz = 856 #default
raioSun=0.117 #star's radius relative to the sun's radius
raioStar=raioSun*696340 #multiplying by the solar radius in km
coeficienteHum=0.65
coeficienteDois=0.28


#create star
estrela_ = Estrela(raio,raioSun,intensidadeMaxima,coeficienteHum,coeficienteDois,tamanhoMatriz)

Nx= estrela_.getNx() #Nx and Ny needed for eclipse plotting
Ny= estrela_.getNy()
dtor = np.pi/180.  

periodo = 6.099 # in days
anguloInclinacao = 89.86  # in hours
ecc = 0
anom = 0 
raioPlanJup = 0.0819 #relative to the radius of jupiter
raioPlanetaRstar = (raioPlanJup*69911)/raioStar #multiplying by the radius of jupiter in km


dec=int(input("Do you want to calculate the semi-orbital axis of the planet using KEPLER'S 3rd LAW? | 1. Yes |  2.No |"))
if dec==1:
    mass=0. #put the mass of the star in relation to the mass of the sun
    semieixoorbital = calSemiEixo(periodo,mass)
    semiEixoRaioStar = ((semieixoorbital/1000)/raioStar)
    #turns into km to do with respect to the radius of the star
else:
    semiEixoUA = Validar('Semi eixo (em UA:)')
    # in units of RStar
    semiEixoRaioStar = ((1.469*(10**8))*semiEixoUA)/raioStar
    #multiplying by the AU (turning into Km) and converting with respect to the radius of the star

latsugerida = calculaLat(semiEixoRaioStar,anguloInclinacao)
print("The suggested latitude for the spot to influence the star's light curve is:", latsugerida)

#spots
count = 0
quantidade = 0 #desired amount of spots, if you want to add, change this variable
#creates quantity-sized arrays to place the smear parameters
fa = [0.]*quantidade #array area spots
fi = [0.]*quantidade #array intensity spots
li = [0.]*quantidade #array longitude spots

while count!=quantidade: # the loop will rotate the number of spots selected by the user
    print('\033[1;35m\n\n══════════════════ Spots Parameters',count+1,'═══════════════════\n\n\033[m')
    r = Validar('Enter the spot radius as a function of the star radius in pixels:')
                
    intensidadeMancha= float(input('Enter the intensity of the spot as a function of the maximum intensity of the star:'))
    fi[count]=intensidadeMancha
    lat=float(input('Latitude of the spot:'))
    longt=float(input('Longitude of the spot::'))
    li[count]=longt

    raioMancha= r*raioStar
    area = np.pi *(raioMancha**2)
    fa[count]= area

    estrela=estrela_.manchas(r,intensidadeMancha,lat,longt) #receives the choice of whether to receive stains or not
    count+=1

#print array of intensity, longitude and area of the spot for testingprint("Intensity:",fi)
print("Longitudes:",li)
print("Areas:",fa)

estrela = estrela_.getEstrela()
#to plot the star
#if you don't want to plot the star, comment lines below
if (quantidade>0): #se manchas foram adicionadas. plotar
    estrela_.Plotar(tamanhoMatriz,estrela)

#creating moon
lua = True #if you don't want moons, change it to False
eclipse= Eclipse(Nx,Ny,raio,estrela)
estrela_.Plotar(tamanhoMatriz,estrela)
eclipse.geraTempoHoras()
tempoHoras=eclipse.getTempoHoras()
#instantiating MOON
rmoon = 0.5 #relative to the radius of the earth
rmoon = rmoon *6371 #multiplying by the R of the land in km
mass = 0.001 #in relation to the mass of the earth
mass = mass * (5.972*(10**24))
massPlaneta = 0.002 # relative to Jupiter's Radius
massPlaneta = massPlaneta * (1.898 *(10**27)) # convert to grams because of the constant G
G = (6.674184*(10**(-11)))
perLua = 0.1 #in days
distancia=((((perLua*24.*3600./2./np.pi)**2)*G*(massPlaneta+mass))**(1./3))/raioStar
distancia = distancia/100
moon = eclipse.criarLua(rmoon,mass,raio,raioStar,tempoHoras,anguloInclinacao,periodo,distancia)
estrela = estrela_.getEstrela()



#eclipse
eclipse.criarEclipse(semiEixoRaioStar, semiEixoUA, raioPlanetaRstar,raioPlanJup,periodo,anguloInclinacao,lua,ecc,anom)


print ("Total Time (Transit):",eclipse.getTempoTransito()) 
tempoTransito=eclipse.getTempoTransito()
curvaLuz=eclipse.getCurvaLuz()
tempoHoras=eclipse.getTempoHoras()

#Light curve plot
pyplot.plot(tempoHoras,curvaLuz)
pyplot.axis([-tempoTransito/2,tempoTransito/2,min(curvaLuz)-0.001,1.001])                       
pyplot.show()
