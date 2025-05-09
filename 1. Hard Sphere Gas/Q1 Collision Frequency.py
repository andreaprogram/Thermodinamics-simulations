
from vpython import *
#Web VPython 3.2
import numpy as np
# Hard-sphere gas.

# Bruce Sherwood

win = 500 # tamny de la finestra on s'executa el programa

Natoms = 500  # fixem N=500

# CONDICIONS TERMODINAMIQUES DEL SISTEMA-----------------------------
L = 1 # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container
mass = 4E-3/6E23 # helium mass (massa molar)
Ratom = 0.03 # wildly exaggerated size of helium atom
k = 1.4E-23 # Boltzmann constant [J/K]
T = 300 # around room temperature [K]
dt = 1E-5
cross_sec= 0.21E-9
n = Natoms / (L**3)  # densitat numèrica
v_avg = np.sqrt(8 * k * T / (np.pi * mass))  # velocitat mitjana

Z = (np.sqrt(2)*cross_sec*(n**2)*v_avg)/2 # freqüència de col·lisions
temps_teo = 1 / Z
print("Temps teòric mig entre col·lisions", temps_teo)

# AJUSTOS VISUALITZACIÓ----------------------------------------------------------
animation = canvas( width=win, height=win, align='left')
animation.range = L
animation.title = 'A "hard-sphere" gas'
s = """  Theoretical and averaged speed distributions (meters/sec).
  Initially all atoms have the same speed, but collisions
  change the speeds of the colliding atoms. One of the atoms is
  marked and leaves a trail so you can follow its path.
  
"""
animation.caption = s

# Dibuix del cub on se situa el sistema 
#12 aristas del cub:
#4 del fons (boxbottom)
#4 del sostre (boxtop)
# 4 verticals (las vert)

d = L/2+Ratom
r = 0.005
boxbottom = curve(color=gray, radius=r)
boxbottom.append([vector(-d,-d,-d), vector(-d,-d,d), vector(d,-d,d), vector(d,-d,-d), vector(-d,-d,-d)])
boxtop = curve(color=gray, radius=r)
boxtop.append([vector(-d,d,-d), vector(-d,d,d), vector(d,d,d), vector(d,d,-d), vector(-d,d,-d)])
vert1 = curve(color=gray, radius=r)
vert2 = curve(color=gray, radius=r)
vert3 = curve(color=gray, radius=r)
vert4 = curve(color=gray, radius=r)
vert1.append([vector(-d,-d,-d), vector(-d,d,-d)])
vert2.append([vector(-d,-d,d), vector(-d,d,d)])
vert3.append([vector(d,-d,d), vector(d,d,d)])
vert4.append([vector(d,-d,-d), vector(d,d,-d)])

# INICIALITZACIÓ DEL PROCÉS-------------------------------------------
Atoms = [] #llista d'atoms
p = [] #llista de moments lineals
apos = [] #llista de posicions
pavg = sqrt(2*mass*1.5*k*T) # average kinetic energy p**2/(2mass) = (3/2)kT

#per cada atom, li assignem posicions inicials randoms (que no surtin de la caixa) i moments inicials randoms (que segueixin Ta equipartició)
for i in range(Natoms):
    x = L*random()-L/2
    y = L*random()-L/2
    z = L*random()-L/2
    if i == 0: #fas que a l'àtom nº 1 li segueixi una traça per on passa
        Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=color.magenta, make_trail=True, retain=100, trail_radius=0.3*Ratom))
    else: Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
    apos.append(vec(x,y,z))
    theta = pi*random()
    phi = 2*pi*random()
    px = pavg*sin(theta)*cos(phi)
    py = pavg*sin(theta)*sin(phi)
    pz = pavg*cos(theta)
    p.append(vector(px,py,pz))



# DEFINICIÓ DE L'HISTOGRAMA DE VELOCITATS (CAL QUE SEGUEIXI UNA GAUSSIANA, TAL COM ASSIGNEM EL MOMENT)-----------------------------------------------------
deltav = 100 # binning for v histogram

def barx(v):
    return int(v/deltav) # index into bars array

nhisto = int(4500/deltav)
histo = []
for i in range(nhisto): histo.append(0.0)
histo[barx(pavg/mass)] = Natoms

gg = graph( width=win, height=0.4*win, xmax=3000, align='left',
    xtitle='speed, m/s', ytitle='Number of atoms', ymax=Natoms*deltav/1000)

theory = gcurve( color=color.blue, width=2 )
dv = 10
for v in range(0,3001+dv,dv):  # theoretical prediction
    theory.plot( v, (deltav/dv)*Natoms*4*pi*((mass/(2*pi*k*T))**1.5) *exp(-0.5*mass*(v**2)/(k*T))*(v**2)*dv ) #GAUSSIANA TEÒRICA (MB)

accum = []
for i in range(int(3000/deltav)): accum.append([deltav*(i+.5),0])
vdist = gvbars(color=color.red, delta=deltav )

def interchange(v1, v2):  # remove from v1 bar, add to v2 bar
    barx1 = barx(v1)
    barx2 = barx(v2)
    if barx1 == barx2:  return
    if barx1 >= len(histo) or barx2 >= len(histo): return
    histo[barx1] -= 1
    histo[barx2] += 1


#RECOMPTE DE COL·LISIONS----------------------------------------------------------------------------------------------------------------------------
def checkCollisions():
    hitlist = []
    r2 = 2*Ratom  #distancia mínima entre 2 àtoms
    r2 *= r2      #quadrat d'aquesta distància
    for i in range(Natoms):  #compara cada parell d'àtoms, fixant l'atom 'i' i comparant-lo amb els 'j' altres
        ai = apos[i]
        for j in range(i) :
            aj = apos[j]
            dr = ai - aj #distancia entre l'atom 'i' i 'j'
            if mag2(dr) < r2: hitlist.append([i,j]) #si la distancia es menor a la distancia minima, es una col·lisio
    return hitlist


# EXECUCIÓ DE L'ALGORISME: una vegada establertes les CI, si dos particules xoquen s'els assigna un nou moment --------------------------------------
nhisto = 0 # number of histogram snapshots to average

colisions = 0
temps_total = 0

while True:
    rate(300)
  
    # Accumulate and average histogram snapshots
    temps_total += dt #contar el temps
    for i in range(len(accum)): accum[i][1] = (nhisto*accum[i][1] + histo[i])/(nhisto+1)
    if nhisto % 10 == 0:
        vdist.data = accum
    nhisto += 1

    # Update all positions
    for i in range(Natoms): Atoms[i].pos = apos[i] = apos[i] + (p[i]/mass)*dt
    
    # Check for collisions
    hitlist = checkCollisions()
    colisions += len(hitlist)
    # If any collisions took place, update momenta of the two atoms
    for ij in hitlist:
        i = ij[0]
        j = ij[1]
        ptot = p[i]+p[j]
        posi = apos[i]
        posj = apos[j]
        vi = p[i]/mass
        vj = p[j]/mass
        vrel = vj-vi
        a = vrel.mag2
        if a == 0: continue;  # exactly same velocities
        rrel = posi-posj
        if rrel.mag > Ratom: continue # one atom went all the way through another
    
        # theta is the angle between vrel and rrel:
        dx = dot(rrel, vrel.hat)       # rrel.mag*cos(theta)
        dy = cross(rrel, vrel.hat).mag # rrel.mag*sin(theta)
        # alpha is the angle of the triangle composed of rrel, path of atom j, and a line
        #   from the center of atom i to the center of atom j where atome j hits atom i:
        alpha = asin(dy/(2*Ratom)) 
        d = (2*Ratom)*cos(alpha)-dx # distance traveled into the atom from first contact
        deltat = d/vrel.mag         # time spent moving from first contact to position inside atom
        
        posi = posi-vi*deltat # back up to contact configuration
        posj = posj-vj*deltat
        mtot = 2*mass
        pcmi = p[i]-ptot*mass/mtot # transform momenta to cm frame
        pcmj = p[j]-ptot*mass/mtot
        rrel = norm(rrel)
        pcmi = pcmi-2*pcmi.dot(rrel)*rrel # bounce in cm frame
        pcmj = pcmj-2*pcmj.dot(rrel)*rrel
        p[i] = pcmi+ptot*mass/mtot # transform momenta back to lab frame
        p[j] = pcmj+ptot*mass/mtot
        apos[i] = posi+(p[i]/mass)*deltat # move forward deltat in time
        apos[j] = posj+(p[j]/mass)*deltat
        interchange(vi.mag, p[i].mag/mass)
        interchange(vj.mag, p[j].mag/mass)


    if colisions % 5000 == 0:    #cada 5000 colisions mirem el temps mig
            print("Nombre de colisions:", colisions)
            print("Temps total:", temps_total)
            print("Temps mig entre colisions:", temps_total / colisions)
            print("Freqüència colisions:", colisions / temps_total)
            temps_total=0
            colisions=0
      
    
    for i in range(Natoms):
        loc = apos[i]
        if abs(loc.x) > L/2: #considerem quan una particula 'sobrepassa' la paret
            if loc.x < 0: p[i].x =  abs(p[i].x) #si loc<0 estarà a la paret esquerra, aleshores cal que p>0
            else: p[i].x =  -abs(p[i].x)  # si loc>0 estara a la paret dreta, aleshores cla que p<0
        
        if abs(loc.y) > L/2:
            if loc.y < 0: p[i].y = abs(p[i].y)
            else: p[i].y =  -abs(p[i].y)
        
        if abs(loc.z) > L/2:
            if loc.z < 0: p[i].z =  abs(p[i].z)
            else: p[i].z =  -abs(p[i].z)
    
