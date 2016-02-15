"""
Werner `95 Dune code

Slabs are moved a specified number of lattice sites,
in the transport direction and is deposited at this site with a probability
that depends upon the presence/absence of sand slabs.

If the slab is not deposited, then it repeatedly is moved to other sites
(specified by thesaltation distance) in the transport direction until
deposition, following which another slab is
chosen randomly for transport. This procedure is repeated to construct the
time evolution of the surface. Slab movement not parallel to the transport
direction originates only from enforcement of the angle of repose.
"""
import time
start_time = time.time()

import matplotlib.pyplot as plt
import numpy as np

#max time for model run (this is real timesteps (i.e., each one is L**2 moves))
Tmax=1;
timesteps=np.arange(1,(Tmax+1));

# define the domain size
L=1000;
Li=L-1

#how many slabs in the domain
#SlabsinD=1000000000;
#How many slabs per cell
each=15

#make the domain and spread them equally (with a random component)
D=np.random.random_integers(0,1,[L,L])+each
np.savetxt("grid.txt", D)
#Saltation Length
Salt=1

#probability of deposition:
#pns= prob. of deposition at a site with no sand slabs
#ps= prob. of deposition at a site with at least one sand slab
pns=0.4;
ps=0.6;

movedslabs=0

for t in timesteps:

    #pick the random cell orders for each single  timestep(which has L**2 moves)
    Picks=np.random.random_integers(0,99,[2,L**2])
    for R,C in np.nditer([Picks[0,:],Picks[1,:],]):

    #if there is a slab at the site, pick it up
        if D[R,C] > 0:
            D[R,C]=D[R,C]-1
            slab=1;
            movedslabs=movedslabs+1
            #move the slab
            while slab == 1:

                #set the initial receiver cell
                hop=R+Salt

                #periodic BCs
                if hop > (L-1):
                    hop=hop-L

                # probablity for depositon
                prob=np.random.random()


                #define the shadow zone
                if hop==0:
                    shadow1=Li
                    shadow2=Li-1
                    shadow3=Li-2
                    shadow4=Li-3
                elif hop==1:
                    shadow1=0
                    shadow2=Li
                    shadow3=Li-1
                    shadow4=Li-2
                elif hop==2:
                    shadow1=1
                    shadow2=0
                    shadow3=Li
                    shadow4=Li-1
                elif hop==3:
                    shadow1=2
                    shadow2=1
                    shadow3=0
                    shadow4=Li
                else:
                    shadow1=hop-1
                    shadow2=hop-2
                    shadow3=hop-3
                    shadow4=hop-4

                #DEPOSITION RULES:
                #shadow Zone
                if D[shadow1,C]-D[hop,C]>=1 or D[shadow2,C]-D[hop,C]>=2 or D[shadow3,C]-D[hop,C]>=3 or D[shadow4,C]-D[hop,C]>=1:
                    D[hop,C]=D[hop,C]+1;
                    slab=0;
                #sandy cell
                elif D[hop,C] > 0 and prob < ps:
                    D[hop,C]=D[hop,C]+1;
                    slab=0;
                #empty cell
                elif D[hop,C]==0 and prob < pns:
                    D[hop,C]=D[hop,C]+1;
                    slab=0;
                #hop again
                else:
                    R=hop

            #Angle of repose: only gradient of '2' is allowed for the lattice
            RE=R
            CE=C
            erosioncheck=1
            while erosioncheck==1:
                if RE==Li:
                    REp=0
                else:
                    REp=RE+1
                if RE==0:
                    REm=Li
                else:
                    REm=RE-1
                if CE ==Li:
                    CEp=0
                else:
                    CEp=CE+1
                if CE ==0:
                    CEm=Li
                else:
                    CEm=CE-1

                N=[D[REp,CE],D[REm,CE],D[RE,CEp],D[RE,CEm]]
                Tallestneighborindex=np.argmax(N)
                Tallestneighbor=N[Tallestneighborindex]
                if Tallestneighbor-D[RE,CE] > 2:
                    #add one to the erosion site
                    D[RE,CE]=D[RE,CE]+1
                    #subtract from tallest neighbor
                    #change the site of interest to the old tallest dune site. (steepest gradient)
                    if Tallestneighborindex == 0:
                        D[REp,CE]=D[REp,CE]-1
                        RE=REp
                    elif Tallestneighborindex == 1:
                        D[REm,C]=D[REm,CE]-1
                        RE=REm
                    elif Tallestneighborindex == 2:
                        D[RE,CEp]=D[RE,CEp]-1
                        CE=CEp
                    elif Tallestneighborindex == 3:
                        D[RE,CEm]=D[RE,CEm]-1
                        CE=CEm
                else:
                    erosioncheck=0

            depcheck=1
            while depcheck==1:
                if hop==Li:
                    hopp=0
                else:
                    hopp=hop+1
                if hop==0:
                    hopm=Li
                else:
                    hopm=hop-1
                if C ==Li:
                    Cp=0
                else:
                    Cp=C+1
                if C ==0:
                    Cm=Li
                else:
                    Cm=C-1

                N=[D[hopp,C],D[hopm,C],D[hop,Cp],D[hop,Cm]]
                Smallestneighborindex=np.argmin(N)
                Smallestneighbor=N[Smallestneighborindex]
                if D[hop,C]-Smallestneighbor > 2:
                    #subtract one from the dep site
                    D[hop,C]=D[hop,C]-1
                    #add 1 to shortest neighbor
                    #and change the site of interest to the old shortest site. (steepest gradient)
                    if Smallestneighborindex == 0:
                        D[hopp,C]=D[hopp,C]+1
                        hop=hopp
                    elif Smallestneighborindex == 1:
                        D[hopm,C]=D[hopm,C]+1
                        hop=hopm
                    elif Smallestneighborindex == 2:
                        D[hop,Cp]=D[hop,Cp]+1
                        C=Cp
                    elif Smallestneighborindex == 3:
                        D[hop,Cm]=D[hop,Cm]+1
                        C=Cm
                else:
                    depcheck=0


#plt.matshow(D)
Di = np.genfromtxt("grid.txt")
#plt.matshow(Di-D)
Dchange=Di-D

#print(np.sum(Di)-np.sum(D))

#print(movedslabs/(L*L))
print(Tmax/(L*L))

print("--- %s seconds ---" % (time.time() - start_time))
