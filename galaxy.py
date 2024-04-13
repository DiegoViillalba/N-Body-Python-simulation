import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import math
from tqdm import tqdm


# Universal constants

G = 1.0 # Universal gravity constant
M = 150.0 # Mass of the center massive object

# Time step
dt = 0.05 # Time step
t = np.arange(0.0, 10, dt) # create a time array from 0..n seconds sampled at 0.05 second steps

# Number of particles
N = 10

# Mass of the black hole in the center
M_bh = 5000

# Velocity of the galaxy
v1_g = [0,0]

# Storage array
A =[]

# number of massive objects
nmo=2

for i in range(0,N+nmo):
    A.append([[],[]])


# Math functions
def dotproduct(u, v):
  return sum((a*b) for a, b in zip(u, v))

def length(u):
  return math.sqrt(dotproduct(u, u))

# Class for vectors and particles
class vector:
    def __init__(self,x,y):
        self.x=x
        self.y=y
    
    # Method to print values
    def __str__(self):
        return f"{self.x} , {self.y}"
    def inner_prod(self, v):
        return self.x*v.x + self.y*v.y
    def length(self):
        return math.sqrt(self.inner_prod(self))
    def distance(self, v):
        return math.sqrt((v.x-self.x)**2+(v.y-self.y)**2)
    

class particle:
    def __init__(self,n,r,v,m):
        # Atributes of the particle
        self.n = n
        self.m = m
        self.r = r
        self.v = v

    def __str__(self):
        return f"Particle {self.n} \nposition -> ({self.r}) \nvelocity -> ({self.v})"
        # Calculation of the force
    
    def R(self):
        R = np.array([self.r.x , self.r.y])
        return R
    
    def position(self):
        x = round(self.r.x +(self.v.x*dt),3)
        y = round(self.r.y +(self.v.y*dt),3)
        return x,y

    def velocity(self,P):
        v_x = round(self.v.x  - (self.acc(P)[0]*dt),3)
        v_y = round(self.v.y  - (self.acc(P)[1]*dt),3)
        return v_x,v_y

        
    def acc(self,P):


        # Acceleration function
        A = np.array([0,0])
        epsilon = 0.5

        for i in range(len(P)):
            if P[i].n == self.n:
                A = A
            else:
                r_rel = P[i].R()-self.R()+epsilon
                m2 = P[i].m
                r2 = length(r_rel)**2
                A = A +((-1*G)*(m2/(r2)))*r_rel

        a_x = round(A[0],3)
        a_y = round(A[1],3)
        return a_x,a_y




######### Main program ########
    
print("- Nbody simulation -\n")
print(f"N = {N}\n")
print("T = 60s \n")
    
P = [] # Particle list

# We define the inicial velocities of the particles to generate stable orbits 
v_orb = []

for n in range(0,N):
    r_dist = round(random.uniform(10,100),3)
    a = r_dist + round(random.uniform(0,30),3)

    vo = round(math.sqrt(2*G*M_bh*((2/(r_dist))-(1/a))),4)


    # if n<int(N/4):
    #     r_dist = -1*r_dist 
    #     vo = -1*vo
    #     inf = [vo,r_dist]
    #     v_orb.append(inf)
    #     #print(n)


    # elif n>int(N/4) and n<int(N/2):
    #     inf = [vo,r_dist]
    #     v_orb.append(inf)



    # elif n>int(N/2) and n<int(3*(N/4)):
    #     vo = -1*vo
    #     inf = [vo,r_dist]
    #     v_orb.append(inf)
    #     #print(n)


    # else:

    #     r_dist = -1*r_dist 
    #     inf = [vo,r_dist]
    #     v_orb.append(inf)


    # Dual ax galaxy

    if n<int(N/2):
        r_dist = -1*r_dist 
        vo = -1*vo
        inf = [vo,r_dist]
        v_orb.append(inf)
        #print(n)

    else:
        inf = [vo,r_dist]
        v_orb.append(inf)



    

#print(v_orb)

for n in range(0,N):
    #P.append(particle(n,vector(round(random.uniform(-5,5),3),round(random.uniform(-5,5))),vector(round(random.uniform(-3,3),3),round(random.uniform(-3,3),3))))
    #P.append(particle(n,vector(round(random.uniform(-50,50),3),round(random.uniform(-50,50))),vector(round(v_orb[n][0]),0),round(random.uniform(1,10),3)))
    P.append(particle(n,vector(v_orb[n][1],0),vector(v1_g[0],v_orb[n][0]+v1_g[1]),round(random.uniform(1,10),3)))

    # if n<int(N/2):
    #     P.append(particle(n,vector(v_orb[n][1],0),vector(v1_g[0],v_orb[n][0]+v1_g[1]),round(random.uniform(1,10),3)))

    # else:
    #     P.append(particle(n,vector(0,v_orb[n][1]),vector(v_orb[n][0]+v1_g[0],v1_g[1]),round(random.uniform(1,10),3)))

    


P.append(particle(N,vector(0,0),vector(v1_g[0],v1_g[1]),M_bh))

# Second black hole
P.append(particle(N,vector(100,-20),vector(0,0),M_bh))
#----- Frame creation -------

#print(P[100])
#
for n in tqdm(t):

    for i in range(0,N+nmo):

        A[i][0].append(P[i].r.x)
        A[i][1].append(P[i].r.y)

        v_x,v_y = P[i].velocity(P)
        P[i].v = vector(v_x,v_y)
        x,y = P[i].position()
        P[i].r = vector(x,y)


#print(P)
#----- Frame calibration -----

# Animation
fig = plt.figure()
fig.patch.set_facecolor('xkcd:black')
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-300, 500), ylim=(-300, 300))
ax.set_aspect('equal')

ax.set_facecolor((0, 0, 0))

ax.grid()

#Animated particles

particles = []
for n in range(N):
    dot, = ax.plot([], [], '*', lw=2, markersize=2)
    particles.append(dot,)

# Black hole
dot,= ax.plot([], [], 'o', lw=2)
particles.append(dot,)

# Black hole
dot,= ax.plot([], [], 'o', lw=2)
particles.append(dot,)
# Time display
    
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes, color='white')




# Animation function


def animate(i):

    # Time text
    time_text.set_text(time_template % (i*dt))

    #Animated particles
    result = [time_text]

    for n in range(N+nmo):
        particles[n].set_data(A[n][0][i],A[n][1][i])
        result.append(particles[n])

    
    return [result[i] for i in range(N+nmo+1)]



ani = animation.FuncAnimation(fig, animate, np.arange(1, len(A[0][1])),
                              interval=25, blit=True)
# ani.save('n_body.mp4', fps=15)
plt.show()

    