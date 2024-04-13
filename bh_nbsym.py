import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import random
import math
from tqdm import tqdm

# Universal constants

G = 1.0 # Universal gravity constant
M = 150.0 # Mass of the center massive object

# Time step
dt = 0.05 # Time step
t = np.arange(0.0, 30, dt) # create a time array from 0..n seconds sampled at 0.05 second steps

# Number of particles
N = 10

# Mass of the black hole in the center
M_bh = 50000

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


# ------- CLASSES ----------

# --- Particle classes

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

    def velocity(self,P,a_x,a_y):
        v_x = round(self.v.x  - (a_x*dt),3)
        v_y = round(self.v.y  - (a_y*dt),3)
        return v_x,v_y


# --- Quad tree class
    
class node:
    def __init__(self,r,s):

        self.r = r # Position vector of the node (numpy array)
        self.s = s # Width of the node
        self.child = False # Initialize the node without childrem
        self.p = []
        self.M = 0
        self.cm = np.array([0,0])

    def divide_quad(self):
        # Divide the quadrant in subquadrants (nw,ne,sw,se)
        r = self.r
        s = self.s

        r_nw = r + np.array([-1*s/4,s/4])
        self.nw = node(r_nw,s/2)#North west wuadrant
        r_ne = r + np.array([s/4,s/4])
        self.ne = node(r_ne,s/2)#North west wuadrant
        r_sw = r + np.array([-1*s/4,-1*s/4])
        self.sw = node(r_sw,s/2)#North west wuadrant
        r_se = r + np.array([s/4,-1*s/4])
        self.se = node(r_se,s/2)#North west wuadrant

        self.child=True

    def contains(self,p):
        if p.r.x <= self.r[0]+ (self.s/2) and p.r.x >= self.r[0] - (self.s/2) and p.r.y <= self.r[1]+ (self.s/2) and p.r.y >= self.r[1] - (self.s/2):
            return True
        else:
            return False

    
    def inser_point(self,p):
        #print("\n\n -----------")
        #print(self.p)
        
        #print(self.M)
        #print("Center of mass")
        #print(self.cm)
        
        if self.contains(p) == False:
            #print("Node skipped")
            return False
        
        else:
            if len(self.p)==0:
                #print("End of node")
                self.p.append(p)
                self.center_of_mass()
                #print(p.n)
                #print(self.p)

                return
            else:
                self.p.append(p)
                self.center_of_mass()
                #print("Inserion in quadrants of node")
                if self.child == False:

                    self.divide_quad()

                    for n in range(0,len(self.p)):
                        self.nw.inser_point(self.p[n])
                        self.ne.inser_point(self.p[n])
                        self.sw.inser_point(self.p[n])
                        self.se.inser_point(self.p[n])              
                    pass

                else:
                    self.nw.inser_point(p)
                    self.ne.inser_point(p)
                    self.sw.inser_point(p)
                    self.se.inser_point(p)  

    def reset(self):
        if self.child==True:
            self.p=[]
            self.nw.reset()
            self.nw = None
            self.ne.reset()
            self.ne = None
            self.sw.reset()
            self.sw = None
            self.se.reset()
            self.se = None
            self.child = False
        else:
            return None
        
    def center_of_mass(self):
        #print("Number of particles in the node")
        #print(self.p)

        if len(self.p)==0:
            return False
        
        M=0
        x=0
        y=0

        for n in range(0,len(self.p)):
            M+=self.p[n].m
        self.M = M

        if self.M ==0 :
            return
        else:
            for n in range(0,len(self.p)):
                x += self.p[n].r.x*self.p[n].m
                y += self.p[n].r.y*self.p[n].m
            x = x/M
            y = y/M
            #print("Center of mass")
            #print([x,y])
            self.cm = np.array([x,y])

    def acceleration(self,p):
        #print(p)
        theta = 5 # Parameter to define the ratio to calculate a full brute force 

        A = np.array([0,0])
        epsilon=10
        if len(self.p)==0:
            #print("No particles in the node")
            a_x=0
            a_y=0
            
            return a_x,a_y
    
        else:
            #print("Calc1")
            #print(self.cm[0])
            d_cm = math.sqrt((self.cm[0]-p.R()[0])**2+(self.cm[1]-p.R()[1])**2)  # Distance of the particle to the center of mass of the node
            #print(d_cm)

            if d_cm==0: # Skipping the case where the body is itself
                a_x=0
                a_y=0
                return a_x,a_y
            
            q = self.s / d_cm
            #print(q)

            if q < theta: # Case of interest, the node is a single body
                #print("node treated as a single body")
                #print(p.R())
                #print(self.cm)
                r_rel = self.cm - p.R() +epsilon
                #print(r_rel)
                m2 = self.M
                r2 = length(r_rel)**2
                A = A +((-1*G)*(m2/(r2)))*r_rel
                a_x = round(A[0],3)
                a_y = round(A[1],3)

                return a_x,a_y
            else:
                #print("Calculus of the child nodes")
                a_x = 0
                a_y = 0

                if self.child==False:
                    r_rel = self.cm - p.R() +epsilon
                    #print(r_rel)
                    m2 = self.M
                    r2 = length(r_rel)**2
                    A = A +((-1*G)*(m2/(r2)))*r_rel
                    a_x = round(A[0],3)
                    a_y = round(A[1],3)

                    return a_x,a_y 


                a_x = self.nw.acceleration(p)[0]+self.ne.acceleration(p)[0]+self.sw.acceleration(p)[0]+self.se.acceleration(p)[0]
                a_y = self.nw.acceleration(p)[1]+self.ne.acceleration(p)[1]+self.sw.acceleration(p)[1]+self.se.acceleration(p)[1]

                return a_x,a_y


#    elif self.child==False:
#             r_rel = self.cm - p.R() 
#             #print(r_rel)
#             m2 = self.M
#             r2 = length(r_rel)**2
#             A = A +((-1*G)*(m2/(r2)))*r_rel
#             a_x = round(A[0],3)
#             a_y = round(A[1],3)

#             return a_x,a_y             


# Class for debugging and drawing the child boxes of the tree

def draw_child(quad):
    if quad.child== False:
        return None
    else:
        plt.gca().add_patch(patches.Rectangle((quad.nw.r[0]-(quad.ne.s)/2,quad.nw.r[1]-(quad.ne.s)/2),quad.nw.s,quad.nw.s,linewidth=1, edgecolor='r', facecolor='none'))
        draw_child(quad.nw)
        plt.gca().add_patch(patches.Rectangle((quad.ne.r[0]-(quad.ne.s)/2,quad.ne.r[1]-(quad.ne.s)/2),quad.ne.s,quad.ne.s,linewidth=1, edgecolor='r', facecolor='none'))
        draw_child(quad.ne)
        plt.gca().add_patch(patches.Rectangle((quad.sw.r[0]-(quad.nw.s)/2,quad.sw.r[1]-(quad.sw.s)/2),quad.sw.s,quad.sw.s,linewidth=1, edgecolor='r', facecolor='none'))
        draw_child(quad.sw)
        plt.gca().add_patch(patches.Rectangle((quad.se.r[0]-(quad.se.s)/2,quad.se.r[1]-(quad.se.s)/2),quad.se.s,quad.se.s,linewidth=1, edgecolor='r', facecolor='none'))
        draw_child(quad.se)



        


######### Main program ########
    
print("- Nbody simulation -\n")
print(f"N = {N}\n")
print("T = 60s \n")
    
P = [] # Particle list

# We define the inicial velocities of the particles to generate stable orbits 
v_orb = []

# for n in range(0,N):
#     r_dist = round(random.uniform(10,300),3)
#     a = r_dist + round(random.uniform(0,30),3)

#     vo = round(math.sqrt(2*G*M_bh*((2/(r_dist))-(1/a))),4)


#     # Dual ax galaxy

#     if n<int(N/2):
#         r_dist = -1*r_dist 
#         vo = -1*vo
#         inf = [vo,r_dist]
#         v_orb.append(inf)
#         #print(n)

#     else:
#         inf = [vo,r_dist]
#         v_orb.append(inf)

for n in range(0,N):
    r_dist = round(random.uniform(50,400),3) -10
    a = r_dist + round(random.uniform(0,10),3)

    vo = round(7*math.sqrt(2*G*M_bh*((2/(r_dist))-(1/a))),4)


    if n<int(N/4):
        r_dist = -1*r_dist 
        vo = -1*vo
        inf = [vo,r_dist]
        v_orb.append(inf)
        #print(n)


    elif n>int(N/4) and n<int(N/2):
        inf = [vo,r_dist]
        v_orb.append(inf)



    elif n>int(N/2) and n<int(3*(N/4)):
        vo = -1*vo
        inf = [vo,r_dist]
        v_orb.append(inf)
        #print(n)


    else:

        r_dist = -1*r_dist 
        inf = [vo,r_dist]
        v_orb.append(inf)

#print(v_orb)

for n in range(0,N):
    #P.append(particle(n,vector(round(random.uniform(-5,5),3),round(random.uniform(-5,5))),vector(round(random.uniform(-3,3),3),round(random.uniform(-3,3),3))))
    #P.append(particle(n,vector(round(random.uniform(-50,50),3),round(random.uniform(-50,50))),vector(round(v_orb[n][0]),0),round(random.uniform(1,10),3)))
    
    #P.append(particle(n,vector(v_orb[n][1],0),vector(v1_g[0],v_orb[n][0]+v1_g[1]),round(random.uniform(1,10),3)))

    if n<int(N/2):
        P.append(particle(n,vector(v_orb[n][1],0),vector(v1_g[0],v_orb[n][0]+v1_g[1]),round(random.uniform(0.1,2),3)))

    else:
        P.append(particle(n,vector(0,v_orb[n][1]),vector(v_orb[n][0]+v1_g[0],v1_g[1]),round(random.uniform(0.11,2),3)))

    


P.append(particle(N,vector(-1,0),vector(0,-7),M_bh))

# Second black hole
P.append(particle(N+1,vector(5000,-2),vector(-10,0),M_bh))


quad_tree = node(np.array([0,-10]),800)


#----- Frame creation -------


k=0 # FOR DEBUGGING
for n in tqdm(t):
    quad_tree.reset()
    #print(quad_tree.child)

    # Creation of the quad tree

    for i in range(0,N+nmo):
        #print(P[i])
        quad_tree.inser_point(P[i])
    

    for i in range(0,N+nmo):
        A[i][0].append(P[i].r.x)
        A[i][1].append(P[i].r.y)

        # calculate the velocities using the barnes hunt optimization algorithm
        #print(P[i])
        #print("Accelerration calc")
        a_x,a_y = quad_tree.acceleration(P[i])
        #print("Particle")
        #print(i)
        #print([a_x,a_y])
        v_x,v_y = P[i].velocity(P,a_x,a_y)
        P[i].v = vector(v_x,v_y)
        x,y = P[i].position()
        P[i].r = vector(x,y)


        


    # --- DEBUGGING   Show quadtree 
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-500, 500), ylim=(-500, 500))
    rect = patches.Rectangle((quad_tree.r[0]-(quad_tree.s)/2,quad_tree.r[1]-(quad_tree.s)/2),quad_tree.s,quad_tree.s,linewidth=1, edgecolor='r', facecolor='none')
    plt.gca().add_patch(rect)
    draw_child(quad_tree)

    
    for n in range(N+nmo):
        ax.plot(A[n][0][k],A[n][1][k], '*', lw=2, markersize=2)

    plt.show()
    fig.clf()
    k+=1

#exit()

k=None

#print(P)
#----- Frame calibration -----

# Animation
fig = plt.figure()
fig.patch.set_facecolor('xkcd:black')
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-500, 500), ylim=(-500, 500))
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
rect = patches.Rectangle((quad_tree.r[0]-(quad_tree.s)/2,quad_tree.r[1]-(quad_tree.s)/2),quad_tree.s,quad_tree.s,linewidth=1, edgecolor='r', facecolor='none')
plt.gca().add_patch(rect)




def animate(i):

    # Time text
    time_text.set_text(time_template % (i*dt))

    #Animated particles
    result = [time_text]

    for n in range(N+nmo):
        particles[n].set_data(A[n][0][i],A[n][1][i])
        result.append(particles[n])
    # print(quad_array[i])

    # result.append(draw_child(quad_array[0]))
    
    return [result[i] for i in range(N+nmo+1)]



ani = animation.FuncAnimation(fig, animate, np.arange(1, len(A[0][1])),
                              interval=25, blit=True)
# ani.save('n_body.mp4', fps=15)
plt.show()
