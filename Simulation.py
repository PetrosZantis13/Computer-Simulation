'''
Computer Simulation Project
Author: Petros Zantis

The following class, called Simulation, is using the animation 'on the fly' approach, 
to simulate and animate the orbits of the inner planets of the Solar System.
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import CelestialBody as CelestialBody
import numpy.linalg as linAlg
    
class Simulation(object):
    
    '''
    Below is the constructor of the class, taking as parameter the time-step value from main().
    The instance variables assigned are: the time-step, the 'allPlanets' list, a 'patches' list to 
    store the patches used in the animation, an 'energies' list to store the total energy of the system 
    throughout the simulation alongside a 'TEtime' list to store the times when an energy is appended in the 
    previous list, and a 'satDists' list to store all the distances of the satellite and Mars, alongside an
    'SDtime' list which stores the corresponding times. Then, the newAcceleration method is called upon the 
    created object, which in turn initialises the accelerations for all the bodies in the simulation.
    '''
    
    def __init__(self, dt) : 
        
        self.dt = dt
        self.allPlanets = self.readFile()
        self.patches = []
        self.energies = []
        self.TEtime = []
        self.satDists = []
        self.SDtime = []
        self.satDistsEarth = []
        #initial represents a 'boolean' to indicate that this is the initial acceleration
        self.newAcceleration(initial = 1) 
        
    '''
    Below is the newAcceleration method, which takes as argument a value named 'initial',
    which as explained above, represents a 'boolean' to indicate when is the initial acceleration
    (i.e. initial=1 ), and when not (i.e. initial=0 ). The functionality of the method is explained
    below.
    ''' 
    
    def newAcceleration(self, initial) :
             
        for pl1 in range(0, len(self.allPlanets)):  # for each planet in the list
            
            # an initialised value for the acceleration; an array of 2 zeros (x and y)
            acc = np.zeros(2) 
            
            for pl2 in range(0, len(self.allPlanets)):  # loop through the rest of the planets
                   
                if(pl1!=pl2):  # to ignore same planet
                     
                    otherMass = self.allPlanets[pl2].mass  # the other planet's mass
                    otherPos = self.allPlanets[pl2].position # the other planet's position
                    
                    # call the gravLaw method and sum the returned acceleration each time
                    acc += self.allPlanets[pl1].gravLaw(otherMass , otherPos)
            
            if(initial==1):  #if its the initial acceleration
               
                #a(t + dt) is assigned the value of the initial acceleration 
                self.allPlanets[pl1].accelerations[2] = acc
                #a(t) is assigned the same value as a(t + dt), which is a rough approximation to be later fixed 
                self.allPlanets[pl1].accelerations[1] = self.allPlanets[pl1].accelerations[2]
                #a(t - dt) is assigned the same value as a(t), which is a rough approximation to be later fixed  
                self.allPlanets[pl1].accelerations[0] = self.allPlanets[pl1].accelerations[1]
                
                initial = 0  #turn the 'boolean' to false, since the initial acceleration is now complete
                
            else:
                #a(t - dt) is assigned the previous value of a(t)
                self.allPlanets[pl1].accelerations[0] = self.allPlanets[pl1].accelerations[1]  
                #a(t) is assigned the previous value of a(t + dt)
                self.allPlanets[pl1].accelerations[1] = self.allPlanets[pl1].accelerations[2] 
                # a(t + dt) now takes the value of the above calculated acceleration
                self.allPlanets[pl1].accelerations[2] = acc                  
                          
    
    '''
    Below is the energy() method, which calculated the kinetic and potential energies of
    each planet in the system, and sums it up to a total energy, before appending it in the
    energies list. More functionality is explained below.
    ''' 
                    
    def energy(self):   
        
        totalEnergy = 0   #initialise the total energy to 0
        
        for pl1 in range(0, len(self.allPlanets)):  # for each planet in the list
            
            pl1Energy = 0    #initialise the planet's total energy to 0
            pl1Pos = self.allPlanets[pl1].position   # the position of the planet
            pl1Vel = self.allPlanets[pl1].velocity   # the velocity of the planet
            pl1Mass = self.allPlanets[pl1].mass      # the mass of the planet
            KE = (1/2)*(pl1Mass)*( linAlg.norm(pl1Vel) )**2      # kinetic energy = (1/2)*(m*v^2)
            PE = 0      #initialise the planet's potential energy to 0 
        
            for pl2 in range(0, len(self.allPlanets)):  #loop through the rest of the planets
                
                if(pl1!=pl2):  #to ignore same planet
                    
                    pl2Mass = self.allPlanets[pl2].mass     # the mass of the other planet
                    pl2Pos = self.allPlanets[pl2].position  # the position of the other planet
                    
                    Gconstant = (6.67408e-11)
                    relativePos = pl2Pos - pl1Pos
                    numerator = (-1/2)*(Gconstant * pl2Mass * pl1Mass)
                    denominator = linAlg.norm( relativePos )
                    PE += numerator/denominator    # potential energy = (-1/2)*(G*M*m)/r
                
            pl1Energy += (KE + PE)      #total energy of each planet = KE + PE 
            totalEnergy += pl1Energy    #total energy of system is the sum of all the planets' energies
               
        self.energies.append(totalEnergy)

    '''
    The orbit() method is used to calculate and print the orbital period of a planet, in addition
    to the ratio of this period to an Earth year.
    ''' 
    def orbit(self, pl, time):
        
        earthYear = 365.25      # an earth year in days
            
        # the name of the planet that completed an orbit
        name = str(self.allPlanets[pl].name)  
        # an index of the previous time the planet had completed an orbit
        prevIdx = len(self.allPlanets[pl].periodTimes) - 1
        #previous crossing time
        prevT = self.allPlanets[pl].periodTimes[prevIdx]
        
        # append the current time in the periodTimes of the planet
        self.allPlanets[pl].periodTimes.append(time)
        
        # subtract the previous time the planet orbited from the current time, to calculate the 'new' period
        # divide by (60*60*24) to convert from seconds to days
        period = (time - prevT)/(60*60*24)
        ratio = period / earthYear
        
        print("The orbital period of " + name +"is: " + str(period) + " days, and thus its ratio to an Earth Year is: " + str(ratio))
        
    '''
    Below is a simple method to calculate the distance between Mars and the satellite from Earth.
    That value is then added in to the satDists list, so that the closest approach can
    be determined at the end of the simulation.
    '''   
    def closestApproach(self, pl):  
    
        satPos = self.allPlanets[pl].position                #the satellite's position
        
        for pl2 in range(0 , len(self.allPlanets)):          #loop through the rest of the planets
            
            if(pl2 == 4):                                    #when Mars is found
                MarsPos = self.allPlanets[pl2].position
                distance = linAlg.norm( satPos - MarsPos )   #calculate the distance between them
                self.satDists.append(distance)
                
            if(pl2 == 3):                                    #when Earth is found
                EarthPos = self.allPlanets[pl2].position
                distance = linAlg.norm( satPos - EarthPos )   #calculate the distance between them
                self.satDistsEarth.append(distance)
                           
    def init(self):   # initialiser for animator
    
        return self.patches
    
    '''
    Below is the animate method which takes as argument the frameNo of the animation.
    That value times the time-step corresponds to the time that has passed after the 
    animation started. The rest of the functionality is explained below.
    '''     
    def animate(self, frameNo):
        
        # the time value is multiplied by 30, to create a faster but also smooth simulation
        time = frameNo * self.dt * 30  

        # a loop that runs the calculations 30 times before drawing a patch.
        # this results to a faster, smoother and more accurate simulation.
        for iter in range(0,30):
            
            for pl in range(0, len(self.allPlanets)):  # for each planet in the list
                
                # the x position before the update
                xPosBefore = self.allPlanets[pl].position[0] 
                # the y position before the update
                yPosBefore = self.allPlanets[pl].position[1]
                
                #call the upadatePosition() method from CelestialBody class, to update position
                self.allPlanets[pl].updatePos()
                
                # increment the time by the time-step
                time += self.dt
                
                # the x position after the update
                xPosAfter = self.allPlanets[pl].position[0]
                # the y position after the update
                yPosAfter = self.allPlanets[pl].position[1]
                
                # change the patch centre according to the new position
                self.patches[pl].center = (xPosAfter, yPosAfter) 
                
                '''
                To check if a planet completes a full orbit, two 'boolean' values are used, quadrant4 
                and quadrant1, each true if and only if the planet is located in the specific quadrant.
                So if the planet was in quadrant4 before the position update and in quadrant1 after the 
                position update, a full orbit is registered.
                ''' 
                quadrant4 = (xPosBefore>0) and (yPosBefore<0)   #x>0 and y<0
                quadrant1 = (xPosAfter>0) and (yPosAfter>0)     #x>0 and y>0
                
                if(quadrant4 and quadrant1):
                    
                    #register the orbit by calling the orbit() method
                    self.orbit(pl, time)
                
                #if the current body is the satellite to mars, call the closestApproach() method
                if( pl == 5):   
                    
                    #append the time value (in days) in the time list
                    self.SDtime.append(time/(60*60*24))
                    self.closestApproach(pl)
        
            #calculate the new acceleration for all the bodies, based on their new positions
            self.newAcceleration(initial=0)
        
            for pl in range(0, len(self.allPlanets)):  # for each planet in the list
                #calculate and update the new velocity, based on the new acceleration calculated above
                self.allPlanets[pl].updateVel()
        
        # call the energy() method which calculates and appends the total energy of the system, in the energies list
        self.energy()  
        #append the time value (in days) in the time list            
        self.TEtime.append(time/(60*60*24))  
      
        return self.patches 
    
    '''
    Below is the display() method, which creates and displays a MatPlotLib figure. Further explanations
    are given below.
    ''' 
        
    def display(self):
        
        #create the figure and axes
        fig = plt.figure()
        ax = plt.axes()
        
        '''
        'radii' is a list containing all the inner planets and the Satellite to Mars radii
        in meters, but exaggerated so as to be visible on the plot, as well as the Sun's but
        way less to again be visible on the plot.
        'colors' is a list containing some chosen colours for the planets (i.e. the patches)
        
        '''
        radii = [14e9,2.4e9,6.1e9,6.4e9,3.4e9,1e9]   
        colors = ['yellow','grey','peru','b','r','k']

        for i in range(0, len(self.allPlanets)):  #loop through all the planets
            
            xPos = self.allPlanets[i].position[0]  #planet's x position
            yPos = self.allPlanets[i].position[1]  #planet's y position

            '''
            The condition below allows for the inner planets to be drawn using the values for radii
            and colours stated in the lists above, but also allows extra planets to be added, using 
            a default radius value of 8e9 meters, and a dark orange colour.
            '''
            if(i < len(radii)):   
                planetRadius = radii[i]
                planetColor = colors[i]
            else:
                planetRadius = 8e9
                planetColor = 'darkorange'
            
            #append the patch in the list of patches and plot it
            self.patches.append(plt.Circle( ( xPos, yPos ) , planetRadius , color=planetColor , animated=True ))
            ax.add_patch(self.patches[i])
             
        ax.axis('scaled')    #scale the axis
        
        # appropriate limits chosen for a 'zoomed in' animation 
        ax.set_xlim(-2.7e11, 2.7e11)  
        ax.set_ylim(-2.7e11, 2.7e11)
        
        # FuncAnimation() method. 'frames' is a value for the total length of the simulation.
        anim = FuncAnimation(fig, self.animate , init_func = self.init, frames = 100000, repeat=False, interval=1, blit= True)
        
        plt.show()
    
    '''
    Below is the energyPlot() method, which simply plots all the values of total energy of
    the system, against the total time, and the satEarthPlot() which plots the distance between
    the satellite and Earth over time.
    Both methods are called after the simulation is ended.
    '''    
    def energyPlot(self):
        
        fileout = open("totalEnergy.txt" , "w")
        
        for i in range(0 , len(self.energies)):
            
            day = int(self.TEtime[i])     
            energy = self.energies[i]
            fileout.write("Total energy of the system at day " + str(day) + " is " + str(energy) + " Joules.\n")
              
        plt.plot(self.TEtime, self.energies)
        ax = plt.gca()
        #the values below are chosen for a zoomed out view, to indicate the consistency of the total energy
        ax.set_ylim( 2 * min(self.energies) , -min(self.energies) )      
        plt.title("Total energy of the system against time")  
        plt.ylabel("Total energy of the system in Joules")
        plt.xlabel("Total time in days")
        plt.show()
        
        fileout.close 
        
    def satEarthPlot(self):

        plt.plot(self.SDtime, self.satDistsEarth)
        plt.title("Satellite distance from Earth over time")  
        plt.ylabel("Satellite distance from Earth in meters")
        plt.xlabel("Total time in days")
        plt.show()
    
    '''
    The readFile() method. First, the file containing all the planet details and other body 
    details is opened in reading mode, and the lines are read. In the file, each planet corresponds 
    to 4 lines of data; its name, its mass, its initial position and initial velocity, in that order. 
    A list called allPlanets is created and returned, which stores all the CelestialBody objects that 
    are created. 
    '''    
    def readFile(self):
        
        filein = open("planetDetails.txt","r")   
        if filein.mode == 'r':
            lines = filein.readlines()
    
        planetNumber = len(lines)/4    # divide by 4 because there are 4 data for each planet
        allPlanets = []                # a list to store all planet objects
        
        for i in range(0, int(planetNumber)):   # loop through all planets in the file
            
            name = lines[4 * i].replace('\n' , ' ')   #replace the line separator with a space, for future 'printing' purposes
            mass = float(lines[4 * i + 1])
            posStr = lines[4 * i + 2].split(',')   # split the position string to x and y directions
            velStr = lines[4 * i + 3].split(',')   # split the velocity string to x and y directions
            
            # convert the above values to floats, and create NumPy arrays to store them as vectors; one for position and one for velocity
            initPos = np.array([float(posStr[0]),float(posStr[1])])
            initVel= np.array([float(velStr[0]),float(velStr[1])])
            
            #create and initialise a CelestialBody object as a planet, and append it in the list
            #(also includes the satellite for the experiment) 
            Planet = CelestialBody.CelestialBody(name, mass , initVel, initPos, self.dt )
            allPlanets.append(Planet)
            
        filein.close()   #close the file with the planet details
        
        return allPlanets
            
    '''
    Below is the main() function, which runs as soon as the program is launched. The method prompts
    for the time-step, to be later used in the implementation of the Beeman algorithm.
    Then, a simulation object called SolarSystem is created by passing the value of the time-step 
    as an argument, and then it is run and displayed using the display() method explained above, in 
    addition to calling the energyPlot() method explained above.
    After the simulation is closed, the closest approach of the satellite from Earth to Mars is
    calculated and printed to the console.
    '''
        
def main():
    
    promptStr1 = "Type in the time-step (dt) value. Note that a small time-step results in a more accurate "
    promptStr2 = "but slow simulation, while a larger one corresponds to a faster but less accurate simulation."
    
    #the time-step to be used in the calculations 
    timestep = float(input(promptStr1 + promptStr2 + "\n" + "Time-step (in seconds) = ")) 

    # create the simulation object, using the specified time-step
    SolarSystem = Simulation(timestep) 
    # run the simulation and display it
    SolarSystem.display()
    
    #find the minimum value within the distances list of the satellite to Mars
    closestApproach = min(SolarSystem.satDists)
    #find the corresponding time of that closest approach
    journeyidx =  SolarSystem.satDists.index(closestApproach)
    journeyTime = SolarSystem.SDtime[journeyidx]
    
    print( "In this simulation, the closest approach of the satellite from Earth to Mars was " + str(closestApproach) +" meters.")
    print( "This was achieved " + str(journeyTime) +" days after the simulation started.")
    
    # call the energyPlot() function when the simulation is closed
    SolarSystem.energyPlot()
    # call the satEarthPlot() function
    SolarSystem.satEarthPlot()

main()

