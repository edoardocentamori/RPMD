import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pickle

#This file is used to visualize the data

# Data coordinate

location='/Users/edoardo/Desktop/simulazione_prova/record/'
pathcoord='3-5-0.001-4'+'.txt'

#5-10-0.01-1  #example of system without thermalization #use time_multiplier=1
#3-5-0.001-4  #example of system with constraint #use time_multiplier=5
#1-10-0.001-18 #example of system with thermalization #use time_multiplier=1

path=location+pathcoord

time_multiplier=5 #speed up the video

# Data extraction

in_file=open(path, 'rb')
A=pickle.load(in_file)
in_file.close() 



Q_n, V_n, E_n, L_n ,dt= A[0], A[1], A[2], A[3], A[4]
a,b,c,d=Q_n.shape
Q_n=Q_n.reshape(a,b*c,d)

# Energy and Angular Momentum plots

fig, ((ax11,ax22),(ax12,ax21)) = plt.subplots(2,2)
x = np.linspace(0,len(Q_n)*dt,len(Q_n)-3)
E_n=E_n[3:]
L_n=L_n[3:]
ax11.scatter(x,E_n,c='yellow',s=0.1)
ax11.set_title('energy')
ax22.scatter(x,np.abs(L_n),c='red',s=0.1)
ax22.set_title('angular momentum')

ax12.scatter(x,E_n,c='yellow',s=0.1)
ax21.scatter(x,np.abs(L_n),c='red',s=0.1)
ax12.set_yscale('log')
ax21.set_yscale('log')
ax12.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax21.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
plt.show()

# time scaling : simply ignore some data to show things faster

S_n=np.zeros((a//time_multiplier,b*c,d))
for i in range(len(Q_n)//time_multiplier):
    S_n[i]=Q_n[i*time_multiplier]
Q_n=S_n
    
steps, N, dim = Q_n.shape


# The animated object

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, numpoints=N):
        self.numpoints = numpoints
        self.stream = self.data_stream()

        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots()
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=5, 
                                          init_func=self.setup_plot, blit=True)
    
        
    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        x, y, s, c = next(self.stream).T
        self.scat = self.ax.scatter(x, y, c=c, s=s, vmin=0, vmax=1,
                                    cmap="jet", edgecolor="k")
        self.ax.axis([-100, 100, -100, 100])
        return self.scat,

    def data_stream(self):
        """ data stream generation """
        x = Q_n.T[0]
        y = Q_n.T[1]
        s = np.ones(self.numpoints).T*0.5
        c = np.random.random(self.numpoints).T
        t=0
        while True:
            t+=1
            t%=steps
            x = Q_n[t].T[0] 
            y = Q_n[t].T[1] 
            yield np.c_[x, y, s, c]

    def update(self, i):
        """Update the scatter plot."""
        data = next(self.stream)

        # Set x and y data...
        self.scat.set_offsets(data[:, :2])
        # Set sizes...
        self.scat.set_sizes(300 * abs(data[:, 2])**1.5 + 100)
        # Set colors..
        self.scat.set_array(data[:, 3])

        return self.scat,


if __name__ == '__main__':
    a = AnimatedScatter()
    plt.show()