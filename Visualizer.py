import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pickle


location='/Users/edoardo/Desktop/simulazione_prova/record/'
pathcoord='2-5-0.001-4'+'.txt'
path=location+pathcoord

#path1 = '/Users/edoardo/Desktop/simulazione_prova/pickle_norm.txt'

time_multiplier=10

in_file=open(path, 'rb')
A=pickle.load(in_file)
Q_n, V_n, E_n, L_n ,dt= A[0], A[1], A[2], A[3], A[4]
a,b,c,d=Q_n.shape
Q_n=Q_n.reshape(a,b*c,d)



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



S_n=np.zeros((a//time_multiplier,b*c,d))
for i in range(len(Q_n)//time_multiplier):
    S_n[i]=Q_n[i*time_multiplier]
Q_n=S_n
    
steps, N, dim = Q_n.shape





'''
location='/Users/edoardo/Desktop/simulazione_prova/pickle_norm.txt'
in_file=open(location, 'rb')
Q_n=pickle.load(in_file).real

steps, N, dim = Q_n.shape
'''
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
        self.ax.axis([-20, 20, -20, 20])
        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,

    def data_stream(self):
        """Generate a random walk (brownian motion). Data is scaled to produce
        a soft "flickering" effect."""
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
            #x += 0.03 * (np.random.random((self.numpoints)) - 0.5)
            #y += 0.03 * (np.random.random((self.numpoints)) - 0.5)
            #s += 0.05 * (np.random.random(self.numpoints) - 0.5)
            #c += 0.02 * (np.random.random(self.numpoints) - 0.5)
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

        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,


if __name__ == '__main__':
    a = AnimatedScatter()
    plt.show()