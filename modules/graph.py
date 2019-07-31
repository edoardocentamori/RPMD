from setting import *
import matplotlib.pyplot as plt


def static_plot(q):
    x = q.swapaxes(0,-1)[0].flatten()
    y = q.swapaxes(0,-1)[1].flatten()
    fig, ax = plt.subplots()
    foo = large_scale
    ax.axis([-foo, foo, -foo, foo])
    ax.scatter(x,y,c='#2300A8')
    ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
    plt.show()