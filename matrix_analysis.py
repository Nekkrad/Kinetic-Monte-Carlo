import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.colors as colors
from matplotlib.animation import FuncAnimation

widths = [4,2]
heights = [4]

datas = pd.read_csv("Data.csv")
time = datas["Time"]
xdata = list(range(0,499))
for i in xdata:
    xdata[i] = xdata[i]


finalmegarray = []
for i in range(1,499):
    df = pd.read_csv(f"{i}.0.csv")
    df = df.drop(df.columns[0], axis=1)
    GA = [l for l in range(0,21)]
    LA = [k for k in range(0,21)]
    array = df.to_numpy()
    finalmegarray.append(array)

my_cmap = matplotlib.pyplot.get_cmap("viridis")
my_cmap.set_bad((0,0,0))
fig = plt.figure(figsize = (11,5))
gs = fig.add_gridspec(1, 2,width_ratios=widths,height_ratios=heights)
ax1 = fig.add_subplot(gs[:,0])
ax2 = fig.add_subplot(gs[:,1])
c = ax1.pcolor(finalmegarray[0], norm=matplotlib.colors.LogNorm(),cmap = my_cmap)
cbar = plt.colorbar(c,ax=ax1)
cbar.set_label('# of chains', rotation = 270,labelpad=15)
ax2.plot(xdata[0],time[0])

# We want to show all ticks...
ax1.set_xticks(GA)
ax1.set_yticks(LA)
# ... and label them with the respective list entries
ax1.set_xticklabels(GA)
ax1.set_yticklabels(LA)
ax1.set_xlabel("GA")
ax1.set_ylabel("LA")
# Rotate the tick labels and set their alignment.
plt.setp(ax1.get_xticklabels(), rotation=45,
         rotation_mode="anchor")
ax1.set_title("LA/GA ratio")
def animate(i):
    #ax.pcolormesh(finalmegarray[i])
    c.set_array(finalmegarray[i].ravel())
    ax2.cla() # clear the previous image
    ax2.plot(xdata[:i], time[:i]) # plot the line
    ax2.set_xlim([xdata[0], xdata[497]]) # fix the x axis
    ax2.set_ylim([time[0],time[497]]) # fix the y axis
    ax2.set_xlabel("Dump")
    ax2.set_ylabel("Time")


ani = FuncAnimation(fig,animate,frames=600,interval=20)
ani.save('animation.mp4')


plt.show()