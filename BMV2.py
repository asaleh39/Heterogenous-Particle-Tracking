from math import sqrt
from scipy.stats import norm
import numpy as np
from pylab import plot, show, grid, axis, xlabel, ylabel, title
import xlwt
import pandas as pd



#%% Brownian Function

def brownian(m, n, delta, dt, pixel):


    x = np.random.randint(1,pixel,[m,1])
    r = np.empty([m,n])
    delta = np.asarray(delta)
    for k in range(m):
        r[k,:] = norm.rvs(size=(1,n), scale=sqrt(2*delta[k]*dt))

    out = x + np.cumsum(r, axis=-1)

    return out

#%% My 2D

pixel = 512
n = 300
dt = 0.02
# Number of realizations to generate.
m = 150
# Create an empty array to store the realizations.
delta = np.random.rand(m,1)

x = brownian(m, n, delta, dt,pixel)
y = brownian(m, n, delta, dt,pixel)
######################################################

# Plot the 2D trajectory.
for k in range(m):
    plot(x[k], y[k])
    #plot(x[k, 0], y[k, 0], 'go')
    #plot(x[k, -1], y[k, -1], 'ro')
title('2D Brownian Motion')
xlabel('x')
ylabel('y')
axis([1,512,1,512])
grid(True)
show()

#%% Shift
xgood = x.copy()
ygood = y.copy()

shift = 0.15

for kidx1 in np.arange(0,m):
#    shift = np.random.random()
    for kidx2 in np.arange(1,n):
        xgood[kidx1,kidx2] = xgood[kidx1,kidx2] + kidx2*shift
        ygood[kidx1,kidx2] = ygood[kidx1,kidx2] + kidx2*shift


for kidx3 in range(m):
    plot(x[kidx3], y[kidx3],xgood[kidx3],ygood[kidx3])
    #plot(x[k, 0], y[k, 0], 'go')
    #plot(x[k, -1], y[k, -1], 'ro')
title('2D Brownian Motion with shift')
xlabel('x')
ylabel('y')
axis([1,512,1,512])
grid(True)
show()




#%% Write Data

#x2 = xgood.T
#y2 = ygood.T
# Use for printing shift data

x2 = x.T
y2 = y.T

numidx = np.arange(m*n).T

xdir = x2[:,0]

for k in range(1,m):
    xdir = np.concatenate((xdir,x2[:,k]))

ydir = y2[:,0]

for k in range(1,m):
    ydir = np.concatenate((ydir,y2[:,k]))
    
stepnum = np.arange(n).T

for k in range(1,m):
    stepnum = np.concatenate((stepnum,np.arange(n).T))
    
particle = 0 * np.ones((n))

for k in range(1,m):
    particle = np.concatenate((particle,k*np.ones(n).T))

wb = xlwt.Workbook()

sheet1 = wb.add_sheet('Sheet 1')

sheet1.write(0, 1, 'particle') 
sheet1.write(0, 2, 'frame') 
sheet1.write(0, 3, 'x') 
sheet1.write(0, 4, 'y') 



for k in range(m*n):
    sheet1.write(k+1,0,numidx[k].item())
    sheet1.write(k+1,1,particle[k].item())
    sheet1.write(k+1,2,stepnum[k].item())
    sheet1.write(k+1,3,xdir[k].item())
    sheet1.write(k+1,4,ydir[k].item())
    
wb.save('bmsimdata.xls')

#%% Read Data
partdata = pd.read_excel(r'C:\Users\abodb\OneDrive - University of North Carolina at Chapel Hill\Documents\Lai Lab Research\Shift Correction\06222020_100nm_G0F_Slide1_5ugml_02.tif.xlsx')
partdata = partdata.to_numpy()

x5 = partdata[:,3].copy()
y5 = partdata[:,4].copy()
partidx = partdata[:,1]

endpart = partdata[-1,1]

for h in np.arange(0,endpart+1):
    k = np.nonzero(partidx == h)
    x10 = x5[k]
    y10 = y5[k]

#^^ Indexes all the particles since they have difference dimensions they cant 
# be in one array. We can still work on each data set and corect for shift
# individually using this for loop and we can also pick with particle to
# work with if we so desire    
