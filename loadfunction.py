from DataPointReader import readerx
from math import pi
from math import cos,sin
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate
import numpy
def coordinates():
    data_points = readerx()
    ca = 0.547
    la = 2.771
    ylist = []
    xlist = []
    zlist = []
    for i in range(1,(len(data_points)+1)):
        for j in range(1,(len(data_points[0])+1)):
            thetaz = (((i-1)*pi)/81)
            thetazplus = (((i)*pi)/81)
            zcoord = (-0.5)*((ca/2)*(1-cos(thetaz))+(ca/2)*(1-cos(thetazplus)))
            thetax = (((j-1)*pi)/41)
            thetaxplus = (((j)*pi)/41)
            xcoord = (-0.5)*((la/2)*(1-cos(thetax))+(la/2)*(1-cos(thetaxplus)))
            #print(zcoord)
            #print(xcoord)
            #print(data_points[i][j])
            #point = [data_points[i][j],zcoord,xcoord]
            ylist.append(data_points[i-1][j-1])
            xlist.append(xcoord)
            zlist.append(zcoord)
    return ylist,xlist,zlist
ylist,xlist,zlist = coordinates()
print(ylist)
fig = pyplot.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.scatter(xlist,zlist,ylist)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
pyplot.show()
