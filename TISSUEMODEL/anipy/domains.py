import numpy as np
import matplotlib
import matplotlib.pyplot as plt
 
fig = plt.figure()
Fn="../Celltypes.dat"
print Fn
Nk=np.loadtxt(Fn)
#Domain number of nodes and variables 
Lx=300    #Nodes in x
Ly=150    #Nodes in y

#print Nk
Z=np.zeros((Lx*Ly))
j=0
for i in Nk:
	Z[j]=i[2]
	j=j+1


print Z
print len(Z)
raw_input()
dz=Z.reshape(Lx,Ly)

raw_input()

plt.imshow(dz,cmap='YlGnBu')
cax = plt.axes([0.80, 0.1, 0.04, 0.8])
plt.colorbar(cax=cax,  ticks=[1, 2])
zmin=0.9
zmax=2.1
plt.clim(zmin,zmax)
fig.suptitle('Tissue Geometry', fontsize=14, fontweight='bold')
plt.show()
fig.savefig('tissue.png', dpi=fig.dpi)
fig.savefig('tissue.svg', dpi=fig.dpi)
