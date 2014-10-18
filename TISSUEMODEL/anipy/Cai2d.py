import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
        comment='Movie support!')
writer = FFMpegWriter(fps=15, metadata=metadata)

fig = plt.figure()
To=6000
Tf=20000
tmk=[]
zi=To
while zi<=Tf:
	tmk.append(zi)
	zi+=5

print tmk
raw_input()

#Initial Setup
Sn='Cai'
Fn="../snaps/"+str(To)+".000000"+Sn+".dat"
print Fn
Nk=np.loadtxt(Fn)

with writer.saving(fig, "CAI2D.mp4", 100):
	plt.imshow(Nk)#, extent = [0,50, 100, 0] )
	cax = plt.axes([0.80, 0.1, 0.04, 0.8])
	zmin=0.0000
	zmax=0.0014
	plt.clim(zmin,zmax)
	plt.colorbar(cax=cax)
	#plt.colorbar()
	fig.suptitle('T= '+str(To), fontsize=14, fontweight='bold')
        writer.grab_frame()
	plt.clf()
	#for i in range(1+To,1+Tf):
	for i in tmk:
		Fn="../snaps/"+str(i)+".000000"+Sn+".dat"
		print Fn
		Nk=np.loadtxt(Fn)
		plt.imshow(Nk)#, extent = [0,50, 100, 0] )
		cax = plt.axes([0.80, 0.1, 0.04, 0.8])
		plt.colorbar(cax=cax)
		plt.clim(zmin,zmax)
		#plt.colorbar(cax=cax)
		#plt.colorbar()
		fig.suptitle('T= '+str(i), fontsize=14, fontweight='bold')
        	writer.grab_frame()
		plt.clf()
