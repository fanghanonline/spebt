import numpy as np
import matplotlib.pyplot as plt
import sys

inFname=sys.argv[1]
data = np.fromfile(inFname,dtype='single')
rows=data.shape[0]/4
myarray=data.reshape((int(rows),4))



# ax.set(xlim=(0,10), ylim=(0, 11), aspect="equal")

# det_indices=myarray[:,3]

points=myarray[np.nonzero(myarray[:,3]==15)][:,0:2]
print(points.shape)
im_template=np.zeros((200,200))
im_template[np.floor(points[:,0]*2).astype(int),np.floor(points[:,1]*2).astype(int)]=1.0
# # print(np.floor(points[:,0]))

fig, ax = plt.subplots(figsize=(13,10))
im=ax.imshow(im_template.T,cmap=plt.get_cmap("gray"),origin="lower")
ax.set(aspect="equal")
ax.set_xticks(np.arange(0,200,20))
ax.set_xticklabels(np.arange(0,100,10),size=20)
ax.set_yticks(np.arange(0,200,20))
ax.set_yticklabels(np.arange(0,100,10),size=20)
cbar=plt.colorbar(im)
cbar.ax.tick_params(labelsize=16) 
plt.tight_layout()
fig.savefig('%s.jpg'%inFname)