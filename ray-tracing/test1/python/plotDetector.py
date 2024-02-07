import sys
import yaml
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
configFname=sys.argv[1]
print(configFname)
with open(configFname, 'r') as file:
    configs = yaml.safe_load(file)
detectors=np.array(configs["detector geometry"])
detectors=detectors[1:-1]
# print(detectors.shape)
# print(detectors)
det_xy=np.array([detectors[:,0],detectors[:,2]]).T
det_inc_xy=np.array([(detectors[:,1]-detectors[:,0]),(detectors[:,3]-detectors[:,2])]).T
# print(det_xy)
# print("increment")
# print(det_inc_xy)
print(detectors.shape)
target=detectors[np.nonzero(detectors[:,6]==15)].flatten()
print(target)
target_xy=[target[0],target[2]]
target_inc_xy=[target[1]-target[0],target[3]-target[2]]
fig,ax= plt.subplots(figsize=(6,10))
target_rect=plt.Rectangle(target_xy,target_inc_xy[0],target_inc_xy[1],color='r')

rect_list = [Rectangle(xy,inc_xy[0],inc_xy[1]) for xy,inc_xy in zip(det_xy,det_inc_xy)]
# rect_list=[Rectangle((1,1),2,2)]
pc = PatchCollection(rect_list,ec='r')
ax.add_collection(pc)
ax.set_xticks(np.arange(1,14,3))
ax.set_xticklabels(np.arange(1,14,3),size=20)
ax.set_yticks(np.arange(0,25,3))
ax.set_yticklabels(np.arange(0,25,3),size=20)
ax.set_xlabel("detector X dimension (mm)",size=20)
ax.set_ylabel("detector Y dimension (mm)",size=20)
ax.plot(1,1)
ax.grid()
ax.set(aspect='equal')
ax.add_patch(target_rect)
plt.tight_layout()
fig.savefig('diagram.png')