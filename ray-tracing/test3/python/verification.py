import sys
import yaml
import numpy as np

configFname = sys.argv[1]
print(configFname)
with open(configFname, "r") as file:
    configs = yaml.safe_load(file)
detectors = np.array(configs["detector geometry"])
pointA = np.array(configs["Debug"]["Point A"])
pointB = np.array(configs["Debug"]["Point B"])
unit_v = (pointB - pointA)/np.linalg.norm(pointB - pointA)
det_edge_x = detectors[:, 1] - detectors[:, 0]
det_edge_y = detectors[:, 3] - detectors[:, 2]
det_edge_z = detectors[:, 5] - detectors[:, 4]

det_face_yz = det_edge_y * det_edge_z
det_face_zx = det_edge_z * det_edge_x
det_face_xy = det_edge_x * det_edge_y
det_faces = np.array([det_face_yz, det_face_zx, det_face_xy]).T
print(det_faces.shape)

print(unit_v)
print(np.matmul(det_faces,unit_v))