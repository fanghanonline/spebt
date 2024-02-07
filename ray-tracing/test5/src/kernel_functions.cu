#include <cassert>
#include <helper.hpp>
#include <kernel_functions.h>

__device__ float center2line_distancesqr_device(float *cuboid, float *pA,
                                                float *pB)
{

  // Vector OA, OB
  // O is the center of the cuboid
  float oa[3], ob[3];
  float ab_sqr = 0;
  for (int i : {0, 1, 2})
  {
    float component_i = 0.5 * (cuboid[2 * i] + cuboid[2 * i + 1]);
    oa[i] = component_i - pA[i];
    ob[i] = component_i - pB[i];
    ab_sqr += (pB[i] - pA[i]) * (pB[i] - pA[i]);
  }
  // OA X OB coss product x, y, z components
  float cpx = oa[1] * ob[2] - oa[2] * ob[1];
  float cpy = oa[2] * ob[0] - oa[0] * ob[2];
  float cpz = oa[0] * ob[1] - oa[1] * ob[0];
  return (cpx * cpx + cpy * cpy + cpz * cpz) / ab_sqr;
}

__device__ float intersection_length_device(float *cuboid, float *pA, float *pB)
{
  float edges[3];
  float t_list[2];
  float vAB[3];
  size_t counter = 0;
  for (int i : {0, 1, 2})
  {
    edges[i] = cuboid[2 * i + 1] - cuboid[2 * i];
    vAB[i] = pB[i] - pA[i];
  }
  float diag_sqr = edges[0] * edges[0] + edges[1] * edges[1] + edges[2] * edges[2];
  if (diag_sqr * 0.25 - center2line_distancesqr_device(cuboid, pA, pB) > 0)
  {
    // Case 1: intersects on face x = x_0 or face x = x_1
    // Note that A_x never equals B_x.
    for (int i : {0, 1})
    {
      float yx =
          (cuboid[i] - pA[0]) / vAB[0] * vAB[1] + pA[1];
      float zx =
          (cuboid[i] - pA[0]) / vAB[0] * vAB[2] + pA[2];
      if ((yx - cuboid[2]) * (yx - cuboid[3]) < 0 and
          (zx - cuboid[4]) * (zx - cuboid[5]) < 0)
      {
        float t = (cuboid[i] - pA[0]) / (vAB[0]);
        if (t < 1)
        {
          t_list[counter] = t;
          counter++;
        }
      }
    }
    // Case 2: intersects on face y = y_0 or face y = y_1
    // Note: we exclude the case when A_y equals B_y
    if (vAB[1] != 0)
    {
      for (int i : {2, 3})
      {
        float xy =
            (cuboid[i] - pA[1]) / vAB[1] * vAB[0] + pA[0];
        float zy =
            (cuboid[i] - pA[1]) / (vAB[1]) * (vAB[2]) + pA[2];
        if ((xy - cuboid[0]) * (xy - cuboid[1]) < 0 and
            (zy - cuboid[4]) * (zy - cuboid[5]) < 0)
        {
          float t = (cuboid[i] - pA[1]) / (vAB[1]);
          if (t < 1)
          {
            t_list[counter] = t;
            counter++;
          }
        }
      }
    }

    // Case 3: intersects on face z = z_0 or face z = z_1
    // Note: we exclude the case when A_z equals B_z
    if (vAB[2] != 0)
    {
      for (int i : {4, 5})
      {
        float xz =
            (cuboid[i] - pA[2]) / (vAB[2]) * (vAB[0]) + pA[0];
        float yz =
            (cuboid[i] - pA[2]) / (vAB[2]) * (vAB[1]) + pA[1];
        if ((xz - cuboid[0]) * (xz - cuboid[1]) < 0 and
            (yz - cuboid[2]) * (yz - cuboid[3]) < 0)
        {
          float t = (cuboid[i] - pA[2]) / (vAB[2]);
          if (t < 1)
          {
            t_list[counter] = t;
            counter++;
          }
        }
      }
    }
  }

  float result = 0.0;
  if (counter > 1)
  {
    result = fabsf(t_list[0] - t_list[1]) * norm3df(vAB[0], vAB[1], vAB[2]);
  }
  return result;
}

__global__ void get_intersection_length_kernel(float *d_objects,
                                               float *d_pointA,
                                               float *d_pointB, int N,
                                               float *d_result)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < N)
  {
    float object[8];
    for (int i = 0; i < 8; i++)
    {
      object[i] = d_objects[tid * 8 + i];
    }
    d_result[tid] = intersection_length_device(object, d_pointA, d_pointB);
  }
}

void set_config_array(YAML::Node config, float *h_config)
{
  const std::vector<float> img_dims = config["image"]["dimension xyz"].as<std::vector<float>>();
  const std::vector<float> img_vpmms = config["image"]["voxel per mm xyz"].as<std::vector<float>>();
  const std::vector<int> img_subdivs = config["image"]["subdivision xyz"].as<std::vector<int>>();
  const std::vector<int> img_subdivs = config["detector"]["crystal n subdivision xyz"].as<std::vector<int>>();
  for (int i : {0, 1, 2})
  {
    h_config[i] = img_dims[i];
    h_config[i + 3] = img_vpmms[i];
  }
}
void get_result_host(YAML::Node config)
{
  const std::vector<std::vector<float>> geomDefinition =
      config["detector geometry"].as<std::vector<std::vector<float>>>();
  std::vector<float> pointA =
      config["Debug"]["Point A"].as<std::vector<float>>();
  std::vector<float> pointB =
      config["Debug"]["Point B"].as<std::vector<float>>();
  // std::vector<float> pointA(3);
  // std::vector<float> pointB(3);

  float h_config[12];

  size_t N_voxels_img = ceil(IMG_DIM_X_ * VXPMM_X_ * IMG_DIM_Y_ * VXPMM_Y_ * IMG_DIM_Z_ * VXPMM_Z_);
  printf("N Image Voxels: %lu\n", N_voxels_img);
  size_t n_objects = geomDefinition.size();
  float *h_objects = (float *)malloc(n_objects * sizeof(float) * 8);
  for (auto object : geomDefinition)
  {
    size_t index = object[6];
    for (size_t i = 0; i < 8; i++)
    {
      h_objects[index * 8 + i] = object[i];
    }
  }

  float *h_result = (float *)malloc(n_objects * sizeof(float));
  float *d_objects, *d_pointA, *d_pointB, *d_result;

  cudaMalloc((void **)&d_result, n_objects * sizeof(float));
  cudaMalloc((void **)&d_objects, n_objects * sizeof(float) * 8);

  cudaMalloc((void **)&d_pointA, 3 * sizeof(float));
  cudaMalloc((void **)&d_pointB, 3 * sizeof(float));
  cudaMemcpy(d_objects, h_objects, n_objects * sizeof(float) * 8,
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_pointA, pointA.data(), 3 * sizeof(float),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_pointB, pointB.data(), 3 * sizeof(float),
             cudaMemcpyHostToDevice);

  get_intersection_length_kernel<<<(n_objects + 255) / 256, 256>>>(
      d_objects, d_pointA, d_pointB, n_objects, d_result);

  cudaMemcpy(h_result, d_result, n_objects * sizeof(float),
             cudaMemcpyDeviceToHost);

  for (auto object : geomDefinition)
  {
    // float cpu_distance = center_to_line_dist_sqr(object, pointA, pointB);
    // std::cout<<cpu_distance<<"\n";
    // float gpu_distance = h_result[(int)object[6]];
    // if (abs(cpu_distance - gpu_distance) < 1e-5)
    //   continue;
    printf("%3.7f\n", h_result[(int)object[6]]);
    // // assert(1 == 1);
    // assert(cpu_distance == gpu_distance);
  }
}
