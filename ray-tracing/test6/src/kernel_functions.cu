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

__device__ void pointA_device(float *pointA, unsigned int index, subdivs_config *config_shared)
{

  // Image N subvoxels x,y,z:         config_shared[0-2],
  // Image subvoxel dimension x,y,z:  config_shared[3-5]
  // Detector N subdivision x,y,z:    config_shared[12-14],
  unsigned int N_img_x = config_shared->img_N_sub_xyz[0];
  unsigned int N_img_y = config_shared->img_N_sub_xyz[1];
  unsigned int id_z = index / (N_img_x * N_img_y);
  unsigned int index_xy = index % (N_img_x * N_img_y);
  unsigned int id_y = index_xy / N_img_x;
  unsigned int id_x = index_xy % N_img_x;
  pointA[0] = (0.5 + id_x) * config_shared->img_dims_sub_xyz[0];
  pointA[1] = (0.5 + id_y) * config_shared->img_dims_sub_xyz[1];
  pointA[2] = (0.5 + id_z) * config_shared->img_dims_sub_xyz[2];
  // printf("Device: %d,%d,%d  %f, %f,%f\n", id_img_x, id_img_y, id_img_z, pointA[0], pointA[1], pointA[2]);
}

__device__ void pointB_device(float *pointB, unsigned int index, subdivs_config *config_shared, float *target_share)
{

  unsigned int id_z = index / (config_shared->det_N_sub_xyz[0] * config_shared->det_N_sub_xyz[1]);
  unsigned int id_xy = index % (config_shared->det_N_sub_xyz[0] * config_shared->det_N_sub_xyz[1]);
  unsigned int id_y = id_xy / config_shared->det_N_sub_xyz[1];
  unsigned int id_x = id_xy % config_shared->det_N_sub_xyz[1];
  float dim_subdiv_x = (target_share[1] - target_share[0]) / (float)config_shared->det_N_sub_xyz[0];
  float dim_subdiv_y = (target_share[3] - target_share[2]) / (float)config_shared->det_N_sub_xyz[1];
  float dim_subdiv_z = (target_share[5] - target_share[4]) / (float)config_shared->det_N_sub_xyz[2];
  pointB[0] = (0.5 + id_x) * dim_subdiv_x;
  pointB[1] = (0.5 + id_y) * dim_subdiv_y;
  pointB[2] = (0.5 + id_z) * dim_subdiv_z;
  // printf("Device: %d,%d,%d  %f, %f, %f PointB: (%f, %f, %f)\n", id_x, id_y, id_z, target_share[1] - target_share[0], target_share[3] - target_share[2], target_share[5] - target_share[4], pointB[0], pointB[1], pointB[2]);
}

__global__ void get_length_kernel(float *d_objects,
                                  subdivs_config *d_config,
                                  unsigned int N_img,
                                  unsigned int N_det,
                                  unsigned int N_objects,
                                  float *d_result)
{
  unsigned int tidx = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int tidy = threadIdx.y + blockIdx.y * blockDim.y;
  if (tidy < N_det and tidx < N_img)
  {
    unsigned int det_N_subs = d_config->det_N_sub_xyz[0] * d_config->det_N_sub_xyz[1] * d_config->det_N_sub_xyz[2];
    unsigned int target_index = tidy / det_N_subs;
    unsigned int subdiv_index = tidy % det_N_subs;
    float pointA[3], pointB[3], *target_share;
    pointA_device(pointA, tidx, d_config);
    target_share = &d_objects[target_index * 8];
    pointB_device(pointB, subdiv_index, d_config, target_share);
    // for (size_t i : {0, 1, 2})
    // {
    //   d_result[tidy * 3 + i] = pointB[i];
    // }

    for (unsigned int i = 0; i < N_objects; i++)
    {
      float *object;
      object = d_objects + i * 8;
      d_result[tidy * N_img + tidx] += intersection_length_device(object, pointA, pointB);
    }

    // printf("Device: ThreadID: (%d, %d) PointB: (%f, %f, %f)\n", tidx, tidy, pointB[0], pointB[1], pointB[2]);
  }
}

void set_config_array(YAML::Node config, subdivs_config *h_config)
{
  const std::vector<float> img_dims = config["image"]["dimension xyz"].as<std::vector<float>>();
  const std::vector<float> img_vpmms = config["image"]["voxel per mm xyz"].as<std::vector<float>>();
  const std::vector<int> img_subdivs = config["image"]["subdivision xyz"].as<std::vector<int>>();
  const std::vector<int> det_subdivs = config["detector"]["crystal n subdivision xyz"].as<std::vector<int>>();
  for (int i : {0, 1, 2})
  {
    // Image N subvoxels x,y,z: 
    h_config->img_N_sub_xyz[i] = ceil(img_dims[i] * img_vpmms[i] * img_subdivs[i]);
    // Image subvoxel dimension x,y,z
    h_config->img_dims_sub_xyz[i] = 1.0 / img_vpmms[i] / img_subdivs[i];
    h_config->det_N_sub_xyz[i] = det_subdivs[i];
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
  subdivs_config h_config, *d_config;
  // std::cout << "h_config size:" << sizeof(subdivs_config) << "\n";
  set_config_array(config, &h_config);
  size_t N_size_img = 1;
  for (int i : {0, 1, 2})
  {
    N_size_img *= h_config.img_N_sub_xyz[i];
  }
  printf("Image N subvoxels: %lu\n", N_size_img);
  std::vector<size_t> sensGeomIndex = config["detector"]["sensitive geometry indices"].as<std::vector<size_t>>();
  size_t N_objects = sensGeomIndex.size();
  float *h_objects = (float *)malloc(N_objects * sizeof(float) * 8);
  for (size_t i = 0; i < N_objects; i++)
  {
    auto geomObj = geomDefinition[sensGeomIndex[i]];
    size_t index = geomObj[6];
    for (size_t j = 0; j < 8; j++)
    {
      h_objects[i * 8 + j] = geomObj[j];
    }
  }

  float *d_objects, *d_result;
  cudaMalloc((void **)&d_config, sizeof(subdivs_config));
  // cudaMalloc((void **)&d_result, N_size_img * sizeof(float) * 3);

  cudaMalloc((void **)&d_objects, N_objects * sizeof(float) * 8);
  cudaMemcpy(d_config, &h_config, sizeof(subdivs_config),
             cudaMemcpyHostToDevice);

  // std::cout << "LINE: " << __LINE__ << "\n";
  // cudaMalloc((void **)&d_pointA, 3 * sizeof(float));
  // cudaMalloc((void **)&d_pointB, 3 * sizeof(float));
  cudaMemcpy(d_objects, h_objects, N_objects * sizeof(float) * 8,
             cudaMemcpyHostToDevice);
  // cudaMemcpy(d_pointA, pointA.data(), 3 * sizeof(float),
  //            cudaMemcpyHostToDevice);
  // cudaMemcpy(d_pointB, pointB.data(), 3 * sizeof(float),
  //            cudaMemcpyHostToDevice);

  // get_intersection_length_kernel<<<(n_objects + 255) / 256, 256>>>(
  //     d_objects, d_pointA, d_pointB, n_objects, d_result);
  unsigned int N_det = 21 * h_config.det_N_sub_xyz[0] * h_config.det_N_sub_xyz[1] * h_config.det_N_sub_xyz[2];
  unsigned int result_size = N_det * N_size_img * sizeof(float);

  float *h_result = (float *)malloc(result_size);
  cudaMalloc((void **)&d_result, result_size);
  dim3 dimBlock(32, 32);
  dim3 dimGrid((N_size_img + 31) / 32, (N_det + 31) / 32);
  get_length_kernel<<<dimGrid, dimBlock>>>(d_objects, d_config, N_size_img, N_det,
                                           N_objects, d_result);
  cudaMemcpy(h_result, d_result, result_size,
             cudaMemcpyDeviceToHost);
  for (size_t i = 0; i < N_det; i++)
  {
    for (size_t j = 0; j < N_size_img; j++)
    {
      printf("%f, ", h_result[i * N_size_img + j]);
    }
    std::cout << "\n";
  }
  // for (auto object : geomDefinition)
  // {
  //   // float cpu_distance = center_to_line_dist_sqr(object, pointA, pointB);
  //   // std::cout<<cpu_distance<<"\n";
  //   // float gpu_distance = h_result[(int)object[6]];
  //   // if (abs(cpu_distance - gpu_distance) < 1e-5)
  //   //   continue;
  //   printf("%3.7f\n", h_result[(int)object[6]]);
  //   // // assert(1 == 1);
  //   // assert(cpu_distance == gpu_distance);
  // }

  cudaFree(d_objects);
  cudaFree(d_result);
  cudaFree(d_config);
}
