#include <cassert>
#include <helper.hpp>
#include <kernel_functions.h>

__device__ float center_to_line_dist_sqr_device(float *cuboid, float *pointA,
                                                float *pointB)
{

  // Vector OA, OB
  // O is the center of the cuboid
  float oa[3], ob[3];
  float ab_sqr = 0;
  for (int i : {0, 1, 2})
  {
    float component_i = 0.5 * (cuboid[2 * i] + cuboid[2 * i + 1]);
    oa[i] = component_i - pointA[i];
    ob[i] = component_i - pointB[i];
    ab_sqr += (pointB[i] - pointA[i]) * (pointB[i] - pointA[i]);
  }
  // OA X OB coss product x, y, z components
  float cpx = oa[1] * ob[2] - oa[2] * ob[1];
  float cpy = oa[2] * ob[0] - oa[0] * ob[2];
  float cpz = oa[0] * ob[1] - oa[1] * ob[0];
  return (cpx * cpx + cpy * cpy + cpz * cpz) / ab_sqr;
}

__device__ float solid_angle_device(float *cuboid, float *unit_v3)
{
  float face_area[3];
  float edges[3];
  for (int i : {0, 1, 2})
  {
    edges[i] = cuboid[2 * i + 1] - cuboid[2 * i];
  }
  face_area[0] = edges[1] * edges[2];
  face_area[1] = edges[0] * edges[2];
  face_area[2] = edges[0] * edges[1];
  return fabsf(unit_v3[0]) * face_area[0] + fabsf(unit_v3[1]) * face_area[1] + fabsf(unit_v3[2]) * face_area[2];
}

__global__ void get_solid_angle_kernel(float *d_objects,
                                       float *d_pointA,
                                       float *d_pointB, int N,
                                       float *d_result)
{
  float unit_v3[3];
  float v3_rnorm = rnorm3df(d_pointB[0] - d_pointA[0], d_pointB[1] - d_pointA[1], d_pointB[2] - d_pointA[2]);
  for (size_t eid : {0, 1, 2})
  {
    unit_v3[eid] = (d_pointB[eid] - d_pointA[eid]) * v3_rnorm;
  }
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < N)
  {
    float object[8];
    for (int i = 0; i < 8; i++)
    {
      object[i] = d_objects[tid * 8 + i];
    }
    d_result[tid] = solid_angle_device(object, unit_v3);
  }
};

void get_solid_angle_host(std::vector<std::vector<float>> geomDefinition,
                          std::vector<float> pointA,
                          std::vector<float> pointB)
{
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
  ;
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

  get_solid_angle_kernel<<<(n_objects + 255) / 256, 256>>>(
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