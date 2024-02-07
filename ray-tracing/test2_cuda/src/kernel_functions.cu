#include <cassert>
#include <helper.hpp>
#include <kernel_functions.h>

__device__ float center_to_line_dist_sqr_dev(float *cuboid, float *pointA,
                                             float *pointB) {

  // Vector OA, OB
  // O is the center of the cuboid
  float oa[3], ob[3];
  float ab_sqr = 0;
  for (int i : {0, 1, 2}) {
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

__global__ void add(int *a, int *b, int *c, int N) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < N) {
    c[tid] = a[tid] + b[tid];
  }
}

__global__ void get_center_to_line_dist_sqr_kernel(float *d_objects,
                                                   float *d_pointA,
                                                   float *d_pointB, int N,
                                                   float *d_result) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < N) {
    float object[8];
    for (int i = 0; i < 8; i++) {
      object[i] = d_objects[tid * 8 + i];
    }
    d_result[tid] = center_to_line_dist_sqr_dev(object, d_pointA, d_pointB);
  }
};

void get_center_to_line_dist_sqr(std::vector<std::vector<float>> geomDefinition,
                                 std::vector<float> pointA,
                                 std::vector<float> pointB) {
  size_t n_objects = geomDefinition.size();
  float *h_objects = (float *)malloc(n_objects * sizeof(float) * 8);
  for (auto object : geomDefinition) {
    size_t index = object[6];
    for (size_t i = 0; i < 8; i++) {
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

  get_center_to_line_dist_sqr_kernel<<<(n_objects + 255) / 256, 256>>>(
      d_objects, d_pointA, d_pointB, n_objects, d_result);

  cudaMemcpy(h_result, d_result, n_objects * sizeof(float),
             cudaMemcpyDeviceToHost);

  for (auto object : geomDefinition) {
    float cpu_distance = center_to_line_dist_sqr(object, pointA, pointB);
    // std::cout<<cpu_distance<<"\n";
    float gpu_distance = h_result[(int)object[6]];
    if (abs(cpu_distance - gpu_distance) < 1e-5)
      continue;
    printf("%3.7f != %3.7f\n", cpu_distance, gpu_distance);
    // // assert(1 == 1);
    // assert(cpu_distance == gpu_distance);
  }
}
void vec_add(int *a, int *b, int *c, int N) {
  int *dev_a, *dev_b, *dev_c;
  cudaMalloc((void **)&dev_a, N * sizeof(int));
  cudaMalloc((void **)&dev_b, N * sizeof(int));
  cudaMalloc((void **)&dev_c, N * sizeof(int));
  for (int i = 0; i < N; i++) {
    a[i] = i;
    b[i] = i;
  }
  cudaMemcpy(dev_a, a, N * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b, b, N * sizeof(int), cudaMemcpyHostToDevice);
  add<<<(N + 255) / 256, 256>>>(dev_a, dev_b, dev_c, N);

  cudaMemcpy(c, dev_c, N * sizeof(int), cudaMemcpyDeviceToHost);

  for (int i = 0; i < N; i++) {
    printf("%d + %d = %d\n", a[i], b[i], c[i]);
  }

  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);
}