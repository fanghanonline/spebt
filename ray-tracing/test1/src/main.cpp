#include "helper.hpp"

inline float vector_length(std::vector<float> v) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
std::vector<float> vector_cross_product(std::vector<float> va,
                                        std::vector<float> vb) {
  return std::vector<float>({va[1] * vb[2] - va[2] * vb[1],
                             va[2] * vb[0] - va[0] * vb[2],
                             va[0] * vb[1] - va[1] * vb[0]});
}

float center_to_line_dist_sqr(std::vector<float> cuboid,
                              std::vector<float> pointA,
                              std::vector<float> pointB) {

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

float intercept(std::vector<float> cuboid, std::vector<float> pA,
                std::vector<float> pB) {
  std::vector<float> t_list({0, 0});
  size_t counter = 0;
  float distance_sqr = center_to_line_dist_sqr(cuboid, pA, pB);
  float diagnal_sqr = 0;
  for (int i : {0, 1, 2}) {
    diagnal_sqr += (cuboid[2 * i + 1] - cuboid[2 * i]) *
                   (cuboid[2 * i + 1] - cuboid[2 * i]);
  }
  // The cuboid center to AB line distance should be smaller than half the
  // cuboid diagonal
  if (distance_sqr < diagnal_sqr * 0.25) {
    // std::cout
    //     << "Valid cuboid intersects with line AB! (Half diagnal - Distance: "
    //     << 0.5 * sqrt(diagnal_sqr) - sqrt(distance_sqr) << ")\n";
    // Case 1: intersects on face x = x_0 or face x = x_1
    // Note that A_x never equals B_x.
    for (int i : {0, 1}) {
      float yx =
          (cuboid[i] - pA[0]) / (pB[0] - pA[0]) * (pB[1] - pA[1]) + pA[1];
      float zx =
          (cuboid[i] - pA[0]) / (pB[0] - pA[0]) * (pB[2] - pA[2]) + pA[2];
      if ((yx - cuboid[2]) * (yx - cuboid[3]) < 0 and
          (zx - cuboid[4]) * (zx - cuboid[5]) < 0) {
        float t = (cuboid[i] - pA[0]) / (pB[0] - pA[0]);
        if (t < 1) {
          t_list[counter] = t;
          counter++;
        }
      }
    }
    // Case 2: intersects on face y = y_0 or face y = y_1
    // Note: we exclude the case when A_y equals B_y
    if (pB[1] - pA[1] != 0) {
      for (int i : {2, 3}) {
        float xy =
            (cuboid[i] - pA[1]) / (pB[1] - pA[1]) * (pB[0] - pA[0]) + pA[0];
        float zy =
            (cuboid[i] - pA[1]) / (pB[1] - pA[1]) * (pB[2] - pA[2]) + pA[2];
        if ((xy - cuboid[0]) * (xy - cuboid[1]) < 0 and
            (zy - cuboid[4]) * (zy - cuboid[5]) < 0) {
          float t = (cuboid[i] - pA[1]) / (pB[1] - pA[1]);
          if (t < 1) {
            t_list[counter] = t;
            counter++;
          }
        }
      }
    }

    // Case 3: intersects on face z = z_0 or face z = z_1
    // Note: we exclude the case when A_z equals B_z
    if (pB[1] - pA[1] != 0) {
      for (int i : {4, 5}) {
        float xz =
            (cuboid[i] - pA[2]) / (pB[2] - pA[2]) * (pB[0] - pA[0]) + pA[0];
        float yz =
            (cuboid[i] - pA[2]) / (pB[2] - pA[2]) * (pB[1] - pA[1]) + pA[1];
        if ((xz - cuboid[0]) * (xz - cuboid[1]) < 0 and
            (yz - cuboid[2]) * (yz - cuboid[3]) < 0) {
          float t = (cuboid[i] - pA[2]) / (pB[2] - pA[2]);
          if (t < 1) {
            t_list[counter] = t;
            counter++;
          }
        }
      }
    }
  }

  float result = 0.0;
  if (counter > 1) {
    // std::cout << "Find " << counter << " intersects!\n";
    float coords[2][3] = {{pA[0], pA[1], pA[2]}, {pA[0], pA[1], pA[2]}};
    for (int i = 0; i < counter; i++) {
      for (int j : {0, 1, 2}) {
        coords[i][j] = t_list[i] * (pB[j] - pA[j]) + pA[j];
      }
    }
    for (int i : {0, 1, 2}) {
      result += (coords[0][i] - coords[1][i]) * (coords[0][i] - coords[1][i]);
    }
  }
  result = sqrt(result);
  return result;
}

int block_counter(std::vector<float> cuboid, std::vector<float> pA,
                  std::vector<float> pB) {
  size_t counter = 0;
  float distance_sqr = center_to_line_dist_sqr(cuboid, pA, pB);
  float diagnal_sqr = 0;
  for (int i : {0, 1, 2}) {
    diagnal_sqr += (cuboid[2 * i + 1] - cuboid[2 * i]) *
                   (cuboid[2 * i + 1] - cuboid[2 * i]);
  }
  // The cuboid center to AB line distance should be smaller than half the
  // cuboid diagonal
  if (distance_sqr < diagnal_sqr * 0.25) {
    // Case 1: intersects on face x = x_0 or face x = x_1
    // Note that A_x never equals B_x.
    for (int i : {0, 1}) {
      float yx =
          (cuboid[i] - pA[0]) / (pB[0] - pA[0]) * (pB[1] - pA[1]) + pA[1];
      float zx =
          (cuboid[i] - pA[0]) / (pB[0] - pA[0]) * (pB[2] - pA[2]) + pA[2];
      if ((yx - cuboid[2]) * (yx - cuboid[3]) < 0 and
          (zx - cuboid[4]) * (zx - cuboid[5]) < 0) {
        float t = (cuboid[i] - pA[0]) / (pB[0] - pA[0]);
        if (t < 1) {
          counter++;
        }
      }
    }
    // Case 2: intersects on face y = y_0 or face y = y_1
    // Note: we exclude the case when A_y equals B_y
    if (pB[1] - pA[1] != 0) {
      for (int i : {2, 3}) {
        float xy =
            (cuboid[i] - pA[1]) / (pB[1] - pA[1]) * (pB[0] - pA[0]) + pA[0];
        float zy =
            (cuboid[i] - pA[1]) / (pB[1] - pA[1]) * (pB[2] - pA[2]) + pA[2];
        if ((xy - cuboid[0]) * (xy - cuboid[1]) < 0 and
            (zy - cuboid[4]) * (zy - cuboid[5]) < 0) {
          float t = (cuboid[i] - pA[1]) / (pB[1] - pA[1]);
          if (t < 1) {
            counter++;
          }
        }
      }
    }

    // Case 3: intersects on face z = z_0 or face z = z_1
    // Note: we exclude the case when A_z equals B_z
    if (pB[1] - pA[1] != 0) {
      for (int i : {4, 5}) {
        float xz =
            (cuboid[i] - pA[2]) / (pB[2] - pA[2]) * (pB[0] - pA[0]) + pA[0];
        float yz =
            (cuboid[i] - pA[2]) / (pB[2] - pA[2]) * (pB[1] - pA[1]) + pA[1];
        if ((xz - cuboid[0]) * (xz - cuboid[1]) < 0 and
            (yz - cuboid[2]) * (yz - cuboid[3]) < 0) {
          float t = (cuboid[i] - pA[2]) / (pB[2] - pA[2]);
          if (t < 1) {
            counter++;
          }
        }
      }
    }
  }
  return counter;
}

vec3_dim image_coordGrid_size(YAML::Node config) {
  // Set up image space coordinate grid
  float image_dimension_x_mm = 180.0;
  float image_dimension_y_mm = 180.0;
  float image_dimension_z_mm = 180.0;
  float image_pixelpermm_x = 1.0;
  float image_pixelpermm_y = 1.0;
  float image_pixelpermm_z = 1.0;
  image_dimension_x_mm = config["image"]["x dimension"].as<float>();
  image_dimension_y_mm = config["image"]["y dimension"].as<float>();
  image_dimension_z_mm = config["image"]["z dimension"].as<float>();

  image_pixelpermm_x = config["image"]["pixel per mm x"].as<float>();
  image_pixelpermm_y = config["image"]["pixel per mm y"].as<float>();
  image_pixelpermm_z = config["image"]["pixel per mm z"].as<float>();

  const size_t image_nPixels_x =
      ceil(image_dimension_x_mm * image_pixelpermm_x);
  const size_t image_nPixels_y =
      ceil(image_dimension_y_mm * image_pixelpermm_y);
  const size_t image_nPixels_z =
      ceil(image_dimension_z_mm * image_pixelpermm_z);
  return vec3_dim(image_nPixels_x, image_nPixels_y, image_nPixels_z);
}

std::vector<std::vector<float>>
coord_transform(std::vector<std::vector<float>> image_coordgrid,
                YAML::Node config) {
  float angle_rad = config["detector"]["detector rotation"].as<float>();
  float x_shift = config["detector"]["detector x-shift"].as<float>();
  float y_shift = 0.5 * config["detector"]["detector y-dimension"].as<float>();
  float trans_x = config["image"]["x dimension"].as<float>() * 0.5;
  float trans_y = config["image"]["y dimension"].as<float>() * 0.5;
  std::vector<std::vector<float>> transformed_coords;
  for (auto input_vec3 : image_coordgrid) {
    // First translational tranformation
    std::vector<float> vec3(3);
    vec3[0] = input_vec3[0] - trans_x;
    vec3[1] = input_vec3[1] - trans_y;
    // Rotational
    // Angle in radians
    vec3[0] = vec3[0] * cos(angle_rad) + vec3[1] * sin(angle_rad);
    vec3[1] = vec3[1] * cos(angle_rad) - vec3[0] * sin(angle_rad);
    // Translational transformation
    vec3[0] = vec3[0] - x_shift;
    // vec3[0] = vec3[0]
    vec3[1] = vec3[1] + y_shift;
    transformed_coords.push_back(vec3);
  }
  return transformed_coords;
}

std::vector<std::vector<float>> image_coordinate_grid(YAML::Node config) {
  // Set up image space coordinate grid
  float image_dimension_x_mm = 180.0;
  float image_dimension_y_mm = 180.0;
  float image_dimension_z_mm = 180.0;
  float image_pixelpermm_x = 1.0;
  float image_pixelpermm_y = 1.0;
  float image_pixelpermm_z = 1.0;
  image_dimension_x_mm = config["image"]["x dimension"].as<float>();
  image_dimension_y_mm = config["image"]["y dimension"].as<float>();
  image_dimension_z_mm = config["image"]["z dimension"].as<float>();

  image_pixelpermm_x = config["image"]["pixel per mm x"].as<float>();
  image_pixelpermm_y = config["image"]["pixel per mm y"].as<float>();
  image_pixelpermm_z = config["image"]["pixel per mm z"].as<float>();

  const size_t image_nPixels_x =
      ceil(image_dimension_x_mm * image_pixelpermm_x);
  const size_t image_nPixels_y =
      ceil(image_dimension_y_mm * image_pixelpermm_y);
  const size_t image_nPixels_z =
      ceil(image_dimension_z_mm * image_pixelpermm_z);
  // float **image_coordgrid = (float **)malloc(image_nPixels_x *
  // image_nPixels_y *
  //                                            image_nPixels_z * sizeof(float
  //                                            *));
  std::vector<std::vector<float>> image_coordgrid(
      image_nPixels_x * image_nPixels_y * image_nPixels_z,
      (std::vector<float>(3)));
  for (size_t i = 0; i < image_nPixels_x; i++) {
    for (size_t j = 0; j < image_nPixels_y; j++) {
      for (size_t k = 0; k < image_nPixels_z; k++) {
        size_t index =
            i * image_nPixels_x * image_nPixels_z + j * image_nPixels_z + k;
        // image_coordgrid[index] = (float *)malloc(3 * sizeof(float));
        image_coordgrid[index][0] = (1.0 * i + 0.5) / image_pixelpermm_x;
        image_coordgrid[index][1] = (1.0 * j + 0.5) / image_pixelpermm_y;
        image_coordgrid[index][2] = (1.0 * k + 0.5) / image_pixelpermm_z;
      }
    }
  }
  // printf("Size of the coordgrid: %lu\n", sizeof(&image_coordgrid));
  return image_coordgrid;
}

int main(int argc, char *argv[]) {
  std::cout << "Hello, from main!\n";
  std::string config_fname = argv[1];
  std::cout << "Using config file: " << config_fname << "\n";
  YAML::Node config = YAML::LoadFile(config_fname);
  //   const std::string username = config["username"].as<std::string>();
  const std::vector<std::vector<float>> geomDefinition =
      config["detector geometry"].as<std::vector<std::vector<float>>>();
  // float **image_coordgrid = image_coordinate_grid(config);
  std::vector<std::vector<float>> image_coordgrid =
      image_coordinate_grid(config);
  std::vector<std::vector<float>> transformed_image_coords =
      coord_transform(image_coordgrid, config);
  printf("Image voxel array size: %lu\n", transformed_image_coords.size());
  std::vector<int> xtal_nsubdivs =
      config["detector"]["crystal n subdivision x,y,z"].as<std::vector<int>>();
  float xtal_total_nsubdiv =
      xtal_nsubdivs[0] * xtal_nsubdivs[1] * xtal_nsubdivs[2];
  // std::vector<std::vector<float>> destination_coords;

  std::vector<std::vector<float>> xtal_subdiv_linearspace(
      xtal_total_nsubdiv, std::vector<float>(3));

  for (int iz = 0; iz < xtal_nsubdivs[2]; iz++) {
    for (int iy = 0; iy < xtal_nsubdivs[1]; iy++) {
      for (int ix = 0; ix < xtal_nsubdivs[0]; ix++) {
        size_t local_index = ix + iy * xtal_nsubdivs[0] +
                             iz * xtal_nsubdivs[0] * xtal_nsubdivs[1];
        xtal_subdiv_linearspace[local_index][0] = ix + 0.5;
        xtal_subdiv_linearspace[local_index][1] = iy + 0.5;
        xtal_subdiv_linearspace[local_index][2] = iz + 0.5;
      }
    }
  }

  // Separate the inner layers and the outmost layer of scintillators
  std::vector<std::vector<float>> inner_layers;
  std::vector<std::vector<float>> outmost_layer;
  float inner_layers_uplim =
      config["detector"]["inner layer x up-lim"].as<float>();
  std::cout << "inner layer x up-lim: " << inner_layers_uplim << "\n";
  for (std::vector<float> scintillator : geomDefinition) {
    if (scintillator[0] > inner_layers_uplim) {
      outmost_layer.push_back(scintillator);
    } else {
      inner_layers.push_back(scintillator);
    }
  }
  printf("inner layer N scintillators:   %lu\noutmost layer N scintillators: "
         "%lu\n",
         inner_layers.size(), outmost_layer.size());

  auto total_n_image_voxels = transformed_image_coords.size();
  std::vector<std::vector<float>> output(
      geomDefinition.size(), std::vector<float>(total_n_image_voxels));
  for (auto detector : geomDefinition) {
    std::vector<float> increments(3);
    for (size_t i = 0; i < 3; i++) {
      increments[i] =
          (detector[2 * i + 1] - detector[2 * i]) / xtal_nsubdivs[i];
    }

    for (size_t i_voxel = 0; i_voxel < total_n_image_voxels; i_voxel++) {
      auto pointA = transformed_image_coords[i_voxel];
      for (auto grid_point : xtal_subdiv_linearspace) {
        int n_blocks = 0;
        std::vector<float> pointB(3);
        for (int i : {0, 1, 2}) {
          pointB[i] = detector[i * 2] + increments[i] * grid_point[i];
        }
        for (std::vector<float> on_path_cuboid : inner_layers) {
          // float line_segment = intercept(on_path_cuboid, pointA,
          // pointB);
          n_blocks += block_counter(on_path_cuboid, pointA, pointB);
          // float line_segment = 0;
          // std::cout << "t: {" << t[0] << " ," << t[1] << "}\n";
        }
        if (n_blocks < 2) {
          output[int(detector[6])][i_voxel]=1.0;
        }
      }
    }
  }
  std::string outfname = "data_" + datime_str() + ".dat";
  std::cout << "Output file name: " << outfname << "\n";
  std::ofstream FILE(outfname, std::ios::out | std::ios::binary);
  if (!FILE) {
    std::cout << "Cannot open file!\n";
    return 1;
  }
  int outsize = output.size();
  printf("Output array size: %d\n", outsize);
  // Now write each vector one by one
  for (auto &v : output) {
    int size = v.size();
    // Store its contents
    FILE.write(reinterpret_cast<const char *>(&v[0]), v.size() * sizeof(float));
  }
  FILE.close();

  //
  //   int counter = 0;
  //   for (std::vector<float> coords : destination_coords) {
  //     counter++;
  //     printf("(%f, %f, %f) ", coords[0], coords[1], coords[2]);
  //     if (counter % xtal_nsubdivs[0] == 0) {
  //       printf("\n");
  //     }
  //   }
  // }

  return EXIT_SUCCESS;
}
