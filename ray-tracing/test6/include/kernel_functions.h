#include <vector>
#include <yaml-cpp/yaml.h>
void vec_add(int *, int *, int *, int);
void get_result_host(YAML::Node);
struct  subdivs_config
{
    unsigned int img_N_sub_xyz[3];
    unsigned int det_N_sub_xyz[3];
    float img_dims_sub_xyz[3];
};