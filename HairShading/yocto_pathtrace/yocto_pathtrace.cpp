//
// Implementation for Yocto/RayTrace.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#include "yocto_pathtrace.h"

#include <yocto/yocto_color.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_parallel.h>
#include <yocto/yocto_sampling.h>
#include <yocto/yocto_shading.h>
#include <yocto/yocto_shape.h>

// hairshading step 1: import library as in the original paper work
//#include <numeric>

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SCENE EVALUATION
// -----------------------------------------------------------------------------
namespace yocto {

//////////////////////////////////////////////////////////////////////////////////////////////////////
// hairshading step 2
// init the base params for geometry of HairBSDF
static const int   pMax = 3;      // segments of scattered light: check eval
const float        h    = 0.0f;   // [-1,1] offset along the curve width
static const float eta  = 1.75f;  // refraction index for hair interior
vec3f sigma = {2.9, 5.1, 10.2};   // absorption coefficient (+=black, -=white)
static const float beta_m = 0.125f;  // longitudinal hair roughness
static const float beta_n = 0.3f;    // azimuthal hair roughness
float sin2kAlpha[3], cos2kAlpha[3];  // angles where scales are offset

// main important params for longitudinal and azimuthal scattering
float gammaO = asin(h);  // angle between surface normal and direction
float v[pMax + 1];       // roughness variation trick (L)
static const float SqrtPiOver8 = 0.626657069f;  // constant used for mapping (A)
float              s           = SqrtPiOver8 *
          (0.265f * beta_n + 1.194f * beta_n * beta_n +
              5.372f * pow(beta_n, 20));  // logistic scale factor from beta(A)

//////////////////////////////////////////////////////////////////////////////////////////////////////

// Check texture size
static vec2i texture_size(const pathtrace_texture* texture) {
  if (!texture->hdr.empty()) {
    return texture->hdr.imsize();
  } else if (!texture->ldr.empty()) {
    return texture->ldr.imsize();
  } else {
    return zero2i;
  }
}

// Evaluate a texture
static vec4f lookup_texture(const pathtrace_texture* texture, const vec2i& ij,
    bool ldr_as_linear = false) {
  if (!texture->hdr.empty()) {
    return texture->hdr[ij];
  } else if (!texture->ldr.empty()) {
    return ldr_as_linear ? byte_to_float(texture->ldr[ij])
                         : srgb_to_rgb(byte_to_float(texture->ldr[ij]));
  } else {
    return {1, 1, 1};
  }
}

// Evaluate a texture
static vec4f eval_texture(const pathtrace_texture* texture, const vec2f& uv,
    bool ldr_as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false) {
  // get texture
  if (!texture) return {1, 1, 1};

  // get yimage width/height
  auto size = texture_size(texture);

  // get coordinates normalized for tiling
  auto s = 0.0f, t = 0.0f;
  if (clamp_to_edge) {
    s = clamp(uv.x, 0.0f, 1.0f) * size.x;
    t = clamp(uv.y, 0.0f, 1.0f) * size.y;
  } else {
    s = fmod(uv.x, 1.0f) * size.x;
    if (s < 0) s += size.x;
    t = fmod(uv.y, 1.0f) * size.y;
    if (t < 0) t += size.y;
  }

  // get yimage coordinates and residuals
  auto i = clamp((int)s, 0, size.x - 1), j = clamp((int)t, 0, size.y - 1);
  auto ii = (i + 1) % size.x, jj = (j + 1) % size.y;
  auto u = s - i, v = t - j;

  if (no_interpolation) return lookup_texture(texture, {i, j}, ldr_as_linear);

  // handle interpolation
  return lookup_texture(texture, {i, j}, ldr_as_linear) * (1 - u) * (1 - v) +
         lookup_texture(texture, {i, jj}, ldr_as_linear) * (1 - u) * v +
         lookup_texture(texture, {ii, j}, ldr_as_linear) * u * (1 - v) +
         lookup_texture(texture, {ii, jj}, ldr_as_linear) * u * v;
}

// Generates a ray from a camera for yimage plane coordinate uv and
// the lens coordinates luv.
static ray3f eval_camera(const pathtrace_camera* camera, const vec2f& image_uv,
    const vec2f& lens_uv) {
  auto q  = vec3f{camera->film.x * (0.5f - image_uv.x),
      camera->film.y * (image_uv.y - 0.5f), camera->lens};
  auto dc = -normalize(q);
  auto e  = vec3f{
      lens_uv.x * camera->aperture / 2, lens_uv.y * camera->aperture / 2, 0};
  auto p = dc * camera->focus / abs(dc.z);
  auto d = normalize(p - e);
  return ray3f{
      transform_point(camera->frame, e), transform_direction(camera->frame, d)};
}

// Samples a camera ray at pixel ij of an image of size size with puv and luv
// as random numbers for pixel and lens respectively
static ray3f sample_camera(const pathtrace_camera* camera, const vec2i& ij,
    const vec2i& size, const vec2f& puv, const vec2f& luv) {
  return eval_camera(camera,
      vec2f{(ij.x + puv.x) / size.x, (ij.y + puv.y) / size.y},
      sample_disk(luv));
}

// Eval position
static vec3f eval_position(
    const pathtrace_instance* instance, int element, const vec2f& uv) {
  auto shape = instance->shape;
  if (!shape->triangles.empty()) {
    auto t = shape->triangles[element];
    return transform_point(
        instance->frame, interpolate_triangle(shape->positions[t.x],
                             shape->positions[t.y], shape->positions[t.z], uv));
  } else if (!shape->lines.empty()) {
    auto l = shape->lines[element];
    return transform_point(instance->frame,
        interpolate_line(shape->positions[l.x], shape->positions[l.y], uv.x));
  } else if (!shape->points.empty()) {
    return transform_point(
        instance->frame, shape->positions[shape->points[element]]);
  } else {
    return zero3f;
  }
}

// Shape element normal.
static vec3f eval_element_normal(
    const pathtrace_instance* instance, int element) {
  auto shape = instance->shape;
  if (!shape->triangles.empty()) {
    auto t = shape->triangles[element];
    return transform_normal(
        instance->frame, triangle_normal(shape->positions[t.x],
                             shape->positions[t.y], shape->positions[t.z]));
  } else if (!shape->lines.empty()) {
    auto l = shape->lines[element];
    return transform_normal(instance->frame,
        line_tangent(shape->positions[l.x], shape->positions[l.y]));
  } else if (!shape->points.empty()) {
    return {0, 0, 1};
  } else {
    return {0, 0, 0};
  }
}

// Eval normal
static vec3f eval_normal(
    const pathtrace_instance* instance, int element, const vec2f& uv) {
  auto shape = instance->shape;
  if (shape->normals.empty()) return eval_element_normal(instance, element);
  if (!shape->triangles.empty()) {
    auto t = shape->triangles[element];
    return transform_normal(
        instance->frame, normalize(interpolate_triangle(shape->normals[t.x],
                             shape->normals[t.y], shape->normals[t.z], uv)));
  } else if (!shape->lines.empty()) {
    auto l = shape->lines[element];
    return transform_normal(instance->frame,
        normalize(
            interpolate_line(shape->normals[l.x], shape->normals[l.y], uv.x)));
  } else if (!shape->points.empty()) {
    return transform_normal(
        instance->frame, normalize(shape->normals[shape->points[element]]));
  } else {
    return zero3f;
  }
}

// Eval texcoord
static vec2f eval_texcoord(
    const pathtrace_instance* instance, int element, const vec2f& uv) {
  auto shape = instance->shape;
  if (shape->texcoords.empty()) return uv;
  if (!shape->triangles.empty()) {
    auto t = shape->triangles[element];
    return interpolate_triangle(shape->texcoords[t.x], shape->texcoords[t.y],
        shape->texcoords[t.z], uv);
  } else if (!shape->lines.empty()) {
    auto l = shape->lines[element];
    return interpolate_line(shape->texcoords[l.x], shape->texcoords[l.y], uv.x);
  } else if (!shape->points.empty()) {
    return shape->texcoords[shape->points[element]];
  } else {
    return zero2f;
  }
}

// Shape element normal.
static pair<vec3f, vec3f> eval_element_tangents(
    const pathtrace_instance* instance, int element) {
  auto shape = instance->shape;
  if (!shape->triangles.empty() && !shape->texcoords.empty()) {
    auto t        = shape->triangles[element];
    auto [tu, tv] = triangle_tangents_fromuv(shape->positions[t.x],
        shape->positions[t.y], shape->positions[t.z], shape->texcoords[t.x],
        shape->texcoords[t.y], shape->texcoords[t.z]);
    return {transform_direction(instance->frame, tu),
        transform_direction(instance->frame, tv)};
  } else {
    return {};
  }
}

static vec3f eval_normalmap(
    const pathtrace_instance* instance, int element, const vec2f& uv) {
  // YOUR CODE GOES HERE ----------------------------SECONDA RICHIESTA
  auto shape      = instance->shape;
  auto normal_tex = instance->material->normal_tex;
  // apply normal mapping
  auto normal   = eval_normal(instance, element, uv);
  auto texcoord = eval_texcoord(instance, element, uv);
  if (normal_tex && !shape->triangles.empty()) {
    auto normalmap = -1 + 2 * eval_texture(normal_tex, texcoord, true);
    auto [tu, tv]  = eval_element_tangents(instance, element);
    auto frame     = frame3f{tu, tv, normal, zero3f};
    frame.x        = orthonormalize(frame.x, frame.z);
    frame.y        = normalize(cross(frame.z, frame.x));
    auto flip_v    = dot(frame.y, tv) < 0;
    normalmap.y *= flip_v ? 1 : -1;  // flip vertical axis
    normal = transform_normal(frame, xyz(normalmap));
  }
  return normal;
}

// Eval shading normal
static vec3f eval_shading_normal(const pathtrace_instance* instance,
    int element, const vec2f& uv, const vec3f& outgoing) {
  auto shape  = instance->shape;
  auto normal = eval_normal(instance, element, uv);
  if (!shape->triangles.empty()) {
    if (instance->material->normal_tex) {
      normal = eval_normalmap(instance, element, uv);
    }
    if (!instance->material->thin) return normal;
    return dot(normal, outgoing) >= 0 ? normal : -normal;
  } else if (!shape->lines.empty()) {
    return orthonormalize(outgoing, normal);
  } else if (!shape->points.empty()) {
    return -outgoing;
  } else {
    return zero3f;
  }
}

// Brdf
struct pathtrace_brdf {
  // brdf lobes
  vec3f diffuse      = {0, 0, 0};
  vec3f specular     = {0, 0, 0};
  vec3f metal        = {0, 0, 0};
  vec3f transmission = {0, 0, 0};
  vec3f refraction   = {0, 0, 0};
  float roughness    = 0;
  float opacity      = 1;
  float ior          = 1;
  vec3f meta         = {0, 0, 0};
  vec3f metak        = {0, 0, 0};
  // weights
  float diffuse_pdf      = 0;
  float specular_pdf     = 0;
  float metal_pdf        = 0;
  float transmission_pdf = 0;
  float refraction_pdf   = 0;
};

// Eval material to obatain emission, brdf and opacity.
static vec3f eval_emission(const pathtrace_instance* instance, int element,
    const vec2f& uv, const vec3f& normal, const vec3f& outgoing) {
  auto material = instance->material;
  auto texcoord = eval_texcoord(instance, element, uv);
  return material->emission *
         xyz(eval_texture(material->emission_tex, texcoord));
}

// Eval material to obatain emission, brdf and opacity.
static pathtrace_brdf eval_brdf(const pathtrace_instance* instance, int element,
    const vec2f& uv, const vec3f& normal, const vec3f& outgoing) {
  auto material = instance->material;
  auto texcoord = eval_texcoord(instance, element, uv);
  auto brdf     = pathtrace_brdf{};
  auto color    = material->color *
               xyz(eval_texture(material->color_tex, texcoord, false));
  auto specular = material->specular *
                  eval_texture(material->specular_tex, texcoord, true).x;
  auto metallic = material->metallic *
                  eval_texture(material->metallic_tex, texcoord, true).x;
  auto roughness = material->roughness *
                   eval_texture(material->roughness_tex, texcoord, true).x;
  auto ior          = material->ior;
  auto transmission = material->transmission *
                      eval_texture(material->emission_tex, texcoord, true).x;
  auto opacity = material->opacity *
                 eval_texture(material->opacity_tex, texcoord, true).x;
  auto thin = material->thin || material->transmission == 0;

  // factors
  auto weight = vec3f{1, 1, 1};
  brdf.metal  = weight * metallic;
  weight *= 1 - metallic;
  brdf.refraction = thin ? zero3f : weight * transmission;
  weight *= 1 - (thin ? 0 : transmission);
  brdf.specular = weight * specular;
  weight *= 1 - specular * fresnel_dielectric(ior, outgoing, normal);
  brdf.transmission = weight * transmission * color;
  weight *= 1 - transmission;
  brdf.diffuse   = weight * color;
  brdf.meta      = reflectivity_to_eta(color);
  brdf.metak     = zero3f;
  brdf.roughness = roughness * roughness;
  brdf.ior       = ior;
  brdf.opacity   = opacity;

  // textures
  if (brdf.diffuse != zero3f || brdf.roughness != 0) {
    brdf.roughness = clamp(brdf.roughness, 0.03f * 0.03f, 1.0f);
  }
  if (brdf.specular == zero3f && brdf.metal == zero3f &&
      brdf.transmission == zero3f && brdf.refraction == zero3f) {
    brdf.roughness = 1;
  }
  if (brdf.opacity > 0.999f) brdf.opacity = 1;

  // weights
  brdf.diffuse_pdf  = max(brdf.diffuse);
  brdf.specular_pdf = max(
      brdf.specular * fresnel_dielectric(brdf.ior, normal, outgoing));
  brdf.metal_pdf = max(
      brdf.metal * fresnel_conductor(brdf.meta, brdf.metak, normal, outgoing));
  brdf.transmission_pdf = max(brdf.transmission);
  brdf.refraction_pdf   = max(brdf.refraction);
  auto pdf_sum = brdf.diffuse_pdf + brdf.specular_pdf + brdf.metal_pdf +
                 brdf.transmission_pdf + brdf.refraction_pdf;
  if (pdf_sum != 0) {
    brdf.diffuse_pdf /= pdf_sum;
    brdf.specular_pdf /= pdf_sum;
    brdf.metal_pdf /= pdf_sum;
    brdf.transmission_pdf /= pdf_sum;
    brdf.refraction_pdf /= pdf_sum;
  }
  return brdf;
}

// check if a brdf is a delta
static bool is_delta(const pathtrace_brdf& brdf) { return brdf.roughness == 0; }

// vsdf
struct pathtrace_vsdf {
  vec3f density    = {0, 0, 0};
  vec3f scatter    = {0, 0, 0};
  float anisotropy = 0;
};

// evaluate volume
static pathtrace_vsdf eval_vsdf(
    const pathtrace_instance* instance, int element, const vec2f& uv) {
  auto material = instance->material;
  // initialize factors
  auto texcoord = eval_texcoord(instance, element, uv);
  auto base     = material->color *
              xyz(eval_texture(material->color_tex, texcoord, false));
  auto transmission = material->transmission *
                      eval_texture(material->emission_tex, texcoord, true).x;
  auto thin = material->thin || material->transmission == 0;
  auto scattering =
      material->scattering *
      xyz(eval_texture(material->scattering_tex, texcoord, false));
  auto scanisotropy = material->scanisotropy;
  auto trdepth      = material->trdepth;

  // factors
  auto vsdf    = pathtrace_vsdf{};
  vsdf.density = (transmission != 0 && !thin)
                     ? -log(clamp(base, 0.0001f, 1.0f)) / trdepth
                     : zero3f;
  vsdf.scatter    = scattering;
  vsdf.anisotropy = scanisotropy;

  return vsdf;
}

// check if we have a volume
static bool has_volume(const pathtrace_instance* instance) {
  return !instance->material->thin && instance->material->transmission != 0;
}

// Evaluate all environment color.
static vec3f eval_environment(const pathtrace_scene* scene, const ray3f& ray) {
  auto emission = zero3f;
  for (auto environment : scene->environments) {
    auto wl       = transform_direction(inverse(environment->frame), ray.d);
    auto texcoord = vec2f{
        atan2(wl.z, wl.x) / (2 * pif), acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
    if (texcoord.x < 0) texcoord.x += 1;
    emission += environment->emission *
                xyz(eval_texture(environment->emission_tex, texcoord));
  }
  return emission;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SHAPE/SCENE BVH
// -----------------------------------------------------------------------------
namespace yocto {

// primitive used to sort bvh entries
struct pathtrace_bvh_primitive {
  bbox3f bbox      = invalidb3f;
  vec3f  center    = zero3f;
  int    primitive = 0;
};

// Splits a BVH node. Returns split position and axis.
static pair<int, int> split_middle(
    vector<pathtrace_bvh_primitive>& primitives, int start, int end) {
  // initialize split axis and position
  auto axis = 0;
  auto mid  = (start + end) / 2;

  // compute primintive bounds and size
  auto cbbox = invalidb3f;
  for (auto i = start; i < end; i++) cbbox = merge(cbbox, primitives[i].center);
  auto csize = cbbox.max - cbbox.min;
  if (csize == zero3f) return {mid, axis};

  // split along largest
  if (csize.x >= csize.y && csize.x >= csize.z) axis = 0;
  if (csize.y >= csize.x && csize.y >= csize.z) axis = 1;
  if (csize.z >= csize.x && csize.z >= csize.y) axis = 2;

  // split the space in the middle along the largest axis
  mid = (int)(std::partition(primitives.data() + start, primitives.data() + end,
                  [axis, middle = center(cbbox)[axis]](auto& primitive) {
                    return primitive.center[axis] < middle;
                  }) -
              primitives.data());

  // if we were not able to split, just break the primitives in half
  if (mid == start || mid == end) {
    // throw runtime_error("bad bvh split");
    mid = (start + end) / 2;
  }

  return {mid, axis};
}

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// Build BVH nodes
static void build_bvh(vector<pathtrace_bvh_node>& nodes,
    vector<pathtrace_bvh_primitive>&              primitives) {
  // prepare to build nodes
  nodes.clear();
  nodes.reserve(primitives.size() * 2);

  // queue up first node
  auto queue = deque<vec3i>{{0, 0, (int)primitives.size()}};
  nodes.emplace_back();

  // create nodes until the queue is empty
  while (!queue.empty()) {
    // grab node to work on
    auto next = queue.front();
    queue.pop_front();
    auto nodeid = next.x, start = next.y, end = next.z;

    // grab node
    auto& node = nodes[nodeid];

    // compute bounds
    node.bbox = invalidb3f;
    for (auto i = start; i < end; i++)
      node.bbox = merge(node.bbox, primitives[i].bbox);

    // split into two children
    if (end - start > bvh_max_prims) {
      // get split
      auto [mid, axis] = split_middle(primitives, start, end);

      // make an internal node
      node.internal = true;
      node.axis     = (byte)axis;
      node.num      = 2;
      node.start    = (int)nodes.size();
      nodes.emplace_back();
      nodes.emplace_back();
      queue.push_back({node.start + 0, start, mid});
      queue.push_back({node.start + 1, mid, end});
    } else {
      // Make a leaf node
      node.internal = false;
      node.num      = (short)(end - start);
      node.start    = start;
    }
  }

  // cleanup
  nodes.shrink_to_fit();
}

static void init_bvh(pathtrace_shape* shape, const pathtrace_params& params) {
  // build primitives
  auto primitives = vector<pathtrace_bvh_primitive>{};
  if (!shape->points.empty()) {
    for (auto idx = 0; idx < shape->points.size(); idx++) {
      auto& p             = shape->points[idx];
      auto& primitive     = primitives.emplace_back();
      primitive.bbox      = point_bounds(shape->positions[p], shape->radius[p]);
      primitive.center    = center(primitive.bbox);
      primitive.primitive = idx;
    }
  } else if (!shape->lines.empty()) {
    for (auto idx = 0; idx < shape->lines.size(); idx++) {
      auto& l         = shape->lines[idx];
      auto& primitive = primitives.emplace_back();
      primitive.bbox = line_bounds(shape->positions[l.x], shape->positions[l.y],
          shape->radius[l.x], shape->radius[l.y]);
      primitive.center    = center(primitive.bbox);
      primitive.primitive = idx;
    }
  } else if (!shape->triangles.empty()) {
    for (auto idx = 0; idx < shape->triangles.size(); idx++) {
      auto& primitive = primitives.emplace_back();
      auto& t         = shape->triangles[idx];
      primitive.bbox  = triangle_bounds(
          shape->positions[t.x], shape->positions[t.y], shape->positions[t.z]);
      primitive.center    = center(primitive.bbox);
      primitive.primitive = idx;
    }
  }

  // build nodes
  if (shape->bvh) delete shape->bvh;
  shape->bvh = new pathtrace_bvh_tree{};
  build_bvh(shape->bvh->nodes, primitives);

  // set bvh primitives
  shape->bvh->primitives.reserve(primitives.size());
  for (auto& primitive : primitives) {
    shape->bvh->primitives.push_back(primitive.primitive);
  }
}

void init_bvh(pathtrace_scene* scene, const pathtrace_params& params,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{0, 1 + (int)scene->shapes.size()};

  // shapes
  for (auto idx = 0; idx < scene->shapes.size(); idx++) {
    if (progress_cb) progress_cb("build shape bvh", progress.x++, progress.y);
    init_bvh(scene->shapes[idx], params);
  }

  // handle progress
  if (progress_cb) progress_cb("build scene bvh", progress.x++, progress.y);

  // instance bboxes
  auto primitives  = vector<pathtrace_bvh_primitive>{};
  auto instance_id = 0;
  for (auto instance : scene->instances) {
    auto& primitive = primitives.emplace_back();
    primitive.bbox  = instance->shape->bvh->nodes.empty()
                         ? invalidb3f
                         : transform_bbox(instance->frame,
                               instance->shape->bvh->nodes[0].bbox);
    primitive.center    = center(primitive.bbox);
    primitive.primitive = instance_id++;
  }

  // build nodes
  if (scene->bvh) delete scene->bvh;
  scene->bvh = new pathtrace_bvh_tree{};
  build_bvh(scene->bvh->nodes, primitives);

  // set bvh primitives
  scene->bvh->primitives.reserve(primitives.size());
  for (auto& primitive : primitives) {
    scene->bvh->primitives.push_back(primitive.primitive);
  }

  // handle progress
  if (progress_cb) progress_cb("build bvh", progress.x++, progress.y);
}

// Intersect ray with a bvh->
static bool intersect_shape_bvh(pathtrace_shape* shape, const ray3f& ray_,
    int& element, vec2f& uv, float& distance, bool find_any) {
  // get bvh and shape pointers for fast access
  auto bvh = shape->bvh;

  // check empty
  if (bvh->nodes.empty()) return false;

  // node stack
  int  node_stack[128];
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // shared variables
  auto hit = false;

  // copy ray to modify it
  auto ray = ray_;

  // prepare ray for fast queries
  auto ray_dinv  = vec3f{1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z};
  auto ray_dsign = vec3i{(ray_dinv.x < 0) ? 1 : 0, (ray_dinv.y < 0) ? 1 : 0,
      (ray_dinv.z < 0) ? 1 : 0};

  // walking stack
  while (node_cur) {
    // grab node
    auto& node = bvh->nodes[node_stack[--node_cur]];

    // intersect bbox
    // if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;
    if (!intersect_bbox(ray, ray_dinv, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // for internal nodes, attempts to proceed along the
      // split axis from smallest to largest nodes
      if (ray_dsign[node.axis]) {
        node_stack[node_cur++] = node.start + 0;
        node_stack[node_cur++] = node.start + 1;
      } else {
        node_stack[node_cur++] = node.start + 1;
        node_stack[node_cur++] = node.start + 0;
      }
    } else if (!shape->points.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& p = shape->points[shape->bvh->primitives[idx]];
        if (intersect_point(
                ray, shape->positions[p], shape->radius[p], uv, distance)) {
          hit      = true;
          element  = shape->bvh->primitives[idx];
          ray.tmax = distance;
        }
      }
    } else if (!shape->lines.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& l = shape->lines[shape->bvh->primitives[idx]];
        if (intersect_line(ray, shape->positions[l.x], shape->positions[l.y],
                shape->radius[l.x], shape->radius[l.y], uv, distance)) {
          hit      = true;
          element  = shape->bvh->primitives[idx];
          ray.tmax = distance;
        }
      }
    } else if (!shape->triangles.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& t = shape->triangles[shape->bvh->primitives[idx]];
        if (intersect_triangle(ray, shape->positions[t.x],
                shape->positions[t.y], shape->positions[t.z], uv, distance)) {
          hit      = true;
          element  = shape->bvh->primitives[idx];
          ray.tmax = distance;
        }
      }
    }

    // check for early exit
    if (find_any && hit) return hit;
  }

  return hit;
}

// Intersect ray with a bvh->
static bool intersect_scene_bvh(const pathtrace_scene* scene, const ray3f& ray_,
    int& instance, int& element, vec2f& uv, float& distance, bool find_any,
    bool non_rigid_frames) {
  // get bvh and scene pointers for fast access
  auto bvh = scene->bvh;

  // check empty
  if (bvh->nodes.empty()) return false;

  // node stack
  int  node_stack[128];
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // shared variables
  auto hit = false;

  // copy ray to modify it
  auto ray = ray_;

  // prepare ray for fast queries
  auto ray_dinv  = vec3f{1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z};
  auto ray_dsign = vec3i{(ray_dinv.x < 0) ? 1 : 0, (ray_dinv.y < 0) ? 1 : 0,
      (ray_dinv.z < 0) ? 1 : 0};

  // walking stack
  while (node_cur) {
    // grab node
    auto& node = bvh->nodes[node_stack[--node_cur]];

    // intersect bbox
    // if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;
    if (!intersect_bbox(ray, ray_dinv, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // for internal nodes, attempts to proceed along the
      // split axis from smallest to largest nodes
      if (ray_dsign[node.axis]) {
        node_stack[node_cur++] = node.start + 0;
        node_stack[node_cur++] = node.start + 1;
      } else {
        node_stack[node_cur++] = node.start + 1;
        node_stack[node_cur++] = node.start + 0;
      }
    } else {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto instance_ = scene->instances[scene->bvh->primitives[idx]];
        auto inv_ray   = transform_ray(
            inverse(instance_->frame, non_rigid_frames), ray);
        if (intersect_shape_bvh(
                instance_->shape, inv_ray, element, uv, distance, find_any)) {
          hit      = true;
          instance = scene->bvh->primitives[idx];
          ray.tmax = distance;
        }
      }
    }

    // check for early exit
    if (find_any && hit) return hit;
  }

  return hit;
}

// Intersect ray with a bvh->
static bool intersect_instance_bvh(const pathtrace_instance* instance,
    const ray3f& ray, int& element, vec2f& uv, float& distance, bool find_any,
    bool non_rigid_frames) {
  auto inv_ray = transform_ray(inverse(instance->frame, non_rigid_frames), ray);
  return intersect_shape_bvh(
      instance->shape, inv_ray, element, uv, distance, find_any);
}

pathtrace_intersection intersect_scene_bvh(const pathtrace_scene* scene,
    const ray3f& ray, bool find_any, bool non_rigid_frames) {
  auto intersection = pathtrace_intersection{};
  intersection.hit  = intersect_scene_bvh(scene, ray, intersection.instance,
      intersection.element, intersection.uv, intersection.distance, find_any,
      non_rigid_frames);
  return intersection;
}
pathtrace_intersection intersect_instance_bvh(
    const pathtrace_instance* instance, const ray3f& ray, bool find_any,
    bool non_rigid_frames) {
  auto intersection = pathtrace_intersection{};
  intersection.hit = intersect_instance_bvh(instance, ray, intersection.element,
      intersection.uv, intersection.distance, find_any, non_rigid_frames);
  return intersection;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate emission
static vec3f eval_emission(
    const vec3f& emission, const vec3f& normal, const vec3f& outgoing) {
  return emission;
}

// Evaluates/sample the BRDF scaled by the cosine of the incoming direction.
static vec3f eval_brdfcos(const pathtrace_brdf& brdf, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (brdf.roughness == 0) return zero3f;

  // accumulate the lobes
  auto brdfcos = zero3f;
  if (brdf.diffuse != zero3f) {
    brdfcos += brdf.diffuse *
               eval_diffuse_reflection(normal, outgoing, incoming);
  }
  if (brdf.specular != zero3f) {
    brdfcos += brdf.specular * eval_microfacet_reflection(brdf.ior,
                                   brdf.roughness, normal, outgoing, incoming);
  }
  if (brdf.metal != zero3f) {
    brdfcos += brdf.metal * eval_microfacet_reflection(brdf.meta, brdf.metak,
                                brdf.roughness, normal, outgoing, incoming);
  }
  if (brdf.transmission != zero3f) {
    brdfcos += brdf.transmission * eval_microfacet_transmission(brdf.ior,
                                       brdf.roughness, normal, outgoing,
                                       incoming);
  }
  if (brdf.refraction != zero3f) {
    brdfcos += brdf.refraction * eval_microfacet_refraction(brdf.ior,
                                     brdf.roughness, normal, outgoing,
                                     incoming);
  }
  return brdfcos;
}

static vec3f eval_delta(const pathtrace_brdf& brdf, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (brdf.roughness != 0) return zero3f;

  auto brdfcos = zero3f;

  if (brdf.specular != zero3f && brdf.refraction == zero3f) {
    brdfcos += brdf.specular *
               eval_delta_reflection(brdf.ior, normal, outgoing, incoming);
  }
  if (brdf.metal != zero3f) {
    brdfcos += brdf.metal * eval_delta_reflection(brdf.meta, brdf.metak, normal,
                                outgoing, incoming);
  }
  if (brdf.transmission != zero3f) {
    brdfcos += brdf.transmission *
               eval_delta_transmission(brdf.ior, normal, outgoing, incoming);
  }
  if (brdf.refraction != zero3f) {
    brdfcos += brdf.refraction *
               eval_delta_refraction(brdf.ior, normal, outgoing, incoming);
  }

  return brdfcos;
}

/*
// Picks a direction based on the BRDF
static vec3f sample_brdfcos(const pathtrace_brdf& brdf, const vec3f& normal,
    const vec3f& outgoing, float rnl, const vec2f& rn) {
  if (brdf.roughness == 0) return zero3f;

  auto cdf = 0.0f;

  if (brdf.diffuse_pdf != 0) {
    cdf += brdf.diffuse_pdf;
    if (rnl < cdf) return sample_diffuse_reflection(normal, outgoing, rn);
  }

  if (brdf.specular_pdf != 0 && brdf.refraction_pdf == 0) {
    cdf += brdf.specular_pdf;
    if (rnl < cdf)
      return sample_microfacet_reflection(
          brdf.ior, brdf.roughness, normal, outgoing, rn);
  }

  if (brdf.metal_pdf != 0) {
    cdf += brdf.metal_pdf;
    if (rnl < cdf)
      return sample_microfacet_reflection(
          brdf.meta, brdf.metak, brdf.roughness, normal, outgoing, rn);
  }

  if (brdf.transmission_pdf != 0) {
    cdf += brdf.transmission_pdf;
    if (rnl < cdf)
      return sample_microfacet_transmission(
          brdf.ior, brdf.roughness, normal, outgoing, rn);
  }

  if (brdf.refraction_pdf != 0) {
    cdf += brdf.refraction_pdf;
    if (rnl < cdf)
      return sample_microfacet_refraction(
          brdf.ior, brdf.roughness, normal, outgoing, rnl, rn);
  }

  return zero3f;
}
*/
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// hairshading step3
// Given the paper work here there is an adaptation of it
inline float Sqr(float v) { return v * v; }

inline float Phi(int p, float gammaO, float gammaT) {
  return 2 * p * gammaT - 2 * gammaO + p * pi;
}

float FrDielectric(float cosThetaI, float etaI, float etaT) {
  cosThetaI = clamp(cosThetaI, -1.0f, 1.0f);
  // Potentially swap indices of refraction
  bool entering = cosThetaI > 0.f;
  if (!entering) {
    std::swap(etaI, etaT);
    cosThetaI = abs(cosThetaI);
  }

  // Compute _cosThetaT_ using Snell's law
  float sinThetaI = sqrt(max((float)0, 1 - cosThetaI * cosThetaI));
  float sinThetaT = etaI / etaT * sinThetaI;

  // Handle total internal reflection
  if (sinThetaT >= 1) return 1;
  float cosThetaT = std::sqrt(std::max((float)0, 1 - sinThetaT * sinThetaT));
  float Rparl     = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                ((etaT * cosThetaI) + (etaI * cosThetaT));
  float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
                ((etaI * cosThetaI) + (etaT * cosThetaT));
  return (Rparl * Rparl + Rperp * Rperp) / 2;
}
static std::array<vec3f, pMax + 1> Ap(
    float cosThetaO, float eta, float h, const vec3f& T) {
  std::array<vec3f, pMax + 1> ap;
  // Compute $p=0$ attenuation at initial cylinder intersection
  float cosGammaO = sqrt(1 - h * h);
  float cosTheta  = cosThetaO * cosGammaO;
  float f         = FrDielectric(cosTheta, 1.f, eta);
  ap[0]           = {f, f, f};

  // Compute $p=1$ attenuation term
  ap[1] = (1 - f) * (1 - f) * T;

  // Compute attenuation terms up to $p=_pMax_$
  for (int p = 2; p < pMax; ++p) ap[p] = ap[p - 1] * T * f;

  // Compute attenuation term accounting for remaining orders of scattering
  // ap[pMax] = ap[pMax - 1] * f * T / (vec3f(1.f) - T * f); //error1
  return ap;
}
inline float I0(float x) {
  float   val   = 0;
  float   x2i   = 1;
  int64_t ifact = 1;
  int     i4    = 1;
  // I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
  for (int i = 0; i < 10; ++i) {
    if (i > 1) ifact *= i;
    val += x2i / (i4 * Sqr(ifact));
    x2i *= x * x;
    i4 *= 4;
  }

  return val;
}

inline float LogI0(float x) {
  if (x > 12)
    return x + 0.5 * (-log(2 * pi) + log(1 / x) + 1 / (8 * x));
  else
    return log(I0(x));
}

inline float LogisticCDF(float x, float s) {
  return 1 / (1 + exp(-x / s));
  //
}

inline float Logistic(float x, float s) {
  x = abs(x);
  return exp(-x / s) / (s * Sqr(1 + exp(-x / s)));
}

inline float TrimmedLogistic(float x, float s, float a, float b) {
  return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
}

static float Mp(float cosThetaI, float cosThetaO, float sinThetaI,
    float sinThetaO, float v) {
  float a  = cosThetaI * cosThetaO / v;
  float b  = sinThetaI * sinThetaO / v;
  float mp = (v <= 0.1f)
                 ? (exp(LogI0(a) - b - 1 / v + 0.6931f + log(1 / (2 * v))))
                 : (exp(-b) * I0(a)) / (sinh(1 / v) * 2 * v);

  return mp;
}
inline float AbsCosTheta(const vec3f& w) { return abs(w.z); }

inline float Np(float phi, int p, float s, float gammaO, float gammaT) {
  float dphi = phi - Phi(p, gammaO, gammaT);
  // Remap _dphi_ to $[-\pi,\pi]$
  while (dphi > pi) dphi -= 2 * pi;
  while (dphi < -pi) dphi += 2 * pi;
  return TrimmedLogistic(dphi, s, -pi, pi);
}

static vec3f f(const pathtrace_brdf& brdf, const vec3f& wo, const vec3f& wi) {
  // Compute hair coordinate system terms related to _wo_
  float sinThetaO = wo.x;
  float cosThetaO = sqrt(1 - Sqr(sinThetaO));
  float phiO      = std::atan2(wo.z, wo.y);

  // Compute hair coordinate system terms related to _wi_
  float sinThetaI = wi.x;
  float cosThetaI = sqrt(1 - Sqr(sinThetaI));
  float phiI      = std::atan2(wi.z, wi.y);

  // Compute $\cos \thetat$ for refracted ray
  float sinThetaT = sinThetaO / eta;
  float cosThetaT = sqrt(1 - Sqr(sinThetaT));

  // Compute $\gammat$ for refracted ray
  float etap      = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
  float sinGammaT = h / etap;
  float cosGammaT = sqrt(1 - Sqr(sinGammaT));
  float gammaT    = asin(sinGammaT);

  // Compute the transmittance _T_ of a single path through the cylinder
  vec3f T = exp(-vec3f(sigma) * (2 * cosGammaT / cosThetaT));

  // Evaluate hair BSDF
  float                       phi = phiI - phiO;
  std::array<vec3f, pMax + 1> ap  = Ap(cosThetaO, eta, h, T);
  vec3f                       fsum;
  for (int p = 0; p < pMax; ++p) {
    // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
    float sinThetaOp, cosThetaOp;
    if (p == 0) {
      sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
      cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
    }

    // Handle remainder of $p$ values for hair scale tilt
    else if (p == 1) {
      sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
      cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
    } else if (p == 2) {
      sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
      cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
    } else {
      sinThetaOp = sinThetaO;
      cosThetaOp = cosThetaO;
    }
    float s      = SqrtPiOver8 * (0.265f * beta_n + 1.194f * beta_n * beta_n +
                                5.372f * pow(beta_n, 20));
    float gammaO = asin(h);
    // Handle out-of-range $\cos \thetao$ from scale adjustment
    cosThetaOp = std::abs(cosThetaOp);
    fsum += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) * ap[p] *
            Np(phi, p, s, gammaO, gammaT);
  }

  // Compute contribution of remaining terms after _pMax_
  fsum += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[0]) * ap[pMax] /
          (2.f * pi);

  if (AbsCosTheta(wi) > 0) fsum /= AbsCosTheta(wi);

  return fsum;
}

std::array<float, pMax + 1> ComputeApPdf(float cosThetaO, float roughness) {
  // Compute array of $A_p$ values for _cosThetaO_
  float sinThetaO = sqrt(1 - cosThetaO * cosThetaO);

  // Compute $\cos \thetat$ for refracted ray
  float sinThetaT = sinThetaO / eta;
  float cosThetaT = sqrt(1 - sinThetaT * sinThetaT);

  // Compute $\gammat$ for refracted ray
  float etap      = sqrt(eta * eta - sinThetaO * sinThetaO) / cosThetaO;
  float sinGammaT = h / etap;
  float cosGammaT = sqrt(1 - sinGammaT * sinGammaT);

  // Compute the transmittance _T_ of a single path through the cylinder
  vec3f                       T  = exp(-sigma * (2 * cosGammaT / cosThetaT));
  std::array<vec3f, pMax + 1> ap = Ap(cosThetaO, eta, h, T);

  // Compute $A_p$ PDF from individual $A_p$ terms
  std::array<float, pMax + 1> apPdf;
  // float sumY = accumulate(ap.begin(), ap.end(), float(0),
  //[](float s, const vec3f& ap) { return s + ap.y; }); // error2
  for (int i = 0; i <= pMax; ++i) apPdf[i] = ap[i].y;  //  / sumY;
  return apPdf;
}

static float SampleTrimmedLogistic(float u, float s, float a, float b) {
  float k = LogisticCDF(b, s) - LogisticCDF(a, s);
  float x = -s * std::log(1 / (u * k + LogisticCDF(a, s)) - 1);
  return clamp(x, a, b);
}

///////////////////////////////////////////////////////////////////////////////
/////////////////////////////  hair brdf 1 ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
static vec3f sample_brdfcos(const pathtrace_brdf& brdf, const vec3f& normal,
    const vec3f& outgoing, float rnl, const vec2f& rn, float* pdf) {
  float sinThetaO = outgoing.x;
  float cosThetaO = sqrt(1 - sinThetaO * sinThetaO);
  float phiO      = std::atan2(outgoing.z, outgoing.y);

  vec2f u[2] = {rn, rn};
  vec3f wi;

  float h = -1 + 2 * rnl;

  // Determine which term $p$ to sample for hair scattering
  std::array<float, pMax + 1> apPdf = ComputeApPdf(cosThetaO, brdf.roughness);
  int                         p;
  for (p = 0; p < pMax; ++p) {
    if (u[0][0] < apPdf[p]) break;
    u[0][0] -= apPdf[p];
  }
  sin2kAlpha[0] = sin(radians(2));
  cos2kAlpha[0] = sqrt(1 - sin2kAlpha[0] * sin2kAlpha[0]);
  for (int i = 1; i < 3; ++i) {
    sin2kAlpha[i] = 2 * cos2kAlpha[i - 1] * sin2kAlpha[i - 1];
    cos2kAlpha[i] = cos2kAlpha[i - 1] * cos2kAlpha[i - 1] -
                    sin2kAlpha[i - 1] * sin2kAlpha[i - 1];
  }
  // Rotate $\sin \thetao$ and $\cos \thetao$ to account for hair scale tilt
  float sinThetaOp, cosThetaOp;
  if (p == 0) {
    sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
    cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
  } else if (p == 1) {
    sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
    cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
  } else if (p == 2) {
    sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
    cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
  } else {
    sinThetaOp = sinThetaO;
    cosThetaOp = cosThetaO;
  }

  // Sample $M_p$ to compute $\thetai$
  v[0] = (0.726f * beta_m + 0.812f * beta_m * beta_m + 3.7f * pow(beta_m, 20)) *
         (0.726f * beta_m + 0.812f * beta_m * beta_m + 3.7f * pow(beta_m, 20));
  v[1]            = .25 * v[0];
  v[2]            = 4 * v[0];
  v[3]            = v[2];
  v[4]            = v[2];
  u[1][0]         = max(u[1][0], float(1e-5));
  float cosTheta  = 1 + v[p] * log(u[1][0] + (1 - u[1][0]) * exp(-2 / v[p]));
  float sinTheta  = sqrt(1 - cosTheta * cosTheta);
  float cosPhi    = cos(2 * pi * u[1][1]);
  float sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
  float cosThetaI = sqrt(1 - sinThetaI * sinThetaI);

  // Compute $\gammat$ for refracted ray
  float etap      = sqrt(eta * eta - sinThetaO * sinThetaO) / cosThetaO;
  float sinGammaT = h / etap;
  float gammaT    = asin(sinGammaT);

  float dphi;
  if (p < pMax)
    dphi = Phi(p, gammaO, gammaT) + SampleTrimmedLogistic(u[0][1], s, -pi, pi);
  else
    dphi = 2 * pi * u[0][1];

  // // Compute _wi_ from sampled hair scattering angles
  float phiI = phiO + dphi;
  wi         = {sinThetaI, cosThetaI * cos(phiI), cosThetaI * sin(phiI)};

  *pdf = 0;
  for (int p = 0; p < pMax; ++p) {
    // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
    float sinThetaOp, cosThetaOp;
    if (p == 0) {
      sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
      cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
    }

    // Handle remainder of $p$ values for hair scale tilt
    else if (p == 1) {
      sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
      cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
    } else if (p == 2) {
      sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
      cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
    } else {
      sinThetaOp = sinThetaO;
      cosThetaOp = cosThetaO;
    }

    // Handle out-of-range $\cos \thetao$ from scale adjustment
    cosThetaOp = abs(cosThetaOp);
    *pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) * apPdf[p] *
            Np(dphi, p, s, gammaO, gammaT);
  }
  *pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) *
          apPdf[pMax] * (1 / (2 * pi));

  return (wi);
}
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

static vec3f sample_delta(const pathtrace_brdf& brdf, const vec3f& normal,
    const vec3f& outgoing, float rnl) {
  if (brdf.roughness != 0) return zero3f;

  // keep a weight sum to pick a lobe
  auto cdf = 0.0f;
  cdf += brdf.diffuse_pdf;

  if (brdf.specular_pdf != 0 && brdf.refraction_pdf == 0) {
    cdf += brdf.specular_pdf;
    if (rnl < cdf) {
      return sample_delta_reflection(brdf.ior, normal, outgoing);
    }
  }

  if (brdf.metal_pdf != 0) {
    cdf += brdf.metal_pdf;
    if (rnl < cdf) {
      return sample_delta_reflection(brdf.meta, brdf.metak, normal, outgoing);
    }
  }

  if (brdf.transmission_pdf != 0) {
    cdf += brdf.transmission_pdf;
    if (rnl < cdf) {
      return sample_delta_transmission(brdf.ior, normal, outgoing);
    }
  }

  if (brdf.refraction_pdf != 0) {
    cdf += brdf.refraction_pdf;
    if (rnl < cdf) {
      return sample_delta_refraction(brdf.ior, normal, outgoing, rnl);
    }
  }

  return zero3f;
}

// Compute the weight for sampling the BRDF
/*
static float sample_brdfcos_pdf(const pathtrace_brdf& brdf, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (brdf.roughness == 0) return 0;

  auto pdf = 0.0f;

  if (brdf.diffuse_pdf != 0) {
    pdf += brdf.diffuse_pdf *
           sample_diffuse_reflection_pdf(normal, outgoing, incoming);
  }

  if (brdf.specular_pdf != 0 && brdf.refraction_pdf == 0) {
    pdf += brdf.specular_pdf * sample_microfacet_reflection_pdf(brdf.ior,
                                   brdf.roughness, normal, outgoing, incoming);
  }

  if (brdf.metal_pdf != 0) {
    pdf += brdf.metal_pdf * sample_microfacet_reflection_pdf(brdf.meta,
                                brdf.metak, brdf.roughness, normal, outgoing,
                                incoming);
  }

  if (brdf.transmission_pdf != 0) {
    pdf += brdf.transmission_pdf * sample_microfacet_transmission_pdf(brdf.ior,
                                       brdf.roughness, normal, outgoing,
                                       incoming);
  }

  if (brdf.refraction_pdf != 0) {
    pdf += brdf.refraction_pdf * sample_microfacet_refraction_pdf(brdf.ior,
                                     brdf.roughness, normal, outgoing,
                                     incoming);
  }

  return pdf;
}
*/

//////////////////////////////////////////////////////////////////////
/////////////////////////7777  hair brdf 2 ////////////////////////////
///////////////////////////////////////////////////////////////////7

static float sample_brdfcos_pdf(const pathtrace_brdf& brdf, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming, const vec2f& rn) {
  // Compute hair coordinate system terms related to _wo_
  if (!brdf.roughness) return 0;
  float sinThetaO = outgoing.x;
  float cosThetaO = sqrt(1 - Sqr(sinThetaO));
  float phiO      = atan2(outgoing.z, outgoing.y);

  // Compute hair coordinate system terms related to _wi_
  float sinThetaI = incoming.x;
  float cosThetaI = sqrt(1 - Sqr(sinThetaI));
  float phiI      = atan2(incoming.z, incoming.y);

  // Compute $\gammat$ for refracted ray
  float etap      = sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
  float sinGammaT = h / etap;
  float gammaT    = asin(sinGammaT);

  // Compute PDF for $A_p$ terms
  std::array<float, pMax + 1> apPdf = ComputeApPdf(cosThetaO, brdf.roughness);

  // Compute PDF sum for hair scattering events
  float phi = phiI - phiO;
  float pdf = 0;
  for (int p = 0; p < pMax; ++p) {
    // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
    float sinThetaOp, cosThetaOp;
    if (p == 0) {
      sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
      cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
    }

    // Handle remainder of $p$ values for hair scale tilt
    else if (p == 1) {
      sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
      cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
    } else if (p == 2) {
      sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
      cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
    } else {
      sinThetaOp = sinThetaO;
      cosThetaOp = cosThetaO;
    }

    // Handle out-of-range $\cos \thetao$ from scale adjustment
    cosThetaOp = std::abs(cosThetaOp);
    pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) * apPdf[p] *
           Np(phi, p, s, gammaO, gammaT);
  }
  pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) * apPdf[pMax] *
         (1 / (2 * pi));
  return pdf;
}

////////////////////////////////////////////////////////////////////////////////////77
///////////////////////////////////////////////////////////////////////////////////////7

static float sample_delta_pdf(const pathtrace_brdf& brdf, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (brdf.roughness != 0) return 0;

  auto pdf = 0.0f;
  if (brdf.specular_pdf != 0 && brdf.refraction_pdf == 0) {
    pdf += brdf.specular_pdf *
           sample_delta_reflection_pdf(brdf.ior, normal, outgoing, incoming);
  }
  if (brdf.metal_pdf != 0) {
    pdf += brdf.metal_pdf * sample_delta_reflection_pdf(brdf.meta, brdf.metak,
                                normal, outgoing, incoming);
  }
  if (brdf.transmission_pdf != 0) {
    pdf += brdf.transmission_pdf *
           sample_delta_transmission_pdf(brdf.ior, normal, outgoing, incoming);
  }
  if (brdf.refraction_pdf != 0) {
    pdf += brdf.refraction_pdf *
           sample_delta_refraction_pdf(brdf.ior, normal, outgoing, incoming);
  }
  return pdf;
}

// Sample lights wrt solid angle
static vec3f sample_lights(const pathtrace_scene* scene, const vec3f& position,
    float rl, float rel, const vec2f& ruv) {
  auto light_id = sample_uniform((int)scene->lights.size(), rl);
  auto light    = scene->lights[light_id];
  if (light->instance) {
    auto element   = sample_discrete_cdf(light->cdf, rel);
    auto uv        = sample_triangle(ruv);
    auto lposition = eval_position(light->instance, element, uv);
    return normalize(lposition - position);
  } else if (light->environment) {
    if (light->environment->emission_tex) {
      auto emission_tex = light->environment->emission_tex;
      auto idx          = sample_discrete_cdf(light->cdf, rel);
      auto size         = texture_size(emission_tex);
      auto uv           = vec2f{
          (idx % size.x + 0.5f) / size.x, (idx / size.x + 0.5f) / size.y};
      return transform_direction(light->environment->frame,
          {cos(uv.x * 2 * pif) * sin(uv.y * pif), cos(uv.y * pif),
              sin(uv.x * 2 * pif) * sin(uv.y * pif)});
    } else {
      return sample_sphere(ruv);
    }
  } else {
    return zero3f;
  }
}

// Sample lights pdf
static float sample_lights_pdf(const pathtrace_scene* scene,
    const vec3f& position, const vec3f& direction) {
  auto pdf = 0.0f;
  for (auto& light : scene->lights) {
    if (light->instance) {
      // check all intersection
      auto lpdf          = 0.0f;
      auto next_position = position;
      for (auto bounce = 0; bounce < 100; bounce++) {
        auto intersection = intersect_instance_bvh(
            light->instance, {next_position, direction});
        if (!intersection.hit) break;
        // accumulate pdf
        auto lposition = eval_position(
            light->instance, intersection.element, intersection.uv);
        auto lnormal = eval_element_normal(
            light->instance, intersection.element);
        // prob triangle * area triangle = area triangle mesh
        auto area = light->cdf.back();
        lpdf += distance_squared(lposition, position) /
                (abs(dot(lnormal, direction)) * area);
        // continue
        next_position = lposition + direction * 1e-3f;
      }
      pdf += lpdf;
    } else if (light->environment) {
      if (light->environment->emission_tex) {
        auto emission_tex = light->environment->emission_tex;
        auto size         = texture_size(emission_tex);
        auto wl           = transform_direction(
            inverse(light->environment->frame), direction);
        auto texcoord = vec2f{atan2(wl.z, wl.x) / (2 * pif),
            acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
        if (texcoord.x < 0) texcoord.x += 1;
        auto i    = clamp((int)(texcoord.x * size.x), 0, size.x - 1);
        auto j    = clamp((int)(texcoord.y * size.y), 0, size.y - 1);
        auto prob = sample_discrete_cdf_pdf(light->cdf, j * size.x + i) /
                    light->cdf.back();
        auto angle = (2 * pif / size.x) * (pif / size.y) *
                     sin(pif * (j + 0.5f) / size.y);
        pdf += prob / angle;
      } else {
        pdf += 1 / (4 * pif);
      }
    }
  }
  pdf *= sample_uniform_pdf((int)scene->lights.size());
  return pdf;
}

static vec3f eval_scattering(
    const pathtrace_vsdf& vsdf, const vec3f& outgoing, const vec3f& incoming) {
  // YOUR CODE GOES HERE -------------------------- TERZA RICHIESTA
  if (vsdf.density == zero3f) return zero3f;
  return vsdf.scatter * vsdf.density *
         eval_phasefunction(vsdf.anisotropy, outgoing, incoming);
}

static vec3f sample_scattering(const pathtrace_vsdf& vsdf,
    const vec3f& outgoing, float rnl, const vec2f& rn) {
  // YOUR CODE GOES HERE -------------------------- TERZA RICHIESTA
  if (vsdf.density == zero3f) return zero3f;
  return sample_phasefunction(vsdf.anisotropy, outgoing, rn);
}

static float sample_scattering_pdf(
    const pathtrace_vsdf& vsdf, const vec3f& outgoing, const vec3f& incoming) {
  // YOUR CODE GOES HERE -------------------------- TERZA RICHIESTA
  if (vsdf.density == zero3f) return 0;
  return sample_phasefunction_pdf(vsdf.anisotropy, outgoing, incoming);
}

static vec4f shade_path(const pathtrace_scene* scene, const ray3f& ray_,
    rng_state& rng, const pathtrace_params& params);

// Path tracing.
static vec4f shade_volpath(const pathtrace_scene* scene, const ray3f& ray_,
    rng_state& rng, const pathtrace_params& params) {
  // YOUR CODE GOES HERE -------------------------------- TERZA RICHIESTA
  // initialize
  auto radiance = zero3f;
  auto weight   = vec3f{1, 1, 1};
  auto ray      = ray_;
  auto vstack   = vector<pathtrace_vsdf>{};
  auto hit      = false;

  // trace  path
  for (auto bounce = 0; bounce < params.bounces; bounce++) {
    // intersect next point
    auto intersection = intersect_scene_bvh(scene, ray);
    if (!intersection.hit) {
      radiance += weight * eval_environment(scene, ray);
      break;
    }

    // handle transmission if inside a volume
    // SAMPLE TRASMITTANCE (1)
    auto in_volume = false;
    if (!vstack.empty()) {
      auto extinction = vstack.back().density;
      auto distance   = sample_transmittance(
          extinction, intersection.distance, rand1f(rng), rand1f(rng));
      weight *= eval_transmittance(extinction, distance) /
                sample_transmittance_pdf(
                    extinction, distance, intersection.distance);
      in_volume             = distance < intersection.distance;
      intersection.distance = distance;
    }
    //

    // if not in volume:  HANDLE SURFACE (2)

    // switch between surface and volume
    if (!in_volume) {
      // prepare shading point
      auto outgoing = -ray.d;
      auto object   = scene->instances[intersection.instance];
      auto element  = intersection.element;
      auto uv       = intersection.uv;
      auto position = eval_position(object, element, uv);
      // now eval shading normal is implemented with evalnormalmap!
      auto normal   = eval_shading_normal(object, element, uv, outgoing);
      auto emission = eval_emission(object, element, uv, normal, outgoing);
      auto brdf     = eval_brdf(object, element, uv, normal, outgoing);

      // handle opacity
      if (brdf.opacity < 1 && rand1f(rng) >= brdf.opacity) {
        ray = {position + ray.d * 1e-2f, ray.d};
        bounce -= 1;
        continue;
      }
      hit = true;

      // accumulate emission
      radiance += weight * eval_emission(emission, normal, outgoing);

      auto incoming = zero3f;  // next dr

      float pdf = 0;  // hairshading insertion

      if (!is_delta(brdf)) {
        if (rand1f(rng) < 0.5f) {
          incoming = sample_brdfcos(
              brdf, normal, outgoing, rand1f(rng), rand2f(rng), &pdf);
        } else {
          incoming = sample_lights(
              scene, position, rand1f(rng), rand1f(rng), rand2f(rng));
        }

        /*weight *= eval_brdfcos(brdf, normal, outgoing, incoming) /
        //(0.5f * sample_brdfcos_pdf(brdf, normal, outgoing, incoming) +
        // 0.5f * sample_lights_pdf(scene, position, incoming));*/
        // changed here for hairshading
        weight *= f(brdf, outgoing, incoming) /
                  (0.5f * (pdf) +
                      0.5f * sample_lights_pdf(scene, position, incoming));
      } else {
        incoming = sample_delta(brdf, normal, outgoing, rand1f(rng));
        weight *= eval_delta(brdf, normal, outgoing, incoming) /
                  sample_delta_pdf(brdf, normal, outgoing, incoming);
      }
      // setup next iteration
      ray = {position, incoming};

      // update volume stack
      if (has_volume(object) &&
          dot(normal, outgoing) * dot(normal, incoming) < 0) {
        if (vstack.empty()) {
          auto volpoint = eval_vsdf(object, element, uv);
          vstack.push_back(volpoint);
        } else {
          vstack.pop_back();
        }
      }

      // else:  HANDLE VOLUME (3)

    } else {
      // prepare shading point
      auto  outgoing = -ray.d;
      auto  position = ray.o + ray.d * intersection.distance;  // eval ray
      auto& vsdf     = vstack.back();

      // handle opacity
      hit = true;

      // accumulate emission //is correct??
      radiance += weight;  // * eval_emission(vsdf, position, outgoing);

      // next direction
      auto incoming = zero3f;
      if (rand1f(rng) < 0.5f) {
        incoming = sample_scattering(vsdf, outgoing, rand1f(rng), rand2f(rng));
      } else {
        incoming = sample_lights(
            scene, position, rand1f(rng), rand1f(rng), rand2f(rng));
      }
      weight *= eval_scattering(vsdf, outgoing, incoming) /
                (0.5f * sample_scattering_pdf(vsdf, outgoing, incoming) +
                    0.5f * sample_lights_pdf(scene, position, incoming));

      // setup next iteration
      ray = {position, incoming};
    }

    // check over weight
    if (weight == zero3f || !isfinite(weight)) break;

    // russian roulette
    if (bounce > 3) {
      auto rr_prob = min((float)0.99, max(weight));
      if (rand1f(rng) >= rr_prob) break;
      weight *= 1 / rr_prob;
    }
  }

  return {radiance.x, radiance.y, radiance.z, hit ? 1.0f : 0.0f};
  // return shade_path(scene, ray_, rng, params);
}

// Path tracing.
static vec4f shade_path(const pathtrace_scene* scene, const ray3f& ray_,
    rng_state& rng, const pathtrace_params& params) {
  // initialize
  auto radiance = zero3f;
  auto weight   = vec3f{1, 1, 1};
  auto ray      = ray_;
  auto hit      = false;

  // trace  path
  for (auto bounce = 0; bounce < params.bounces; bounce++) {
    // intersect next point
    auto intersection = intersect_scene_bvh(scene, ray);
    if (!intersection.hit) {
      radiance += weight * eval_environment(scene, ray);
      break;
    }

    // prepare shading point
    auto outgoing = -ray.d;
    auto instance = scene->instances[intersection.instance];
    auto element  = intersection.element;
    auto uv       = intersection.uv;
    auto position = eval_position(instance, element, uv);
    auto normal   = eval_shading_normal(instance, element, uv, outgoing);
    auto emission = eval_emission(instance, element, uv, normal, outgoing);
    auto brdf     = eval_brdf(instance, element, uv, normal, outgoing);

    // handle opacity
    if (brdf.opacity < 1 && rand1f(rng) >= brdf.opacity) {
      ray = {position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }
    hit = true;

    // accumulate emission
    radiance += weight * eval_emission(emission, normal, outgoing);

    // next direction
    auto incoming = zero3f;

    float pdf = 0;  // hair insertion

    if (!is_delta(brdf)) {
      if (rand1f(rng) < 0.5f) {
        incoming = sample_brdfcos(
            brdf, normal, outgoing, rand1f(rng), rand2f(rng), &pdf);
      } else {
        incoming = sample_lights(
            scene, position, rand1f(rng), rand1f(rng), rand2f(rng));
      }
      /*weight *= eval_brdfcos(brdf, normal, outgoing, incoming) /
                (0.5f * sample_brdfcos_pdf(brdf, normal, outgoing, incoming) +
                    0.5f * sample_lights_pdf(scene, position, incoming));*/
      // changed here for hairshading
      weight *=
          f(brdf, outgoing, incoming) /
          (0.5f * (pdf) + 0.5f * sample_lights_pdf(scene, position, incoming));
    } else {
      incoming = sample_delta(brdf, normal, outgoing, rand1f(rng));
      weight *= eval_delta(brdf, normal, outgoing, incoming) /
                sample_delta_pdf(brdf, normal, outgoing, incoming);
    }

    // setup next iteration
    ray = {position, incoming};

    // check weight
    if (weight == zero3f || !isfinite(weight)) break;

    // russian roulette
    if (bounce > 3) {
      auto rr_prob = min((float)0.99, max(weight));
      if (rand1f(rng) >= rr_prob) break;
      weight *= 1 / rr_prob;
    }
  }

  return {radiance.x, radiance.y, radiance.z, hit ? 1.0f : 0.0f};
}

// Recursive path tracing.
static vec4f shade_naive(const pathtrace_scene* scene, const ray3f& ray_,
    rng_state& rng, const pathtrace_params& params) {
  // initialize
  auto radiance = zero3f;
  auto weight   = vec3f{1, 1, 1};
  auto ray      = ray_;
  auto hit      = false;

  // trace  path
  for (auto bounce = 0; bounce < params.bounces; bounce++) {
    // intersect next point
    auto intersection = intersect_scene_bvh(scene, ray);
    if (!intersection.hit) {
      radiance += weight * eval_environment(scene, ray);
      break;
    }

    // prepare shading point
    auto outgoing = -ray.d;
    auto instance = scene->instances[intersection.instance];
    auto element  = intersection.element;
    auto uv       = intersection.uv;
    auto position = eval_position(instance, element, uv);
    auto normal   = eval_shading_normal(instance, element, uv, outgoing);
    auto emission = eval_emission(instance, element, uv, normal, outgoing);
    auto brdf     = eval_brdf(instance, element, uv, normal, outgoing);

    // handle opacity
    if (brdf.opacity < 1 && rand1f(rng) >= brdf.opacity) {
      ray = {position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }
    hit = true;

    // accumulate emission
    radiance += weight * eval_emission(emission, normal, outgoing);

    // next direction
    auto incoming = zero3f;

    if (!is_delta(brdf)) {
      incoming = sample_brdfcos(
          brdf, normal, outgoing, rand1f(rng), rand2f(rng), NULL);
      weight *= eval_brdfcos(brdf, normal, outgoing, incoming) /
                sample_brdfcos_pdf(brdf, normal, outgoing, incoming,
                    rand2f(rng));  // changed here for hairshading

    } else {
      incoming = sample_delta(brdf, normal, outgoing, rand1f(rng));
      weight *= eval_delta(brdf, normal, outgoing, incoming) /
                sample_delta_pdf(brdf, normal, outgoing, incoming);
    }

    // setup next iteration
    ray = {position, incoming};

    // check weight
    if (weight == zero3f || !isfinite(weight)) break;

    // russian roulette
    if (bounce > 3) {
      auto rr_prob = min((float)0.99, max(weight));
      if (rand1f(rng) >= rr_prob) break;
      weight *= 1 / rr_prob;
    }
  }

  return {radiance.x, radiance.y, radiance.z, hit ? 1.0f : 0.0f};
}

// Eyelight for quick previewing.
static vec4f shade_eyelight(const pathtrace_scene* scene, const ray3f& ray_,
    rng_state& rng, const pathtrace_params& params) {
  // initialize
  auto radiance = zero3f;
  auto weight   = vec3f{1, 1, 1};
  auto ray      = ray_;
  auto hit      = false;

  // trace  path
  for (auto bounce = 0; bounce < max(params.bounces, 4); bounce++) {
    // intersect next point
    auto intersection = intersect_scene_bvh(scene, ray);
    if (!intersection.hit) {
      radiance += weight * eval_environment(scene, ray);
      break;
    }

    // prepare shading point
    auto outgoing = -ray.d;
    auto instance = scene->instances[intersection.instance];
    auto element  = intersection.element;
    auto uv       = intersection.uv;
    auto position = eval_position(instance, element, uv);
    auto normal   = eval_shading_normal(instance, element, uv, outgoing);
    auto emission = eval_emission(instance, element, uv, normal, outgoing);
    auto brdf     = eval_brdf(instance, element, uv, normal, outgoing);

    // handle opacity
    if (brdf.opacity < 1 && rand1f(rng) >= brdf.opacity) {
      ray = {position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }
    hit = true;

    // accumulate emission
    radiance += weight * eval_emission(emission, normal, outgoing);

    // brdf * light
    auto incoming = outgoing;
    radiance += weight * pif * eval_brdfcos(brdf, normal, outgoing, incoming);

    // continue path
    if (!is_delta(brdf)) break;
    incoming = sample_delta(brdf, normal, outgoing, rand1f(rng));
    weight *= eval_delta(brdf, normal, outgoing, incoming) /
              sample_delta_pdf(brdf, normal, outgoing, incoming);
    if (weight == zero3f || !isfinite(weight)) break;

    // setup next iteration
    ray = {position, incoming};
  }

  return {radiance.x, radiance.y, radiance.z, hit ? 1.0f : 0.0f};
}

// Trace a single ray from the camera using the given algorithm.
using pathtrace_shader_func = vec4f (*)(const pathtrace_scene* scene,
    const ray3f& ray, rng_state& rng, const pathtrace_params& params);
static pathtrace_shader_func get_shader(const pathtrace_params& params) {
  switch (params.shader) {
    case pathtrace_shader_type::naive: return shade_naive;
    case pathtrace_shader_type::path: return shade_path;
    case pathtrace_shader_type::eyelight: return shade_eyelight;
    case pathtrace_shader_type::volpath: return shade_volpath;
    default: {
      throw std::runtime_error("sampler unknown");
      return nullptr;
    }
  }
}

// Trace a block of samples
void render_sample(pathtrace_state* state, const pathtrace_scene* scene,
    const pathtrace_camera* camera, const vec2i& ij,
    const pathtrace_params& params) {
  auto shader = get_shader(params);
  auto ray    = sample_camera(camera, ij, state->render.imsize(),
      rand2f(state->rngs[ij]), rand2f(state->rngs[ij]));
  auto shaded = shader(scene, ray, state->rngs[ij], params);
  if (!isfinite(xyz(shaded))) shaded = {shaded.x, shaded.y, shaded.z, 1};
  if (max(xyz(shaded)) > params.clamp) {
    auto scale = params.clamp / max(xyz(shaded));
    shaded = {shaded.x * scale, shaded.y * scale, shaded.z * scale, shaded.w};
  }
  state->accumulation[ij] += shaded;
  state->samples[ij] += 1;
  state->render[ij] = state->accumulation[ij] / state->samples[ij];
}

// Forward declaration
pathtrace_light* add_light(pathtrace_scene* scene);

// Init trace lights
void init_lights(pathtrace_scene* scene, const pathtrace_params& params,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{0, 1};
  if (progress_cb) progress_cb("build light", progress.x++, progress.y);

  for (auto light : scene->lights) delete light;
  scene->lights.clear();

  for (auto instance : scene->instances) {
    if (instance->material->emission == zero3f) continue;
    auto shape = instance->shape;
    if (shape->triangles.empty()) continue;
    if (progress_cb) progress_cb("build light", progress.x++, ++progress.y);
    auto light      = add_light(scene);
    light->instance = instance;
    light->cdf      = vector<float>(shape->triangles.size());
    for (auto idx = 0; idx < light->cdf.size(); idx++) {
      auto& t         = shape->triangles[idx];
      light->cdf[idx] = triangle_area(
          shape->positions[t.x], shape->positions[t.y], shape->positions[t.z]);
      if (idx) light->cdf[idx] += light->cdf[idx - 1];
    }
  }
  for (auto environment : scene->environments) {
    if (environment->emission == zero3f) continue;
    if (progress_cb) progress_cb("build light", progress.x++, ++progress.y);
    auto light         = add_light(scene);
    light->environment = environment;
    if (environment->emission_tex) {
      auto texture = environment->emission_tex;
      auto size    = texture_size(texture);
      light->cdf   = vector<float>(size.x * size.y);
      for (auto i = 0; i < light->cdf.size(); i++) {
        auto ij       = vec2i{i % size.x, i / size.x};
        auto th       = (ij.y + 0.5f) * pif / size.y;
        auto value    = lookup_texture(texture, ij);
        light->cdf[i] = max(value) * sin(th);
        if (i) light->cdf[i] += light->cdf[i - 1];
      }
    }
  }

  // handle progress
  if (progress_cb) progress_cb("build light", progress.x++, progress.y);
}

// Init a sequence of random number generators.
void init_state(pathtrace_state* state, const pathtrace_scene* scene,
    const pathtrace_camera* camera, const pathtrace_params& params) {
  auto image_size =
      (camera->film.x > camera->film.y)
          ? vec2i{params.resolution,
                (int)round(params.resolution * camera->film.y / camera->film.x)}
          : vec2i{
                (int)round(params.resolution * camera->film.x / camera->film.y),
                params.resolution};
  state->render.assign(image_size, zero4f);
  state->accumulation.assign(image_size, zero4f);
  state->samples.assign(image_size, 0);
  state->rngs.assign(image_size, {});
  auto init_rng = make_rng(1301081);
  for (auto& rng : state->rngs) {
    rng = make_rng(params.seed, rand1i(init_rng, 1 << 31) / 2 + 1);
  }
}

// Progressively compute an image by calling render_samples multiple times.
void render_samples(pathtrace_state* state, const pathtrace_scene* scene,
    const pathtrace_camera* camera, const pathtrace_params& params) {
  if (params.noparallel) {
    for (auto j = 0; j < state->render.imsize().y; j++) {
      for (auto i = 0; i < state->render.imsize().x; i++) {
        render_sample(state, scene, camera, {i, j}, params);
      }
    }
  } else {
    parallel_for(state->render.imsize().x, state->render.imsize().y,
        [state, scene, camera, &params](int i, int j) {
          render_sample(state, scene, camera, {i, j}, params);
        });
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SCENE TESSELATION
// -----------------------------------------------------------------------------
namespace yocto {

// perform one level of subdivision and modify
template <typename T>
static void tesselate_catmullclark(
    std::vector<vec4i>& quads, std::vector<T>& vert, bool lock_boundary) {
  // YOUR CODE GOES HERE ---------------------- Prima richiesta
  //<init edges>
  // construct edge map and get edges and boundary
  auto emap     = make_edge_map(quads);
  auto edges    = get_edges(emap);
  auto boundary = get_boundary(emap);
  // initialize number of elements
  auto nv = (int)vert.size();      // nverts
  auto ne = (int)edges.size();     // nedges
  auto nb = (int)boundary.size();  // nboundary
  auto nf = (int)quads.size();     // nquads

  // create vertices --> phase 1
  auto tverts = vector<T>();
  for (auto v : vert)  // vertices [0, nv)
    tverts.push_back(v);
  for (auto e : edges)  // edge vertices [nv, nv+ne)
    tverts.push_back((vert[e.x] + vert[e.y]) / 2);
  for (auto q : quads)  // face vertices [nv+ne, nv+ne+nf)
    if (q.z != q.w)
      // quads
      tverts.push_back((vert[q.x] + vert[q.y] + vert[q.z] + vert[q.w]) / 4);
    else
      // triangles
      tverts.push_back((vert[q.x] + vert[q.y] + vert[q.z]) / 3);

  //<create faces>
  // Create faces by splitting each face into smaller ones
  auto tquads = std::vector<vec4i>();
  for (auto i = 0; i < quads.size(); i++) {
    auto q = quads[i];
    if (q.z != q.w) {  // If quad
      auto first_quad  = vec4i{q.x, nv + edge_index(emap, {q.x, q.y}),
          nv + ne + i, nv + edge_index(emap, {q.w, q.x})};  // Bottom left quad
      auto second_quad = vec4i{q.y, nv + edge_index(emap, {q.y, q.z}),
          nv + ne + i, nv + edge_index(emap, {q.x, q.y})};
      auto third_quad  = vec4i{q.z, nv + edge_index(emap, {q.z, q.w}),
          nv + ne + i, nv + edge_index(emap, {q.y, q.z})};
      auto fourth_quad = vec4i{q.w, nv + edge_index(emap, {q.w, q.x}),
          nv + ne + i, nv + edge_index(emap, {q.z, q.w})};
      tquads.push_back(first_quad);
      tquads.push_back(second_quad);
      tquads.push_back(third_quad);
      tquads.push_back(fourth_quad);
    } else {  // If triangle
      auto first_quad  = vec4i{q.x, nv + edge_index(emap, {q.x, q.y}),
          nv + ne + i, nv + edge_index(emap, {q.z, q.x})};  // Bottom left quad
      auto second_quad = vec4i{q.y, nv + edge_index(emap, {q.y, q.z}),
          nv + ne + i, nv + edge_index(emap, {q.x, q.y})};
      auto third_quad  = vec4i{q.z, nv + edge_index(emap, {q.z, q.x}),
          nv + ne + i, nv + edge_index(emap, {q.y, q.z})};
      tquads.push_back(first_quad);
      tquads.push_back(second_quad);
      tquads.push_back(third_quad);
    }
  }
  //<setup boundary>
  // Setup boundary arrays for collecting indices of boundary sides
  auto tboundary = vector<vec2i>();
  for (auto e : boundary) {
    // auto e = boundary[i];
    tboundary.push_back({e.x, nv + edge_index(emap, e)});
    tboundary.push_back({nv + edge_index(emap, e), e.y});
  }
  auto tcrease_edges = std::vector<vec2i>();
  auto tcrease_verts = std::vector<int>();  // List of vertices that need to
                                            // stay fixed (not movable)
  if (lock_boundary) {
    for (auto& b : tboundary) {
      tcrease_verts.push_back(b.x);
      tcrease_verts.push_back(b.y);
    }
  } else {
    for (auto& b : tboundary) tcrease_edges.push_back(b);
  }

  // Store an explicit distinction of vertices that can be moved and those not
  // 0 = Vertex is a crease vertex ==> Not movable
  // 1 = Vertex belonging to a boundary ==>
  // 2 = Vertex belonging to a face ==> Movable define vertex valence
  auto tverts_val = vector<int>(tverts.size(), 2);
  for (auto& e : tboundary) {
    tverts_val[e.x] = (lock_boundary) ? 0 : 1;
    tverts_val[e.y] = (lock_boundary) ? 0 : 1;
  }

  // averaging --> phase 2
  auto avert  = std::vector<T>(tverts.size(), T());
  auto acount = std::vector<int>(tverts.size(), 0);
  for (auto vertex_index : tcrease_verts) {
    if (tverts_val[vertex_index] != 0) continue;
    avert[vertex_index] += tverts[vertex_index];
    acount[vertex_index] += 1;
  }
  for (auto& edge : tcrease_edges) {
    auto centroid = (tverts[edge.x] + tverts[edge.y]) / 2;
    // The centroid correspond to the midpoint for edges
    for (auto edge_index : {edge.x, edge.y}) {
      if (tverts_val[edge_index] != 1) continue;
      // Update vertex position
      avert[edge_index] += centroid;
      acount[edge_index] += 1;
    }
  }
  for (auto& quad : tquads) {
    auto centroid =
        (tverts[quad.x] + tverts[quad.y] + tverts[quad.z] + tverts[quad.w]) / 4;
    for (auto vertex_index : {quad.x, quad.y, quad.z,
             quad.w}) {  // Loop over the vertices of the quad
      if (tverts_val[vertex_index] != 2)
        continue;  // Don't apply correction if the vertex is an unmovable one
      // Update vertex position
      avert[vertex_index] += centroid;
      acount[vertex_index] += 1;
    }
  }
  for (auto i = 0; i < tverts.size(); i++) {
    // Smooth vertices based on their neighbors (number of adjacent vertices)
    avert[i] /= (float)acount[i];
  }

  for (auto i = 0; i < tverts.size(); i++) {
    if (tverts_val[i] != 2) continue;
    // Correct the vertex position if the numb of adjacent faces is less than 4
    avert[i] = tverts[i] + (avert[i] - tverts[i]) * (4 / (float)acount[i]);
  }
  tverts = avert;

  swap(tquads, quads);
  swap(tverts, vert);
}
// end

void tesselate_shape(pathtrace_shape* shape) {
  if (shape->subdivisions != 0) {
    if (!shape->quadspos.empty()) {
      for (auto level = 0; level < shape->subdivisions; level++)
        tesselate_catmullclark(shape->quadspos, shape->positions, false);
      for (auto level = 0; level < shape->subdivisions; level++)
        tesselate_catmullclark(shape->quadstexcoord, shape->texcoords, true);
      if (shape->smooth) {
        shape->normals   = compute_normals(shape->quadspos, shape->positions);
        shape->quadsnorm = shape->quadspos;
      } else {
        shape->normals   = {};
        shape->quadsnorm = {};
      }
    } else {
      throw std::runtime_error("not supported yet");
    }
    shape->subdivisions = 0;
  }

  if (!shape->quadspos.empty()) {
    auto [quads, positions, normals, texcoords] = split_facevarying(
        shape->quadspos, shape->quadsnorm, shape->quadstexcoord,
        shape->positions, shape->normals, shape->texcoords);
    shape->positions     = positions;
    shape->normals       = normals;
    shape->texcoords     = texcoords;
    shape->triangles     = quads_to_triangles(quads);
    shape->quadspos      = {};
    shape->quadsnorm     = {};
    shape->quadstexcoord = {};
    shape->points        = {};
    shape->lines         = {};
    shape->radius        = {};
  }

  if (shape->displacement != 0 && shape->displacement_tex != nullptr &&
      !shape->triangles.empty()) {
    if (shape->normals.empty())
      shape->normals = compute_normals(shape->triangles, shape->positions);
    for (auto idx = 0; idx < shape->positions.size(); idx++) {
      auto disp = mean(xyz(
          eval_texture(shape->displacement_tex, shape->texcoords[idx], true)));
      if (!shape->displacement_tex->ldr.empty()) disp -= 0.5f;
      shape->positions[idx] += shape->normals[idx] * shape->displacement * disp;
    }
    if (shape->smooth) {
      shape->normals = compute_normals(shape->triangles, shape->positions);
    } else {
      shape->normals = {};
    }

    shape->displacement     = 0;
    shape->displacement_tex = nullptr;
  }
}

void tesselate_shapes(pathtrace_scene* scene, const pathtrace_params& params,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{0, (int)scene->shapes.size()};

  // tesselate shapes
  for (auto shape : scene->shapes) {
    if (progress_cb) progress_cb("tesselate shape", progress.x++, progress.y);
    tesselate_shape(shape);
  }

  // done
  if (progress_cb) progress_cb("tesselate shape", progress.x++, progress.y);
}

// Initialize subdivision surfaces
void init_subdivs(pathtrace_scene* scene, const pathtrace_params& params,
    progress_callback progress_cb) {
  tesselate_shapes(scene, params, progress_cb);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// cleanup
pathtrace_shape::~pathtrace_shape() {
  if (bvh) delete bvh;
}

// cleanup
pathtrace_scene::~pathtrace_scene() {
  if (bvh) delete bvh;
  for (auto camera : cameras) delete camera;
  for (auto instance : instances) delete instance;
  for (auto shape : shapes) delete shape;
  for (auto material : materials) delete material;
  for (auto texture : textures) delete texture;
  for (auto environment : environments) delete environment;
  for (auto light : lights) delete light;
}

// Add element
pathtrace_camera* add_camera(pathtrace_scene* scene) {
  return scene->cameras.emplace_back(new pathtrace_camera{});
}
pathtrace_texture* add_texture(pathtrace_scene* scene) {
  return scene->textures.emplace_back(new pathtrace_texture{});
}
pathtrace_shape* add_shape(pathtrace_scene* scene) {
  return scene->shapes.emplace_back(new pathtrace_shape{});
}
pathtrace_material* add_material(pathtrace_scene* scene) {
  return scene->materials.emplace_back(new pathtrace_material{});
}
pathtrace_instance* add_instance(pathtrace_scene* scene) {
  return scene->instances.emplace_back(new pathtrace_instance{});
}
pathtrace_environment* add_environment(pathtrace_scene* scene) {
  return scene->environments.emplace_back(new pathtrace_environment{});
}
pathtrace_light* add_light(pathtrace_scene* scene) {
  return scene->lights.emplace_back(new pathtrace_light{});
}

// Set cameras
void set_frame(pathtrace_camera* camera, const frame3f& frame) {
  camera->frame = frame;
}
void set_lens(pathtrace_camera* camera, float lens, float aspect, float film) {
  camera->lens = lens;
  camera->film = aspect >= 1 ? vec2f{film, film / aspect}
                             : vec2f{film * aspect, film};
}
void set_focus(pathtrace_camera* camera, float aperture, float focus) {
  camera->aperture = aperture;
  camera->focus    = focus;
}

// Add texture
void set_texture(pathtrace_texture* texture, const image<vec4b>& img) {
  texture->ldr = img;
  texture->hdr = {};
}
void set_texture(pathtrace_texture* texture, const image<vec4f>& img) {
  texture->ldr = {};
  texture->hdr = img;
}

// Add shape
void set_points(pathtrace_shape* shape, const vector<int>& points) {
  shape->points = points;
}
void set_lines(pathtrace_shape* shape, const vector<vec2i>& lines) {
  shape->lines = lines;
}
void set_triangles(pathtrace_shape* shape, const vector<vec3i>& triangles) {
  shape->triangles = triangles;
}
void set_fvquads(pathtrace_shape* shape, const vector<vec4i>& quadspos,
    const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord) {
  shape->quadspos      = quadspos;
  shape->quadsnorm     = quadsnorm;
  shape->quadstexcoord = quadstexcoord;
}
void set_positions(pathtrace_shape* shape, const vector<vec3f>& positions) {
  shape->positions = positions;
}
void set_normals(pathtrace_shape* shape, const vector<vec3f>& normals) {
  shape->normals = normals;
}
void set_texcoords(pathtrace_shape* shape, const vector<vec2f>& texcoords) {
  shape->texcoords = texcoords;
}
void set_radius(pathtrace_shape* shape, const vector<float>& radius) {
  shape->radius = radius;
}
void set_subdivision(
    pathtrace_shape* shape, int subdivisions, bool catmullclark, bool smooth) {
  shape->subdivisions = subdivisions;
  shape->catmullclark = catmullclark;
  shape->smooth       = smooth;
}
void set_displacement(pathtrace_shape* shape, float displacement,
    pathtrace_texture* displacement_tex) {
  shape->displacement     = displacement;
  shape->displacement_tex = displacement_tex;
}
// Add instance
void set_frame(pathtrace_instance* instance, const frame3f& frame) {
  instance->frame = frame;
}
void set_shape(pathtrace_instance* instance, pathtrace_shape* shape) {
  instance->shape = shape;
}
void set_material(pathtrace_instance* instance, pathtrace_material* material) {
  instance->material = material;
}

// Add material
void set_emission(pathtrace_material* material, const vec3f& emission,
    pathtrace_texture* emission_tex) {
  material->emission     = emission;
  material->emission_tex = emission_tex;
}
void set_color(pathtrace_material* material, const vec3f& color,
    pathtrace_texture* color_tex) {
  material->color     = color;
  material->color_tex = color_tex;
}
void set_specular(pathtrace_material* material, float specular,
    pathtrace_texture* specular_tex) {
  material->specular     = specular;
  material->specular_tex = specular_tex;
}
void set_metallic(pathtrace_material* material, float metallic,
    pathtrace_texture* metallic_tex) {
  material->metallic     = metallic;
  material->metallic_tex = metallic_tex;
}
void set_ior(pathtrace_material* material, float ior) { material->ior = ior; }
void set_transmission(pathtrace_material* material, float transmission,
    bool thin, float trdepth, pathtrace_texture* transmission_tex) {
  material->transmission     = transmission;
  material->thin             = thin;
  material->trdepth          = trdepth;
  material->transmission_tex = transmission_tex;
}
void set_thin(pathtrace_material* material, bool thin) {
  material->thin = thin;
}
void set_roughness(pathtrace_material* material, float roughness,
    pathtrace_texture* roughness_tex) {
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
}
void set_opacity(pathtrace_material* material, float opacity,
    pathtrace_texture* opacity_tex) {
  material->opacity     = opacity;
  material->opacity_tex = opacity_tex;
}
void set_scattering(pathtrace_material* material, const vec3f& scattering,
    float scanisotropy, pathtrace_texture* scattering_tex) {
  material->scattering     = scattering;
  material->scanisotropy   = scanisotropy;
  material->scattering_tex = scattering_tex;
}
void set_normalmap(
    pathtrace_material* material, pathtrace_texture* normal_tex) {
  material->normal_tex = normal_tex;
}

// Add environment
void set_frame(pathtrace_environment* environment, const frame3f& frame) {
  environment->frame = frame;
}
void set_emission(pathtrace_environment* environment, const vec3f& emission,
    pathtrace_texture* emission_tex) {
  environment->emission     = emission;
  environment->emission_tex = emission_tex;
}

}  // namespace yocto
