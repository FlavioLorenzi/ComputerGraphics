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

// import libs
#include "yocto_raytrace.h"

#include <yocto/yocto_color.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_parallel.h>
#include <yocto/yocto_shading.h>

// new libs
#include <yocto_gui/yocto_shade.h>

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SCENE EVALUATION
// -----------------------------------------------------------------------------
namespace yocto {

// Check texture size
static vec2i texture_size(const raytrace_texture* texture) {
  if (!texture->hdr.empty()) {
    return texture->hdr.imsize();
  } else if (!texture->ldr.empty()) {
    return texture->ldr.imsize();
  } else {
    return zero2i;
  }
}

// Evaluate a texture
static vec4f lookup_texture(const raytrace_texture* texture, const vec2i& ij,
    bool ldr_as_linear = false) {
  if (!texture->hdr.empty()) {
    return texture->hdr[ij];
  } else if (!texture->ldr.empty()) {
    return ldr_as_linear ? byte_to_float(texture->ldr[ij])
                         : srgb_to_rgb(byte_to_float(texture->ldr[ij]));
  } else {
    return {1, 1, 1, 1};
  }
}

// Evaluate a texture
static vec4f eval_texture(const raytrace_texture* texture, const vec2f& uv,
    bool ldr_as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false) {
  // get texture
  if (!texture) return {1, 1, 1};

  // get yimg::image width/height
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

  // get yimg::image coordinates and residuals
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

// Generates a ray from a camera for yimg::image plane coordinate uv and
// the lens coordinates luv
static ray3f eval_camera(const raytrace_camera* camera, const vec2f& image_uv) {
  // YOUR CODE GOES HERE ----------------------- code as in the slides
  // pinhole camera con raggio del tipo r(u,v) = { o, (q-o / |q-o|) }
  auto q = vec3f{camera->film.x * (0.5f - image_uv.x),
      camera->film.y * (image_uv.y - 0.5f), camera->lens};
  auto e = vec3f{0};
  // abbiamo -q perche cambia di segno
  auto d = normalize(-q - e);
  return ray3f{
      transform_point(camera->frame, e), transform_direction(camera->frame, d)};
}

// Eval position
static vec3f eval_position(
    const raytrace_shape* shape, int element, const vec2f& uv) {
  if (!shape->triangles.empty()) {
    auto t = shape->triangles[element];
    return interpolate_triangle(shape->positions[t.x], shape->positions[t.y],
        shape->positions[t.z], uv);
  } else if (!shape->lines.empty()) {
    auto l = shape->lines[element];
    return interpolate_line(shape->positions[l.x], shape->positions[l.y], uv.x);
  } else if (!shape->points.empty()) {
    return shape->positions[shape->points[element]];
  } else {
    return zero3f;
  }
}

// Shape element normal.
// se mancano le coordinate...
static vec3f eval_element_normal(const raytrace_shape* shape, int element) {
  // YOUR CODE GOES HERE -----------------------
  // se ho retta torno la tangente a questa
  // se ho triangolo torno la normale a questo
  if (!shape->triangles.empty()) {
    auto triangle = shape->triangles[element];
    return triangle_normal(shape->positions[triangle.x],
        shape->positions[triangle.y], shape->positions[triangle.z]);
  } else {
    auto line = shape->lines[element];
    return line_tangent(shape->positions[line.x], shape->positions[line.y]);
  }
}

// Eval normal
static vec3f eval_normal(
    const raytrace_shape* shape, int element, const vec2f& uv) {
  // YOUR CODE GOES HERE -----------------utile per terza e quinta richiesta
  if (shape->normals.empty()) return eval_element_normal(shape, element);

  if (!shape->triangles.empty()) {
    auto t = shape->triangles[element];
    return normalize(interpolate_triangle(
        shape->normals[t.x], shape->normals[t.y], shape->normals[t.z], uv));
  } else if (!shape->lines.empty()) {
    auto l = shape->lines[element];
    return normalize(
        interpolate_line(shape->normals[l.x], shape->normals[l.y], uv.x));
  } else {
    return {0, 0, 1};
  }
}

// namespace yocto

// Eval texcoord
static vec2f eval_texcoord(
    const raytrace_shape* shape, int element, const vec2f& uv) {
  // YOUR CODE GOES HERE ----------------------- identico a eval normal con txc
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

// Evaluate all environment color.
static vec3f eval_environment(const raytrace_scene* scene, const ray3f& ray) {
  // YOUR CODE GOES HERE ----------------------- come nelle slide:
  auto radiance = zero3f;
  for (auto environment : scene->environments) {
    auto envradiance = environment->emission;

    if (environment->emission_tex) {
      auto local_dir = transform_direction(inverse(environment->frame), ray.d);
      auto texcoord  = vec2f{atan2(local_dir.z, local_dir.x) / (2 * pif),
          acos(clamp(local_dir.y, -1.0f, 1.0f)) / pif};
      if (texcoord.x < 0) texcoord.x += 1;
      envradiance *= xyz(eval_texture(environment->emission_tex, texcoord));
    }
    radiance += envradiance;
  }
  return radiance;
}
}  // namespace yocto
// namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SHAPE/SCENE BVH
// -----------------------------------------------------------------------------
namespace yocto {

// primitive used to sort bvh entries
struct raytrace_bvh_primitive {
  bbox3f bbox      = invalidb3f;
  vec3f  center    = zero3f;
  int    primitive = 0;
};

// Splits a BVH node. Returns split position and axis.
static pair<int, int> split_middle(
    vector<raytrace_bvh_primitive>& primitives, int start, int end) {
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
static void build_bvh(vector<raytrace_bvh_node>& nodes,
    vector<raytrace_bvh_primitive>&              primitives) {
  // prepare to build nodes
  nodes.clear();
  nodes.reserve(primitives.size() * 2);

  // queue up first node
  auto queue = std::deque<vec3i>{{0, 0, (int)primitives.size()}};
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
      node.axis     = axis;
      node.num      = 2;
      node.start    = (int)nodes.size();
      nodes.emplace_back();
      nodes.emplace_back();
      queue.push_back({node.start + 0, start, mid});
      queue.push_back({node.start + 1, mid, end});
    } else {
      // Make a leaf node
      node.internal = false;
      node.num      = end - start;
      node.start    = start;
    }
  }

  // cleanup
  nodes.shrink_to_fit();
}

static void init_bvh(raytrace_shape* shape, const raytrace_params& params) {
  // build primitives
  auto primitives = vector<raytrace_bvh_primitive>{};
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
  shape->bvh = new raytrace_bvh_tree{};
  build_bvh(shape->bvh->nodes, primitives);

  // set bvh primitives
  shape->bvh->primitives.reserve(primitives.size());
  for (auto& primitive : primitives) {
    shape->bvh->primitives.push_back(primitive.primitive);
  }
}

void init_bvh(raytrace_scene* scene, const raytrace_params& params,
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
  auto primitives = vector<raytrace_bvh_primitive>{};
  auto object_id  = 0;
  for (auto instance : scene->instances) {
    auto& primitive = primitives.emplace_back();
    primitive.bbox  = instance->shape->bvh->nodes.empty()
                         ? invalidb3f
                         : transform_bbox(instance->frame,
                               instance->shape->bvh->nodes[0].bbox);
    primitive.center    = center(primitive.bbox);
    primitive.primitive = object_id++;
  }

  // build nodes
  if (scene->bvh) delete scene->bvh;
  scene->bvh = new raytrace_bvh_tree{};
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
static bool intersect_shape_bvh(raytrace_shape* shape, const ray3f& ray_,
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
static bool intersect_scene_bvh(const raytrace_scene* scene, const ray3f& ray_,
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
static bool intersect_instance_bvh(const raytrace_instance* instance,
    const ray3f& ray, int& element, vec2f& uv, float& distance, bool find_any,
    bool non_rigid_frames) {
  auto inv_ray = transform_ray(inverse(instance->frame, non_rigid_frames), ray);
  return intersect_shape_bvh(
      instance->shape, inv_ray, element, uv, distance, find_any);
}

raytrace_intersection intersect_scene_bvh(const raytrace_scene* scene,
    const ray3f& ray, bool find_any, bool non_rigid_frames) {
  auto intersection = raytrace_intersection{};
  intersection.hit  = intersect_scene_bvh(scene, ray, intersection.instance,
      intersection.element, intersection.uv, intersection.distance, find_any,
      non_rigid_frames);
  return intersection;
}
raytrace_intersection intersect_instance_bvh(const raytrace_instance* instance,
    const ray3f& ray, bool find_any, bool non_rigid_frames) {
  auto intersection = raytrace_intersection{};
  intersection.hit = intersect_instance_bvh(instance, ray, intersection.element,
      intersection.uv, intersection.distance, find_any, non_rigid_frames);
  return intersection;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace yocto {

///////////////////////////////////////////////////
////////////7 EXTRA POINT 2 ////////////////////////
///////////////////////////////////////////////////
////// my own shader

static vec4f shade_cartoon(const raytrace_scene* scene, const ray3f& ray,
    int bounce, rng_state& rng, const raytrace_params& params) {
  auto isec = intersect_scene_bvh(scene, ray);

  auto ee_x = eval_environment(scene, ray).x;
  auto ee_y = eval_environment(scene, ray).y;
  auto ee_z = eval_environment(scene, ray).z;
  if (!isec.hit) return {ee_x, ee_y, ee_z, 1};

  // can be useful to implemetnt
  auto instance = scene->instances[isec.instance];
  auto shape    = instance->shape;
  auto material = instance->material;

  auto position = transform_point(
      instance->frame, eval_position(shape, isec.element, isec.uv));

  auto normal = transform_direction(
      instance->frame, eval_normal(shape, isec.element, isec.uv));
  auto outgoing = -ray.d;
  ///

  if (!instance->shape->lines.empty()) {
    // con la linea ci serve la tangente!
    normal = orthonormalize(outgoing, normal);
  } else if (!instance->shape->triangles.empty()) {
    // se la normale guarda verso l esterno, cambia segno
    if (dot(outgoing, normal) < 0) normal = -normal;
  }
  // WE TRY TO APPROXIMATE COLORS

  // calcolo illuminazione dalla luce direzionale e del colore
  // auto incoming = sample_hemisphere(normal, rand2f(rng));
  auto incoming = vec3f{1, 2, 5};  // prendo un punto fisso

  auto WorldSpaceLP = vec3f{1, 10, 10};  // vettore che punta all'opposto

  auto texcoord = eval_texcoord(shape, isec.element, isec.uv);
  auto color    = material->color *
               xyz(eval_texture(material->color_tex, texcoord));

  auto NdotL = dot(WorldSpaceLP, normal);

  // SHADOWS ! ! !
  // se prima di colpire la sorgente interseca il pavimento
  // quel punto sara in ombra!
  auto  isecx = intersect_scene_bvh(scene, {position, incoming});
  float shadow;
  if (isecx.hit)
    shadow = 0;
  else
    shadow = 1;

  // NB: è una funzione simile alla sigmoide
  // con input : real number x, left edge, right edge
  // ritorna 0 se x è minoreuguale al left edge, 1 altrimenti
  // ritorna invece un'interpolazione smoothata tra 0 e 1
  // float lightIntensity = smoothstep(0, 0.01, NdotL);

  float lightIntensity = smoothstep(0, 0.01, NdotL * shadow);

  // Ambient Light
  auto AmbientColor = vec4f{0.4, 0.4, 0.4, 1};

  auto light = lightIntensity * LightColor;

  // Specular Reflection
  auto  SpecularColor = vec4f{0.9, 0.9, 0.9, 1};  // titns the reflection
  float Glossiness    = 32;                       // size of the reflection

  auto  halfVector        = normalize(WorldSpaceLP + outgoing);  // halfway
  float NdotH             = dot(normal, halfVector);
  float specularIntensity = pow(
      NdotH * lightIntensity, Glossiness * Glossiness);

  // once again use smoothstep to toonify the reflextion
  float specularIntensitySmooth = smoothstep(0.005, 0.01, specularIntensity);
  auto  specular                = specularIntensitySmooth * SpecularColor;

  // Rim light: simulate the reflection of the light
  auto rimDot = 1 - dot(outgoing, normal);

  auto  RimColor  = vec4f{1, 1, 1, 1};
  float RimAmount = 0.716;  // Tra 0 e 1

  float RimThreshold = 0.1;  // Tra 0 e 1

  float rimIntensity = rimDot * pow(NdotL, RimThreshold);
  rimIntensity = smoothstep(RimAmount - 0.01, RimAmount + 0.01, rimIntensity);
  auto rim     = rimIntensity * RimColor;

  auto radiance = color * xyz(AmbientColor + light + specular + rim);

  return vec4f{radiance.x, radiance.y, radiance.z, 1};
}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Raytrace renderer.
static vec4f shade_raytrace(const raytrace_scene* scene, const ray3f& ray,
    int bounce, rng_state& rng, const raytrace_params& params) {
  // YOUR CODE GOES HERE -----------------------
  // intersect next point: stesso setup dell'eyelight
  auto isec = intersect_scene_bvh(scene, ray);
  auto ee_x = eval_environment(scene, ray).x;
  auto ee_y = eval_environment(scene, ray).y;
  auto ee_z = eval_environment(scene, ray).z;
  if (!isec.hit) return {ee_x, ee_y, ee_z, 1};

  // shading point, eval geometry and normal corrextions
  auto instance = scene->instances[isec.instance];
  auto shape    = instance->shape;
  auto material = instance->material;
  auto position = transform_point(
      instance->frame, eval_position(shape, isec.element, isec.uv));
  auto normal = transform_direction(
      instance->frame, eval_normal(shape, isec.element, isec.uv));
  auto outgoing = -ray.d;

  // dot prod
  auto NdotO = dot(outgoing, normal);

  ////handling normals and lines
  if (!instance->shape->lines.empty()) {
    // con la linea ci serve la tangente!
    normal = orthonormalize(outgoing, normal);
  } else if (!instance->shape->triangles.empty()) {
    // se la normale guarda verso l esterno, cambia segno
    if (NdotO < 0) normal = -normal;
  }

  auto texcoord = eval_texcoord(shape, isec.element, isec.uv);

  // Materials and textures: each one is a f. scalar (tranne color e emission)

  auto color = material->color *
               xyz(eval_texture(material->color_tex, texcoord));
  auto emission = material->emission *
                  xyz(eval_texture(material->emission_tex, texcoord));
  auto specular = material->specular *
                  (eval_texture(material->specular_tex, texcoord)).x;
  auto metallic = material->metallic *
                  (eval_texture(material->metallic_tex, texcoord)).x;
  auto roughness = material->roughness *
                   (eval_texture(material->roughness_tex, texcoord)).x;

  // conv
  roughness = roughness * roughness;

  auto transmission = material->transmission *
                      (eval_texture(material->transmission_tex, texcoord)).x;

  // Handling opacity
  auto opacity = material->opacity;
  if (material->opacity_tex)
    // take the first coord as above
    opacity = opacity * (eval_texture(material->opacity_tex, texcoord)).x;

  // Check on opacity
  if (rand1f(rng) > opacity) {
    auto incoming = -outgoing;
    return shade_raytrace(scene, {position, incoming}, bounce + 1, rng, params);
  }

  // setup color and exit if enough bounces are done
  // accumulate emission
  auto radiance = emission;
  if (bounce >= params.bounces) {
    return {radiance.x, radiance.y, radiance.z, 1};
  }

  // Handling materials now!
  //// final loop: shader per una varietà di materiali

  if (transmission) {
    // Polished dieletric:
    // Scatter light both REFLECTING and TRANSMITTING
    /* Vogliamo solo un raggio di luce continuo e scegliamo casualmente la
     direzione in cui andare in base al termine fresnel ---> non aggiornare i
     pesi Fresnel. Il coefficiente di riflessione è piccolo pari a K_s = 0.04
     */ //NB: se il materiale è fino ed è un polished dieletric -> refraction

    // Evaluate Fresnel term
    auto fs = fresnel_schlick(vec3f{0.04, 0.04, 0.04}, normal, outgoing);

    if (material->thin) {  // check per l'EXTRA POINT 1 ////////////////////NB//

      // Random
      if (rand1f(rng) < mean(fs)) {
        auto incoming = reflect(outgoing, normal);
        radiance += xyz(shade_raytrace(
            scene, {position, incoming}, bounce + 1, rng, params));
      } else {
        auto incoming = -outgoing;
        radiance += color * xyz(shade_raytrace(scene, {position, incoming},
                                bounce + 1, rng, params));
      }
    } else {
      // take index of refraction and its inverse: Apply Refraction! !
      auto ior = reflectivity_to_eta(vec3f{0.04, 0.04, 0.04});
      auto eta = 1.f / mean(ior);

      if (rand1f(rng) < mean(fs)) {
        auto incoming = reflect(outgoing, normal);
        radiance += xyz(shade_raytrace(
            scene, {position, incoming}, bounce + 1, rng, params));
      } else {
        // inverti eta se il prodotto e' minore di zero
        auto incoming = NdotO < 0 ? refract(outgoing, normal, 1 / eta)
                                  : refract(outgoing, normal, eta);

        radiance += color * xyz(shade_raytrace(scene, {position, incoming},
                                bounce + 1, rng, params));
      }
    }

  } else if (metallic && !roughness) {
    // Polished metal:
    // Diffondono la luce come specchi (riflettono). Il colore
    // della superficie cambia ad incidenza normale seguendo le leggi di
    // Fresnel

    auto incoming = reflect(outgoing, normal);
    auto fs       = fresnel_schlick(color, normal, outgoing);
    radiance += fs * xyz(shade_raytrace(
                         scene, {position, incoming}, bounce + 1, rng, params));
  } else if (metallic && roughness) {
    // Rough metal:
    // Perché bisettrice (o halfway) e non normale? Questo perché la normale
    // effettiva del microfacet è la metà h tra l'ingresso e l'uscita, ovvero
    //  h = o+1/|o+i|

    // Incoming dir
    auto incoming = sample_hemisphere(normal, rand2f(rng));

    // Bisettrice
    auto halfway = normalize(outgoing + incoming);

    // Fresnel term
    auto fs = fresnel_schlick(color, halfway, outgoing);

    // Ogni micro sfaccettatura diffonde la luce come uno specchio metallico.
    // Pertanto, la distribuzione rappresenta il rapporto dei microfacets
    // orientati lungo la bisettrice h
    auto distribution = microfacet_distribution(roughness, normal, halfway);

    // Il termine shadowing cattura il rapporto tra microfacet visibile da un
    // particolare angolo in entrata e in uscita
    auto geometric = microfacet_shadowing(
        roughness, normal, halfway, outgoing, incoming);

    // den di ogni termine sumation
    //    4*|normal * outgoing|*|normal * incoming|
    auto denominator = 4 * abs(dot(normal, incoming)) *
                       abs(dot(normal, outgoing));
    // Summation term
    auto sum_term = fs * distribution * geometric / denominator;

    // Computing lighting (recursively)
    radiance += (2 * pif) * abs(dot(normal, incoming)) * sum_term *
                xyz(shade_raytrace(
                    scene, {position, incoming}, bounce + 1, rng, params));
  } else if (specular) {
    // Plastic:
    // La plastica è modellata come una superficie opaca rivestita con uno
    // strato dielettrico. Approssimato come somma di un contributo diffuso e
    // speculare:
    //  lo Specular layer riflette la luce (K_s = 0.04 is good)
    //  e la Diffuse Component è ponderata di uno meno il contributo di quella
    //  speculare per conservare l'energia

    // Sempre gli stessi steps!!!
    auto incoming = sample_hemisphere(normal, rand2f(rng));

    auto halfway = normalize(outgoing + incoming);

    auto fs = fresnel_schlick(vec3f{0.04}, halfway, outgoing);

    auto distribution = microfacet_distribution(roughness, normal, halfway);

    auto geometric = microfacet_shadowing(
        roughness, normal, halfway, outgoing, incoming);

    auto denominator = 4 * abs(dot(normal, incoming)) *
                       abs(dot(normal, outgoing));

    auto sum_term = fs * distribution * geometric / denominator;

    radiance += (2 * pif) * (color / pif * (1 - fs) + sum_term) *
                abs(dot(normal, incoming)) *
                xyz(shade_raytrace(
                    scene, {position, incoming}, bounce + 1, rng, params));

  } else {
    // Diffuse light
    auto incoming = sample_hemisphere(normal, rand2f(rng));
    // Errori di ombreggiatura: aggiungiamo un piccolo epsilon se il raggio
    // interseca la forma e il punto di ombreggiatura (già fatto prima)
    radiance += (2 * pif) * (color / pif) * abs(dot(normal, incoming)) *
                xyz(shade_raytrace(
                    scene, {position, incoming}, bounce + 1, rng, params));
  }

  // auto dp = color * abs(dot(normal, incoming));
  return vec4f{radiance.x, radiance.y, radiance.z, 1.0f};
}  // namespace yocto

//////////////////////////////////////////////
/////////////////////////////////////////////////

// Eyelight for quick previewing.
static vec4f shade_eyelight(const raytrace_scene* scene, const ray3f& ray,
    int bounce, rng_state& rng, const raytrace_params& params) {
  // YOUR CODE GOES HERE ----------------------- quinta richiesta
  // molto simile a shade normal,  solo con l'aggiunta di incoming
  // interseco punti successivi
  // intersect next point: stesso setup dell'eyelight
  auto isec = intersect_scene_bvh(scene, ray);
  if (!isec.hit) {
    auto ee = eval_environment(scene, ray);
    return vec4f{ee.x, ee.y, ee.z, 1};
  }

  // shading point
  auto instance = scene->instances[isec.instance];
  auto shape    = instance->shape;
  auto normal   = transform_direction(
      instance->frame, eval_normal(shape, isec.element, isec.uv));

  auto dir      = -ray.d;
  auto texcoord = eval_texcoord(shape, isec.element, isec.uv);
  auto position = transform_point(
      instance->frame, eval_position(shape, isec.element, isec.uv));

  // return the color and the absolute value between the dot prod.
  // of the normal and the outgoing direction
  auto color = instance->material->color;
  auto sc    = color * abs(dot(normal, dir));  // shortcut
  return {sc.x, sc.y, sc.z, 1};
}

static vec4f shade_normal(const raytrace_scene* scene, const ray3f& ray,
    int bounce, rng_state& rng, const raytrace_params& params) {
  // YOUR CODE GOES HERE ----------------------- terza richiesta
  // esattamente come per lo shade color, ma stavolta uso eval_normal
  // (implementato su)
  // interseco punti successivi
  auto isec = intersect_scene_bvh(scene, ray);
  if (!isec.hit) return zero4f;
  // shading point
  auto instance = scene->instances[isec.instance];
  auto shape    = instance->shape;
  auto normal   = transform_direction(
      instance->frame, eval_normal(shape, isec.element, isec.uv));

  return vec4f{normal.x * 0.5f + 0.5f, normal.y * 0.5f + 0.5f,
      normal.z * 0.5f + 0.5f, 1};
}

static vec4f shade_texcoord(const raytrace_scene* scene, const ray3f& ray,
    int bounce, rng_state& rng, const raytrace_params& params) {
  // YOUR CODE GOES HERE ----------------------- quarta richiesta
  // simile a shade normal, ma con instances e return diverso con fmod
  // intersezione con next point
  auto isec = intersect_scene_bvh(scene, ray);
  if (!isec.hit) return zero4f;
  // preparazione shading point
  auto instance = scene->instances[isec.instance];
  auto shape    = instance->shape;
  auto texcoord = eval_texcoord(shape, isec.element, isec.uv);
  // note: usare fmode() per forzarli in un range 0-1
  return {fmod(texcoord.x, 1), fmod(texcoord.y, 1), 0, 1};
}

static vec4f shade_color(const raytrace_scene* scene, const ray3f& ray,
    int bounce, rng_state& rng, const raytrace_params& params) {
  // YOUR CODE GOES HERE ----------------------- seconda richiesta
  // intersezione con prossimo punto
  auto isec = intersect_scene_bvh(scene, ray);
  if (!isec.hit) return zero4f;

  // preparo il punto di shading
  auto instance = scene->instances[isec.instance];
  return {instance->material->color.x, instance->material->color.y,
      instance->material->color.z, 1};  // return color
}

// Trace a single ray from the camera using the given algorithm.
using raytrace_shader_func = vec4f (*)(const raytrace_scene* scene,
    const ray3f& ray, int bounce, rng_state& rng,
    const raytrace_params& params);
static raytrace_shader_func get_shader(const raytrace_params& params) {
  switch (params.shader) {
    case raytrace_shader_type::raytrace: return shade_raytrace;
    case raytrace_shader_type::eyelight: return shade_eyelight;
    case raytrace_shader_type::normal: return shade_normal;
    case raytrace_shader_type::texcoord: return shade_texcoord;
    case raytrace_shader_type::color: return shade_color;
    case raytrace_shader_type::cartoon: return shade_cartoon;  // extra point 2
    default: {
      throw std::runtime_error("sampler unknown");
      return nullptr;
    }
  }
}

// MAIN LOOP

// Trace a block of samples
void render_sample(raytrace_state* state, const raytrace_scene* scene,
    const raytrace_camera* camera, const vec2i& ij,
    const raytrace_params& params) {
  // YOUR CODE GOES HERE ----------------------- prima richiesta
  return render_samples(state, scene, camera, params);
}

// Init a sequence of random number generators.
void init_state(raytrace_state* state, const raytrace_scene* scene,
    const raytrace_camera* camera, const raytrace_params& params) {
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

//
//
/////////////////////////////////////////////////////////////////////////////////////

// Progressively compute an image by calling trace_samples multiple times.
void render_samples(raytrace_state* state, const raytrace_scene* scene,
    const raytrace_camera* camera, const raytrace_params& params) {
  // YOUR CODE GOES HERE
  auto shader = get_shader(params);

  if (params.noparallel) {
    // YOUR CODE GOES HERE ----------------------- Prima richiesta main loop
    // loop over image pixels
    for (auto j = 0; j < state->render.imsize().y; j++) {
      for (auto i = 0; i < state->render.imsize().x; i++) {
        // grab the pixel
        auto puv = rand2f(state->rngs[i, j]);
        // get pixel uv from rng
        auto p  = vec2f{(float)i, (float)j} + puv;
        auto uv = vec2f{
            p.x / state->render.imsize().x, p.y / state->render.imsize().y};
        // get camera ray : 2 ways available
        // auto ray = camera_ray(camera->frame, camera->lens, camera->film,
        // uv);
        auto ray = eval_camera(camera, uv);
        // call shader img
        auto color = shader(scene, ray, 0, (state->rngs[i, j]), params);
        // clamp to max value
        if (length(xyz(color)) > params.clamp)
          color = normalize(color) * params.clamp;

        // update state accumulation, samples and render
        state->accumulation[{i, j}] += color;
        state->samples[{i, j}] += 1;
        state->render[{i, j}] = (state->accumulation[{i, j}]) /
                                (state->samples[{i, j}]);
      };
    }
  } else {
    // YOUR CODE GOES HERE -----------------------
    parallel_for(state->render.imsize().x, state->render.imsize().y,
        [state, scene, camera, shader, &params](int i, int j) {
          // Image dimension
          // auto sy = state->render.imsize().y;
          // auto sx = state->render.imsize().x;
          auto sy = state->render.imsize().y;
          auto sx = state->render.imsize().x;

          auto ij_v = vec2i{i, j};

          // we get uv pixel from rng
          auto puv = rand2f(state->rngs[ij_v]);
          auto d   = vec2f{(float)i, (float)j} + puv;

          auto uv = vec2f{d.x / sx, d.y / sy};

          // we get the camera ray
          auto ray = eval_camera(camera, uv);

          // We get the proper shader
          auto color = shader(scene, ray, 0, state->rngs[i, j], params);

          // Here we clamp to max value
          if (length(xyz(color)) > params.clamp) {
            color = normalize(color) * params.clamp;
          }

          // update state accumulation, samples and render

          state->accumulation[ij_v] += color;
          state->samples[ij_v] += 1;
          state->render[ij_v] = (state->accumulation[ij_v]) /
                                (state->samples[ij_v]);
        });
  }
}
///////////////////////////////////////////////////////////////////////////////////

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// cleanup
raytrace_shape::~raytrace_shape() {
  if (bvh) delete bvh;
}

// cleanup
raytrace_scene::~raytrace_scene() {
  if (bvh) delete bvh;
  for (auto camera : cameras) delete camera;
  for (auto instance : instances) delete instance;
  for (auto shape : shapes) delete shape;
  for (auto material : materials) delete material;
  for (auto texture : textures) delete texture;
  for (auto environment : environments) delete environment;
}

// Add element
raytrace_camera* add_camera(raytrace_scene* scene) {
  return scene->cameras.emplace_back(new raytrace_camera{});
}
raytrace_texture* add_texture(raytrace_scene* scene) {
  return scene->textures.emplace_back(new raytrace_texture{});
}
raytrace_shape* add_shape(raytrace_scene* scene) {
  return scene->shapes.emplace_back(new raytrace_shape{});
}
raytrace_material* add_material(raytrace_scene* scene) {
  return scene->materials.emplace_back(new raytrace_material{});
}
raytrace_instance* add_instance(raytrace_scene* scene) {
  return scene->instances.emplace_back(new raytrace_instance{});
}
raytrace_environment* add_environment(raytrace_scene* scene) {
  return scene->environments.emplace_back(new raytrace_environment{});
}

// Set cameras
void set_frame(raytrace_camera* camera, const frame3f& frame) {
  camera->frame = frame;
}
void set_lens(raytrace_camera* camera, float lens, float aspect, float film) {
  camera->lens = lens;
  camera->film = aspect >= 1 ? vec2f{film, film / aspect}
                             : vec2f{film * aspect, film};
}
void set_focus(raytrace_camera* camera, float aperture, float focus) {
  camera->aperture = aperture;
  camera->focus    = focus;
}

// Add texture
void set_texture(raytrace_texture* texture, const image<vec4b>& img) {
  texture->ldr = img;
  texture->hdr = {};
}
void set_texture(raytrace_texture* texture, const image<vec4f>& img) {
  texture->ldr = {};
  texture->hdr = img;
}

// Add shape
void set_points(raytrace_shape* shape, const vector<int>& points) {
  shape->points = points;
}
void set_lines(raytrace_shape* shape, const vector<vec2i>& lines) {
  shape->lines = lines;
}
void set_triangles(raytrace_shape* shape, const vector<vec3i>& triangles) {
  shape->triangles = triangles;
}
void set_positions(raytrace_shape* shape, const vector<vec3f>& positions) {
  shape->positions = positions;
}
void set_normals(raytrace_shape* shape, const vector<vec3f>& normals) {
  shape->normals = normals;
}
void set_texcoords(raytrace_shape* shape, const vector<vec2f>& texcoords) {
  shape->texcoords = texcoords;
}
void set_radius(raytrace_shape* shape, const vector<float>& radius) {
  shape->radius = radius;
}

// Add instance
void set_frame(raytrace_instance* instance, const frame3f& frame) {
  instance->frame = frame;
}
void set_shape(raytrace_instance* instance, raytrace_shape* shape) {
  instance->shape = shape;
}
void set_material(raytrace_instance* instance, raytrace_material* material) {
  instance->material = material;
}

// Add material
void set_emission(raytrace_material* material, const vec3f& emission,
    raytrace_texture* emission_tex) {
  material->emission     = emission;
  material->emission_tex = emission_tex;
}
void set_color(raytrace_material* material, const vec3f& color,
    raytrace_texture* color_tex) {
  material->color     = color;
  material->color_tex = color_tex;
}
void set_specular(raytrace_material* material, float specular,
    raytrace_texture* specular_tex) {
  material->specular     = specular;
  material->specular_tex = specular_tex;
}
void set_metallic(raytrace_material* material, float metallic,
    raytrace_texture* metallic_tex) {
  material->metallic     = metallic;
  material->metallic_tex = metallic_tex;
}
void set_ior(raytrace_material* material, float ior) { material->ior = ior; }
void set_transmission(raytrace_material* material, float transmission,
    bool thin, float trdepth, raytrace_texture* transmission_tex) {
  material->transmission     = transmission;
  material->thin             = thin;
  material->trdepth          = trdepth;
  material->transmission_tex = transmission_tex;
}
void set_thin(raytrace_material* material, bool thin) { material->thin = thin; }
void set_roughness(raytrace_material* material, float roughness,
    raytrace_texture* roughness_tex) {
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
}
void set_opacity(
    raytrace_material* material, float opacity, raytrace_texture* opacity_tex) {
  material->opacity     = opacity;
  material->opacity_tex = opacity_tex;
}
void set_scattering(raytrace_material* material, const vec3f& scattering,
    float scanisotropy, raytrace_texture* scattering_tex) {
  material->scattering     = scattering;
  material->scanisotropy   = scanisotropy;
  material->scattering_tex = scattering_tex;
}

// Add environment
void set_frame(raytrace_environment* environment, const frame3f& frame) {
  environment->frame = frame;
}
void set_emission(raytrace_environment* environment, const vec3f& emission,
    raytrace_texture* emission_tex) {
  environment->emission     = emission;
  environment->emission_tex = emission_tex;
}

}  // namespace yocto
