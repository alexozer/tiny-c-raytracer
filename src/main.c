#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "export.h"

#define IMG_WIDTH 1000
#define IMG_HEIGHT 1000
#define EPSILON 0.00001

typedef struct {
    double x;
    double y;
    double z;
} vec3;

typedef struct {
    vec3 pt;
    vec3 dir;
} ray;

vec3 vec3_scale(vec3 vec, double scale) {
    return (vec3){vec.x * scale, vec.y * scale, vec.z * scale};
}

double vec3_dot(vec3 a, vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double vec3_norm_sq(vec3 vec) {
    return vec3_dot(vec, vec);
}

double vec3_norm(vec3 vec) {
    return sqrt(vec3_norm_sq(vec));
}

vec3 vec3_normalize(vec3 vec) {
    return vec3_scale(vec, 1 / vec3_norm(vec));
}

vec3 vec3_add(vec3 a, vec3 b) {
    return (vec3){a.x + b.x, a.y + b.y, a.z + b.z};
}

vec3 vec3_sub(vec3 a, vec3 b) {
    return vec3_add(a, vec3_scale(b, -1));
}

bool vec3_eq(vec3 a, vec3 b) {
    return vec3_norm_sq(vec3_sub(a, b)) < EPSILON;
}

const vec3 sphere_center = {0.2, 0.2, -3};
const double sphere_radius = 1;
const vec3 light_pos = {6, 6, 4};
const double light_intensity = 1;

bool intersect_ray_sphere(ray ray, vec3 sphere_center, double sphere_radius, vec3 *near_isec,
        vec3 *far_isec) {
    // Define formula vars
    double a = ray.dir.x;
    double b = ray.dir.y;
    double c = ray.dir.z;

    double d = ray.pt.x;
    double e = ray.pt.y;
    double f = ray.pt.z;

    double j = sphere_center.x;
    double k = sphere_center.y;
    double l = sphere_center.z;
    double r = sphere_radius;

    // Compute the coefficients for the quadratic
    double A = a * a + b * b + c * c;
    double B = 2 * (a * (d - j) + b * (e - k) + c * (f - l));
    double C = (d - j) * (d - j) + (e - k) * (e - k) + (f - l) * (f - l) - r * r;

    double det = B * B - 4 * A * C;
    if (det < 0) return false;

    double s1 = (-B + sqrt(det)) / (2 * A);
    double s2 = (-B - sqrt(det)) / (2 * A);

    double near_s, far_s;
    if (fabs(s1) < fabs(s2)) {
        near_s = s1;
        far_s = s2;
    } else {
        near_s = s2;
        far_s = s1;
    }

    if (near_isec) {
        near_isec->x = d + near_s * a;
        near_isec->y = e + near_s * b;
        near_isec->z = f + near_s * c;
    }
    if (far_isec) {
        far_isec->x = d + far_s * a;
        far_isec->y = e + far_s * b;
        far_isec->z = f + far_s * c;
    }

    return true;
}

int main() {
    unsigned num_pixels = IMG_HEIGHT * IMG_WIDTH;
    double *render_buffer = calloc(num_pixels, sizeof(double));
    uint8_t *render_byte_buffer = malloc(num_pixels * 4);

    for (int x = 0; x < IMG_WIDTH; x++) {
        for (int y = 0; y < IMG_HEIGHT; y++) {
            // Compute world coordinate of pixel in frustum
            // This is also the vec3 of the ray going out of the frustum
            double x_world = x / (double)IMG_WIDTH - 0.5;
            double y_world = (IMG_HEIGHT - y) / (double)IMG_HEIGHT - 0.5;
            double z_world = -1;

            // Compute the incoming ray vec3
            ray incoming_ray = {
                    {0, 0, 0},
                    {x_world, y_world, z_world},
            };

            // Compute the near-intersection with the sphere
            vec3 cam_sphere_intersect;
            if (!intersect_ray_sphere(incoming_ray, sphere_center, sphere_radius,
                    &cam_sphere_intersect, NULL)) continue;

            // Compute intersect-light ray
            vec3 light_sphere_ray = vec3_scale(
                    vec3_normalize(vec3_sub(cam_sphere_intersect, light_pos)),
                    light_intensity);

            // Check if it's occluded by the other side of the sphere
            ray intersect_light_ray = {light_pos, light_sphere_ray};
            vec3 near_light_sphere_intersect;
            if (!intersect_ray_sphere(intersect_light_ray, sphere_center, sphere_radius,
                    &near_light_sphere_intersect, NULL)) continue;
            if (!vec3_eq(near_light_sphere_intersect, cam_sphere_intersect)) continue;

            // Evaluate the surface scattering distribution function (matte)
            vec3 sphere_intersect_norm = vec3_normalize(vec3_sub(cam_sphere_intersect,
                    sphere_center));
            double refl_intensity = -vec3_dot(light_sphere_ray, sphere_intersect_norm);

            render_buffer[y * IMG_HEIGHT + x] = refl_intensity;
        }
    }

    double max_intensity = 0;
    for (int i = 0; i < num_pixels; i++) {
        if (render_buffer[i] > max_intensity) max_intensity = render_buffer[i];
    }

    if (max_intensity > EPSILON) {
        for (int pixel = 0; pixel < num_pixels; pixel++) {
            uint8_t intensity = (uint8_t)(render_buffer[pixel] * 255 / max_intensity);
            render_byte_buffer[pixel * 4 + 0] = intensity;
            render_byte_buffer[pixel * 4 + 1] = intensity;
            render_byte_buffer[pixel * 4 + 2] = intensity;
            render_byte_buffer[pixel * 4 + 3] = 255;
        }
    }

    write_png_file("test.png", IMG_WIDTH, IMG_HEIGHT, render_byte_buffer);

    free(render_buffer);
    free(render_byte_buffer);
}
