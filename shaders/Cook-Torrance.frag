#version 330
#define PI 3.1415926538

uniform vec4 u_color;
uniform vec3 u_cam_pos; // view direction
uniform vec3 u_light_pos; // light direction
uniform vec3 u_light_intensity;

in vec4 v_position;
in vec4 v_normal; // surface normal
in vec4 v_tangent;

out vec4 out_color;


void main() {
 // ambient math
  vec3 ambient_light = vec3(1.0);
  float k_a = 0.3;

  // diffuse math
  float r = distance(u_light_pos, vec3(v_position));
  vec4 I_c = vec4(u_light_intensity, 0.0) / pow(r, 2);
  float s = 0.25; // constant ratio that balances out the diffuse and specular contributions to our reflectance
  float n_dot_l = dot(vec3(normalize(v_normal)), vec3(normalize(u_light_pos)));
  
  // D_blinn (specular math)
  vec3 h = normalize(u_cam_pos + u_light_pos);
  float alpha = 0.5; // also interpreted as roughness
  float first = 1 / (PI * pow(alpha, 2));
  float second = dot(h, vec3(normalize(v_normal)));
  float third = (2 / pow(alpha, 2)) - 2;
  float D_blinn = pow(first * second, third);

  // Geometric Attenuation (specular math)
  float h_dot_n = dot(h, vec3(normalize(v_normal)));
  float n_dot_v = dot(vec3(normalize(v_normal)), vec3(normalize(v_position)));
  float v_dot_h = dot(vec3(normalize(v_position)), h);
  float G = min(1, min(2 * h_dot_n * n_dot_v / v_dot_h, 2 * h_dot_n * n_dot_l / v_dot_h));

  // Fresnel (specular math)
  float n = 2; // index of refraction
  float F_o = pow(n - 1, 2) / pow(n + 1, 2);
  float F = F_o + (1 - F_o) * pow(1 - v_dot_h, 5);
  // specular math
  float r_s = D_blinn * G * F / (4 * n_dot_l * dot(vec3(normalize(v_normal)), vec3(normalize(v_normal))));
  
  
  // ambient lighting addition
  out_color += k_a * vec4(ambient_light, 0.0);

  // diffuse lighting addition
  out_color += I_c * max(0, n_dot_l) * s * u_color;

  // specular lighting addition
  out_color += I_c * max(0, n_dot_l) * (1-s) * r_s;


  out_color.a = 1;
}
