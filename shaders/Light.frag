#version 330

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;

out vec4 out_color;

void main() {
  float r = distance(u_light_pos, vec3(v_position));
  // Specular Reflection
  vec3 h = normalize(u_cam_pos + u_light_pos);
  float k_s = 0.5;
  int p = 100;
  out_color = k_s * vec4(u_light_intensity, 0.0) / pow(r, 2) * pow(max(0.0, dot(vec3(normalize(v_normal)), h)), p);
  // Ambient Lighting
  vec3 ambient_light = vec3(1.0);
  float k_a = 0.1;
  out_color += k_a * vec4(ambient_light, 0.0);
  // Diffuse Lighting
  out_color += u_color * vec4(u_light_intensity, 0.0) / pow(r, 2) * max(0.0, dot(vec3(normalize(v_normal)), normalize(u_light_pos)));
  out_color.a = 1;
}
