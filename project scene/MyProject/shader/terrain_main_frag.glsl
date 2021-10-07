#version 330 core

struct TerrainRegion
{
    float min;
    float max;
	float ymin;
	float ymax;
};

in struct fragment_data
{
    vec3 position;
    vec3 normal;
    vec3 color;
    vec2 uv;

	vec3 eye;
} fragment;

layout(location=0) out vec4 FragColor;

uniform TerrainRegion region1;
uniform TerrainRegion region2;
uniform TerrainRegion region3;


uniform sampler2D region1ColorMap;
uniform sampler2D region2ColorMap;
uniform sampler2D region3ColorMap;


uniform vec3 light = vec3(1.0, 1.0, 1.0);

uniform vec3 color = vec3(1.0, 1.0, 1.0); 
uniform float alpha = 1.0f;
uniform float Ka = 0.4;
uniform float Kd = 0.8; 
uniform float Ks = 0.4f;
uniform float specular_exp = 64.0;


vec4 GenerateTerrainColor()
{
    
	vec2 uv_image = vec2(fragment.uv.x, 1-fragment.uv.y);
	
	vec4 terrainColor = vec4(0.0, 0.0, 0.0, 1.0);
    float height = fragment.position.z;
    float width = fragment.position.y;
	float regionMin = 0.0;
    float regionMax = 0.0;
    float regionRange = 0.0;
    float regionWeight = 0.0;
    
	float regionWidthMin = 0.0;
    float regionWidthMax = 0.0;
    float regionWidthRange = 0.0;
	float regionWeightWidth = 0.0;
	
    // Terrain region 1.
    regionMin = region1.min;
    regionMax = region1.max;
    regionRange = regionMax - regionMin;
    regionWeight = (regionRange - abs(height - regionMax)) / regionRange;
    regionWeight = max(0.0, regionWeight);
	
	regionWidthMin = region1.ymin;
    regionWidthMax = region1.ymax;
    regionWidthRange = regionWidthMax - regionWidthMin;
    regionWeightWidth = (regionWidthRange - abs(width - regionWidthMax)) / regionWidthRange;
    regionWeightWidth = max(0.0, regionWeightWidth);
	
    terrainColor += max(0,regionWeightWidth)* texture(region1ColorMap, uv_image);
	
	
    // Terrain region 2.
    regionMin = region2.min;
    regionMax = region2.max;
    regionRange = regionMax - regionMin;
    regionWeight = (regionRange - abs(height - regionMax)) / regionRange;
    regionWeight = max(0.0, regionWeight);
	
	regionWidthMin = region2.ymin;
    regionWidthMax = region2.ymax;
    regionWidthRange = regionWidthMax - regionWidthMin;
	regionWidthRange = regionWidthMax - regionWidthMin;
    regionWeightWidth = (regionWidthRange - abs(width - regionWidthMax)) / regionWidthRange;
    regionWeightWidth = max(0.0, regionWeightWidth);
	
    terrainColor += max(0,regionWeightWidth) * texture(region2ColorMap, uv_image);

	  // Terrain region 3.
    regionMin = region3.min;
    regionMax = region3.max;
    regionRange = regionMax - regionMin;
    regionWeight = (regionRange - abs(height - regionMax)) / regionRange;
    regionWeight = max(0.0, regionWeight);
    
	regionWidthMin = region3.ymin;
    regionWidthMax = region3.ymax;
    regionWidthRange = regionWidthMax - regionWidthMin;
	regionWidthRange = regionWidthMax - regionWidthMin;
    regionWeightWidth = (regionWidthRange - abs(width - regionWidthMax)) / regionWidthRange;
    regionWeightWidth = max(0.0, regionWeightWidth);
	
	terrainColor += max(0,regionWeightWidth) * texture(region3ColorMap, uv_image);

    return terrainColor;
}

void main()
{   
    vec3 N = normalize(fragment.normal);
	if (gl_FrontFacing == false) {
		N = -N;
	}
	vec3 L = normalize(light-fragment.position);

	float diffuse = max(dot(N,L),0.0);
	float specular = 0.0;
	if(diffuse>0.0){
		vec3 R = reflect(-L,N);
		vec3 V = normalize(fragment.eye-fragment.position);
		specular = pow( max(dot(R,V),0.0), specular_exp );
	}

	vec4 color_image_texture = GenerateTerrainColor();
	
	
	
	
	
	float alpha_blend = fragment.uv.y;
	
	
	
	vec3 color_object  = fragment.color * color*color_image_texture.rgb;
	vec3 color_shading = (Ka + Kd * diffuse) * color_object + Ks * specular * vec3(1.0, 1.0, 1.0);
	
	FragColor = vec4(color_shading, alpha * color_image_texture.a);

}
