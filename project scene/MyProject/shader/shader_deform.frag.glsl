#version 330 core

in struct fragment_data
{
    vec3 position;
    vec3 normal;
    vec3 color;
    vec2 uv;
	vec3 eye;
	vec4 gl_Position;
	vec3 ToCamVec;
} fragment;

layout(location=0) out vec4 FragColor;

uniform sampler2D image_texture;  //reflection
uniform sampler2D image_texture_2; // refraction
uniform sampler2D image_texture_3; // dudv
uniform sampler2D image_texture_4; //normal map
uniform sampler2D image_texture_5; //depth map

uniform vec3 light = vec3(1.0, 1.0, 1.0);

uniform vec3 color = vec3(0.0, 0.0, 1.0); // Unifor color of the object

uniform float alpha = 1.0f; // alpha coefficient

uniform float Ka = 1.0; // Ambient coefficient
uniform float Kd = 0.0; // Diffuse coefficient
uniform float Ks = 1.0f;// Specular coefficient

uniform float specular_exp = 64.0; // Specular exponent
uniform bool use_texture = true;
uniform bool texture_inverse_y = false;

uniform float Movefactor;

vec3 normal;

const float Wavestrength=0.04;
void main()
{	
	// finding texture coordinates
	vec2 ndc= (fragment.gl_Position.xy/(fragment.gl_Position.w*2.0)) +0.5 ;
	vec2 refracTexCoords= vec2(ndc.x,ndc.y);
	vec2 reflecTexCoords= vec2(ndc.x,-ndc.y);

	// converting depth [0,1] into water depth_distance
	float near= 0.1f;
	float far= 1000.0f;
	float depth= texture(image_texture_5,refracTexCoords).r;
	float floor_dist=2.0 * near * far / (far + near - (2.0 * depth - 1.0) * (far - near));

	depth= gl_FragCoord.z;
	float water_dist=2.0 * near * far / (far + near - (2.0 * depth - 1.0) * (far - near));

	float water_depth= floor_dist-water_dist;

	//using image_texture_3: dudv map to calculate distortion on texture coordinates
	vec2 distortion1= Wavestrength*(texture(image_texture_3,vec2((fragment.uv.x+Movefactor), fragment.uv.y)).rg*2.0-1.0);
	vec2 distortion2= Wavestrength*(texture(image_texture_3,vec2((-fragment.uv.x+Movefactor), (fragment.uv.y+Movefactor))).rg*2.0-1.0);
	
	vec2 distortion= (distortion1+distortion2)*clamp(water_depth/20.0,0.0,1.0);
	
	reflecTexCoords+=distortion;
	reflecTexCoords.x=clamp(reflecTexCoords.x,0.001,0.999);
	reflecTexCoords.y=clamp(reflecTexCoords.y,-0.999,-0.001);
	refracTexCoords+=distortion;
	refracTexCoords=clamp(refracTexCoords,0.001,0.999);

	//calculating normals via corresponding normal map
	vec4 normalMapColour=texture(image_texture_4,distortion);
	normal= vec3(normalMapColour.r*2.0-1.0, normalMapColour.b,normalMapColour.g*2.0-1.0);

	//viewVec is the direction from fragment to camera, we are calculating the refractiv factor
	vec3 viewVec=normalize(fragment.ToCamVec);
	float refractivFactor= dot(viewVec,normal);
	refractivFactor=pow(refractivFactor,2);

	
	//calculating specular light using normal
	vec3 N = normalize(normal);
	if (gl_FrontFacing == false) {
		N = -N;
	}
	vec3 L = normalize(light-fragment.position);

	float diffuse = max(dot(N,L),0.0);
	float specular = 0.0;
	if(diffuse>0.0){
		vec3 R = reflect(-L,N);
		vec3 V = normalize(fragment.eye-fragment.position);
		specular = pow( max(dot(R,V),0.0), specular_exp )*clamp(water_depth/5.0,0.0,1.0);
	}

	//computing colors via textures and textures coordinates
	vec4 color_image_texture = texture(image_texture, reflecTexCoords);
	vec4 color_image_texture2 = texture(image_texture_2, refracTexCoords);


	// Finally mix the color of the two textures
	vec4 color_image_texture_ = mix(color_image_texture,color_image_texture2,1-refractivFactor);
	

	
	vec3 color_object  = (color*color_image_texture_.rgb)/color;
	vec3 color_shading = (Ka + Kd * diffuse) * color_object + Ks * specular * vec3(1.0, 1.0, 1.0); 
	
	


	
	FragColor = vec4(color_shading, alpha);// *clamp(water_depth/5.0,0.0,1.0)/alpha)*color_image_texture.a/alpha);
	
	
}
