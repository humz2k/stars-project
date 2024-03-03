#version 330

in vec3 fragPosition;
in vec2 fragTexCoord;
in vec4 fragColor;
in vec3 fragNormal;

uniform sampler2D texture0;
uniform sampler2D texture1;
uniform vec4 colDiffuse;
uniform vec3 camPos;

uniform float R;

out vec4 finalColor;

void main(){

    float dist = distance(fragPosition,vec3(0.0,0.0,0.0)) / R;

    finalColor = vec4(1.0,1.0-dist,0.0,1.0);
}