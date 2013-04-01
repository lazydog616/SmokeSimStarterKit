#pragma once
#ifndef PARTICLE_H
#define PARTICLE_H


class Particle {
public:
	float s; // -1 for negative, 1 for positive
	vec3 pos;
	float r;


	Particle(float px, float py, float pz) {
		pos = vec3(px, py, pz);
	};

	Particle(vec3 p) {
		pos = p;
	};
};

#endif