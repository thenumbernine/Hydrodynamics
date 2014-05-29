#pragma once

struct IHydro {
	virtual void update() = 0;
	virtual void draw() = 0;
	virtual void resize(int width, int height) = 0;
	virtual void pan(int dx, int dy) = 0;
	virtual void zoom(int dz) = 0;
};

