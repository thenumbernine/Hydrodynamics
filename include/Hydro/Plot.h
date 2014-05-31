#pragma once

#include "Quat.h"
#include "TensorMath/Vector.h"
#include <OpenGL/gl.h>
#include <algorithm>

template<int rank>
struct Plot {
	Quat viewAngle;
	float dist;

	Plot() : dist(2.) {}

	template<typename Hydro>
	void draw(Hydro &hydro) {
		typedef typename Hydro::CellGrid CellGrid;
		typedef typename Hydro::Cell Cell;
		typedef typename Hydro::StateVector StateVector;
		typedef typename Hydro::IVector IVector;
		typedef typename Hydro::Real Real;
		
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(0,0,-dist);
		Quat angleAxis = viewAngle.toAngleAxis();
		glRotatef(angleAxis(3) * 180. / M_PI, angleAxis(0), angleAxis(1), angleAxis(2));
		
		glEnable(GL_TEXTURE_1D);
		glDisable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_ONE, GL_ONE);
		glBegin(GL_POINTS);
		for (typename CellGrid::value_type &v : hydro.cells) {
			IVector index = v.first;
			Cell &cell = v.second;
			bool edge = false;
			for (int side = 0; side < rank; ++side) {
				if (index(side) < hydro.nghost-1 || index(side) >= hydro.size(side) - hydro.nghost) {
					edge = true;
					break;
				}
			}
			if (!edge) {
				//color by state value, neglect height or use it for coordinates
				glTexCoord1f(2. * cell.state(0));	//color by density
				glVertex3d(cell.x(0), cell.x(1), cell.x(2));
			}
		}
		glEnd();
		glDisable(GL_TEXTURE_1D);
	}

	static void resize(int width, int height) {
		const float zNear = .01;
		const float zFar = 10;
		float aspectRatio = (float)width / (float)height;
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glFrustum(-aspectRatio * zNear, aspectRatio * zNear, -zNear, zNear, zNear, zFar);
	}

	void pan(int dx, int dy) {
		float magn = sqrt(dx * dx + dy * dy);
		float fdx = (float)dx / magn;
		float fdy = (float)dy / magn;
		Quat rotation = Quat(fdy, fdx, 0, magn * M_PI / 180.).fromAngleAxis();
		viewAngle = rotation * viewAngle;
		viewAngle /= Quat::length(viewAngle);
	}

	void zoom(int dz) {
		dist *= (float)exp((float)dz * -.03f);
	}
};

//1D case
template<>
struct Plot<1> {
	::Vector<float, 2> viewPos;
	float viewZoom;

	Plot() : viewZoom(1.) {}

	template<typename Hydro>
	void draw(Hydro &hydro) {
		typedef typename Hydro::CellGrid CellGrid;
		typedef typename Hydro::Cell Cell;
		typedef typename Hydro::StateVector StateVector;
		typedef typename Hydro::IVector IVector;
		glPushMatrix();
		glTranslatef(-viewPos(0), -viewPos(1), 0);
		glScalef(viewZoom, viewZoom, viewZoom);
#if 0	//show piecewise step functions - in anticipation of getting PPM method working	
		for (int state = 0; state < 3; ++state) {
			::Vector<float,3> color;
			color(state) = 1;
			glColor3fv(color.v);
			for (int i = 0; i < hydro.size(0); ++i) {
				Cell &cell = hydro.cells(IVector(i));
				StateVector primitivesLeft = hydro.equationOfState->getPrimitives(cell.stateLeft(0));
				StateVector primitives = hydro.equationOfState->getPrimitives(cell.state);
				StateVector primitivesRight = hydro.equationOfState->getPrimitives(cell.stateRight(0));
				glBegin(GL_LINE_STRIP);
				glVertex2d(cell.interfaces(0).x(0), primitivesLeft(state));
				glVertex2d(cell.x(0), primitives(state));
				glVertex2d(hydro.cells(i+1).second.interfaces(0).x(0), primitivesRight(state));
				glEnd();
			}
		}
#endif
#if 1	//good ol fashioned graph
		for (int state = 0; state < 3; ++state) {
			::Vector<float,3> color;
			color(state) = 1;
			glColor3fv(color.v);
			glBegin(GL_LINE_STRIP);
			for (typename CellGrid::value_type &v : hydro.cells) {
				Cell &cell = v.second;
				StateVector primitives = hydro.equationOfState->getPrimitives(cell.state);
				glVertex2d(cell.x(0), primitives(state));
			}
			glEnd();
		}
#endif
		glPopMatrix();
	}

	static void resize(int width, int height) {
		float aspectRatio = (float)width / (float)height;
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(-aspectRatio, aspectRatio, -1., 3., -1., 1.);
		glMatrixMode(GL_MODELVIEW);
	}

	void pan(int dx, int dy) {
		viewPos(0) -= (float)dx * 0.01f;
		viewPos(1) += (float)dy * 0.01f;
	}
	
	void zoom(int dz) {
		float scale = (float)exp((float)dz * -.03f); 
		viewPos *= scale; 
		viewZoom *= scale;
	}
};

//2D case
template<>
struct Plot<2> {
	::Vector<float, 2> viewPos;
	float viewZoom;

	Plot() : viewZoom(1.) {}

	template<typename Hydro>
	void draw(Hydro &hydro) {
		typedef typename Hydro::CellGrid CellGrid;
		typedef typename Hydro::Cell Cell;
		typedef typename Hydro::IVector IVector;
		typedef typename Hydro::Real Real;
		glPushMatrix();
		glTranslatef(-viewPos(0), -viewPos(1), 0);
		glScalef(viewZoom, viewZoom, viewZoom);
		glEnable(GL_TEXTURE_1D);
		for (int y = hydro.nghost-1; y < hydro.size(1)-1; ++y) {
			glBegin(GL_TRIANGLE_STRIP);
			for (int x = hydro.nghost-1; x < hydro.size(0); ++x) {
			
				for (int offset = 0; offset < 2; ++offset) {
					int index = x + hydro.cells.size(0) * (y + offset);
					Cell &cell = hydro.cells.v[index].second;

					//color by state value, neglect height or use it for coordinates
					glTexCoord1f(2. * cell.state(0));
					glVertex2d(cell.x(0), cell.x(1));
				}
			}
			glEnd();
		}
		glDisable(GL_TEXTURE_1D);
		glPopMatrix();
	}

	static void resize(int width, int height) {
		float aspectRatio = (float)width / (float)height;
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(-aspectRatio, aspectRatio, -1., 1., -1., 1.);
		glMatrixMode(GL_MODELVIEW);
	}

	void pan(int dx, int dy) {
		viewPos(0) -= (float)dx * 0.01f;
		viewPos(1) += (float)dy * 0.01f;
	}
	
	void zoom(int dz) {
		float scale = (float)exp((float)dz * -.03f); 
		viewPos *= scale; 
		viewZoom *= scale;
	}
};



