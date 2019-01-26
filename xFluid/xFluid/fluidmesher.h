#pragma once

#include "array3d.h"
#include "material.h"
#include "trianglemesh.h"
#include "directeddistancefield.h"

class FluidMesher
{
public:
	FluidMesher();

	//void setParticleRadius(double r);


	/**
	Generate fluid surface mesh.
	INPUT:
	field, the directed distance field object. The points with the same distance value as threshold will be regarded as surface point.
	mesh, the result mesh object.
	*/
	void generateMesh(DirectedDistanceField& field, TriangleMesh& mesh);


private:
	/**
	For a certain surface cell, generate surface mesh. The METHOD we use is *Marching Cubes*.
	INPUT:
	g, grid cell index.
	mesh, the result mesh.
	edgeU, edgeV, edgeW, the structure that record which vertex has been push into mesh vertics. As each vertex may be shared by multi triangles.
	*/
	void _generateSurfaceCellMesh(/*DirectedDistanceField & field, */const GridIndex3d& g, TriangleMesh& mesh, Array3d<int>& edgeU, Array3d<int>& edgeV, Array3d<int>& edgeW);

	void _smoothMesh(TriangleMesh& mesh);
private:
	//double _particleRadius;

	DirectedDistanceField * _cur_field = nullptr;

	static const int _edgeTable[256];
	static const int _triTable[256][16];
	static const int _edge2vertice[12][2];								// table storing the relationship between edge index and the index of the two vertices of the edge.
	static const int _edge2dir[12];										// record the direction of each edge
};