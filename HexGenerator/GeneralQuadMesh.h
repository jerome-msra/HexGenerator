/*! \file QuadMesh.h
*   \brief Quad Mesh
*   \author Jerome Jiang
*   \date   documented on 12/4/2015
*
*   Mesh for Quad
*/

#ifndef _GENERAL_QUAD_MESH_H_
#define _GENERAL_QUAD_MESH_H_

#include "Mesh\QuadBaseMesh.h"

namespace MeshLib
{
	template<typename V, typename E, typename F, typename H>
	class CGeneralQuadMesh : public CQuadMesh < V, E, F, H >
	{
	public:
		typedef V CVertex;
		typedef E CEdge;
		typedef F CFace;
		typedef H CHalfEdge;

		typedef CBoundary<V, E, F, H> CBoundary;
		typedef CLoop<V, E, F, H> CLoop;

		typedef MeshVertexIterator<V, E, F, H>			MeshVertexIterator;
		typedef MeshFaceIterator<V, E, F, H>			MeshFaceIterator;
		typedef MeshEdgeIterator<V, E, F, H>			MeshEdgeIterator;
		typedef VertexVertexIterator<V, E, F, H>		VertexVertexIterator;
		typedef FaceVertexIterator<V, E, F, H>			FaceVertexIterator;
		typedef VertexEdgeIterator<V, E, F, H>			VertexEdgeIterator;
		typedef VertexFaceIterator<V, E, F, H>			VertexFaceIterator;
		typedef VertexOutHalfedgeIterator<V, E, F, H>	VertexOutHalfedgeIterator;
		typedef VertexInHalfedgeIterator<V, E, F, H>	VertexInHalfedgeIterator;
		typedef FaceHalfedgeIterator<V, E, F, H>		FaceHalfedgeIterator;
	};

	typedef CGeneralQuadMesh<CQuadVertex, CQuadEdge, CQuadFace, CQuadHalfEdge> CGQMesh;
};

#endif