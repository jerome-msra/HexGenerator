#ifndef _TET_GENERAL_TMESH_H_
#define _TET_GENERAL_TMESH_H_

#include "TMeshLibHeaders.h"

namespace MeshLib
{
	namespace TMeshLib
	{
		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		class CTetGeneralTMesh : public CTMesh<TV, V, HE, TE, E, HF, F, T>
		{

			typedef TV	CTVertex;
			typedef V	CVertex;
			typedef HE	CHalfEdge;
			typedef TE	CTEdge;
			typedef E	CEdge;
			typedef HF	CHalfFace;
			typedef F	CFace;
			typedef T	CTet;

		public:
			CTetGeneralTMesh();
			~CTetGeneralTMesh();

			typedef TMeshTetIterator<TV, V, HE, TE, E, HF, F, T> MeshTetIterator;
			typedef TetHalfFaceIterator<TV, V, HE, TE, E, HF, F, T> TetHFIterator;
			typedef TMeshEdgeIterator<TV, V, HE, TE, E, HF, F, T> MeshEdgeIterator;
			typedef TMeshFaceIterator<TV, V, HE, TE, E, HF, F, T> MeshFaceIterator;
			typedef HalfFaceVertexIterator<TV, V, HE, TE, E, HF, F, T> HalfFaceVertexIterator;
			typedef FaceVertexIterator<TV, V, HE, TE, E, HF, F, T> FaceVertexIterator;
			typedef TMeshVertexIterator<TV, V, HE, TE, E, HF, F, T> MeshVertexIterator;
			typedef TVertexVertexIterator<TV, V, HE, TE, E, HF, F, T> VertexVertexIterator;

		private:

		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		CTetGeneralTMesh<TV, V, HE, TE, E, HF, F, T>::CTetGeneralTMesh()
		{
		}

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		CTetGeneralTMesh<TV, V, HE, TE, E, HF, F, T>::~CTetGeneralTMesh()
		{
		}

		typedef CTetGeneralTMesh<CTVertex, CVertex, CHalfEdge, CTEdge, CEdge, CHalfFace, CFace, CTet> CTGTMesh;
	};
};

#endif