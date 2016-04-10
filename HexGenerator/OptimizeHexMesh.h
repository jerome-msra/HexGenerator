#ifndef _OPTIMIZE_HEX_MESH_H_
#define _OPTIMIZE_HEX_MESH_H_

#include "HMeshLibHeaders.h"

namespace MeshLib
{
	namespace HMeshLib
	{

		class CHexOptimizeVertex : public CVertex
		{
		public:
			bool & diskBoundary() { return m_diskBoundary; };
			bool & boundary() { return m_boundary; };
			bool & disk() { return m_disk; }
			//std::vector<CHexOptimizeVertex *> & sameVertex() { return m_sameVertex; };

			void _from_string();

		protected:
			bool m_diskBoundary = false;
			bool m_boundary = false;
			//std::vector<CHexOptimizeVertex *> m_sameVertex;
			bool m_disk = false;
		};

		inline void CHexOptimizeVertex::_from_string()
		{
			CParser parser(m_string);

			for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
			{
				CToken * token = *iter;
				
				if (token->m_key == "disk")
				{
					m_disk = true;
				}
			}
		};


		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		class OptimizeHexMesh : public CHMesh < HXV, V, HE, HXE, E, HF, F, HX >
		{
		public:
			typedef HXV CHVertex;
			typedef V	CVertex;
			typedef HE	CHalfEdge;
			typedef HXE CHEdge;
			typedef E	CEdge;
			typedef HF	CHalfFace;
			typedef HX	CHex;

			typedef HMeshHexIterator<HXV, V, HE, HXE, E, HF, F, HX> MeshTetIterator;
			typedef HexHalfFaceIterator<HXV, V, HE, HXE, E, HF, F, HX> TetHFIterator;
			typedef HMeshEdgeIterator<HXV, V, HE, HXE, E, HF, F, HX> MeshEdgeIterator;
			typedef HMeshFaceIterator<HXV, V, HE, HXE, E, HF, F, HX> MeshFaceIterator;
			typedef HalfFaceVertexIterator<HXV, V, HE, HXE, E, HF, F, HX> HalfFaceVertexIterator;
			typedef FaceVertexIterator<HXV, V, HE, HXE, E, HF, F, HX> FaceVertexIterator;
			typedef HMeshVertexIterator<HXV, V, HE, HXE, E, HF, F, HX> MeshVertexIterator;
			typedef HVertexVertexIterator<HXV, V, HE, HXE, E, HF, F, HX> VertexVertexIterator;

		public:
			OptimizeHexMesh() { ; };
			~OptimizeHexMesh() { ; };

			void labelBoundaryVertices();

		private:

		};

		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		void OptimizeHexMesh< HXV, V, HE, HXE, E, HF, F, HX >::labelBoundaryVertices()
		{
			for (std::list<CHalfFace*>::iterator hfIter = m_pHalfFaces.begin(); hfIter != m_pHalfFaces.end(); hfIter++)
			{
				CHalfFace * pHF = *hfIter;
				CHalfFace * pHFD = HalfFaceDual(pHF);
				if (pHFD == nullptr)
				{
					CHalfEdge * pHalfEdge = HalfFaceHalfEdge(pHF);
					do
					{
						CVertex * pV0 = HalfEdgeTarget(pHalfEdge);
						pV0->boundary() = true;
						pHalfEdge = HalfEdgeNext(pHalfEdge);

					} while (pHalfEdge != HalfFaceHalfEdge(pHF));
				}
			}
		};

		typedef OptimizeHexMesh<CHVertex, CHexOptimizeVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex> COHMesh;


	}; // namespace HMeshLib
}; // namespace MeshLib

#endif