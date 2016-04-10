#ifndef _OPTIMIZE_TET_MESH_H_
#define _OPTIMIZE_TET_MESH_H_

#include "TMeshLibHeaders.h"
#include "Parser\parser.h"
#include "Parser\traits_io.h"

namespace MeshLib
{
	namespace TMeshLib
	{

		class CTetOptimizeVertex : public CVertex
		{
		public:
			bool & diskBoundary() { return m_diskBoundary; };
			bool & boundary() { return m_boundary; };
			bool & disk() { return m_disk; }
			//std::vector<CTetOptimizeVertex *> & sameVertex() { return m_sameVertex; };

			void _from_string();

		protected:
			bool m_diskBoundary = false;
			bool m_boundary = false;
			//std::vector<CTetOptimizeVertex *> m_sameVertex;
			bool m_disk = false;
		};

		inline void CTetOptimizeVertex::_from_string()
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
		class OptimizeTetMesh : public CTMesh < HXV, V, HE, HXE, E, HF, F, HX >
		{
		public:
			typedef HXV CHVertex;
			typedef V	CVertex;
			typedef HE	CHalfEdge;
			typedef HXE CHEdge;
			typedef E	CEdge;
			typedef HF	CHalfFace;
			typedef HX	CTet;

			typedef TMeshTetIterator<HXV, V, HE, HXE, E, HF, F, HX> MeshTetIterator;
			typedef TetHalfFaceIterator<HXV, V, HE, HXE, E, HF, F, HX> TetHFIterator;
			typedef TMeshEdgeIterator<HXV, V, HE, HXE, E, HF, F, HX> MeshEdgeIterator;
			typedef TMeshFaceIterator<HXV, V, HE, HXE, E, HF, F, HX> MeshFaceIterator;
			typedef HalfFaceVertexIterator<HXV, V, HE, HXE, E, HF, F, HX> HalfFaceVertexIterator;
			typedef FaceVertexIterator<HXV, V, HE, HXE, E, HF, F, HX> FaceVertexIterator;
			typedef TMeshVertexIterator<HXV, V, HE, HXE, E, HF, F, HX> MeshVertexIterator;
			typedef TVertexVertexIterator<HXV, V, HE, HXE, E, HF, F, HX> VertexVertexIterator;

		public:
			OptimizeTetMesh() { ; };
			~OptimizeTetMesh() { ; };

			void labelBoundaryVertices();

		private:

		};

		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		void OptimizeTetMesh< HXV, V, HE, HXE, E, HF, F, HX >::labelBoundaryVertices()
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

		typedef OptimizeTetMesh<CTVertex, CTetOptimizeVertex, CHalfEdge, CTEdge, CEdge, CHalfFace, CFace, CTet> COTMesh;


	}; // namespace HMeshLib
}; // namespace MeshLib

#endif