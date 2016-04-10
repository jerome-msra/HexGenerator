#ifndef _MERGE_HEX_MESH_H_
#define _MERGE_HEX_MESH_H_

#include "HMeshLibHeaders.h"

#include "Parser\parser.h"
#include "Parser\traits_io.h"

#include <queue>

namespace MeshLib
{
	namespace HMeshLib
	{

		class CMergeHexVertex : public CVertex
		{
		public:
			bool & boundary() { return m_boundary; };
			bool & zero() { return m_zero; };
			bool & onDisk() { return m_onDisk; };
			bool & visited() { return m_visited; };
			int & icoord() { return m_iCoord; };
			int & jcoord() { return m_jCoord; };

			int & zeroDistance() { return m_zeroDistance; };

			CPoint & cylinder() { return m_cylinder; };

			void _from_string();

			void _to_string();

		protected:
			bool m_boundary = false;
			bool m_zero = false;
			bool m_onDisk = false;
			bool m_visited = false;
			CPoint m_cylinder;
			int m_iCoord = 0;
			int m_jCoord = 0;

			int m_zeroDistance = 0;
		};

		inline void CMergeHexVertex::_from_string()
		{
			CParser parser(m_string);

			for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
			{
				CToken * token = *iter;
				if (token->m_key == "icoord")
				{
					std::string line = strutil::trim(token->m_value, "()");
					m_iCoord = strutil::parseString<int>(line);
				}
				else if (token->m_key == "jcoord")
				{
					std::string line = strutil::trim(token->m_value, "()");
					m_jCoord = strutil::parseString<int>(line);
				}
				else if (token->m_key == "zero") // manually marked
				{
					m_zero = true;
				}
				else if (token->m_key == "disk")
				{
					m_onDisk = true;
				}
				else if (token->m_key == "cylinder")
				{
					token->m_value >> m_cylinder;
				}
			}
		};
		
		inline void CMergeHexVertex::_to_string()
		{
			std::string line;
			std::stringstream iss(line);
			iss << "zerodist=(" << m_zeroDistance << ")";
			if (m_string.length() > 0)
			{
				m_string += " ";
				m_string += iss.str();
			}
		};

		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		class CMergeHexMesh : public CHMesh < HXV, V, HE, HXE, E, HF, F, HX >
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
			CMergeHexMesh() { ; };
			~CMergeHexMesh() { ; };

			std::vector<CVertex*> & baseZero() { return m_baseZero; };

			void labelBoundaryVertices();

			void labelBaseZeros();

		private:
			std::vector<CVertex*> m_baseZero;
		};

		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		void CMergeHexMesh< HXV, V, HE, HXE, E, HF, F, HX >::labelBoundaryVertices()
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

		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		void CMergeHexMesh< HXV, V, HE, HXE, E, HF, F, HX >::labelBaseZeros()
		{
			std::queue<CVertex*> myqueue;
			for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
			{
				CVertex * pV = *vIter;
				if (pV->onDisk())
				{
					if (pV->zero())
					{
						pV->visited() = true;
						m_baseZero.push_back(pV);
						myqueue.push(pV);
					}
				}
			}

			while (! myqueue.empty())
			{
				CVertex * pV = myqueue.front();
				myqueue.pop();
				for (VertexVertexIterator vvIter(this, pV); !vvIter.end(); vvIter++)
				{
					CVertex * pVV = *vvIter;
					if (pVV->onDisk() && pVV->icoord() == pV->icoord() && pVV->visited() == false)
					{
						pVV->visited() = true;
						m_baseZero.push_back(pVV);
						myqueue.push(pVV);
					}
				}
			}

			// compute the distance to the correspondent base zero points
			for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
			{
				CVertex * pV = *vIter;
				if (pV->onDisk())
				{
					for (int i = 0; i < m_baseZero.size(); i++)
					{
						CVertex * pZero = m_baseZero[i];
						if (pV->cylinder()[2] == pZero->cylinder()[2] && pV->jcoord() == pZero->jcoord())
						{
							pV->zeroDistance() = pV->icoord() - pZero->icoord();
						}
					}
				}
			}
		};

		typedef CMergeHexMesh<CHVertex, CMergeHexVertex, CHalfEdge, CHEdge, CEdge, CHalfFace, CFace, CHex> CMHMesh;

	}; // namespace HMeshLib
}; // namespace MeshLib

#endif