#ifndef _HEX_GENERATE_H_MESH_H_
#define _HEX_GENERATE_H_MESH_H_

#include "HMeshLibHeaders.h"

namespace MeshLib
{
	namespace HMeshLib
	{

		class CHexGenerateVertex : public CVertex
		{
		public:
			bool & boundary() { return m_boundary; };
			bool & disk() { return m_onDisk; };
			CPoint & normal() { return m_normal; };
			CPoint & cylinder() { return m_cylinderCoord; };
			bool & onCircle() { return m_onCircle; };
			int & icoord() { return m_iCoord; };
			int & jcoord() { return m_jCoord; };

			void _from_string();

			void _to_string();

		protected:
			bool m_onDisk = false;
			bool m_boundary = false;
			CPoint m_normal;
			CPoint m_cylinderCoord;
			bool m_onCircle = false;
			int m_iCoord = 0;
			int m_jCoord = 0;
		};
		
		inline void CHexGenerateVertex::_from_string()
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
				else if (token->m_key == "cylinder")
				{
					token->m_value >> m_cylinderCoord;
				}
				else if (token->m_key == "disk")
				{
					m_onDisk = true;
				}
			}
		};

		inline void CHexGenerateVertex::_to_string()
		{
			m_string = "";

			std::string line;
			std::stringstream iss(line);
			
			if (m_onDisk)
			{
				iss << "disk ";
			}

			iss << "icoord=(" << m_iCoord << ") ";
			iss << "jcoord=(" << m_jCoord << ") ";
			iss << "cylinder=(" << m_cylinderCoord << ")";

			if (m_string.size() > 0)
			{
				m_string += " ";
			}

			m_string += iss.str();
		};

		class CHexGenerateHexHalfFace : public CHalfFace
		{
		public:
			CPoint & normal() { return m_normal; };

		protected:
			CPoint m_normal;
		};

		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		class HexGenerateHMesh : public CHMesh < HXV, V, HE, HXE, E, HF, F, HX >
		{
		public:
			typedef HXV CHVertex;
			typedef V	CVertex;
			typedef HE	CHalfEdge;
			typedef HXE CHEdge;
			typedef E	CEdge;
			typedef HF	CHalfFace;
			typedef HX	CHex;

			typedef HMeshHexIterator<HXV, V, HE, HXE, E, HF, F, HX> MeshHexIterator;
			typedef HexHalfFaceIterator<HXV, V, HE, HXE, E, HF, F, HX> HexHFIterator;
			typedef HMeshEdgeIterator<HXV, V, HE, HXE, E, HF, F, HX> MeshEdgeIterator;
			typedef HMeshFaceIterator<HXV, V, HE, HXE, E, HF, F, HX> MeshFaceIterator;
			typedef HalfFaceVertexIterator<HXV, V, HE, HXE, E, HF, F, HX> HalfFaceVertexIterator;
			typedef FaceVertexIterator<HXV, V, HE, HXE, E, HF, F, HX> FaceVertexIterator;
			typedef HMeshVertexIterator<HXV, V, HE, HXE, E, HF, F, HX> MeshVertexIterator;
			typedef HVertexVertexIterator<HXV, V, HE, HXE, E, HF, F, HX> VertexVertexIterator;

		public:
			HexGenerateHMesh();
			~HexGenerateHMesh();

			void generatePositionMap();

			void labelBoundaryVertices();

			void _halfface_normal();

			V * createVertex(int id);

			V * findVertex(CPoint p);

			void constructAllHexes(std::vector<std::vector<int>> hexList);

		private:

			std::map<CPoint, int> m_mapPositionToId;
		};

		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		HexGenerateHMesh< HXV, V, HE, HXE, E, HF, F, HX >::HexGenerateHMesh()
		{
		};

		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		HexGenerateHMesh< HXV, V, HE, HXE, E, HF, F, HX >::~HexGenerateHMesh()
		{
		};

		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		void HexGenerateHMesh<HXV, V, HE, HXE, E, HF, F, HX>::_halfface_normal()
		{
			for (std::list<CHalfFace*>::iterator hiter = m_pHalfFaces.begin(); hiter != m_pHalfFaces.end(); hiter++)
			{
				CHalfFace *pHF = *hiter;
				CHalfEdge *pHE = HalfFaceHalfEdge(pHF);
				std::vector<CPoint> ps;
				for (int k = 0; k < 4; k++)
				{
					ps.push_back(HalfEdgeTarget(pHE)->position());
					pHE = HalfEdgeNext(pHE);
				}
				CPoint n0 = (ps[1] - ps[0]) ^ (ps[2] - ps[0]);
				CPoint n1 = (ps[2] - ps[0]) ^ (ps[3] - ps[0]);
				CPoint n = (n0 + n1) / 2;
				n = n / n.norm();
				pHF->normal() = n;
			}
		};


		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		V * HexGenerateHMesh<HXV, V, HE, HXE, E, HF, F, HX>::findVertex(CPoint p)
		{
			//for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
			//{
			//	CVertex * pV = *vIter;
			//	CPoint pt = pV->position();
			//	if ((pt - p).norm() == 0)
			//	{
			//		return pV;
			//	}
			//}

			auto it = m_mapPositionToId.find(p);
			if (it != m_mapPositionToId.end())
			{
				int id = it->second;
				CVertex * pV = idVertex(id);
				return pV;
			}
			else
			{
				return nullptr;
			}
		}

		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		V * HexGenerateHMesh< HXV, V, HE, HXE, E, HF, F, HX >::createVertex(int id)
		{
			int vid = id;
			CVertex * v = new CVertex();
			v->id() = vid;
			m_pVertices.push_back(v);
			m_nVertices++;
			m_map_Vertices.insert(std::pair<int, CVertex *>(vid, v));
			return v;
		};

		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		void HexGenerateHMesh< HXV, V, HE, HXE, E, HF, F, HX >::generatePositionMap()
		{
			for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
			{
				CVertex * pV = *vIter;
				CPoint p = pV->position();
				m_mapPositionToId[p] = pV->id();
			}
		};

		template<typename HXV, typename V, typename HE, typename HXE, typename E, typename HF, typename F, typename HX>
		void HexGenerateHMesh< HXV, V, HE, HXE, E, HF, F, HX >::labelBoundaryVertices()
		{
			for (std::list<HF*>::iterator hfIter = m_pHalfFaces.begin(); hfIter != m_pHalfFaces.end(); hfIter++)
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
		void HexGenerateHMesh< HXV, V, HE, HXE, E, HF, F, HX >::constructAllHexes(std::vector<std::vector<int>> hexList)
		{
			m_nHexs = 0;

			for (std::vector<std::vector<int>>::iterator hIter = hexList.begin(); hIter != hexList.end(); hIter++)
			{
				int vId[8];
				int hId = (int)(hIter - hexList.begin());
				std::vector<int> hexVertId = *hIter;
				for (size_t k = 0; k < 8; k++)
				{
					vId[k] = hexVertId[k];
				}

				CHex * pHex = new CHex();
				m_pHexs.push_back(pHex);
				m_map_Hexs.insert(std::pair<int, CHex *>(hId, pHex));
				m_nHexs++;
				_construct_hex(pHex, hId, vId);

			}

			_construct_faces();
			_construct_edges();

			m_nEdges = (int)m_pEdges.size();

			for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
			{
				CVertex * pV = *vIter;
				if (pV->id() > m_maxVertexId)
				{
					m_maxVertexId = pV->id();
				}
			}
		};

		typedef HexGenerateHMesh<CHVertex, CHexGenerateVertex, CHalfEdge, CHEdge, CEdge, CHexGenerateHexHalfFace, CFace, CHex> CHGHMesh;
	};
};
#endif