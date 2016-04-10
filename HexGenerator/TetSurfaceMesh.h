#ifndef _TET_SURFACE_MESH_H_
#define _TET_SURFACE_MESH_H_

#include <map>
#include <vector>
#include <list>

#include "Mesh\Vertex.h"
#include "Mesh\Edge.h"
#include "Mesh\Face.h"
#include "Mesh\HalfEdge.h"

#include "Mesh\BaseMesh.h"
#include "Mesh\iterators.h"
#include "Parser\parser.h"
#include "Parser\traits_io.h"

#include "Mesh\boundary.h"

namespace MeshLib
{
	class CSurfaceVertex : public CVertex
	{
	public:

		CSurfaceVertex() : CVertex() {};

		~CSurfaceVertex() {};

		CSurfaceVertex(const CSurfaceVertex & v)
		{
			m_id = v.m_id;
			m_father = v.m_father;
			m_rgb = v.m_rgb;
			m_u = v.m_u;
			m_huv = v.m_huv;
			m_group = v.m_group;
			m_point = v.m_point;
			m_z = v.m_z;
			m_uv = v.m_uv;
		};

		int & group() {
			return m_group;
		};

		CPoint & rgb() { return m_rgb; };

		/*! \brief return the reference to the id on the tet mesh*/
		int & father() { return m_father; };

		int & idx() { return m_idx; };

		double & u() { return m_u; };

		CPoint2 & huv() { return m_huv; };
		
		double & z() { return m_z; };
		
		bool & fixed() { return m_fixed; };

		bool & diskBoundary() { return m_diskBoundary; };

		/*! \brief Save vertex uv coordinates to the vertex string */
		void _to_string();

		void _from_string();


	protected:
		/*! \brief the id of the correspond vertex on the tet mesh*/
		int m_father;
		int m_group;
		int m_idx;
		CPoint m_rgb;
		// for laplacian solver
		double m_u;
		// result of disk harmonic mapping
		CPoint2 m_huv;
		double m_z;
		bool m_fixed;
		bool m_boundary;

		bool m_diskBoundary;
	};

	inline void CSurfaceVertex::_from_string()
	{
		CParser parser(m_string);

		for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
		{
			CToken * token = *iter;
			if (token->m_key == "father")
			{
				std::string line = strutil::trim(token->m_value, "()");
				m_father = strutil::parseString<int>(line);
			}
			else if (token->m_key == "uv")
			{
				token->m_value >> m_uv;
			}
			else if (token->m_key == "rgb")
			{
				token->m_value >> m_rgb;
			}
			else if (token->m_key == "group")
			{
				std::string line = strutil::trim(token->m_value, "()");
				m_group = strutil::parseString<int>(line);
			}
		}
	};

	// Save vertex uv coordinates to the vertex string 
	inline void CSurfaceVertex::_to_string()
	{
		m_string = "";

		std::string line;
		std::stringstream iss(line);

		iss << "group=(" << m_group << ")";
		iss << " father=(" << m_father << ")";
		iss << " rgb=(" << m_rgb[0] << " " << m_rgb[1] << " " << m_rgb[2] << ")";

		if (m_string.size() > 0)
		{
			m_string += " ";
		}

		m_string += iss.str();
	};

	class CSurfaceEdge :public CEdge
	{
	public:
		double & length() {
			return m_length;
		};

		/*! \brief edge weight*/
		double & weight() { return m_weight; };

	protected:
		double m_length = 0;
		double m_weight = 0;
	};

	class CSurfaceHalfEdge : public CHalfEdge
	{
	public:

		/*! \brief corner angle*/
		double & angle() { return m_angle; };

	protected:
		double m_angle = 0;

	};


	template<typename V, typename E, typename F, typename H>
	class CTetSurfaceMesh : public CBaseMesh < V, E, F, H >
	{
	public:
		typedef CLoop<V, E, F, H> CLoop;
		typedef CLoopSegment<V, E, F, H> CSegment;
		typedef CBoundary<V, E, F, H> CBoundary;
	public:
		typedef V CVertex;
		typedef E CEdge;
		typedef F CFace;
		typedef H CHalfEdge;

		typedef MeshVertexIterator<V, E, F, H> MeshVertexIterator;
		typedef MeshEdgeIterator<V, E, F, H> MeshEdgeIterator;
		typedef VertexVertexIterator<V, E, F, H> VertexVertexIterator;
		typedef VertexEdgeIterator<V, E, F, H> VertexEdgeIterator;
		typedef VertexInHalfedgeIterator<V, E, F, H> VertexInHalfedgeIterator;
		typedef VertexFaceIterator<V, E, F, H> VertexFaceIterator;
		typedef FaceVertexIterator<V, E, F, H> FaceVertexIterator;

		//for choosing edge sign
		typedef MeshFaceIterator<V, E, F, H> MeshFaceIterator;
		typedef FaceEdgeIterator<V, E, F, H> FaceEdgeIterator;
		typedef FaceHalfedgeIterator<V, E, F, H> FaceHalfedgeIterator;

	public:

		void generatePointMap();

		void labelBoundaryVertices();

		V * findVertex(CPoint p, int g);

		bool existVertex(int id);

	private:

		std::map<CPoint, int> m_mapPointToId;
	};

	template<typename V, typename E, typename F, typename H>
	bool CTetSurfaceMesh<V, E, F, H>::existVertex(int id)
	{
		std::map<int, CVertex*>::iterator it = m_map_vert.find(id);
		return it != m_map_vert.end();
	};

	template<typename V, typename E, typename F, typename H>
	void CTetSurfaceMesh<V, E, F, H>::generatePointMap()
	{
		for (std::list<CVertex*>::iterator vIter = m_verts.begin(); vIter != m_verts.end(); vIter++)
		{
			CVertex * pV = *vIter;
			CPoint pt = pV->point();
			m_mapPointToId[pt] = pV->id();
		}
	};

	template<typename V, typename E, typename F, typename H>
	void CTetSurfaceMesh<V, E, F, H>::labelBoundaryVertices()
	{
		for (std::list<E*>::iterator eiter = m_edges.begin(); eiter != m_edges.end(); ++eiter)
		{
			CEdge *     edge = *eiter;
			CHalfEdge * he[2];

			he[0] = edgeHalfedge(edge, 0);
			he[1] = edgeHalfedge(edge, 1);

			assert(he[0] != NULL);

			if (he[1] != NULL)
			{
				assert(he[0]->target() == he[1]->source() && he[0]->source() == he[1]->target());

				if (he[0]->target()->id() < he[0]->source()->id())
				{
					edge->halfedge(0) = he[1];
					edge->halfedge(1) = he[0];
				}

				assert(edgeVertex1(edge)->id() < edgeVertex2(edge)->id());

				he[0]->vertex()->boundary() = false;
				he[0]->he_prev()->vertex()->boundary() = false;
			}
			else
			{
				he[0]->vertex()->boundary() = true;
				he[0]->he_prev()->vertex()->boundary() = true;
			}

		}
	};

	template<typename V, typename E, typename F, typename H>
	V * CTetSurfaceMesh<V, E, F, H>::findVertex(CPoint p, int g)
	{
		//for (std::list<CVertex*>::iterator vIter = m_verts.begin(); vIter != m_verts.end(); vIter++)
		//{
		//	CVertex * pV = *vIter;
		//	CPoint pt = pV->point();
		//	if ((pt - p).norm() == 0 && pV->group() == g)
		//	{
		//		return pV;
		//	}
		//}

		auto it = m_mapPointToId.find(p);
		if (it != m_mapPointToId.end())
		{
			int id = it->second;
			CVertex * pV = idVertex(id);
			if (pV->group() == g)
			{
				return pV;
			}
		}

		return nullptr;
	};

	typedef CTetSurfaceMesh<CSurfaceVertex, CSurfaceEdge, CFace, CSurfaceHalfEdge> CTSMesh;
};


#endif