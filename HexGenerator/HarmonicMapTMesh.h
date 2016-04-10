#ifndef _HARMONIC_MAP_TMESH_H_
#define _HARMONIC_MAP_TMESH_H_

#include "TMeshLibHeaders.h"

namespace MeshLib
{
	namespace TMeshLib
	{

		class CHarmonicVertex : public CVertex
		{
		public:

			bool & fixed() { return m_fixed; };

			bool & boundary() { return m_boundary; };

			double & u() { return m_u; };

			int & idx() { return m_idx; };

			CPoint & huv() { return m_huv; };

			CPoint & coord3d() { return m_3dCoord; };

			void _to_string();

		protected:
			int m_idx = 0;
			bool m_fixed = false;
			bool m_boundary = false;
			double m_u;
			CPoint m_huv;
			CPoint m_3dCoord;
		};


		inline void CHarmonicVertex::_to_string()
		{
			m_string = "";

			std::string line;
			std::stringstream iss(line);

			iss << "coord3d=(" << m_3dCoord[0] << " " << m_3dCoord[1] << " " << m_3dCoord[2] << ")";

			if (m_string.size() > 0)
			{
				m_string += " ";
			}

			m_string += iss.str();
		}

		class CHarmonicEdge : public CEdge
		{
		public:
			double & weight() { return m_weight; };

			double & length() { return m_length; };
		protected:
			double m_weight = 0.0;
			double m_length = 1.0;
		};

		class CHarmonicHalfFace : public CHalfFace
		{
		public:
			CPoint & normal() { return m_normal; };

		protected:
			CPoint	 m_normal;
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		class CHarmonicMapTMesh : public CTMesh<TV, V, HE, TE, E, HF, F, T>
		{
		public:

			typedef TV	CTVertex;
			typedef V	CVertex;
			typedef HE	CHalfEdge;
			typedef TE	CTEdge;
			typedef E	CEdge;
			typedef HF	CHalfFace;
			typedef F	CFace;
			typedef T	CTet;

			typedef TMeshTetIterator<TV, V, HE, TE, E, HF, F, T> MeshTetIterator;
			typedef TetHalfFaceIterator<TV, V, HE, TE, E, HF, F, T> TetHFIterator;
			typedef TMeshEdgeIterator<TV, V, HE, TE, E, HF, F, T> MeshEdgeIterator;
			typedef TMeshFaceIterator<TV, V, HE, TE, E, HF, F, T> MeshFaceIterator;
			typedef HalfFaceVertexIterator<TV, V, HE, TE, E, HF, F, T> HalfFaceVertexIterator;
			typedef FaceVertexIterator<TV, V, HE, TE, E, HF, F, T> FaceVertexIterator;
			typedef TMeshVertexIterator<TV, V, HE, TE, E, HF, F, T> MeshVertexIterator;
			typedef TVertexVertexIterator<TV, V, HE, TE, E, HF, F, T> VertexVertexIterator;

			typedef EdgeTEdgeIterator<TV, V, HE, TE, E, HF, F, T> EdgeTEdgeIterator;
			typedef TVertexEdgeIterator<TV, V, HE, TE, E, HF, F, T> VertexEdgeIterator;
		public:
			CHarmonicMapTMesh();
			~CHarmonicMapTMesh();

			void _labelBoundaryVertices();
			void _computeHalfFaceNormal();
		private:

		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		CHarmonicMapTMesh<TV, V, HE, TE, E, HF, F, T>::CHarmonicMapTMesh()
		{
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		CHarmonicMapTMesh<TV, V, HE, TE, E, HF, F, T>::~CHarmonicMapTMesh()
		{
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		void CHarmonicMapTMesh<TV, V, HE, TE, E, HF, F, T>::_labelBoundaryVertices()
		{
			for (std::list<CHalfFace*>::iterator hfIter = m_pHalfFaces.begin(); hfIter != m_pHalfFaces.end(); hfIter++)
			{
				CHalfFace * pHF = *hfIter;
				// boundary
				CHalfFace * pHFDual = HalfFaceDual(pHF);
				if (pHFDual == nullptr)
				{
					for (HalfFaceVertexIterator hfvIter(this, pHF); !hfvIter.end(); hfvIter++)
					{
						CVertex * pV = *hfvIter;
						pV->boundary() = true;
					}
				}
			}
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		void CHarmonicMapTMesh<TV, V, HE, TE, E, HF, F, T>::_computeHalfFaceNormal()
		{
			for (std::list<CHalfFace*>::iterator hiter = m_pHalfFaces.begin(); hiter != m_pHalfFaces.end(); hiter++)
			{
				CHalfFace *pHF = *hiter;
				CHalfEdge *pHE = HalfFaceHalfEdge(pHF);
				std::vector<CPoint> ps;
				for (int k = 0; k < 3; k++)
				{
					ps.push_back(HalfEdgeTarget(pHE)->position());
					pHE = HalfEdgeNext(pHE);
				}
				CPoint n = (ps[1] - ps[0]) ^ (ps[2] - ps[0]);
				n = n / n.norm();
				pHF->normal() = n;
			}
		};

		typedef CHarmonicMapTMesh<CTVertex, CHarmonicVertex, CHalfEdge, CTEdge, CHarmonicEdge, CHarmonicHalfFace, CFace, CTet> CHTMesh;


	};
};

#endif