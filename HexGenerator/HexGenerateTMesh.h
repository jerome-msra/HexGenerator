#ifndef _HEX_GENERATE_T_MESH_H_
#define _HEX_GENERATE_T_MESH_H_

#include "TMeshLibHeaders.h"

#include "Parser\parser.h"
#include "Parser\traits_io.h"

namespace MeshLib
{
	namespace TMeshLib
	{
		class CHexGenerateTetVertex : public CVertex
		{
		public:
			CPoint & coord3d() { return m_3dCoord; };

			bool & zero() { return m_zero; };

			void _from_string();

		protected:
			CPoint m_3dCoord;
			bool m_zero = false;
		};

		// read in traits
		inline void CHexGenerateTetVertex::_from_string()
		{
			CParser parser(m_string);

			for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
			{
				CToken * token = *iter;
				if (token->m_key == "coord3d")
				{
					token->m_value >> m_3dCoord;
				}
				else if (token->m_key == "zero")
				{
					m_zero = true;
				}
			}
		};

		class CHexGenerateTetHalfFace : public CHalfFace
		{
		public:
			CPoint & normal() { return m_normal; };

		protected:
			CPoint m_normal;
		};

		class CHexGenerateTetFace : public CFace
		{
		public:
			double & area() { return m_area; };

		protected:
			double m_area = 0;
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		class CHexGenerateTMesh : public CTMesh < TV, V, HE, TE, E, HF, F, T >
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

		public:
			CHexGenerateTMesh();
			~CHexGenerateTMesh();

			typedef TMeshTetIterator<TV, V, HE, TE, E, HF, F, T> MeshTetIterator;
			typedef TetHalfFaceIterator<TV, V, HE, TE, E, HF, F, T> TetHFIterator;
			typedef TMeshEdgeIterator<TV, V, HE, TE, E, HF, F, T> MeshEdgeIterator;
			typedef TMeshFaceIterator<TV, V, HE, TE, E, HF, F, T> MeshFaceIterator;
			typedef HalfFaceVertexIterator<TV, V, HE, TE, E, HF, F, T> HalfFaceVertexIterator;
			typedef FaceVertexIterator<TV, V, HE, TE, E, HF, F, T> FaceVertexIterator;
			typedef TMeshVertexIterator<TV, V, HE, TE, E, HF, F, T> MeshVertexIterator;
			typedef TVertexVertexIterator<TV, V, HE, TE, E, HF, F, T> VertexVertexIterator;

		public:

			/*! \brief Compute the normal of half face*/
			void _halfface_normal();

			T * _locateVertex(CPoint p, CPoint & localCoord);

			void _locateVertexTest();

			double minHeight() {
				if (m_heightComputed)
				{
					return m_minHeight;
				}
				else
				{
					_findHeight();
					return m_minHeight;
				}
			};

			double maxHeight() {
				if (m_heightComputed)
				{
					return m_maxHeight;
				}
				else
				{
					_findHeight();
					return m_maxHeight;
				}
			};

			double height();

			double middle();

			CPoint baryCoord(CFace * face, CPoint p);

			CPoint pullbackPosition(CFace * face, CPoint baryCoord);

			CPoint onCylinderPosition(CFace * face, CPoint baryCoord);

		protected:
			double m_minHeight = 0.0;
			double m_maxHeight = 0.0;

			bool m_heightComputed = false;

		private:
			void _findHeight();
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		CHexGenerateTMesh<TV, V, HE, TE, E, HF, F, T>::CHexGenerateTMesh()
		{
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		CHexGenerateTMesh<TV, V, HE, TE, E, HF, F, T>::~CHexGenerateTMesh()
		{
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		T * CHexGenerateTMesh<TV, V, HE, TE, E, HF, F, T>::_locateVertex(CPoint p, CPoint & localCoord)
		{
			for (std::list<CTet*>::iterator tIter = m_pTets.begin(); tIter != m_pTets.end(); tIter++)
			{
				CTet * pT = *tIter;

				Eigen::Matrix4d matTest;

				for (int i = 0; i < 4; i++)
				{
					CVertex * pV = TetVertex(pT, i);
					CPoint p = pV->position();
					for (int j = 0; j < 3; j++)
					{
						matTest(j, i) = p[j];
					}
					matTest(3, i) = 1.0;
				}

				//std::cout << matTest.determinant() << std::endl;

				// non diffeomorphic tet
				if (matTest.determinant() <= 0)
				{
					continue;
				}

				bool inside = true;
				for (int i = 0; i < 4; i++)
				{
					Eigen::Matrix4d mat;
					for (int j = 0; j < 4; j++)
					{
						CPoint pos;
						if (j == i)
						{
							pos = p;
						}
						else
						{
							CVertex * pTV = TetVertex(pT, j);
							pos = pTV->position();
						}
						for (int k = 0; k < 3; k++)
						{
							mat(k, j) = pos[k];
						}
					}

					for (int j = 0; j < 4; j++)
					{
						mat(3, j) = 1.0;
					}

					if (mat.determinant() < 0)
					{
						inside = false;
						break;
					}
				}

				if (inside)
				{
					// compute local coordinate

					CPoint pos[4];
					for (int i = 0; i < 4; i++)
					{
						CVertex * pV = TetVertex(pT, i);
						pos[i] = pV->position();
					}

					// local frame
					CPoint base[3];
					for (int b = 0; b < 3; b++)
					{
						base[b] = pos[b + 1] - pos[0];
					}

					CPoint localPt = p - pos[0];

					Eigen::Matrix3d A;
					for (int i = 0; i < 3; i++)
					{
						for (int j = 0; j < 3; j++)
						{
							A(i, j) = base[j][i];
						}
					}
					Eigen::Vector3d b;
					for (int i = 0; i < 3; i++)
					{
						b[i] = localPt[i];
					}

					Eigen::Vector3d x;
					x = A.fullPivLu().solve(b);
					for (int i = 0; i < 3; i++)
					{
						localCoord[i] = x(i);
					}

					return pT;
				}
			}

			return nullptr;
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		void CHexGenerateTMesh<TV, V, HE, TE, E, HF, F, T>::_locateVertexTest()
		{
			CPoint p[4];
			p[0] = CPoint(0, 0, 1);
			p[1] = CPoint(-1, 0, -0.5);
			p[2] = CPoint(1, -2, -0.5);
			p[3] = CPoint(1, 2, -0.5);

			Eigen::Matrix4d matTest;

			for (int i = 0; i < 4; i++)
			{
				CPoint pt = p[i];
				for (int j = 0; j < 3; j++)
				{
					matTest(j, i) = pt[j];
				}
				matTest(3, i) = 1.0;
			}

			std::cout << matTest.determinant() << std::endl;

			CPoint zero(0, 0, 0);

			for (int i = 0; i < 4; i++)
			{
				Eigen::Matrix4d mat;
				for (int j = 0; j < 4; j++)
				{
					CPoint pos;
					if (j == i)
					{
						pos = zero;
					}
					else
					{
						pos = p[j];
					}
					for (int k = 0; k < 3; k++)
					{
						mat(k, j) = pos[k];
					}
				}

				for (int j = 0; j < 4; j++)
				{
					mat(3, j) = 1.0;
				}

				std::cout << mat << std::endl;
				std::cout << mat.determinant() << std::endl;
			}
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		void CHexGenerateTMesh<TV, V, HE, TE, E, HF, F, T>::_halfface_normal()
		{
			for (std::list<HF*>::iterator hiter = m_pHalfFaces.begin(); hiter != m_pHalfFaces.end(); hiter++)
			{
				HF *pHF = *hiter;
				HE *pHE = HalfFaceHalfEdge(pHF);
				std::vector<CPoint> ps;
				for (int k = 0; k < 3; k++)
				{
					ps.push_back(HalfEdgeTarget(pHE)->position());
					pHE = HalfEdgeNext(pHE);
				}
				CPoint n = (ps[1] - ps[0]) ^ (ps[2] - ps[0]);
				F * face = HalfFaceFace(pHF);
				face->area() = n.norm() / 2.0;
				n = n / n.norm();
				pHF->normal() = n;
			}
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		CPoint CHexGenerateTMesh<TV, V, HE, TE, E, HF, F, T>::onCylinderPosition(CFace * face, CPoint baryCoord)
		{
			CPoint q[3];
			int i = 0;

			CHalfFace * pHF = FaceLeftHalfFace(face);
			CHalfEdge * pHalfEdge = HalfFaceHalfEdge(pHF);
			do
			{
				CVertex * pV0 = HalfEdgeTarget(pHalfEdge);
				q[i++] = pV0->position();
				pHalfEdge = HalfEdgeNext(pHalfEdge);

			} while (pHalfEdge != HalfFaceHalfEdge(pHF));

			CPoint p;

			for (i = 0; i < 3; i++)
			{
				p += q[i] * baryCoord[i];
			}

			return p;
		};


		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		CPoint CHexGenerateTMesh<TV, V, HE, TE, E, HF, F, T>::pullbackPosition(CFace * face, CPoint baryCoord)
		{
			CPoint q[3];
			int i = 0;

			CHalfFace * pHF = FaceLeftHalfFace(face);
			CHalfEdge * pHalfEdge = HalfFaceHalfEdge(pHF);
			do
			{
				CVertex * pV0 = HalfEdgeTarget(pHalfEdge);
				q[i++] = pV0->coord3d();
				pHalfEdge = HalfEdgeNext(pHalfEdge);

			} while (pHalfEdge != HalfFaceHalfEdge(pHF));

			CPoint p;

			for (i = 0; i < 3; i++)
			{
				p += q[i] * baryCoord[i];
			}

			return p;
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		CPoint CHexGenerateTMesh<TV, V, HE, TE, E, HF, F, T>::baryCoord(CFace * face, CPoint p)
		{
			CPoint q[3];
			int i = 0;

			CHalfFace * pHF = FaceLeftHalfFace(face);
			CPoint normal = pHF->normal();
			CHalfEdge * pHalfEdge = HalfFaceHalfEdge(pHF);
			do
			{
				CVertex * pV0 = HalfEdgeTarget(pHalfEdge);
				q[i++] = pV0->position();
				pHalfEdge = HalfEdgeNext(pHalfEdge);

			} while (pHalfEdge != HalfFaceHalfEdge(pHF));


			double s = face->area();

			CPoint bary;

			for (int i = 0; i < 3; i++)
			{
				bary[i] = ((q[(i + 1) % 3] - p) ^ (q[(i + 2) % 3] - p)) * normal / (2.0 * s);
			}

			//bary[2] = 1 - (bary[0] + bary[1]);

			return bary;
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		double CHexGenerateTMesh<TV, V, HE, TE, E, HF, F, T>::height()
		{
			if (m_heightComputed)
			{
				return m_maxHeight - m_minHeight;
			}
			else
			{
				_findHeight();
				return m_maxHeight - m_minHeight;
			}
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		double CHexGenerateTMesh<TV, V, HE, TE, E, HF, F, T>::middle()
		{
			if (m_heightComputed)
			{
				return (m_maxHeight + m_minHeight) / 2.0;
			}
			else
			{
				_findHeight();
				return (m_maxHeight + m_minHeight) / 2.0;
			}
		};

		template<typename TV, typename V, typename HE, typename TE, typename E, typename HF, typename F, typename T>
		void CHexGenerateTMesh<TV, V, HE, TE, E, HF, F, T>::_findHeight()
		{
			double maxH = -1e30;
			double minH = 1e30;

			for (std::list<CVertex*>::iterator vIter = m_pVertices.begin(); vIter != m_pVertices.end(); vIter++)
			{
				CVertex * pV = *vIter;
				CPoint p = pV->position();
				maxH = (p[2] > maxH) ? p[2] : maxH;
				minH = (p[2] < minH) ? p[2] : minH;
			}

			m_minHeight = minH;
			m_maxHeight = maxH;

			m_heightComputed = true;
		};


		typedef CHexGenerateTMesh<CTVertex, CHexGenerateTetVertex, CHalfEdge, CTEdge, CEdge, CHexGenerateTetHalfFace, CHexGenerateTetFace, CTet> CHGTMesh;
	};
};

#endif