/*!
*      \file TetToHexGenerator.h
*      \brief Generate the Hexahedron mesh from Tetrahedron Mesh. Input is the tet mesh mapped onto the standard cylinder.
*
*	   \author Jerome Jiang
*      \date 12/18/2015
*
*/


#ifndef _TET_TO_HEX_GENERATOR_H_
#define _TET_TO_HEX_GENERATOR_H_

#include "HexGenerateHMesh.h"
#include "HexGenerateTMesh.h"

#include <Eigen\Dense>

//#include <algorithm>

namespace MeshLib
{
	namespace HMeshLib
	{
		template<typename T, typename H>
		class CTetHexGenerator
		{
		public:
			CTetHexGenerator(T * pTetMesh);
			~CTetHexGenerator();

			/*! \brief Generate the Hexahedron mesh from Tet mesh*/
			void _genHex();

			/*! \brief Return the hex mesh*/
			H * _hexMesh() { return m_pHexMesh; };

			void _writeHex(const char * output);

			/*! \brief Return the number of hexes along the tunnel loop*/
			int & _numLayer() { return m_numLayer; };

			/*! \brief Return the number of hexes along the handle loop*/
			int & _xQuad() { return m_numXQuad; };

			/*! \brief Return the number of hexes along the handle loop*/
			int & _yQuad() { return m_numYQuad; };

			bool & _zero() { return m_zero; };

		protected:

			/*! \brief Generate the standard Hexahedron Mesh with the geometry of cylinder
			*/
			void __genStdHexMesh();

			void __rotateToZeroDirection();

			/*! \brief Locate all the hex vertices in the tet cylinder
			*/
			void __locateVertices();

			/*! \brief stick the surface of standard hex mesh which is now a cube to standard cylinder
			*/
			void __locateSurfaceViaNormal();

			/*! \brief find the intersection between the ray starting from p with direction of np and the surface of tet mesh
			*/
			CPoint __intersection(CPoint startPoint, CPoint normal);

		private:

			void ___detectTetOrientation();

			CPoint __zeroPointDirection();

			Eigen::Matrix3d ___convertToRotationMatrix(CPoint axis, double angle);

		private:
			T * m_pTetMesh;
			H * m_pHexMesh;

			int m_numLayer = 40;
			int m_numXQuad = 7;
			int m_numYQuad = 7;

			bool m_zero = false;

			double m_height;

			double m_middle;

			CPoint m_zeroDirection;
		};


		template<typename T, typename H>
		CTetHexGenerator<T, H>::CTetHexGenerator(T * pTetMesh)
		{
			m_pTetMesh = pTetMesh;
			m_height = m_pTetMesh->height();
			m_pTetMesh->_halfface_normal();


			m_pHexMesh = new H();
		};

		template<typename T, typename H>
		CTetHexGenerator<T, H>::~CTetHexGenerator()
		{
		};

		template<typename T, typename H>
		CPoint CTetHexGenerator<T, H>::__zeroPointDirection()
		{
			std::vector<CPoint> topZero;
			std::vector<CPoint> bottomZero;
			for (T::MeshVertexIterator vIter(m_pTetMesh); !vIter.end(); vIter++)
			{
				T::CVertex * pV = *vIter;
				if (pV->zero())
				{
					CPoint p = pV->position();
					if (p[2] - m_pTetMesh->middle() > 0)
					{
						topZero.push_back(p);
					}
					else
					{
						bottomZero.push_back(p);
					}
				}
			}

			assert(topZero.size() > 0 && topZero.size() % 2 == 0);
			return topZero[1] - topZero[0];
		}

		template<typename T, typename H>
		void CTetHexGenerator<T, H>::__genStdHexMesh()
		{
			//// number of layer of hexes
			//m_numLayer = numLayer;

			//// number of quad in each side, such that the number of hexes in each layer is m_numXQuad^2
			//m_numXQuad = numQuad;
			//m_numYQuad = 7;
			std::vector<std::vector<int>> layerHex;
			//layerHex.resize(m_numLayer);

			double unitHeight = m_height / (double)m_numLayer;

			double zlbase = m_pTetMesh->minHeight();

			double unitXQuad = (sqrt(2) / 2.0) / (double)(m_numXQuad);

			double unitYQuad = (sqrt(2) / 2.0) / (double)(m_numYQuad);

			double basePoint = 0.5 - (sqrt(2) / 4.0);

			int numHexPerLayer = m_numXQuad * m_numYQuad;

			int numVerticesPerLayer = (m_numXQuad + 1) * (m_numYQuad + 1);

			// for layer - z
			int vertId = 0;
			int hexId = 0;

			for (size_t l = 0; l < m_numLayer + 1; l++)
			{
				// from the bottom to the top
				double zl = zlbase + l * unitHeight;
				int vertIdStart = vertId;

				//for x
				for (size_t i = 0; i < m_numXQuad + 1; i++)
				{
					// for y
					for (size_t j = 0; j < m_numYQuad + 1; j++)
					{
						CPoint p;
						p[0] = basePoint + i * unitXQuad;
						p[1] = basePoint + j * unitYQuad;
						p[2] = zl;
						
						bool onCircle = false;
						
						if ((i == 0 && j != 0) || (i == m_numXQuad && j != m_numYQuad))
						{
							double newXX = (0.5 * 0.5) - (p[1] - 0.5) * (p[1] - 0.5);
							double newX0 = sqrt(newXX) + 0.5;
							double newX1 = -sqrt(newXX) + 0.5;
							if (fabs(newX0 - p[0]) < fabs(newX1 - p[0]))
							{
								p[0] = newX0;
							}
							else
							{
								p[0] = newX1;
							}
							onCircle = true;
						}

						if ((i != 0 && j == 0) || (i != m_numXQuad && j == m_numYQuad))
						{
							double newXX = (0.5 * 0.5) - (p[0] - 0.5) * (p[0] - 0.5);
							double newX0 = sqrt(newXX) + 0.5;
							double newX1 = -sqrt(newXX) + 0.5;
							if (fabs(newX0 - p[1]) < fabs(newX1 - p[1]))
							{
								p[1] = newX0;
							}
							else
							{
								p[1] = newX1;
							}
							onCircle = true;
						}

						if ((i == 0 && j == 0) || (i == m_numXQuad && j == m_numYQuad))
						{
							onCircle = true;
						}

						H::CVertex * hexVert = m_pHexMesh->createVertex(vertId);
						hexVert->disk() = (zl == zlbase || zl == m_height);
						hexVert->icoord() = (int)i;
						hexVert->jcoord() = (int)j;
						hexVert->position() = p;
						if (zl == m_height)
						{
							hexVert->position()[2] -= 0.0001;
						}
						else if (zl == zlbase)
						{
							hexVert->position()[2] += 0.0001;
						}
						hexVert->onCircle() = onCircle;
						vertId++;
					}
				}

				int vertIdEnd = vertId;

				// create list of tet
				if (l != m_numLayer)
				{
					for (int hi = 0; hi < m_numXQuad; hi++)
					{
						for (int hj = 0; hj < m_numYQuad; hj++)
						{
							std::vector<int> hexVertList;
							hexVertList.resize(8);
							hexVertList[0] = vertIdStart + hi * (m_numYQuad + 1) + hj;
							hexVertList[1] = hexVertList[0] + (m_numYQuad + 1);
							hexVertList[2] = hexVertList[1] + 1;
							hexVertList[3] = hexVertList[0] + 1;

							hexVertList[4] = hexVertList[3] + numVerticesPerLayer;
							hexVertList[5] = hexVertList[2] + numVerticesPerLayer;
							hexVertList[6] = hexVertList[1] + numVerticesPerLayer;
							hexVertList[7] = hexVertList[0] + numVerticesPerLayer;

							layerHex.push_back(hexVertList);
						}
					}
				}
			}

			m_pHexMesh->constructAllHexes(layerHex);

			m_pHexMesh->labelBoundaryVertices();

			m_pHexMesh->_halfface_normal();
		};

		template<typename T, typename H>
		CPoint CTetHexGenerator<T, H>::__intersection(CPoint startPoint, CPoint normal)
		{
			CPoint result;

			std::vector<double> lambdaList;
			std::vector<CPoint> baryCoordList;
			std::vector<T::CFace*> faceList;
			double minBary = 1e+30;
			CPoint minBaryCoord;
			for (T::MeshFaceIterator fIter(m_pTetMesh); !fIter.end(); fIter++)
			{
				T::CFace * pF = *fIter;
				T::CHalfFace * pFLeft = m_pTetMesh->FaceLeftHalfFace(pF);
				T::CHalfFace * pFRight = m_pTetMesh->FaceRightHalfFace(pF);
				CPoint planeNormal;
				CPoint planePoint;

				if (pFRight != NULL)
				{
					continue;
				}

				planeNormal = pFLeft->normal() * (-1);
				T::CHalfEdge * he = m_pTetMesh->HalfFaceHalfEdge(pFLeft);
				T::CVertex * vert = m_pTetMesh->HalfEdgeTarget(he);
				/*if (vert->id() == 367 || vert->id() == 850 || vert->id() == 866)
				{
				std::cout << "This Face" << std::endl;
				}*/
				planePoint = vert->position();

				double parameter = planeNormal * planePoint;

				double lambda = (parameter - planeNormal * startPoint) / (planeNormal * normal);

				CPoint intersection = startPoint + normal * lambda;

				CPoint baryCoord = m_pTetMesh->baryCoord(pF, intersection);

				if ((fabs(baryCoord[0]) < 0.3 || baryCoord[0] >= 0) && (fabs(baryCoord[1]) < 0.3 || baryCoord[1] >= 0) && (fabs(baryCoord[2]) < 0.3 || baryCoord[2] >= 0))
				//if (baryCoord[0] >= 0 && baryCoord[1] >= 0 && baryCoord[2] >= 0)
				{
					lambdaList.push_back(lambda);
					baryCoordList.push_back(baryCoord);
					faceList.push_back(pF);
				}
				else if (minBary > fabs(std::min(std::min(baryCoord[0], baryCoord[1]), baryCoord[2])))
				{
					minBary = fabs(std::min(std::min(baryCoord[0], baryCoord[1]), baryCoord[2]));
					minBaryCoord = baryCoord;
				}

			}

			if (lambdaList.size() == 0)
			{
				std::cout << "[GenHex]minBary " << minBary << std::endl;
				std::cout << "[GenHex]minBaryCoord " << minBaryCoord << std::endl;
				return CPoint(0, 0, 0);
			}

			double minLambda = 1e+30;
			size_t minLambdaIdx = 0;
			for (size_t i = 0; i < lambdaList.size(); i++)
			{
				double currentLambda = lambdaList[i];
				if (fabs(currentLambda) < fabs(minLambda))
				{
					minLambda = currentLambda;
					minLambdaIdx = i;
				}
			}

			CPoint baryCoord = baryCoordList[minLambdaIdx];
			T::CFace * pF = faceList[minLambdaIdx];

			result = m_pTetMesh->pullbackPosition(pF, baryCoord);
			return result;
		};

		template<typename T, typename H>
		void CTetHexGenerator<T, H>::___detectTetOrientation()
		{
			for (T::MeshTetIterator tIter(m_pTetMesh); !tIter.end(); tIter++)
			{
				T::CTet * pT = *tIter;
				T::CVertex * vertList[4];
				for (int i = 0; i < 4; i++)
				{
					T::CVertex * pV = m_pTetMesh->TetVertex(pT, i);
					vertList[i] = pV;
				}

				Eigen::Matrix4d mat;

				for (int i = 0; i < 4; i++)
				{
					T::CVertex * pV = vertList[i];
					CPoint p = pV->position();
					for (int j = 0; j < 3; j++)
					{
						mat(j, i) = p[j];
					}
					mat(3, i) = 1.0;
				}

				double det = mat.determinant();

				if (det <= 0)
				{
					std::cout << det << std::endl;
				}
			}
		};

		template<typename T, typename H>
		void CTetHexGenerator<T, H>::__locateVertices()
		{

			//___detectTetOrientation();

			for (H::MeshVertexIterator hvIter(m_pHexMesh); !hvIter.end(); hvIter++)
			{
				H::CVertex * pHV = *hvIter;
				if (pHV->onCircle())
				{
					continue;
				}
				CPoint p = pHV->position();
				CPoint localCoord;

				std::cout << "\r";
				std::cout << "[GenHex]Locating Vertex " << pHV->id() << "/" << m_pHexMesh->numVertices() << "...";
				std::cout << p;
				T::CTet * pTet = m_pTetMesh->_locateVertex(p, localCoord);
				if (pTet != nullptr)
				{
					CPoint newPos;

					CPoint pos[4];
					for (int i = 0; i < 4; i++)
					{
						T::CVertex * pV = m_pTetMesh->TetVertex(pTet, i);
						pos[i] = pV->coord3d();
					}

					CPoint base[3];
					for (int b = 0; b < 3; b++)
					{
						base[b] = pos[b + 1] - pos[0];
					}

					Eigen::Matrix3d A;
					for (int i = 0; i < 3; i++)
					{
						for (int j = 0; j < 3; j++)
						{
							A(i, j) = base[j][i];
						}
					}

					Eigen::Vector3d x;
					for (int i = 0; i < 3; i++)
					{
						x[i] = localCoord[i];
					}

					Eigen::Vector3d p = A * x;

					newPos[0] = p(0) + pos[0][0];
					newPos[1] = p(1) + pos[0][1];
					newPos[2] = p(2) + pos[0][2];

					pHV->cylinder() = pHV->position();
					pHV->position() = newPos;

				}
				else
				{
					std::cout << std::endl;
					std::cout << "[GenHex]Locating Failed for vertex " << pHV->id() << std::endl;
				}
			}
		};

		template<typename T, typename H>
		void CTetHexGenerator<T, H>::__locateSurfaceViaNormal()
		{
			for (H::MeshVertexIterator vIter(m_pHexMesh); !vIter.end(); vIter++)
			{
				H::CVertex * pV = *vIter;

				if (!pV->onCircle())
				{
					continue;
				}
				std::list<H::CHalfFace *> * vertexHalfFaces = m_pHexMesh->VertexHalfFaceList(pV);
				CPoint n(0, 0, 0);

				for (std::list<H::CHalfFace*>::iterator hfIter = vertexHalfFaces->begin(); hfIter != vertexHalfFaces->end(); hfIter++)
				{
					H::CHalfFace * pHF = *hfIter;
					H::CHalfFace * pHFD = m_pHexMesh->HalfFaceDual(pHF);
					if (pHFD != nullptr)
					{
						continue;
					}
					CPoint normal = pHF->normal();
					n += normal;
				}

				n /= (double)(vertexHalfFaces->size());
				pV->normal() = (n * (-1)) / n.norm();
			}

			int numNoIntersection = 0;

			std::cout << std::endl;

			for (H::MeshVertexIterator vIter(m_pHexMesh); !vIter.end(); vIter++)
			{
				H::CVertex * pV = *vIter;

				if (!pV->onCircle())
				{
					continue;
				}
				CPoint p = pV->position();
				CPoint np = pV->normal();
				
				std::cout << "\r";
				std::cout << "[GenHex]Locating Surface Vertex " << pV->id() << "/" << m_pHexMesh->numVertices() << "...";

				CPoint newPos = __intersection(p, np);
				if (newPos == CPoint(0, 0, 0))
				{
					std::cout << pV->position() << std::endl;
					numNoIntersection++;;
				}
				pV->cylinder() = pV->position();
				pV->position() = newPos;
			}

			std::cout << "[GenHex] # of No Intersections " << numNoIntersection << std::endl;
		};

		template<typename T, typename H>
		Eigen::Matrix3d CTetHexGenerator<T, H>::___convertToRotationMatrix(CPoint axis, double angle)
		{
			Eigen::Matrix3d rotation;

			rotation(0, 0) = cos(angle) + axis[0] * axis[0] * (1 - cos(angle));
			rotation(0, 1) = axis[0] * axis[1] * (1 - cos(angle)) - axis[2] * sin(angle);
			rotation(0, 2) = axis[0] * axis[2] * (1 - cos(angle)) + axis[1] * sin(angle);

			rotation(1, 0) = axis[0] * axis[1] * (1 - cos(angle)) + axis[2] * sin(angle);
			rotation(1, 1) = cos(angle) + axis[1] * axis[1] * (1 - cos(angle));
			rotation(1, 2) = axis[2] * axis[1] * (1 - cos(angle)) - axis[0] * sin(angle);

			rotation(2, 0) = axis[0] * axis[1] * (1 - cos(angle)) - axis[1] * sin(angle);
			rotation(2, 1) = axis[2] * axis[1] * (1 - cos(angle)) + axis[0] * sin(angle);
			rotation(2, 2) = cos(angle) + axis[2] * axis[2] * (1 - cos(angle));

			return rotation;
		};

		template<typename T, typename H>
		void CTetHexGenerator<T, H>::__rotateToZeroDirection()
		{
			std::cout << "[GenHex]Rotate to zero direction..." << std::endl;
			CPoint2 y_axis(0, 1);
			CPoint2 zero2D(m_zeroDirection[0], m_zeroDirection[1]);
			double angle = acos((y_axis * zero2D) / (y_axis.norm() * zero2D.norm()));
			angle = angle;
			CPoint rotationAxis(0, 0, 1);
			rotationAxis = rotationAxis / rotationAxis.norm();
			Eigen::Matrix3d rotationMatrix = ___convertToRotationMatrix(rotationAxis, angle);

			for (H::MeshVertexIterator vIter(m_pHexMesh); !vIter.end(); vIter++)
			{
				H::CVertex * pV = *vIter;
				CPoint p = pV->position();
				p = p - CPoint(0.5, 0.5, 0);
				Eigen::Vector3d pos;
				for (int k = 0; k < 3; k++)
				{
					pos(k) = p[k];
				}

				Eigen::Vector3d newPos = rotationMatrix * pos;
				CPoint newP;
				for (int k = 0; k < 3; k++)
				{
					newP[k] = newPos(k);
				}
				pV->position() = newP + CPoint(0.5, 0.5, 0);
			}

		};

		template<typename T, typename H>
		void CTetHexGenerator<T, H>::_genHex()
		{
			
			__genStdHexMesh();
			//m_pHexMesh->_write_hm("standardHexTest.hm");
			if (m_zero)
			{
				m_zeroDirection = __zeroPointDirection();
				__rotateToZeroDirection();
			}
			m_pHexMesh->_write_hm("standardHexRotateTest.hm");
			__locateVertices();
			__locateSurfaceViaNormal();
		};

		template<typename T, typename H>
		void CTetHexGenerator<T, H>::_writeHex(const char * output)
		{
			m_pHexMesh->_write_hm(output);
		}
	}
}

#endif