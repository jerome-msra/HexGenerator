/*!
*      \file SurfaceCylinderMapper.h
*      \brief Mapper to map the cylinder to standard cylinder. Determine the boundary conditions for volume harmonic map.
*
*	   \author Jerome Jiang
*      \date 12/14/2015
*	   \modified 12/16/2015
*
*/


#ifndef _SURFACE_CYLINDER_MAPPER_H_
#define _SURFACE_CYLINDER_MAPPER_H_

#include "TetSurfaceMesh.h"
#include "HarmonicMapTMesh.h"

#include "DiskHarmonicMapperNoFixingBoundary.h"

#include "conformal\DiskHarmonicMapper.h"

namespace MeshLib
{
	template<typename M, typename T>
	class CSurfaceCylinderMapper
	{
	public:
		CSurfaceCylinderMapper(M * pMesh);
		~CSurfaceCylinderMapper();

		void _connectTet(T * tetMesh);

		void _map(M * domain, double height = 1.0);

		void _write_m(const char * output);

	protected:

		void __mapDisk();

		void __mapCylinder();

		void __mergeDiskCylinder();

		void __flipFace();

	private:
		// input meshes
		M * m_pMesh;
		M * m_pDomain;
		T * m_tetMesh;

		// result meshes
		/*! \brief mapped to disk meshes */
		std::vector<M*> m_diskMeshes;
		/*! \brief cylinder mesh */
		M * m_cylinderMesh;
		/*! \brief The closed cylinder mesh with disk and cylinder*/
		M * m_closedCylinderMesh;

		double m_height;
	};

	template<typename M, typename T>
	CSurfaceCylinderMapper<M, T>::CSurfaceCylinderMapper(M * pMesh)
	{
		m_pMesh = pMesh;
		int idx = 1;
		for (M::MeshVertexIterator vIter(m_pMesh); !vIter.end(); vIter++)
		{
			M::CVertex * pV = *vIter;
			pV->idx() = idx;
			idx++;
		}
	};

	template<typename M, typename T>
	CSurfaceCylinderMapper<M, T>::~CSurfaceCylinderMapper()
	{
	};

	template<typename M, typename T>
	void CSurfaceCylinderMapper<M, T>::_connectTet(T * tetMesh)
	{
		m_tetMesh = tetMesh;
		for (M::MeshVertexIterator vIter(m_pMesh); !vIter.end(); vIter++)
		{
			M::CVertex * pV = *vIter;
			CPoint pt = pV->point();
			bool found = false;
			double minDist = 1e+30;
			T::CVertex * minPV = nullptr;

			for (T::MeshVertexIterator tvIter(m_tetMesh); !tvIter.end(); tvIter++)
			{
				T::CVertex * ptV = *tvIter;
				CPoint ptvp = ptV->position();
				if (minDist > (pt - ptvp).norm())
				{
					minDist = (pt - ptvp).norm();
					minPV = ptV;
				}
				if (pt[0] == ptvp[0] && pt[1] == ptvp[1] && pt[2] == ptvp[2])
				{
					pV->father() = ptV->id();
					found = true;
					break;
				}
			}

			if (!found && minDist < 1e-6)
			{
				pV->father() = minPV->id();
			}
		}
	};

	template<typename M, typename T>
	void CSurfaceCylinderMapper<M, T>::__flipFace()
	{
		m_closedCylinderMesh = new M();

		for (M::MeshVertexIterator vIter(m_cylinderMesh); !vIter.end(); vIter++)
		{
			M::CVertex * pV = *vIter;
			M::CVertex * newVert = m_closedCylinderMesh->createVertex(pV->id());
			*newVert = (*pV);
		}

		for (M::MeshFaceIterator fIter(m_cylinderMesh); !fIter.end(); fIter++)
		{
			M::CFace * pF = *fIter;
			std::vector<M::CVertex *> vertList;
			for (M::FaceVertexIterator fvIter(pF); !fvIter.end(); fvIter++)
			{
				M::CVertex * pV = *fvIter;
				M::CVertex * pCV = m_closedCylinderMesh->idVertex(pV->id());
				vertList.push_back(pCV);
			}
			std::reverse(vertList.begin(), vertList.end());
			m_closedCylinderMesh->createFace(vertList, pF->id());
		}

		m_closedCylinderMesh->write_m("flippedCylinder.m");
	};

	template<typename M, typename T>
	void CSurfaceCylinderMapper<M, T>::__mapDisk()
	{
		for (int g = 1; g < 3; g++)
		{
			M * diskMesh = new M();

			for (M::MeshVertexIterator vIter(m_pMesh); !vIter.end(); vIter++)
			{
				M::CVertex * pV = *vIter;
				if (pV->group() == g)
				{
					pV->z() = g - 1;

					M::CVertex * diskVertex = diskMesh->createVertex(pV->id());
					diskVertex->point() = pV->point();
					diskVertex->group() = pV->group();
					diskVertex->father() = pV->father();
					diskVertex->rgb() = pV->rgb();
				}
			}

			for (M::MeshFaceIterator fIter(m_pMesh); !fIter.end(); fIter++)
			{
				M::CFace * pF = *fIter;
				std::vector<M::CVertex*> vertexList;
				for (M::FaceVertexIterator fvIter(pF); !fvIter.end(); fvIter++)
				{
					M::CVertex * pV = *fvIter;
					if (pV->group() == g)
					{
						vertexList.push_back(diskMesh->idVertex(pV->id()));
					}
				}
				if (vertexList.size() == 3)
				{
					diskMesh->createFace(vertexList, pF->id());
				}
			}

			diskMesh->labelBoundary();

			Holomorphy::CDiskHarmonicMapper<M> * diskMapper = new Holomorphy::CDiskHarmonicMapper<M>(diskMesh);

			diskMapper->_map();

			// label boundary vertices on the original cylinder mesh
			M::CBoundary * boundary = new M::CBoundary(diskMesh);
			std::vector<M::CLoop*> boundaryLoops = boundary->loops();
			M::CLoop * loop = boundaryLoops[0];
			std::list<M::CHalfEdge*> halfedges = loop->halfedges();
			for (std::list<M::CHalfEdge*>::iterator hIter = halfedges.begin(); hIter != halfedges.end(); hIter++)
			{
				M::CHalfEdge * pH = *hIter;
				M::CVertex * pV = diskMesh->halfedgeTarget(pH);
				M::CVertex * orinVert = m_pMesh->idVertex(pV->id());
				orinVert->diskBoundary() = true;
				orinVert->huv() = pV->huv();
			}


			for (M::MeshVertexIterator vIter(diskMesh); !vIter.end(); vIter++)
			{
				M::CVertex * pV = *vIter;
				CPoint2 phuv = pV->huv();
				pV->point()[0] = phuv[0];
				pV->point()[1] = phuv[1];
				pV->point()[2] = m_height * (g - 1);
			}

			std::string debugFile = "disk";
			std::stringstream iss("");
			iss << g << ".m";

			debugFile += iss.str();

			diskMesh->write_m(debugFile.c_str());

			m_diskMeshes.push_back(diskMesh);
		}
	};

	template<typename M, typename T>
	void CSurfaceCylinderMapper<M, T>::__mapCylinder()
	{
		// dijkstra shortest path
		std::vector<std::vector<M::CVertex*>> diskBoundaryVertices;
		diskBoundaryVertices.resize(2);

		for (M::MeshVertexIterator vIter(m_pMesh); !vIter.end(); vIter++)
		{
			M::CVertex * pV = *vIter;
			if (pV->diskBoundary())
			{
				diskBoundaryVertices[pV->group() - 1].push_back(pV);
			}
		}

		M::CBoundary * domainBoundary = new M::CBoundary(m_pDomain);
		std::vector<M::CLoop*> domainLoops = domainBoundary->loops();

		std::vector<std::vector<M::CHalfEdge*>> diskBoundaries;
		diskBoundaries.resize(2);
		for (std::vector<M::CLoop*>::iterator lIter = domainLoops.begin(); lIter != domainLoops.end(); lIter++)
		{
			M::CLoop * loop = *lIter;
			std::list<M::CHalfEdge*> halfedges = loop->halfedges();
			for (std::list<M::CHalfEdge*>::iterator hIter = halfedges.begin(); hIter != halfedges.end(); hIter++)
			{
				M::CHalfEdge * pH = *hIter;
				M::CVertex * pVS = m_pDomain->halfedgeSource(pH);
				M::CVertex * pVT = m_pDomain->halfedgeTarget(pH);
				M::CVertex * originVertS = m_pMesh->findVertex(pVS->point());
				M::CVertex * originVertT = m_pMesh->findVertex(pVT->point());
				if (originVertS != nullptr)
				{
					int g = originVertS->group();
					diskBoundaries[g - 1].push_back(pH);
				}
				else if (originVertT != nullptr)
				{
					int g = originVertT->group();
					diskBoundaries[g - 1].push_back(pH);
				}
			}
		}

		int numvertMesh = m_pMesh->numVertices();
		// group 0
		for (M::MeshVertexIterator vIter(m_pMesh); !vIter.end(); vIter++)
		{
			M::CVertex * pV = *vIter;
			if (pV->group() != 0 && (!pV->diskBoundary()))
			{
				continue;
			}

			M::CVertex * domainVert = m_pDomain->findVertex(pV->point());
			CPoint2 uv = domainVert->uv();

			assert(domainVert != nullptr);

			std::vector<CPoint2> huvCandidates;
			std::vector<double> domainDists;
			huvCandidates.resize(2);
			domainDists.resize(2);

			for (int g = 1; g < 3; g++)
			{
				std::vector<M::CHalfEdge*> halfedges = diskBoundaries[g - 1];
				bool found = false;
				for (std::vector<M::CHalfEdge*>::iterator hIter = halfedges.begin(); hIter != halfedges.end(); hIter++)
				{
					M::CHalfEdge * pH = *hIter;
					M::CVertex * pV0 = m_pDomain->halfedgeSource(pH);
					M::CVertex * pV1 = m_pDomain->halfedgeTarget(pH);
					CPoint2 uv0 = pV0->uv();
					CPoint2 uv1 = pV1->uv();

					if ((uv0[0] >= uv[0] && uv1[0] <= uv[0]) || (uv0[0] <= uv[0] && uv1[0] >= uv[0]))
					{
						M::CVertex * originV0 = m_pMesh->findVertex(pV0->point());
						M::CVertex * originV1 = m_pMesh->findVertex(pV1->point());
						found = true;
						if (originV0 != nullptr && originV1 != nullptr)
						{
							assert(originV0->group() == originV1->group());

							CPoint2 huv0 = originV0->huv();
							CPoint2 huv1 = originV1->huv();

							double lambda = (uv[0] - uv0[0]) / (uv1[0] - uv0[0]);
							CPoint2  huv = huv0 + (huv1 - huv0) * lambda;

							huvCandidates[originV0->group() - 1] = huv;
							domainDists[originV0->group() - 1] = fabs(uv[1] - (uv0[1] + uv1[1]) / 2.0);
						}
						else if (originV0 != nullptr)
						{
							CPoint2 huv0 = originV0->huv();
							huvCandidates[originV0->group() - 1] = huv0;
							domainDists[originV0->group() - 1] = fabs(uv[1] - (uv0[1] + uv1[1]) / 2.0);
						}
						else if (originV1 != nullptr)
						{
							CPoint2 huv1 = originV1->huv();
							huvCandidates[originV1->group() - 1] = huv1;
							domainDists[originV1->group() - 1] = fabs(uv[1] - (uv0[1] + uv1[1]) / 2.0);
						}

						break;
					}
				}

				// find the nearest
				// The u coordinates are consistent in the middle. just use the largest or smallest
				if (!found)
				{
					std::map<double, M::CVertex*> uvMap;
					for (std::vector<M::CHalfEdge*>::iterator hIter = halfedges.begin(); hIter != halfedges.end(); hIter++)
					{
						M::CHalfEdge * pH = *hIter;
						M::CVertex * pTV = m_pDomain->halfedgeTarget(pH);
						CPoint2 uv = pTV->uv();
						if (m_pMesh->findVertex(pTV->point()) != nullptr)
						{
							uvMap[uv[0]] = pTV;
						}
					}

					std::map<double, M::CVertex*>::iterator first = uvMap.begin();
					std::map<double, M::CVertex*>::iterator last = uvMap.end();
					last--;

					M::CVertex * firstVert = first->second;
					M::CVertex * endVert = last->second;

					M::CVertex * vert = (fabs(uv[0] - firstVert->uv()[0]) < fabs(uv[0] - endVert->uv()[0])) ? firstVert : endVert;
					M::CVertex * originVert = m_pMesh->findVertex(vert->point());
					huvCandidates[originVert->group() - 1] = originVert->huv();
					domainDists[originVert->group() - 1] = fabs(uv[1] - vert->uv()[1]);
				}
			}

			double z = domainDists[0] / (domainDists[0] + domainDists[1]);
			CPoint2 xy = huvCandidates[0];
			pV->huv() = xy;
			pV->z() = z * m_height;
		}

		m_cylinderMesh = new M();

		for (M::MeshVertexIterator vIter(m_pMesh); !vIter.end(); vIter++)
		{
			M::CVertex * pV = *vIter;
			if (pV->group() == 0 || pV->diskBoundary())
			{
				M::CVertex * cylinderVert = m_cylinderMesh->createVertex(pV->id());
				cylinderVert->point() = pV->point();
				cylinderVert->group() = pV->group();
				cylinderVert->father() = pV->father();
				cylinderVert->rgb() = pV->rgb();
				cylinderVert->huv() = pV->huv();
				cylinderVert->z() = pV->z();
			}
		}

		for (M::MeshFaceIterator fIter(m_pMesh); !fIter.end(); fIter++)
		{
			M::CFace * pF = *fIter;
			std::vector<M::CVertex*> vertexList;
			for (M::FaceVertexIterator fvIter(pF); !fvIter.end(); fvIter++)
			{
				M::CVertex * pV = *fvIter;
				if (pV->group() == 0 || pV->diskBoundary())
				{
					vertexList.push_back(m_cylinderMesh->idVertex(pV->id()));
				}
			}
			if (vertexList.size() == 3)
			{
				m_cylinderMesh->createFace(vertexList, pF->id());
			}
		}

		for (M::MeshVertexIterator vIter(m_cylinderMesh); !vIter.end(); vIter++)
		{
			M::CVertex * pV = *vIter;
			pV->point()[0] = pV->huv()[0];
			pV->point()[1] = pV->huv()[1];
			pV->point()[2] = pV->z();
		}

		m_cylinderMesh->labelBoundary();

		m_cylinderMesh->write_m("cylinderTest.m");
	};

	template<typename M, typename T>
	void CSurfaceCylinderMapper<M, T>::__mergeDiskCylinder()
	{

		// disk 2 seems to need to be remapped to disk, using boundary conditions from the cylinder
		M * disk2 = m_diskMeshes[1];
		for (M::MeshVertexIterator vIter(disk2); !vIter.end(); vIter++)
		{
			M::CVertex * pV = *vIter;
			if (pV->boundary())
			{
				M::CVertex * cpV = m_cylinderMesh->idVertex(pV->id());
				assert(cpV != nullptr);
				pV->huv()[0] = cpV->point()[0];
				pV->huv()[1] = cpV->point()[1];
				pV->fixed() = true;
			}
		}

		Holomorphy::CDiskHarmonicMapperNoFixingBoundary<M> * diskMapperNoFixing = new Holomorphy::CDiskHarmonicMapperNoFixingBoundary<M>(disk2);

		diskMapperNoFixing->_mapNoFixingBoundary();

		for (M::MeshVertexIterator vIter(disk2); !vIter.end(); vIter++)
		{
			M::CVertex * pV = *vIter;
			CPoint2 phuv = pV->huv();
			pV->point()[0] = phuv[0];
			pV->point()[1] = phuv[1];
		}

		disk2->write_m("disktest2.m");

		for (size_t i = 0; i < m_diskMeshes.size(); i++)
		{
			M * diskMesh = m_diskMeshes[i];
			for (M::MeshVertexIterator vIter(diskMesh); !vIter.end(); vIter++)
			{
				M::CVertex * pV = *vIter;

				bool exist = m_cylinderMesh->existVertex(pV->id());
				if (!exist)
				{
					M::CVertex * newCVert = m_cylinderMesh->createVertex(pV->id());
					newCVert->point() = pV->point();
					newCVert->group() = pV->group();
					newCVert->father() = pV->father();
					newCVert->rgb() = pV->rgb();
				}
			}

			for (M::MeshFaceIterator fIter(diskMesh); !fIter.end(); fIter++)
			{
				M::CFace * pF = *fIter;
				std::vector<M::CVertex *> vertList;
				for (M::FaceVertexIterator fvIter(pF); !fvIter.end(); fvIter++)
				{
					M::CVertex * pV = *fvIter;

					M::CVertex * pCVert = m_cylinderMesh->idVertex(pV->id());
					if (pCVert == nullptr)
					{
						std::cout << pV->id() << std::endl;
					}
					vertList.push_back(pCVert);
				}

				m_cylinderMesh->createFace(vertList, pF->id());
			}

			std::string debugFile = "closeCylinder";
			std::stringstream iss("");
			iss << i << ".m";

			debugFile += iss.str();

			m_cylinderMesh->write_m(debugFile.c_str());
		}
		
	};


	template<typename M, typename T>
	void CSurfaceCylinderMapper<M, T>::_map(M * domain, double height)
	{
		m_height = height;

		// step 1 map two disk (group 1 and 2)
		std::cout << "[map_surface]Mapping Disks..." << std::endl;
		__mapDisk();

		m_pDomain = domain;
		std::cout << "[map_surface]Mapping Cylinder..." << std::endl;
		__mapCylinder();

		std::cout << "[map_surface]Mergeing Cylinder with Disks..." << std::endl;
		__mergeDiskCylinder();

		std::cout << "[map_surface]Flipping Face Orientation..." << std::endl;
		__flipFace();
	};

	template<typename M, typename T>
	void CSurfaceCylinderMapper<M, T>::_write_m(const char * output)
	{
		m_closedCylinderMesh->write_m(output);
	}

};

#endif