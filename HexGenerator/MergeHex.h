#ifndef _MERGE_HEX_H_
#define _MERGE_HEX_H_

#include "MergeHexMesh.h"

namespace MeshLib
{
	namespace HMeshLib
	{
		
		template<typename H>
		class MergeHex
		{
		public:
			MergeHex();
			~MergeHex();

			void _addMesh(H * hex);

			void _merge();

		private:
	
			void _mergeMesh(H * hex);

			typename H::CVertex * __findCorrespondDiskVertex(typename H::CVertex * vert);

		private:
			std::vector<H*> m_hexlist;
			H * m_hex;
		};

		template<typename H>
		MergeHex<H>::MergeHex()
		{
		};

		template<typename H>
		MergeHex<H>::~MergeHex()
		{
		};

		template<typename H>
		void MergeHex<H>::_addMesh(H * hex)
		{
			hex->labelBoundaryVertices();
			hex->labelBaseZeros();
			m_hexlist.push_back(hex);
		};

		template<typename H>
		typename H::CVertex * MergeHex<H>::__findCorrespondDiskVertex(typename H::CVertex * vert)
		{
			std::vector<H::CVertex*> candidates;

			for (H::MeshVertexIterator vIter(m_hex); !vIter.end(); vIter++)
			{
				H::CVertex * pV = *vIter;
				if (pV->onDisk())
				{
					if (abs(pV->zeroDistance()) == abs(vert->zeroDistance()) && pV->jcoord() == vert->jcoord())
					{
						candidates.push_back(pV);
					}
				}
			}

			double minDist = INT_MAX;
			H::CVertex * minDistVert = nullptr;
			CPoint refp = vert->position();
			for (std::vector<H::CVertex*>::iterator vIter = candidates.begin(); vIter != candidates.end(); vIter++)
			{
				H::CVertex * pV = *vIter;
				CPoint p = pV->position();
				if ((p - refp).norm() < minDist)
				{
					minDist = (p - refp).norm();
					minDistVert = pV;
				}
			}

			return minDistVert;
		};

		template<typename H>
		void MergeHex<H>::_mergeMesh(H * hex)
		{

			for (H::MeshVertexIterator vIter(hex); !vIter.end(); vIter++)
			{
				H::CVertex * pV = *vIter;
				if (!pV->onDisk())
				{
					continue;
				}

				H::CVertex * nearestVert = nullptr;

				nearestVert = __findCorrespondDiskVertex(pV);
				if (nearestVert != nullptr)
				{
					CPoint p1 = pV->position();
					CPoint p2 = nearestVert->position();
					pV->position() = (p1 + p2) / 2;
					nearestVert->position() = (p1 + p2) / 2;
				}
				else
				{
					std::cout << "[MergeHex]:Can't find nearest correspondent vertex " << pV->id() << "..." << std::endl;
				}
			}
		};

		template<typename H>
		void MergeHex<H>::_merge()
		{
			m_hex = m_hexlist[0];

			for (size_t i = 1; i < m_hexlist.size(); i++)
			{
				_mergeMesh(m_hexlist[i]);
			}

		};
	
	}; // namespace HMeshLib
}; // namespace MeshLib

#endif