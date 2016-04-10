#ifndef _OPTIMIZE_HEX_H_
#define _OPTIMIZE_HEX_H_

#include "OptimizeHexMesh.h"

namespace MeshLib
{
	namespace HMeshLib
	{

		template<typename H>
		class OptimizeHex
		{
		public:
			OptimizeHex();
			~OptimizeHex();

			void _addHexMesh(H * hex)
			{
				hex->labelBoundaryVertices();
				m_hexMeshList.push_back(hex);
			};

			void _smoothSurface(int iterations);
			
			void _smoothInterior(int iterations);

			void _writeHex(const char * output);

		private:

			void __labelDiskBoundary();


		private:
			std::vector<H*> m_hexMeshList;
		};

		template<typename H>
		OptimizeHex<H>::OptimizeHex()
		{
		};

		template<typename H>
		OptimizeHex<H>::~OptimizeHex()
		{
		};

		template<typename H>
		void OptimizeHex<H>::_smoothSurface(int interations)
		{
			__labelDiskBoundary();

			for (size_t r = 0; r < interations; r++)
			{
				for (size_t i = 0; i < m_hexMeshList.size(); i++)
				{
					H * hmesh = m_hexMeshList[i];
					for (H::MeshVertexIterator vIter(hmesh); !vIter.end(); vIter++)
					{
						H::CVertex * pV = *vIter;
						if ((!pV->boundary()) || pV->disk() || pV->diskBoundary())
						{
							continue;
						}
						CPoint sumP;
						int sum = 0;
						for (H::VertexVertexIterator vvIter(hmesh, pV); !vvIter.end(); vvIter++)
						{
							H::CVertex * pVV = *vvIter;
							if (pVV->boundary())
							{
								sumP += pVV->position();
								sum++;
							}
						}
						pV->position() = sumP / sum;
					}
				}
			}
		};

		template<typename H>
		void OptimizeHex<H>::_smoothInterior(int interations)
		{
			for (size_t r = 0; r < interations; r++)
			{
				for (size_t i = 0; i < m_hexMeshList.size(); i++)
				{
					H * hmesh = m_hexMeshList[i];
					for (H::MeshVertexIterator vIter(hmesh); !vIter.end(); vIter++)
					{
						H::CVertex * pV = *vIter;
						if (pV->boundary())
						{
							continue;
						}
						CPoint sumP;
						int sum = 0;
						for (H::VertexVertexIterator vvIter(hmesh, pV); !vvIter.end(); vvIter++)
						{
							H::CVertex * pVV = *vvIter;
							if (pVV->boundary())
							{
								sumP += pVV->position();
								sum++;
							}
						}
						pV->position() = sumP / sum;
					}
				}
			}
		};

		template<typename H>
		void OptimizeHex<H>::__labelDiskBoundary()
		{
			for (size_t i = 0; i < m_hexMeshList.size(); i++)
			{
				H * hmesh0 = m_hexMeshList[i];
				for (size_t j = i + 1; j < m_hexMeshList.size(); j++)
				{
					H * hmesh1 = m_hexMeshList[j];

					for (H::MeshVertexIterator vIter0(hmesh0); !vIter0.end(); vIter0++)
					{
						H::CVertex * pV0 = *vIter0;
						CPoint p0 = pV0->position();
						for (H::MeshVertexIterator vIter1(hmesh1); !vIter1.end(); vIter1++)
						{
							H::CVertex * pV1 = *vIter1;
							CPoint p1 = pV1->position();
							if (p0 == p1)
							{
								pV0->diskBoundary() = true;
								//pV0->sameVertex().push_back(pV1);
								pV1->diskBoundary() = true;
								//pV1->sameVertex().push_back(pV0);
							}
						}
					}
				}
			}
		};

		template<typename H>
		void OptimizeHex<H>::_writeHex(const char * output)
		{
			std::string outputName(output);
			size_t index = outputName.find_last_of('.');
			std::string rawName = outputName.substr(0, index);

			for (size_t i = 0; i < m_hexMeshList.size(); i++)
			{
				std::string outputname = rawName + "_" + std::to_string(i) + ".hm";
				H * hmesh = m_hexMeshList[i];
				hmesh->_write_hm(outputname.c_str());
			}
		};


	}; // namespace HMeshLib
}; // namespace MeshLib


#endif