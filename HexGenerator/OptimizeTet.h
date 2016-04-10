#ifndef _OPTIMIZE_Tet_H_
#define _OPTIMIZE_Tet_H_

#include "OptimizeTetMesh.h"

namespace MeshLib
{
	namespace TMeshLib
	{

		template<typename T>
		class OptimizeTet
		{
		public:
			OptimizeTet();
			~OptimizeTet();

			void _addTetMesh(T * Tet)
			{
				Tet->labelBoundaryVertices();
				m_TetMeshList.push_back(Tet);
			};

			void _smoothSurface(int iterations);

			void _smoothInterior(int iterations);

			void _writeTet(const char * output);

		private:

			void __labelDiskBoundary();


		private:
			std::vector<T*> m_TetMeshList;
		};

		template<typename T>
		OptimizeTet<T>::OptimizeTet()
		{
		};

		template<typename T>
		OptimizeTet<T>::~OptimizeTet()
		{
		};

		template<typename T>
		void OptimizeTet<T>::_smoothSurface(int interations)
		{
			__labelDiskBoundary();

			for (size_t r = 0; r < interations; r++)
			{
				for (size_t i = 0; i < m_TetMeshList.size(); i++)
				{
					T * tmesh = m_TetMeshList[i];
					for (T::MeshVertexIterator vIter(tmesh); !vIter.end(); vIter++)
					{
						T::CVertex * pV = *vIter;
						if ((!pV->boundary()) || pV->disk() || pV->diskBoundary())
						{
							continue;
						}
						CPoint sumP;
						int sum = 0;
						for (T::VertexVertexIterator vvIter(tmesh, pV); !vvIter.end(); vvIter++)
						{
							T::CVertex * pVV = *vvIter;
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

		template<typename T>
		void OptimizeTet<T>::_smoothInterior(int interations)
		{
			for (size_t r = 0; r < interations; r++)
			{
				for (size_t i = 0; i < m_TetMeshList.size(); i++)
				{
					T * tmesh = m_TetMeshList[i];
					for (T::MeshVertexIterator vIter(tmesh); !vIter.end(); vIter++)
					{
						T::CVertex * pV = *vIter;
						if (pV->boundary())
						{
							continue;
						}
						CPoint sumP;
						int sum = 0;
						for (T::VertexVertexIterator vvIter(tmesh, pV); !vvIter.end(); vvIter++)
						{
							T::CVertex * pVV = *vvIter;
							sumP += pVV->position();
							sum++;
						}
						pV->position() = sumP / sum;
					}
				}
			}
		};

		template<typename T>
		void OptimizeTet<T>::__labelDiskBoundary()
		{
			for (size_t i = 0; i < m_TetMeshList.size(); i++)
			{
				T * tmesh0 = m_TetMeshList[i];
				for (size_t j = i + 1; j < m_TetMeshList.size(); j++)
				{
					T * tmesh1 = m_TetMeshList[j];

					for (T::MeshVertexIterator vIter0(tmesh0); !vIter0.end(); vIter0++)
					{
						T::CVertex * pV0 = *vIter0;
						CPoint p0 = pV0->position();
						for (T::MeshVertexIterator vIter1(tmesh1); !vIter1.end(); vIter1++)
						{
							T::CVertex * pV1 = *vIter1;
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

		template<typename T>
		void OptimizeTet<T>::_writeTet(const char * output)
		{
			std::string outputName(output);
			size_t index = outputName.find_last_of('.');
			std::string rawName = outputName.substr(0, index);

			for (size_t i = 0; i < m_TetMeshList.size(); i++)
			{
				std::string outputname = rawName + "_" + std::to_string(i) + ".t";
				T * tmesh = m_TetMeshList[i];
				tmesh->_write_t(outputname.c_str());
			}
		};


	}; // namespace tmeshLib
}; // namespace MeshLib


#endif