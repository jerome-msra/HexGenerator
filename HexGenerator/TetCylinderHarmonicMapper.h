#ifndef _TET_CYLINDER_HARMONIC_MAPPER_H_
#define _TET_CYLINDER_HARMONIC_MAPPER_H_

#include "Operator\TOperator.h"
#include "Laplace\TLaplace.h"

#include <ppl.h>

using namespace concurrency;

namespace MeshLib
{
	namespace TMeshLib
	{
		template<typename M, typename T>
		class CTetCylinderHarmonicMapper
		{
		public:
			CTetCylinderHarmonicMapper(M * pSurface, T * pTet);
			~CTetCylinderHarmonicMapper();

			void _map();

			void _iterativeMap(double stepLength = 0.5, double threshold = 1e-6);

			void _write_t(const char * output);

		protected:

			void __setBoundary();

			void __updateHuv(typename T::CVertex * pV, const double stepLength);

			CPoint __computeLaplace(typename T::CVertex * pV);

			double __computeHarmonicEnergy();

		private:
			M * m_pSurface;
			T * m_pTet;
		};

		template<typename M, typename T>
		CTetCylinderHarmonicMapper<M, T>::CTetCylinderHarmonicMapper(M * pSurface, T * pTet)
		{
			m_pSurface = pSurface;
			m_pTet = pTet;

			m_pTet->_labelBoundaryVertices();

			m_pTet->_computeHalfFaceNormal();

			for (T::MeshVertexIterator vIter(m_pTet); !vIter.end(); vIter++)
			{
				T::CVertex * pV = *vIter;
				pV->fixed() = pV->boundary();
			}
		};

		template<typename M, typename T>
		CTetCylinderHarmonicMapper<M, T>::~CTetCylinderHarmonicMapper()
		{
		};

		template<typename M, typename T>
		void CTetCylinderHarmonicMapper<M, T>::__setBoundary()
		{
			for (M::MeshVertexIterator vIter(m_pSurface); !vIter.end(); vIter++)
			{
				M::CVertex * pV = *vIter;
				T::CVertex * pTetVert = m_pTet->idVertex(pV->father());
				pTetVert->huv() = pV->point();
			}
		};

		template<typename M, typename T>
		void CTetCylinderHarmonicMapper<M, T>::_map()
		{
			__setBoundary();
			CTOperator<T> * op = new CTOperator<T>(m_pTet);
			op->_embedding_2_metric();
			op->_metric_2_laplace();

			for (int k = 0; k < 3; k++)
			{
				CTLaplace<T> L(m_pTet);

				for (T::MeshVertexIterator vIter(m_pTet); !vIter.end(); vIter++)
				{
					T::CVertex * pV = *vIter;
					if (!pV->fixed())
					{
						continue;
					}
					pV->u() = pV->huv()[k];
				}

				L.solve();

				for (T::MeshVertexIterator vIter(m_pTet); !vIter.end(); vIter++)
				{
					T::CVertex * pV = *vIter;
					if (pV->fixed())
					{
						continue;
					}
					pV->huv()[k] = pV->u();
				}
			}

			for (T::MeshVertexIterator vIter(m_pTet); !vIter.end(); vIter++)
			{
				T::CVertex * pV = *vIter;
				pV->coord3d() = pV->position();
				pV->position() = pV->huv();
			}
		};

		template<typename M, typename T>
		inline double CTetCylinderHarmonicMapper<M, T>::__computeHarmonicEnergy()
		{
			double harmonic_energy = 0.;

			//serial algorithm
			for (T::MeshEdgeIterator eiter(m_pTet); !eiter.end(); ++eiter)
			{
				T::CEdge * pE = *eiter;
				T::CVertex * pV = m_pTet->EdgeVertex1(pE);
				T::CVertex * pW = m_pTet->EdgeVertex2(pE);
				harmonic_energy += pE->weight() * ((pW->huv() - pV->huv()) * (pW->huv() - pV->huv()));
			}

			return harmonic_energy;
		};


		template<typename M, typename T>
		inline CPoint CTetCylinderHarmonicMapper<M, T>::__computeLaplace(typename T::CVertex * pV)
		{
			CPoint lap(0, 0, 0);
			for (T::VertexEdgeIterator veiter(m_pTet, pV); !veiter.end(); ++veiter)
			{
				T::CEdge * pE = *veiter;
				T::CVertex * pW = m_pTet->EdgeVertex1(pE) == pV ? m_pTet->EdgeVertex2(pE) : m_pTet->EdgeVertex1(pE);
				assert(pW != NULL && pW != pV);
				lap += (pW->huv() - pV->huv()) * pE->weight();
			}

			return lap;
		};

		template<typename M, typename T>
		inline void CTetCylinderHarmonicMapper<M, T>::__updateHuv(typename T::CVertex * pV, const double stepLength)
		{
			if (!pV->boundary())
			{
				CPoint & lap_h = __computeLaplace(pV);

				CPoint & normal_component = pV->huv() * (lap_h * pV->huv());
				CPoint & tangent_component = lap_h - normal_component;

				pV->huv() += tangent_component * stepLength;
			}
		};

		template<typename M, typename T>
		void CTetCylinderHarmonicMapper<M, T>::_iterativeMap(double stepLength, double threshold)
		{
			__setBoundary();

			CTOperator<T> * op = new CTOperator<T>(m_pTet);
			op->_embedding_2_metric();
			op->_metric_2_laplace();
			
			double energy[2];
			energy[1] = __computeHarmonicEnergy();
			size_t iter_count = 0;
			do
			{
				//parallel algorithm
				std::list<T::CVertex *> & vertices = m_pTet->vertices();
				parallel_for_each(std::begin(vertices), std::end(vertices), [&](T::CVertex * pV) {
					__updateHuv(pV, stepLength);
				});

				energy[0] = energy[1];
				energy[1] = __computeHarmonicEnergy();
				if (iter_count % 100 == 0)
					std::cout << iter_count << "\terror:" << fabs(energy[1] - energy[0]) << std::endl;
				++iter_count;

			} while (fabs(energy[1] - energy[0]) > threshold);
		};

		template<typename M, typename T>
		void CTetCylinderHarmonicMapper<M, T>::_write_t(const char * output)
		{
			m_pTet->_write_t(output);
		};
	}
};

#endif