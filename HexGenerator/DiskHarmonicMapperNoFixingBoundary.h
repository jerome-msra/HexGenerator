#ifndef _DISK_HARMONIC_MAPPER_NO_FIXING_BOUNDARY_H_
#define _DISK_HARMONIC_MAPPER_NO_FIXING_BOUNDARY_H_

#include "conformal\DiskHarmonicMapper.h"

namespace MeshLib
{
	namespace Holomorphy
	{
		template<typename M>
		class CDiskHarmonicMapperNoFixingBoundary : public CDiskHarmonicMapper < M >
		{
		public:

			CDiskHarmonicMapperNoFixingBoundary(M * pMesh) : CDiskHarmonicMapper(pMesh) {};

			~CDiskHarmonicMapperNoFixingBoundary(){};

			void _mapNoFixingBoundary();
		};

		template <typename M>
		void CDiskHarmonicMapperNoFixingBoundary<M>::_mapNoFixingBoundary()
		{
			for (int k = 0; k < 2; k++)
			{
				CLaplace<M> L(m_pMesh);

				for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
				{
					M::CVertex * pV = *viter;
					if (!pV->fixed()) continue;
					pV->u() = pV->huv()[k];
				}

				L.solve();

				for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
				{
					M::CVertex * pV = *viter;
					if (pV->fixed()) continue;
					pV->huv()[k] = pV->u();
				}
			}
		};
	};
};

#endif