#include <algorithm>
#include <iomanip>
#include <iostream>
#include <random>
using namespace std;

//constexpr unsigned int N = 196560;
//constexpr unsigned char D = 24;

constexpr unsigned short HG = 1 << 11;

#include "golay.h"
#include "signs.h"

__constant__ unsigned int  CGOLAY[HG    ];
__constant__ unsigned char CSIGNS[1 << 6];


__global__ void evaluate(R* const out, const RVP xy, const int frame) {
	const auto tid = threadIdx.x;
	const auto id = 512*blockIdx.x + tid;

    const auto y = id/W;
    const auto x = id % W;

    if (y >= H || x >= W)
        return;

    const R xf = (x - XO)/RES;
    const R yf = (y - YO)/RES;

	__shared__ unsigned int  golay[HG    ];
	__shared__ unsigned char signs[1 << 6];

	golay[tid        ] = CGOLAY[tid        ];
	golay[tid +   512] = CGOLAY[tid +   512];
	golay[tid + 2*512] = CGOLAY[tid + 2*512];
	golay[tid + 3*512] = CGOLAY[tid + 3*512];

	if (tid < 1 << 6)
		signs[tid] = CSIGNS[tid];

	__syncthreads();

	R coords[D];

	#pragma unroll
	for (unsigned char d = 0; d < D; d++)
		  coords[d] = xf*xy[0] [d] + yf*xy[1] [d];

	R s = 0;

	for (unsigned char i = 1; i < D; i++)
		for (unsigned char j = 0; j < D; j++)
			if (j < i) {
				s += cosf(4*(coords[i] + coords[j] ));
				s += cosf(4*(coords[i] - coords[j] ));
			}

	for (unsigned char sign:  signs)
		for (unsigned short i = 0; i < 759; i++) {
			unsigned int octad = golay[i];
			unsigned char l = __ffs(octad) - 1;
			octad -= 1 << l;

			R a = coords[l];

			#pragma unroll
			for (signed char d = 1; d < 8; d++) {
				l = __ffs(octad) - 1;
				octad -= 1 << l;

				a += (((sign >> d) & 1) ? -1 : 1)*coords[l];
			}

			s += cosf(2*a);
		}

	for (const unsigned int cw:  golay) {
		R a = 0;

		#pragma unroll
		for (unsigned char i = 0; i < D; i++)
			a += (((cw >> i) & 1) ? 1 : -1)*coords[i];

		#pragma unroll
		for (unsigned char t = 0; t < D; t++)
			s += cosf((a + (((cw >> t) & 1) ? -4 : 4)*coords[t] ));
	}

	out[id] = 2*s;

	//out[id] = (maxE*((2*x + 30*frame) % W))/W;
}

R* deviceOutput;

void initCUDA() {
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

	cudaMemcpyToSymbol(CGOLAY, GOLAY, sizeof(GOLAY));
	cudaMemcpyToSymbol(CSIGNS, SIGNS, sizeof(SIGNS));

	cudaMalloc(&deviceOutput, sizeof(R)*P);
}

void evaluateFrame(R* const data, const RVP& xy) {
	constexpr dim3 GRID_SIZE((P + 511)/512);
	constexpr dim3 BLOCK_SIZE(512);

	static int frame = 0;

	evaluate <<<GRID_SIZE, BLOCK_SIZE>>> (deviceOutput, xy, frame++);

     cudaError_t err = cudaGetLastError();

     if ( err != cudaSuccess )
     {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));       

        // Possibly: exit(-1) if program cannot continue....
     }

	cudaMemcpy(data, deviceOutput, sizeof(R)*P, cudaMemcpyDeviceToHost);

	cout << data[0] << endl;
}
