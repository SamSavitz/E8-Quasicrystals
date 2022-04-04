#include <algorithm>
#include <array>
#include <cassert>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

using namespace std;
//using namespace chrono;

#include <tbb/parallel_for.h>

using Z = ptrdiff_t;
using R = float;

constexpr Z D = 24;
constexpr Z N = 196560;

constexpr R minE = -342.961657023412635055311941451;
constexpr R maxE = N;

constexpr R KNEEP = 69;//100;
constexpr R KNEEM = 2./3;

#include "codec.h"
#include "color.h"

using  ZV = array <Z , D>;
using  RV = array <R , D>;
using RVP = array <RV, 2>;

#include "leechData.h"

constexpr Z W = 1980;
constexpr Z H = 1080;
constexpr Z P = H*W;

constexpr R RES = 13;
constexpr const char* const FN = "medium5.mp4";

constexpr R XO = W/2;
constexpr R YO = H/2;

#include "leech.cu"

constexpr Z FS = 360;
constexpr R EPSILON = 1e-9;
constexpr Z SEED = 3;

ZV Vs[N/2];

#include "vectors.h"

void initLeech();

int main(int argc, const char* argv[]) {
    RGBEncoder encoder(FN, W, H, 30, 1);

    preprocess();
    initLeech();
    initCUDA();

    stringstream ss;

	//high_resolution_clock::time_point ot;
	//Z oy;

	//Z yp = -1;

    R* const data = new R[P];
	   RGB* const rgbFrame = new RGB   [P  ];
    //YUV420* const yuvFrame = new YUV420[P/4];

    for (Z n = 0; n < XY.size(); ++n)
        for (Z f = 0; f < FS; ++f) {
            const Z ff = FS*n + f;

            cerr << n << ' ' << f << ' ' << ff << endl;

            RVP xy = interpolate(XY[n], XY[n + 1], f);

            evaluateFrame(data, xy);

            for (int i = 0; i < P; i++)
                rgbFrame[i] = color(data[i] );

            /*
            tbb::parallel_for(
                //execution::par_unseq,
                tbb::blocked_range <int>  (0, H),
                [&] (tbb::blocked_range <int> r) {

            	for (Z y = r.begin(); y < r.end(); ++y) {
                    cout << y << endl;

            		const R yf = (y - YO)/RES;

            		for (Z x = 0; x < W; ++x) {
            			const R xf = (x - XO)/RES;

                        const RV z = mix(xf, xy[0], yf, xy[1] );

                        R sum = 0;

            			for (const ZV& r :  Vs)
            				sum += cosf(dot(r, z));

                        sum *= 2;

                        rgbFrame[W*y + x] = color(sum);
            		}
            	}}
            );
            */

            //rgb2yuv420(W, H, rgbFrame, yuvFrame);

            encoder.writeFrame(rgbFrame);
        }

    delete [] data, rgbFrame;//, yuvFrame;
}


Z golay[1 << 12];
Z ternary[759];

void initLeech() {
	Z count = 0, count1 = 0, count2 = 0;

	for (Z i = 0; i < 1 << 24; i++) {
		const bool parity1 = __builtin_popcount(i & 0x00000f) % 2;
		const bool parity2 = __builtin_popcount(i & 0x0000f0) % 2;
		const bool parity3 = __builtin_popcount(i & 0x000f00) % 2;
		const bool parity4 = __builtin_popcount(i & 0x00f000) % 2;
		const bool parity5 = __builtin_popcount(i & 0x0f0000) % 2;
		const bool parity6 = __builtin_popcount(i & 0xf00000) % 2;
		const bool parity7 = __builtin_popcount(i & 0x111111) % 2;

		constexpr Z W = 2;
		constexpr Z B = 3;

		constexpr Z WORDS[8] = {0, 1, W, B, B, W, 1, 0};
		Z words[6];

		words[0] = WORDS[(i & 0x00000e) >> 1];
		words[1] = WORDS[(i & 0x0000e0) >> 5];
		words[2] = WORDS[(i & 0x000e00) >> 9];
		words[3] = WORDS[(i & 0x00e000) >> 13];
		words[4] = WORDS[(i & 0x0e0000) >> 17];
		words[5] = WORDS[(i & 0xe00000) >> 21];

		constexpr Z HEXES[64] [6] = {
			{0, 0, 0, 0, 0, 0},


			{1, 1, W, W, B, B},
			{W, W, B, B, 1, 1},
			{B, B, 1, 1, W, W},

			{1, 1, B, B, W, W},
			{W, W, 1, 1, B, B},
			{B, B, W, W, 1, 1},


			{0, 0, 1, 1, 1, 1},
			{0, 0, W, W, W, W},
			{0, 0, B, B, B, B},

			{1, 1, 0, 0, 1, 1},
			{W, W, 0, 0, W, W},
			{B, B, 0, 0, B, B},

			{1, 1, 1, 1, 0, 0},
			{W, W, W, W, 0, 0},
			{B, B, B, B, 0, 0},


			{W, B, W, B, W, B},
			{B, 1, B, 1, B, 1},
			{1, W, 1, W, 1, W},

			{W, B, B, W, B, W},
			{B, 1, 1, B, 1, B},
			{1, W, W, 1, W, 1},

			{B, W, W, B, B, W},
			{1, B, B, 1, 1, B},
			{W, 1, 1, W, W, 1},

			{B, W, B, W, W, B},
			{1, B, 1, B, B, 1},
			{W, 1, W, 1, 1, W},



			{0, 1, 0, 1, W, B},
			{0, W, 0, W, B, 1},
			{0, B, 0, B, 1, W},

			{0, 1, W, B, 0, 1},
			{0, W, B, 1, 0, W},
			{0, B, 1, W, 0, B},

			{W, B, 0, 1, 0, 1},
			{B, 1, 0, W, 0, W},
			{1, W, 0, B, 0, B},


			{0, 1, 1, 0, B, W},
			{0, W, W, 0, 1, B},
			{0, B, B, 0, W, 1},

			{0, 1, B, W, 1, 0},
			{0, W, 1, B, W, 0},
			{0, B, W, 1, B, 0},

			{W, B, 1, 0, 1, 0},
			{B, 1, W, 0, W, 0},
			{1, W, B, 0, B, 0},


			{1, 0, 0, 1, B, W},
			{W, 0, 0, W, 1, B},
			{B, 0, 0, B, W, 1},

			{1, 0, W, B, 1, 0},
			{W, 0, B, 1, W, 0},
			{B, 0, 1, W, B, 0},

			{B, W, 0, 1, 1, 0},
			{1, B, 0, W, W, 0},
			{W, 1, 0, B, B, 0},


			{1, 0, 1, 0, W, B},
			{W, 0, W, 0, B, 1},
			{B, 0, B, 0, 1, W},

			{1, 0, B, W, 0, 1},
			{W, 0, 1, B, 0, W},
			{B, 0, W, 1, 0, B},

			{B, W, 1, 0, 0, 1},
			{1, B, W, 0, 0, W},
			{W, 1, B, 0, 0, B},
		};

		if (!(
				(parity1 && parity2 && parity3 && parity4 && parity5 && parity6 && parity7)
			||
				(!parity1 && !parity2 && !parity3 && !parity4 && !parity5 && !parity6 && !parity7)))

			goto fail;

		for (Z j = 0; j < 64; j++) {
			for (Z k = 0; k < 6; k++)
				if (words[k] != HEXES[j] [k])  goto loop;

			goto win;

			loop:;
		}

		goto fail;

		win:

		golay[count1++] = i;

		if (__builtin_popcount(i) == 8)
			ternary[count2++] = i;

    		fail:;
	}

	for (Z i = 0; i < 24; i++)
		for (Z j = 0; j < i; j++) {
			Vs[count] [j] = 4;
			Vs[count++] [i] = 4;

			Vs[count] [j] = 4;
			Vs[count++] [i] = -4;
		}

	for (Z i = 0; i < 759; i++) {
		for (Z signs = 0; signs < 256; signs += 2) {
			count2 = 0;

			if (!__builtin_parity(signs) % 2) {
				for (Z j = 0; j < 24; j++)
					if (ternary[i] & (1 << j))
						Vs[count] [j] = signs & (1 << count2++) ? -2 : 2;

				count++;
			}
		}
	}

	for (Z k = 0; k < 4096; k += 1) {
		for (Z i = 0; i < 24 && count != 98280; i++) {
			for (Z j = 0; j < 24; j++)
				Vs[count] [j] = 1;

			Vs[count] [i] = -3;

			for (Z l = 0; l < 24; l++)
				if (golay[k] & (1 << l))
					Vs[count] [l] = -Vs[count] [l];

			if (Vs[count] [0] > 0)
			   count++;
		}
	}
}
