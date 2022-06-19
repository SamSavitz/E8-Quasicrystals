#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>
//#include <execution>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

using namespace std;
using namespace chrono;

#include <tbb/parallel_for.h>

using Z = ptrdiff_t;
using R = float;

constexpr Z D = 8;
constexpr Z N = 240;

constexpr R SKEW = 15;
constexpr R PAUSE = 1.15;

constexpr R minE = -N/SKEW;
constexpr R maxE =  N     ;

constexpr R KNEEM = 4;
constexpr R KNEEP = 15.5;

constexpr R TANHM = 0.5;
constexpr R TANHP = 1.2;

#include "codec.h"
#include "color.h"

using ZV = array <Z, D>;
using RV = array <R, D>;
using RVP = array <RV, 2>;

#include "data.h"

constexpr Z W = 1920/2;
constexpr Z H = 1920/2;
constexpr Z P = H*W;

constexpr R RES1 = 24;
constexpr R RES2 = 7;

constexpr R XO = W/2;
constexpr R YO = H/2;

constexpr R FPS  = 30;

constexpr R BT = 0.3 ;
constexpr R ST = 1   ;
constexpr R ZT = 10.5;
constexpr R TT = 0.75;
constexpr R RT = 12  ;

constexpr R SECS = ST + ZT + TT + RT;
constexpr Z FS = FPS*SECS;

constexpr R EPSILON = 1e-9;

constexpr Z SEED = 4;

#include "vectors.h"

int main(int argc, const char* argv[]) {
    RGBEncoder encoder("out.mp4", W, H, FPS, 1);

    preprocess();

    stringstream ss;

	high_resolution_clock::time_point ot;
	Z oy;

	Z yp = -1;

	   RGB* const rgbFrame = new RGB   [P  ];
    //YUV420* const yuvFrame = new YUV420[P/4];

    for (Z n = 0; n < XY.size(); ++n) {
        const Z sf = n ? 0 : -FPS*BT;
        const Z ef = FPS*(n == XY.size() - 1 ? ST + ZT + TT + BT : SECS);

        for (Z f = sf; f < ef; ++f) {
            const Z ff = FS*n + f + FPS*BT;

            cerr << n << ' ' << f << ' ' << ff << '/' << ef << endl;

            cout << ff << '\t' << n << '\t' << f << '\t';

            RVP xy = interpolate(
                XY[n],
                XY[n == XY.size() - 1 ? n : n + 1],
                   n == XY.size() - 1 && f >= FPS*(ST + ZT + TT) ? FPS*(ST + ZT + TT) : f
            );

            tbb::parallel_for(
                //execution::par_unseq,
                tbb::blocked_range <int> (0, H),
                [&] (tbb::blocked_range <int> r) {

                	for (Z y = r.begin(); y < r.end(); ++y) {
                    //for (Z y = 0; y < H; ++y) {
                		const R yf = (y - YO);

                		for (Z x = 0; x < W; ++x) {
                			const R xf = (x - XO);

                            const RV z = mix(xf, xy[0], yf, xy[1] );

                            R sum = 0;

                			for (const ZV& r :  ROOTS)
                				sum += cos(dot(r, z));

                            sum *= 2;

                            rgbFrame[W*y + x] = color(sum);
                		}
                    //}
                    }
                }
            );

            //rgb2yuv420(W, H, rgbFrame, yuvFrame);

            encoder.writeFrame(rgbFrame);
        }
    }

    delete [] rgbFrame;//, yuvFrame;
}
