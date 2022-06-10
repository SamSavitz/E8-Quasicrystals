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

constexpr R minE = -N/SKEW;
constexpr R maxE =  N     ;

#include "codec.h"
#include "color.h"

using ZV = array <Z, D>;
using RV = array <R, D>;
using RVP = array <RV, 2>;

#include "data.h"

constexpr Z W = 1920;
constexpr Z H = 1080;
constexpr Z P = H*W;

constexpr R RES = 55;

constexpr R XO = W/2;
constexpr R YO = H/2;

constexpr Z FS = 360;

constexpr R EPSILON = 1e-9;

constexpr Z SEED = 2;

#include "vectors.h"

int main(int argc, const char* argv[]) {
    RGBEncoder encoder("out.mp4", W, H, 30, 1);

    preprocess();

    stringstream ss;

	high_resolution_clock::time_point ot;
	Z oy;

	Z yp = -1;

	   RGB* const rgbFrame = new RGB   [P  ];
    //YUV420* const yuvFrame = new YUV420[P/4];

    for (Z n = 0; n < XY.size(); ++n)
        for (Z f = 0; f < FS; ++f) {
            const Z ff = FS*n + f;

            cerr << n << ' ' << f << ' ' << ff << endl;

            RVP xy = interpolate(XY[n], XY[n + 1], f);

            tbb::parallel_for(
                //execution::par_unseq,
                tbb::blocked_range <int>  (0, H),
                [&] (tbb::blocked_range <int> r) {

                	for (Z y = r.begin(); y < r.end(); ++y) {
                		const R yf = (y - YO)/RES;

                		for (Z x = 0; x < W; ++x) {
                			const R xf = (x - XO)/RES;

                            const RV z = mix(xf, xy[0], yf, xy[1] );

                            R sum = 0;

                			for (const ZV& r :  ROOTS)
                				sum += cos(dot(r, z));

                            sum *= 2;

                            rgbFrame[W*y + x] = color(sum);
                		}
                	}
                }
            );

            //rgb2yuv420(W, H, rgbFrame, yuvFrame);

            encoder.writeFrame(rgbFrame);
        }

    delete [] rgbFrame;//, yuvFrame;
}
