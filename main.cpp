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

constexpr R RES = 10;

constexpr R XO = H/2;
constexpr R YO = H/2;

constexpr Z FS = 60;

constexpr R EPSILON = 1e-9;

constexpr Z SEED = 1;


template <class V1, class V2>
inline R dot(const V1& v, const V2& w) {
    R s = 0;

    for (unsigned char i = 0; i < D; i++)
        s += v[i]*w[i];

    return s;
}

inline RV mix(const R a, const RV& v, const R b, const RV& w) {
    RV r;

    for (unsigned char i = 0; i < D; i++)
        r[i] = a*v[i] + b*w[i];

    return r;
}

inline void rotate(const R c, const R s, RVP& v) {
    const RV t    = mix( c, v[0], s, v[1] );
             v[1] = mix(-s, v[0], c, v[1] );
             v[0] = t;
}

inline void rotate(const R t, RVP& v) {
    rotate(cos(t), sin(t), v);
}

inline void rotate(RVP& v) {
    static default_random_engine gen(SEED);
    static uniform_real_distribution <R> dist(0.0, 2*M_PI);
    rotate(dist(gen), v);
}

inline void preprocess() {
    rotate(XY[0] );

    for (unsigned int n = 0; n < XY.size() - 1; n++) {
        RVP& a = XY[n    ];
        RVP& b = XY[n + 1];

        const R sxx = dot(a[0], b[0] );
        const R sxy = dot(a[0], b[1] );
        const R syx = dot(a[1], b[0] );
        const R syy = dot(a[1], b[1] );

        const R s = sqrt((sxx + syy)*(sxx + syy) + (sxy - syx)*(sxy - syx));

        if (s < EPSILON)
            rotate(b);
        else {
            const R ct = (sxx + syy)/s;
            const R st = (sxy - syx)/s;

            rotate(ct, st, b);
        }

    }
}

inline void normalize(RV& v) {
    R s = 0;

    for (unsigned char i = 0; i < D; i++)
        s += v[i]*v[i];

    s = 1/sqrt(s);

    for (unsigned char i = 0; i < D; i++)
        v[i] *= s;
}

int main(int argc, const char* argv[]) {
    Encoder encoder("out.mp4", W, H, 10, 0.5);

    preprocess();

    stringstream ss;

	high_resolution_clock::time_point ot;
	Z oy;

	Z yp = -1;
	RGB* const frame = new RGB[P];

    for (Z n = 0; n < XY.size(); ++n) {
    	const Z np = n + 1;

        for (Z f = 0; f < FS; ++f) {
            R t = (R) f/FS;
            const Z ff = FS*n + f;

            cerr << np << ' ' << f << ' ' << ff << endl;

            t = 2*(t - 0.5);
            t = (t*(t*t*(t*t*(21 - 5*t*t) - 35) + 35) + 16)/32;

            RV X2 = mix(1 - t, XY[n] [0], t, XY[np] [0] );
            RV Y2 = mix(1 - t, XY[n] [1], t, XY[np] [1] );

            if (f) {
                normalize(X2);
                normalize(Y2);
            }

            const R a = dot(X2, Y2);

            const RV X = mix((sqrt(1 - a) + sqrt(1 + a))/2/sqrt(1 - a*a), X2, (sqrt(1 - a) - sqrt(1 + a))/2/sqrt(1 - a*a), Y2);
            const RV Y = mix((sqrt(1 - a) - sqrt(1 + a))/2/sqrt(1 - a*a), X2, (sqrt(1 - a) + sqrt(1 + a))/2/sqrt(1 - a*a), Y2);

            tbb::parallel_for(
                //execution::par_unseq,
                tbb::blocked_range <int>  (0, H),
                [&] (tbb::blocked_range <int> r) {

                	for (Z y = r.begin(); y < r.end(); ++y) {
                		const R yf = (y - YO)/RES;

                		for (Z x = 0; x < W; ++x) {
                			const R xf = (x - XO)/RES;

                            const RV z = mix(xf, X, yf, Y);

                            R sum = 0;

                			for (const ZV& r :  ROOTS)
                				sum += cos(dot(r, z));

                            sum *= 2;

                            frame[W*y + x] = color(sum);
                		}
                	}
                }
            );

            encoder.writeFrame(frame);
        }
    }
}
