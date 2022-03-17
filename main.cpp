#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

using namespace std;
using namespace chrono;

#include "BMP/EasyBMP.h"

using Z = ptrdiff_t;
using R = double;

constexpr Z D = 8;
constexpr Z N = 240;

using ZV = array <Z, D>;
using RV = array <R, D>;
using RVP = array <RV, 2>;

#include "data.h"

constexpr R SKEW = 15;

constexpr Z W = 1920;
constexpr Z H = 1200;
constexpr Z P = H*W;

constexpr R KNEEP = 4;
constexpr R KNEEM = 1;

constexpr R RES = 10;

constexpr R XO = H/2;
constexpr R YO = H/2;

constexpr R RB = 72187./212655;
constexpr R G = 250000./357579;
constexpr R JUMP = 0.2;

constexpr Z FS = 360;

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

inline unsigned char GAMMA(const R i) {
	R o;

	if (i <= 0.0031308)
		o = 12.92*i;
	else
		o = 1.055*pow(i, 5./12) - 0.055;

	o *= 256;

	if (Z(o) == 256)
		return 255;
	else
		return o;
}


int main(int argc, const char* argv[]) {
	assert(argc == 3);

	const Z n = atoi(argv[1] );
	Z np = n + 1;

    const Z f = atoi(argv[2] );
    R t = (R) f/FS;
    const Z ff = FS*n + f;

    cout << np << ' ' << f << ' ' << ff << endl;

    //if (t < 0.02 || t > 0.98)
    //    t = round(t);
    //else {
    //    t = (t - 0.02)/0.96;
        t = 2*(t - 0.5);
        t = (t*(t*t*(t*t*(21 - 5*t*t) - 35) + 35) + 16)/32;
    //}

    preprocess();

    RV X2 = mix(1 - t, XY[n] [0], t, XY[n + 1] [0] );
    RV Y2 = mix(1 - t, XY[n] [1], t, XY[n + 1] [1] );

    if (f) {
        normalize(X2);
        normalize(Y2);
    }

    const R a = dot(X2, Y2);

    const RV X = mix((sqrt(1 - a) + sqrt(1 + a))/2/sqrt(1 - a*a), X2, (sqrt(1 - a) - sqrt(1 + a))/2/sqrt(1 - a*a), Y2);
    const RV Y = mix((sqrt(1 - a) - sqrt(1 + a))/2/sqrt(1 - a*a), X2, (sqrt(1 - a) + sqrt(1 + a))/2/sqrt(1 - a*a), Y2);

    stringstream ss;

	high_resolution_clock::time_point ot;
	Z oy;

	Z yp = -1;
	R* vs = new R[P];

	for (Z y = 0; y < H; y++) {
		const Z nyp = round(100.*y/H);

		/*if (nyp > yp) {
			yp = nyp;

			ss.str("");

			ss << "C,\t";

			if (np < 10)
				ss << ' ';

			ss << np << ",\t";

			if (yp < 100) {
				ss << ' ';

				if (yp < 10)
					ss << ' ';
			}

			ss << yp << '%';

			const high_resolution_clock::time_point nt = high_resolution_clock::now();

			if (y) {
				const duration <R, nano> ts = nt - ot;
				const Z dy = y - oy;
				const Z ry = H - y;
				const Z rm = round(ry*ts.count()/(1e9*dy));

				ss << ",\t";

				if (rm > 3600) {
					if (rm/3600 < 10)
						ss << ' ';

					ss << rm/3600;

					if ((rm/60) % 60 < 10)
						ss << '0';

					ss << (rm/60) % 60 << ':';

					if (rm % 60 < 10)
						ss << '0';

					ss << rm % 60;
				} else if (rm > 60) {
					ss << " 0:";

					if (rm/60 < 10)
						ss << '0';

					ss << rm/60 << ':';

					if (rm % 60 < 10)
						ss << '0';

					ss << rm % 60;
				} else {
					ss << " 0:00:";

					if (rm < 10)
						ss << '0';

					ss << rm;
				}
			}

			ot = nt;
			oy = y;

			ss << endl;
			cerr << ss.str();
			cerr.flush();
		}*/

		const R yf = (y - YO)/RES;

		for (Z x = 0; x < W; x++) {
			const R xf = (x - XO)/RES;

            const RV z = mix(xf, X, yf, Y);

            R sum = 0;

			for (const ZV& r :  ROOTS)
				sum += cos(dot(r, z));

            sum *= 2;

            vs[W*y + x] = sum;
		}
	}

	R minE, maxE;

	maxE = N;
	minE = -N/SKEW;

	R pKneePF, posPF, nKneePF, negPF;

	pKneePF = KNEEP/maxE;
	posPF = 1/tanh(pKneePF*maxE);
	nKneePF = -KNEEM/minE;
	negPF = 1/tanh(nKneePF*minE);

	BMP out;
        out.SetSize(W, H);
        out.SetBitDepth(24);
        out.SetDPI(0, 0);

	yp = -5;

    for (Z y = 0; y < H; y++) {
		const Z nyp = round(100.*y/H);

		/*if (nyp >= yp + 5) {
			yp = nyp;

			ss.str("");
			ss << "O,\t";

			if (np < 10)
				ss << ' ';

			ss << np << ",\t";

			if (yp < 100) {
				ss << ' ';

				if (yp < 10)
					ss << ' ';
			}

			ss << yp << '%' << endl;
			cerr << ss.str();
			cerr.flush();
		}*/

		const R yf = (y - YO)/RES;

		for (Z x = 0; x < W; x++) {
			const R xf = (x - XO)/RES;
			R e = vs[W*y + x];

			RGBApixel& pixel = *out(x, y);

			if (minE > 0) {
				e = tanh(pKneePF*(e - minE));
				e *= posPF;

				pixel.Red = GAMMA(e);
				pixel.Green = GAMMA((1 - G)*e);
				pixel.Blue = GAMMA(e);
			} else
				if (e > 0) {
					e = tanh(pKneePF*e);
					e *= posPF;

					pixel.Red = GAMMA((1 - RB*JUMP)*e + RB*JUMP);
					pixel.Green = GAMMA((1 - G)*e);
					pixel.Blue = GAMMA(e);
				} else {
					e = tanh(nKneePF*e);
					e *= negPF;

					pixel.Red = 0;
					pixel.Green = GAMMA(G*e);
					pixel.Blue = GAMMA(JUMP*(1 - e));
				}
		}
	}

	ss.str("");

	ss << ff << ".bmp";

    out.WriteToFile(ss.str().c_str());
}
