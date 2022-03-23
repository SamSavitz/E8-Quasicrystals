/* cubehelix() is based on dllu's work:
 * https://github.com/dllu/puppybrot/blob/master/cubehelix.cpp
 */

#include <cmath>
#include <sstream>

inline unsigned char GAMMA(R i) {
	if (i <= 0.0031308)
		i *= 12.92;
	else
		i = 1.055000*pow(i, 5./12) - 0.055001;

    return 256*i;
}

inline RGB barney(R e) {
    constexpr R RB = 72187./212655;
    constexpr R G = 250000./357579;
    constexpr R JUMP = 0.2;

    constexpr R KNEEP = 4;
    constexpr R KNEEM = 1;

	constexpr R pKneePF =  KNEEP/maxE;
	constexpr R posPF = 1/tanh(pKneePF*maxE);
	constexpr R nKneePF = -KNEEM/minE;
	constexpr R negPF = 1/tanh(nKneePF*minE);

    RGB p;

	if (e > 0) {
		e = tanh(pKneePF*e);
		e *= posPF;

		p.r = GAMMA((1 - RB*JUMP)*e + RB*JUMP);
		p.g = GAMMA((1 - G)*e);
		p.b = GAMMA(e);
	} else {
		e = tanh(nKneePF*e);
		e *= negPF;

		p.r = 0;
		p.g = GAMMA(G*e);
		p.b = GAMMA(JUMP*(1 - e));
	}

    return p;
}

/* Based on dllu's work:
 * https://github.com/dllu/puppybrot/blob/master/cubehelix.cpp
 */
inline RGB cubehelix(const R i) {
    constexpr R start = 0.5;
    constexpr R rotations = -2.5;
    constexpr R hue = 1.0;
    constexpr R gamma = 1.0;

    const double phi = 2*M_PI*(start/3 + rotations*i);
    const double lg = pow(i, gamma);
    const double alpha = hue*lg*(1 - lg)/2;

    const double cphi = cos(phi);
    const double sphi = sin(phi);

    const double r = lg + alpha*(-0.14861*cphi + 1.78277*sphi);
    const double g = lg + alpha*(-0.29227*cphi - 0.90649*sphi);
    const double b = lg + alpha*( 1.97294*cphi               );

    return RGB(255*r, 255*g, 255*b);
}

inline RGB hsl(R h, const R s, const R l) {
    const R c = (1 - fabs(2*l - 1))*s;
    h /= 60;
    const R x = c*(1 - fabs(fmod(h, 2) - 1));

    RGB r;

    switch ((Z) floor(h)) {
                case 0:
            r = RGB(c, x, 0);
        break;  case 1:
            r = RGB(x, c, 0);
        break;  case 2:
            r = RGB(0, c, x);
        break;  case 3:
            r = RGB(0, x, c);
        break;  case 4:
            r = RGB(x, 0, c);
        break;  case 5:
            r = RGB(c, 0, x);
        break;
    }

    const R m = l - c/2;

    r.r += m;
    r.g += m;
    r.b += m;

    r.r *= 255;
    r.g *= 255;
    r.b *= 255;

    return r;
}

inline RGB color(R i) {
    if (i < 0) {
        i = (minE - i)/minE;
        return cubehelix(1 - i);
    } else {
        constexpr R KNEEP = 4;
    	constexpr R pKneePF = KNEEP/maxE;
    	constexpr R posPF = 1/tanh(pKneePF*maxE);

		i = (posPF*tanh(pKneePF*i) + i/maxE)/2;

        if      (i < 1./8)
            return hsl(     480*i + 240      , 4*i + 1./2, 5*i/2      );
        else if (i < 1./2)
            return hsl(fmod(480*i + 240, 360), 1         , (2*i + 1)/4);
        else
            return hsl(fmod(480*i + 240, 360), 1         , 1./2       );
    }
}
