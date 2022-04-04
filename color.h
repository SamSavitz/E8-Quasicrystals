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
	const R posPF = 1/tanh(pKneePF*maxE);
	constexpr R nKneePF = -KNEEM/minE;
	const R negPF = 1/tanh(nKneePF*minE);

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

    const R phi = 2*M_PI*(start/3 + rotations*i);
    const R lg = pow(i, gamma);
    const R alpha = hue*lg*(1 - lg)/2;

    const R cphi = cos(phi);
    const R sphi = sin(phi);

    const R r = lg + alpha*(-0.14861*cphi + 1.78277*sphi);
    const R g = lg + alpha*(-0.29227*cphi - 0.90649*sphi);
    const R b = lg + alpha*( 1.97294*cphi               );

    return RGB((unsigned char) (255*r), (unsigned char) (255*g), (unsigned char) (255*b));
}

inline RGB hsl(R h, const R s, const R l) {
    const R c = (1 - fabs(2*l - 1))*s;
    h /= 60;
    const R x = c*(1 - fabs(fmod(h, 2) - 1));

    R r, g, b;

    switch ((Z) floor(h)) {
                case 0:
            r = c;
            g = x;
            b = 0;
        break;  case 1:
            r = x;
            g = c;
            b = 0;
        break;  case 2:
            r = 0;
            g = c;
            b = x;
        break;  case 3:
            r = 0;
            g = x;
            b = c;
        break;  case 4:
            r = x;
            g = 0;
            b = c;
        break;  case 5:
            r = c;
            g = 0;
            b = x;
        break;
    }

    const R m = l - c/2;

    r += m;
    g += m;
    b += m;

    r *= 255.999;
    g *= 255.999;
    b *= 255.999;

    return RGB(r, g, b);
}

inline RGB color(R i) {
    if (i < 0) {
    	constexpr R nKneePF = -KNEEM/minE;
    	const R negPF = 1/tanh(nKneePF*minE);

        i = negPF*tanh(nKneePF*i);

        //i = (minE - i)/minE;
        return cubehelix(i);
    } else {
    	constexpr R pKneePF = KNEEP/maxE;
    	const R posPF = 1/tanh(pKneePF*maxE);

	i = (1.5*posPF*tanh(pKneePF*i) + i/maxE)/2.5;

        //i = posPF*tanh(pKneePF*i);

        if      (i < 1./8)
            return hsl(     3840*i/7 + 240      , 2*i        , 5*i/2      );
        else if (i < 1./2)
            return hsl(fmod(3840*i/7 + 240, 360), 2*i        , (2*i + 1)/4);
        else if (i < 11./12)
            return hsl(fmod(3840*i/7 + 240, 360), 1          , 1./2       );
        else if (i < 17./18)
            return hsl(fmod(3840*i/7 + 240, 360), 1          , -18*i + 17 );
        else
            return RGB(0, 0, 0);
    }
}

// 11/12*x + b = 1/2, 17/18*x + b = 0
