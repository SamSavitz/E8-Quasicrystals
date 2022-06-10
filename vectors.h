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
    cout.precision(18);

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

inline R sigmoid(R t) {
    t = 2*t - 1;
    return (t*(t*t*(t*t*(21 - 5*t*t) - 35) + 35) + 16)/32;
}

inline R zoom(const R t) {
    return 0.0000000000015625*(129140163*t*t*t*t*t*t*(t - 1)*(t - 1)*(t - 1)*(t - 1)*(159190625 + 27*t*(-35247500 + t*(83659250 + 9*t*(4321179*t - 10160300)))));
}

inline RVP interpolate(const RVP& v1, const RVP& v2, const Z f) {
    R s = (R) f/FPS;
    R r = 0, z = 0;
    bool trans = false;

    if (s >= ST) {
        s -= ST;

        if (s < ZT)
            z = zoom(s/ZT);
        else {
            s -= ZT;

            if (s >= TT) {
                r = (s - TT)/RT;
                trans = true;
            }
        }
    }

    r = sigmoid(r);
    z = 1/(RES1 + (RES2 - RES1)*z);

    cout << trans << '\t' << r << '\t' << z << '\n';

    //cerr << f << ' ' << r << ' ' << z << endl;

    RV X2 = mix(1 - r, v1[0], r, v2[0] );
    RV Y2 = mix(1 - r, v1[1], r, v2[1] );

    if (f) {
        normalize(X2);
        normalize(Y2);
    }

    const R a = dot(X2, Y2);

    const RV X = mix(z*(sqrt(1 - a) + sqrt(1 + a))/2/sqrt(1 - a*a), X2, z*(sqrt(1 - a) - sqrt(1 + a))/2/sqrt(1 - a*a), Y2);
    const RV Y = mix(z*(sqrt(1 - a) - sqrt(1 + a))/2/sqrt(1 - a*a), X2, z*(sqrt(1 - a) + sqrt(1 + a))/2/sqrt(1 - a*a), Y2);

    for (Z i = 0; i < D - 1; ++i)
        cout << X[i] << '\t';

    cout << X[D - 1] << '\n';

    for (Z i = 0; i < D - 1; ++i)
        cout << Y[i] << '\t';

    cout << Y[D - 1] << endl;

    return RVP{X, Y};
}
