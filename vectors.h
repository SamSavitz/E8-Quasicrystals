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

inline RVP interpolate(const RVP& v1, const RVP& v2, const Z f) {
    R t = (R) f/FS;

    t = 2*(t - 0.5);
    t = (t*(t*t*(t*t*(21 - 5*t*t) - 35) + 35) + 16)/32;

    RV X2 = mix(1 - t, v1[0], t, v2[0] );
    RV Y2 = mix(1 - t, v1[1], t, v2[1] );

    if (f) {
        normalize(X2);
        normalize(Y2);
    }

    const R a = dot(X2, Y2);

    const RV X = mix((sqrt(1 - a) + sqrt(1 + a))/2/sqrt(1 - a*a), X2, (sqrt(1 - a) - sqrt(1 + a))/2/sqrt(1 - a*a), Y2);
    const RV Y = mix((sqrt(1 - a) - sqrt(1 + a))/2/sqrt(1 - a*a), X2, (sqrt(1 - a) + sqrt(1 + a))/2/sqrt(1 - a*a), Y2);

    return RVP{X, Y};
}
