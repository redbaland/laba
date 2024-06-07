#pragma once
#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

inline int NPS2 = 7;
inline int NPS3 = 16;
inline int N = 26;
inline int NUMBER_SAM = 0;
inline double B = 8.0 / 1000.0;
inline double DO = 630.0 / 1000.0;
inline double DI = DO - 2 * B;
inline double E = 0.2 / 1000 / DI;
inline double VISCOSITY = 15 * pow(10, (-6));
inline double PMAX = 6.3 * pow(10, 6);
inline double DENSITY = 850;
inline double PY = 10 * pow(10, 3);
inline double HK = 35;
inline double SUPPORT = 40;
inline double HMIN = 40;
inline double HMAX = PMAX / (9.81 * DENSITY);
inline double QLAST = 0;
inline double SAMANGLE = 0;
inline double HY = PY / (9.81 * DENSITY);
inline double XSAM = 0;
inline double HSAM = 0;
inline std::vector<double> NASOS2;
inline std::vector<double> NASOS4;
inline std::vector<double> PMAXMASS;
inline std::vector<double> PRESSURE;
inline std::vector<double> HEADMAX(N, 0);
inline std::vector<double> CAVITATION(N, 0);
inline std::vector<double> HEADX(N, 0);
inline std::vector<double> HEAD(N, 0);
inline std::vector<double> Q_FROM_PUMP2{1000.0, 1250.0, 1500.0, 1750.0, 2000.0};
inline std::vector<double> Q_FROM_PUMP1{1000.0, 1250.0, 1500.0, 1750.0, 2000.0};
inline std::vector<double> HEAD_FROM_PUMP2{260.0, 253.0, 245.0, 235.0, 222.0};
inline std::vector<double> HEAD_FROM_PUMP1{243.0, 237.0, 233.0, 225.0, 218.0};
inline std::vector<std::vector<double>> FINAL_SAMOTEK{{0, 0, 0, 0, 0},
                                               {0, 0, 0, 0, 0},
                                               {0, 0, 0, 0, 0},
                                               {0, 0, 0, 0, 0},
                                               {0, 0, 0, 0, 0}};
inline std::vector<double> PROFILE = {150.2, 153.3, 171.6, 96,    125.1, 155,   158,
                               134.5, 134.5, 112.7, 83.1,  102.2, 151.4, 109.5,
                               148.4, 83.2,  180.1, 180.1, 191.8, 177,   77.1,
                               161.3, 203.3, 118.3, 171.5, 203.1};
inline std::vector<double> LENGTH = {
    0,      20000,  40000,  60000,  80000,  100000, 120000, 140000, 140000,
    160000, 180000, 200000, 220000, 240000, 260000, 280000, 300000, 300000,
    320000, 340000, 360000, 380000, 400000, 420000, 440000, 460000};

