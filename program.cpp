#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <matplot/matplot.h>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

std::vector<double> PROFILE = {150.2, 153.3, 171.6, 96,    125.1, 155,   158,
                               134.5, 134.5, 112.7, 83.1,  102.2, 151.4, 109.5,
                               148.4, 83.2,  180.1, 180.1, 191.8, 177,   77.1,
                               161.3, 203.3, 118.3, 171.5, 203.1};
std::vector<double> LENGTH = {
    0,      20000,  40000,  60000,  80000,  100000, 120000, 140000, 140000,
    160000, 180000, 200000, 220000, 240000, 260000, 280000, 300000, 300000,
    320000, 340000, 360000, 380000, 400000, 420000, 440000, 460000};
double B = 8.0 / 1000.0;
double DO = 630.0 / 1000.0;
double DI = DO - 2 * B;
double E = 0.2 / 1000 / DI;
double VISCOSITY = 15 * pow(10, (-6));
int NPS2 = 7;
int NPS3 = 16;
int N = 26;
double PMAX = 6.3 * pow(10, 6);
double DENSITY = 850;
double PY = 10 * pow(10, 3);
double HK = 35;
double SUPPORT = 40;
double HMIN = 40;
double HMAX = PMAX / (9.81 * DENSITY);
double QLAST = 0;
int NUMBER_SAM = 0;
double SAMANGLE = 0;
double HY = PY / (9.81 * DENSITY);
double XSAM = 0;
double HSAM = 0;
std::vector<double> PMAXMASS;
std::vector<double> PRESSURE;
std::vector<double> HEADMAX(PROFILE);
std::vector<double> CAVITATION(N, 0);
std::vector<double> HEADX(N, 0);
std::vector<double> HEAD(N, 0);
std::vector<double> Q_FROM_PUMP2{1000.0, 1250.0, 1500.0, 1750.0, 2000.0};
std::vector<double> Q_FROM_PUMP1{1000.0, 1250.0, 1500.0, 1750.0, 2000.0};
std::vector<double> HEAD_FROM_PUMP2{260.0, 253.0, 245.0, 235.0, 222.0};
std::vector<double> HEAD_FROM_PUMP1{243.0, 237.0, 233.0, 225.0, 218.0};
std::vector<double> NASOS2;
std::vector<double> NASOS4;
std::vector<std::vector<double>> FINAL_SAMOTEK{{0, 0, 0, 0, 0},
                                               {0, 0, 0, 0, 0},
                                               {0, 0, 0, 0, 0},
                                               {0, 0, 0, 0, 0},
                                               {0, 0, 0, 0, 0}};

void CountGlobals();
void CalculateConsumption();
double CalculateDeviation(double);
void CalculateGravityFlowAreas();
void CalculatePressure();
void PlotEverything();
std::vector<double> SolveLlss(std::vector<double>, std::vector<double>);

int main() {
  CountGlobals();
  CalculateConsumption();
  CalculateGravityFlowAreas();
  CalculatePressure();
  PlotEverything();
  return 0;
}

void CountGlobals() {
  for (int i = 0; i < N; i++) {
    HEADMAX[i] += HMAX;
    CAVITATION[i] = HY + PROFILE[i];
  }
  CAVITATION[N - 1] = 0;
  HEADX[N - 1] = HK + PROFILE[N - 1];
  HEAD[N - 1] = HK + PROFILE[N - 1];
  NASOS2 = SolveLlss(Q_FROM_PUMP1, HEAD_FROM_PUMP1);
  NASOS4 = SolveLlss(Q_FROM_PUMP2, HEAD_FROM_PUMP2);
}

std::vector<double> SolveLlss(std::vector<double> Q_FROM_PUMP,
                              std::vector<double> HEAD_FROM_PUMP) {
  std::vector<double> res;
  Eigen::Matrix<double, 5, 2> A;
  Eigen::Matrix<double, 5, 1> b;
  for (int i = 0; i < 5; i++) {
    A.coeffRef(i, 0) = pow(Q_FROM_PUMP[i], 2);
    A.coeffRef(i, 1) = 1;
    b[i] = HEAD_FROM_PUMP[i];
  }
  auto Res = A.colPivHouseholderQr().solve(b);
  res.resize(Res.size());
  for (int i = 0; i < Res.size(); i++)
    res[i] = Res[i];
  res[0] *= (-1);
  return res;
}

void PlotEverything() {
  using namespace matplot;
  tiledlayout(1, 2);
  auto ax1 = nexttile();
  plot(LENGTH, PROFILE, LENGTH, HEADMAX, LENGTH, HEAD);
  ylabel(ax1, "Напор, м");
  xlabel(ax1 "Расстояние, км");
  title(ax1, "График распределения напоров");
  auto ax2 = nexttile();
  plot(ax2, LENGTH, PRESSURE, LENGTH, PMAXMASS);
  ylabel(ax2, "Давление, МПа");
  xlabel(ax2, "Расстояние, км");
  title(ax2, "График распределения давлений");
  ::matplot::legend(ax1, {"Высотные отметки Z(X)",
                          "Несущая способность Hmax(X)", "Напор H(X)"});
  ::matplot::legend(ax2, {"Давление P(X)", "Несущая способность Pmax(X)"});
  show();
}

double Bisection(double a, double b, double (*function)(double)) {
  if (function(a) * function(b) >= 0) {
    std::cout << "Метод бисекции не применим" << std::endl;
    return 0.0;
  }
  int iter_count = 0;
  while ((b - a) / 2 > 0 && iter_count < 500) {
    double c = (a + b) / 2;
    if (function(c) == 0)
      return c;
    else if (function(c) * function(a) < 0)
      b = c;
    else
      a = c;
    iter_count++;
  }
  return (a + b) / 2;
}

double Resistance(double Q) {
  double Re = (Q * 4) / (M_PI * VISCOSITY * DI);
  double hydraulicResistance = 0;
  if (Re <= 2320)
    hydraulicResistance = 64 / Re;
  if (2320 <= Re <= 10 * pow(E, (-1)))
    hydraulicResistance = 0.3164 / pow(Re, 0.25);
  if (10 * pow(E, (-1)) < Re < 500 * pow(E, (-1)))
    hydraulicResistance = 0.11 * pow((E + (68 / Re)), 0.25);
  if (Re >= 500 * pow(E, (-1)))
    hydraulicResistance = 0.11 * pow(E, 0.25);
  return hydraulicResistance;
}

double Fi(double f) {
  double alpha = fabs(atan(SAMANGLE));
  double a1 = pow((f - sin(f)), (5.0 / 3.0)) / pow(f, (2.0 / 3.0));
  double a2 = 0.2419 * QLAST / (pow(DI, (8.0 / 3.0)) * sqrt(sin(alpha)));
  double deviation = a1 - a2;
  return deviation;
}

double PumpsHead(double Q, std::vector<double> coeff_aprocs) {
  double diff_head = coeff_aprocs[1] - coeff_aprocs[0] * pow((Q * 3600), 2);
  return diff_head;
}

double HeadLoses(double distance, double Q) {
  double lymbda = Resistance(Q);
  double lose_h =
      lymbda * distance * 8 * pow(Q, 2) / (9.81 * pow(DI, 5) * pow(M_PI, 2));
  return lose_h;
}

double CalculateDeviation(double Q) {
  int no_iter = 0;
  HEAD = std::vector<double>(N, 0);
  HEAD[N - 1] = HK + PROFILE[N - 1];
  std::vector<double> L(N, 0);
  for (int i = 0; i < N; i++)
    L[i] = LENGTH[i];
  std::vector<double> Z(N, 0);
  for (int i = 0; i < N; i++)
    Z[i] = PROFILE[i];
  std::vector<double> Kav = CAVITATION;
  std::vector<double> MAX = HEADMAX;
  for (int i = N - 2; i >= 0; i--) {
    HEAD[i] = HEAD[i + 1] + HeadLoses(L[i + 1] - L[i], Q);
    if (i == NPS2)
      HEAD[i] = HEAD[i + 1] - 2 * PumpsHead(Q, NASOS4);
    if (i == NPS3)
      HEAD[i] = HEAD[i + 1] - 2 * PumpsHead(Q, NASOS2);
    if (HEAD[i] < 0)
      no_iter = 1;
    if (HEAD[i] < CAVITATION[i] && no_iter == 0) {
      if (i == 0 || i == NPS2 || i == NPS3) {
        std::cout << "Нарушено условие подпора" << std::endl;
      } else {
        XSAM = L[i + 1] +
               (L[i + 1] - L[i]) * (HEAD[i + 1] - Z[i + 1] - HY) /
                   ((Z[i + 1] - Z[i]) + (L[i + 1] - L[i]) * HeadLoses(1, Q));
        HSAM = (Z[i + 1] - Z[i]) / (L[i + 1] - L[i]) * (XSAM - L[i]) + Z[i];
        L.insert(L.begin() + (i + 1), XSAM);
        Z.insert(Z.begin() + (i + 1), HSAM);
        HEAD.insert(HEAD.begin() + (i + 1), HY + Z[i + 1]);
        MAX.insert(MAX.begin() + (i + 1), HMAX + Z[i + 1]);
        Kav.insert(Kav.begin() + (i + 1), HY + Z[i + 1]);
        HEAD[i] = HY + Z[i];
      }
    }
  }
  double a1 = HEAD[0];
  double a2 = SUPPORT + 2 * PumpsHead(Q, NASOS2) + PROFILE[0];
  return a2 - a1;
}

void CalculateConsumption() {
  QLAST = Bisection(0.2, 1.2, &CalculateDeviation);
  std::cout << "Расход [м3/c] =" << QLAST << std::endl;
}

void CalculateGravityFlowAreas() {
  for (int i = N - 2; i >= 0; i--) {
    HEAD[i] = HEAD[i + 1] + HeadLoses(LENGTH[i + 1] - LENGTH[i], QLAST);
    if (i == NPS2)
      HEAD[i] = HEAD[i + 1] - 2 * PumpsHead(QLAST, NASOS4);
    if (i == NPS3)
      HEAD[i] = HEAD[i + 1] - 2 * PumpsHead(QLAST, NASOS2);
    if (HEAD[i] < CAVITATION[i]) {
      NUMBER_SAM++;
      XSAM = LENGTH[i + 1] +
             (LENGTH[i + 1] - LENGTH[i]) * (HEAD[i + 1] - PROFILE[i + 1] - HY) /
                 ((PROFILE[i + 1] - PROFILE[i]) +
                  (LENGTH[i + 1] - LENGTH[i]) * HeadLoses(1, QLAST));

      HSAM = (PROFILE[i + 1] - PROFILE[i]) / (LENGTH[i + 1] - LENGTH[i]) *
                 (XSAM - LENGTH[i]) +
             PROFILE[i];
      LENGTH.insert(LENGTH.begin() + (i + 1), XSAM);
      PROFILE.insert(PROFILE.begin() + (i + 1), HSAM);
      HEAD.insert(HEAD.begin() + (i + 1), HY + PROFILE[i + 1]);
      HEADMAX.insert(HEADMAX.begin() + (i + 1), HMAX + PROFILE[i + 1]);
      CAVITATION.insert(CAVITATION.begin() + (i + 1), HY + PROFILE[i + 1]);
      HEAD[i] = HY + PROFILE[i];
      SAMANGLE = (PROFILE[i] - HSAM) / (LENGTH[i] - XSAM);
      double f = Bisection(0.1, 2 * M_PI, &Fi);
      double part = ((f - sin(f)) / (2 * M_PI) * 100);
      FINAL_SAMOTEK[NUMBER_SAM - 1] = {LENGTH[i], XSAM, PROFILE[i], HSAM, part};
    }
  }
  if (NUMBER_SAM > 0) {
    std::cout << "Количество самотечных участков = " << NUMBER_SAM << std::endl;
    for (int i = 0; i < 5; i++) {
      std::cout << "(координата конца, координата начала, высотка конца, "
                   "высотка начала, степень заполнения):"
                << i + 1 << " = " << FINAL_SAMOTEK[i][0] << FINAL_SAMOTEK[i][1]
                << FINAL_SAMOTEK[i][2] << FINAL_SAMOTEK[i][3]
                << FINAL_SAMOTEK[i][4] << std::endl;
    }
  } else {
    std::cout << "Самотечные участки отсутствуют" << std::endl;
  }
}

void CalculatePressure() {
  PMAXMASS = std::vector<double>(HEAD.size(), PMAX / 1000000);
  PRESSURE = std::vector<double>(HEAD.size(), 0);
  for (int i = (HEAD.size() - 1); i >= 0; i--)
    PRESSURE[i] = (HEAD[i] - PROFILE[i]) * (DENSITY * 9.81) / 1000000;
}
