#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

// TODO: uppercase globals
std::vector<double> PROFILE = {150.2, 153.3, 171.6, 96,    125.1, 155,   158,
                               134.5, 134.5, 112.7, 83.1,  102.2, 151.4, 109.5,
                               148.4, 83.2,  180.1, 180.1, 191.8, 177,   77.1,
                               161.3, 203.3, 118.3, 171.5, 203.1};
std::vector<double> DLINA = {
    0,      20000,  40000,  60000,  80000,  100000, 120000, 140000, 140000,
    160000, 180000, 200000, 220000, 240000, 260000, 280000, 300000, 300000,
    320000, 340000, 360000, 380000, 400000, 420000, 440000, 460000};
double b = 8.0 / 1000.0;
double Do = 630.0 / 1000.0;
double Di = Do - 2 * b;
double E = 0.2 / 1000 / Di;
double viscosity = 15 * pow(10, (-6));
int NPS2 = 7;
int NPS3 = 16;
int n = 26;
double Pmax = 6.3 * pow(10, 6);
double plotnost = 850;
double Py = 10 * pow(10, 3);
double Hk = 35;
double Podpor = 40;
double Hmin = 40;
double Hmax = Pmax / (9.81 * plotnost);
double Q_last = 0;
int NUMBER_SAM = 0;
double sam_angle = 0;
double Hy = Py / (9.81 * plotnost);
double x_sam = 0;
double h_sam = 0;
std::vector<double> Pmax_mass;
std::vector<double> Pressure;
std::vector<double> Head_max(PROFILE);
std::vector<double> kavitacia(n, 0);
std::vector<double> Head_x(n, 0);
std::vector<double> Head(n, 0);
std::vector<double> Q_FROM_PUMP2{1000, 1250, 1500, 1750, 2000};
std::vector<double> Q_FROM_PUMP1{1000, 1250, 1500, 1750, 2000};
std::vector<double> HEAD_FROM_PUMP2{260, 253, 245, 235, 222};
std::vector<double> HEAD_FROM_PUMP1{243, 237, 233, 225, 218};
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
double Bisection(double, double, double (*)(double));
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
  for (int i = 0; i < n; i++) {
    Head_max[i] += Hmax;
    kavitacia[i] = Hy + PROFILE[i];
  }
  kavitacia[n - 1] = 0;
  Head_x[n - 1] = Hk + PROFILE[n - 1];
  Head[n - 1] = Hk + PROFILE[n - 1];
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
// TODO: plots
void PlotEverything() { std::cout << "Here will be plots" << std::endl; }

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
  double Re = (Q * 4) / (M_PI * viscosity * Di);
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
  double alpha = fabs(atan(sam_angle));
  double a1 = pow((f - sin(f)), (5.0 / 3.0)) / pow(f, (2.0 / 3.0));
  double a2 = 0.2419 * Q_last / (pow(Di, (8.0 / 3.0)) * sqrt(sin(alpha)));
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
      lymbda * distance * 8 * pow(Q, 2) / (9.81 * pow(Di, 5) * pow(M_PI, 2));
  return lose_h;
}

double CalculateDeviation(double Q) {
  int no_iter = 0;
  Head = std::vector<double>(n, 0);
  Head[n - 1] = Hk + PROFILE[n - 1];
  std::vector<double> L(n, 0);
  for (int i = 0; i < n; i++)
    L[i] = DLINA[i];
  std::vector<double> Z(n, 0);
  for (int i = 0; i < n; i++)
    Z[i] = PROFILE[i];
  std::vector<double> Kav = kavitacia;
  std::vector<double> MAX = Head_max;
  for (int i = n - 2; i >= 0; i--) {
    Head[i] = Head[i + 1] + HeadLoses(L[i + 1] - L[i], Q);
    if (i == NPS2)
      Head[i] = Head[i + 1] - 2 * PumpsHead(Q, NASOS4);
    if (i == NPS3)
      Head[i] = Head[i + 1] - 2 * PumpsHead(Q, NASOS2);
    if (Head[i] < 0)
      no_iter = 1;
    if (Head[i] < kavitacia[i] && no_iter == 0) {
      if (i == 0 || i == NPS2 || i == NPS3) {
        std::cout << "Нарушено условие подпора" << std::endl;
      } else {
        x_sam = L[i + 1] +
                (L[i + 1] - L[i]) * (Head[i + 1] - Z[i + 1] - Hy) /
                    ((Z[i + 1] - Z[i]) + (L[i + 1] - L[i]) * HeadLoses(1, Q));
        h_sam = (Z[i + 1] - Z[i]) / (L[i + 1] - L[i]) * (x_sam - L[i]) + Z[i];
        L.insert(L.begin() + (i + 1), x_sam);
        Z.insert(Z.begin() + (i + 1), h_sam);
        Head.insert(Head.begin() + (i + 1), Hy + Z[i + 1]);
        MAX.insert(MAX.begin() + (i + 1), Hmax + Z[i + 1]);
        Kav.insert(Kav.begin() + (i + 1), Hy + Z[i + 1]);
        Head[i] = Hy + Z[i];
      }
    }
  }
  double a1 = Head[0];
  double a2 = Podpor + 2 * PumpsHead(Q, NASOS2) + PROFILE[0];
  return a2 - a1;
}

void CalculateConsumption() {
  Q_last = Bisection(0.2, 1.2, &CalculateDeviation);
  std::cout << "Расход [м3/c] = " << Q_last << std::endl;
}

void CalculateGravityFlowAreas() {
  for (int i = n - 2; i >= 0; i--) {
    Head[i] = Head[i + 1] + HeadLoses(DLINA[i + 1] - DLINA[i], Q_last);
    if (i == NPS2)
      Head[i] = Head[i + 1] - 2 * PumpsHead(Q_last, NASOS4);
    if (i == NPS3)
      Head[i] = Head[i + 1] - 2 * PumpsHead(Q_last, NASOS2);
    if (Head[i] < kavitacia[i]) {
      NUMBER_SAM++;
      x_sam = DLINA[i + 1] +
              (DLINA[i + 1] - DLINA[i]) * (Head[i + 1] - PROFILE[i + 1] - Hy) /
                  ((PROFILE[i + 1] - PROFILE[i]) +
                   (DLINA[i + 1] - DLINA[i]) * HeadLoses(1, Q_last));

      h_sam = (PROFILE[i + 1] - PROFILE[i]) / (DLINA[i + 1] - DLINA[i]) *
                  (x_sam - DLINA[i]) +
              PROFILE[i];
      DLINA.insert(DLINA.begin() + (i + 1), x_sam);
      PROFILE.insert(PROFILE.begin() + (i + 1), h_sam);
      Head.insert(Head.begin() + (i + 1), Hy + PROFILE[i + 1]);
      Head_max.insert(Head_max.begin() + (i + 1), Hmax + PROFILE[i + 1]);
      kavitacia.insert(kavitacia.begin() + (i + 1), Hy + PROFILE[i + 1]);
      Head[i] = Hy + PROFILE[i];
      sam_angle = (PROFILE[i] - h_sam) / (DLINA[i] - x_sam);
      double f = Bisection(0.1, 2 * M_PI, &Fi);
      double part = ((f - sin(f)) / (2 * M_PI) * 100);
      FINAL_SAMOTEK[NUMBER_SAM - 1] = {DLINA[i], x_sam, PROFILE[i], h_sam,
                                       part};
    }
  }
  if (NUMBER_SAM > 0) {
    std::cout << "Количество самотечных участков = " << NUMBER_SAM << std::endl;
    for (int i = 0; i < 5; i++) {
      std::cout << "(координата конца, координата начала, высотка конца, "
                   "высотка начала, степень заполнения): "
                << i + 1 << " = " << FINAL_SAMOTEK[i][0] << FINAL_SAMOTEK[i][1]
                << FINAL_SAMOTEK[i][2] << FINAL_SAMOTEK[i][3]
                << FINAL_SAMOTEK[i][4] << std::endl;
    }
  } else {
    std::cout << "Самотечные участки отсутствуют" << std::endl;
  }
}

void CalculatePressure() {
  Pmax_mass = std::vector<double>(Head.size(), Pmax / 1000000);
  Pressure = std::vector<double>(Head.size(), 0);
  for (int i = Head.size() - 1; i >= 0; i--)
    Pressure[i] = (Head[i] - PROFILE[i]) * (plotnost * 9.81) / 1000000;
}
