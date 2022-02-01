#include "car.h"

namespace {
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 5.991464547107981;

/******************************************************************************
 *                       Code generated with sympy 1.8                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3763381325127121296) {
   out_3763381325127121296[0] = delta_x[0] + nom_x[0];
   out_3763381325127121296[1] = delta_x[1] + nom_x[1];
   out_3763381325127121296[2] = delta_x[2] + nom_x[2];
   out_3763381325127121296[3] = delta_x[3] + nom_x[3];
   out_3763381325127121296[4] = delta_x[4] + nom_x[4];
   out_3763381325127121296[5] = delta_x[5] + nom_x[5];
   out_3763381325127121296[6] = delta_x[6] + nom_x[6];
   out_3763381325127121296[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2374270180543176146) {
   out_2374270180543176146[0] = -nom_x[0] + true_x[0];
   out_2374270180543176146[1] = -nom_x[1] + true_x[1];
   out_2374270180543176146[2] = -nom_x[2] + true_x[2];
   out_2374270180543176146[3] = -nom_x[3] + true_x[3];
   out_2374270180543176146[4] = -nom_x[4] + true_x[4];
   out_2374270180543176146[5] = -nom_x[5] + true_x[5];
   out_2374270180543176146[6] = -nom_x[6] + true_x[6];
   out_2374270180543176146[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_5587839573095323498) {
   out_5587839573095323498[0] = 1.0;
   out_5587839573095323498[1] = 0.0;
   out_5587839573095323498[2] = 0.0;
   out_5587839573095323498[3] = 0.0;
   out_5587839573095323498[4] = 0.0;
   out_5587839573095323498[5] = 0.0;
   out_5587839573095323498[6] = 0.0;
   out_5587839573095323498[7] = 0.0;
   out_5587839573095323498[8] = 0.0;
   out_5587839573095323498[9] = 1.0;
   out_5587839573095323498[10] = 0.0;
   out_5587839573095323498[11] = 0.0;
   out_5587839573095323498[12] = 0.0;
   out_5587839573095323498[13] = 0.0;
   out_5587839573095323498[14] = 0.0;
   out_5587839573095323498[15] = 0.0;
   out_5587839573095323498[16] = 0.0;
   out_5587839573095323498[17] = 0.0;
   out_5587839573095323498[18] = 1.0;
   out_5587839573095323498[19] = 0.0;
   out_5587839573095323498[20] = 0.0;
   out_5587839573095323498[21] = 0.0;
   out_5587839573095323498[22] = 0.0;
   out_5587839573095323498[23] = 0.0;
   out_5587839573095323498[24] = 0.0;
   out_5587839573095323498[25] = 0.0;
   out_5587839573095323498[26] = 0.0;
   out_5587839573095323498[27] = 1.0;
   out_5587839573095323498[28] = 0.0;
   out_5587839573095323498[29] = 0.0;
   out_5587839573095323498[30] = 0.0;
   out_5587839573095323498[31] = 0.0;
   out_5587839573095323498[32] = 0.0;
   out_5587839573095323498[33] = 0.0;
   out_5587839573095323498[34] = 0.0;
   out_5587839573095323498[35] = 0.0;
   out_5587839573095323498[36] = 1.0;
   out_5587839573095323498[37] = 0.0;
   out_5587839573095323498[38] = 0.0;
   out_5587839573095323498[39] = 0.0;
   out_5587839573095323498[40] = 0.0;
   out_5587839573095323498[41] = 0.0;
   out_5587839573095323498[42] = 0.0;
   out_5587839573095323498[43] = 0.0;
   out_5587839573095323498[44] = 0.0;
   out_5587839573095323498[45] = 1.0;
   out_5587839573095323498[46] = 0.0;
   out_5587839573095323498[47] = 0.0;
   out_5587839573095323498[48] = 0.0;
   out_5587839573095323498[49] = 0.0;
   out_5587839573095323498[50] = 0.0;
   out_5587839573095323498[51] = 0.0;
   out_5587839573095323498[52] = 0.0;
   out_5587839573095323498[53] = 0.0;
   out_5587839573095323498[54] = 1.0;
   out_5587839573095323498[55] = 0.0;
   out_5587839573095323498[56] = 0.0;
   out_5587839573095323498[57] = 0.0;
   out_5587839573095323498[58] = 0.0;
   out_5587839573095323498[59] = 0.0;
   out_5587839573095323498[60] = 0.0;
   out_5587839573095323498[61] = 0.0;
   out_5587839573095323498[62] = 0.0;
   out_5587839573095323498[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_4671392523931821141) {
   out_4671392523931821141[0] = state[0];
   out_4671392523931821141[1] = state[1];
   out_4671392523931821141[2] = state[2];
   out_4671392523931821141[3] = state[3];
   out_4671392523931821141[4] = state[4];
   out_4671392523931821141[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_4671392523931821141[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_4671392523931821141[7] = state[7];
}
void F_fun(double *state, double dt, double *out_1408835289975819796) {
   out_1408835289975819796[0] = 1;
   out_1408835289975819796[1] = 0;
   out_1408835289975819796[2] = 0;
   out_1408835289975819796[3] = 0;
   out_1408835289975819796[4] = 0;
   out_1408835289975819796[5] = 0;
   out_1408835289975819796[6] = 0;
   out_1408835289975819796[7] = 0;
   out_1408835289975819796[8] = 0;
   out_1408835289975819796[9] = 1;
   out_1408835289975819796[10] = 0;
   out_1408835289975819796[11] = 0;
   out_1408835289975819796[12] = 0;
   out_1408835289975819796[13] = 0;
   out_1408835289975819796[14] = 0;
   out_1408835289975819796[15] = 0;
   out_1408835289975819796[16] = 0;
   out_1408835289975819796[17] = 0;
   out_1408835289975819796[18] = 1;
   out_1408835289975819796[19] = 0;
   out_1408835289975819796[20] = 0;
   out_1408835289975819796[21] = 0;
   out_1408835289975819796[22] = 0;
   out_1408835289975819796[23] = 0;
   out_1408835289975819796[24] = 0;
   out_1408835289975819796[25] = 0;
   out_1408835289975819796[26] = 0;
   out_1408835289975819796[27] = 1;
   out_1408835289975819796[28] = 0;
   out_1408835289975819796[29] = 0;
   out_1408835289975819796[30] = 0;
   out_1408835289975819796[31] = 0;
   out_1408835289975819796[32] = 0;
   out_1408835289975819796[33] = 0;
   out_1408835289975819796[34] = 0;
   out_1408835289975819796[35] = 0;
   out_1408835289975819796[36] = 1;
   out_1408835289975819796[37] = 0;
   out_1408835289975819796[38] = 0;
   out_1408835289975819796[39] = 0;
   out_1408835289975819796[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_1408835289975819796[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_1408835289975819796[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1408835289975819796[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1408835289975819796[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_1408835289975819796[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_1408835289975819796[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_1408835289975819796[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_1408835289975819796[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_1408835289975819796[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_1408835289975819796[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1408835289975819796[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1408835289975819796[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_1408835289975819796[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_1408835289975819796[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_1408835289975819796[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1408835289975819796[56] = 0;
   out_1408835289975819796[57] = 0;
   out_1408835289975819796[58] = 0;
   out_1408835289975819796[59] = 0;
   out_1408835289975819796[60] = 0;
   out_1408835289975819796[61] = 0;
   out_1408835289975819796[62] = 0;
   out_1408835289975819796[63] = 1;
}
void h_25(double *state, double *unused, double *out_8167740060048081459) {
   out_8167740060048081459[0] = state[6];
}
void H_25(double *state, double *unused, double *out_136662046690260990) {
   out_136662046690260990[0] = 0;
   out_136662046690260990[1] = 0;
   out_136662046690260990[2] = 0;
   out_136662046690260990[3] = 0;
   out_136662046690260990[4] = 0;
   out_136662046690260990[5] = 0;
   out_136662046690260990[6] = 1;
   out_136662046690260990[7] = 0;
}
void h_24(double *state, double *unused, double *out_8028196242872299114) {
   out_8028196242872299114[0] = state[4];
   out_8028196242872299114[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4331524086921821580) {
   out_4331524086921821580[0] = 0;
   out_4331524086921821580[1] = 0;
   out_4331524086921821580[2] = 0;
   out_4331524086921821580[3] = 0;
   out_4331524086921821580[4] = 1;
   out_4331524086921821580[5] = 0;
   out_4331524086921821580[6] = 0;
   out_4331524086921821580[7] = 0;
   out_4331524086921821580[8] = 0;
   out_4331524086921821580[9] = 0;
   out_4331524086921821580[10] = 0;
   out_4331524086921821580[11] = 0;
   out_4331524086921821580[12] = 0;
   out_4331524086921821580[13] = 1;
   out_4331524086921821580[14] = 0;
   out_4331524086921821580[15] = 0;
}
void h_30(double *state, double *unused, double *out_472137416670460425) {
   out_472137416670460425[0] = state[4];
}
void H_30(double *state, double *unused, double *out_8969507209450629326) {
   out_8969507209450629326[0] = 0;
   out_8969507209450629326[1] = 0;
   out_8969507209450629326[2] = 0;
   out_8969507209450629326[3] = 0;
   out_8969507209450629326[4] = 1;
   out_8969507209450629326[5] = 0;
   out_8969507209450629326[6] = 0;
   out_8969507209450629326[7] = 0;
}
void h_26(double *state, double *unused, double *out_4661942715490933705) {
   out_4661942715490933705[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5603131469989754690) {
   out_5603131469989754690[0] = 0;
   out_5603131469989754690[1] = 0;
   out_5603131469989754690[2] = 0;
   out_5603131469989754690[3] = 0;
   out_5603131469989754690[4] = 0;
   out_5603131469989754690[5] = 0;
   out_5603131469989754690[6] = 0;
   out_5603131469989754690[7] = 1;
}
void h_27(double *state, double *unused, double *out_5277330746000624512) {
   out_5277330746000624512[0] = state[3];
}
void H_27(double *state, double *unused, double *out_7681925221614004014) {
   out_7681925221614004014[0] = 0;
   out_7681925221614004014[1] = 0;
   out_7681925221614004014[2] = 0;
   out_7681925221614004014[3] = 1;
   out_7681925221614004014[4] = 0;
   out_7681925221614004014[5] = 0;
   out_7681925221614004014[6] = 0;
   out_7681925221614004014[7] = 0;
}
void h_29(double *state, double *unused, double *out_1415647269572122860) {
   out_1415647269572122860[0] = state[1];
}
void H_29(double *state, double *unused, double *out_5145708440321797373) {
   out_5145708440321797373[0] = 0;
   out_5145708440321797373[1] = 1;
   out_5145708440321797373[2] = 0;
   out_5145708440321797373[3] = 0;
   out_5145708440321797373[4] = 0;
   out_5145708440321797373[5] = 0;
   out_5145708440321797373[6] = 0;
   out_5145708440321797373[7] = 0;
}
void h_28(double *state, double *unused, double *out_4685965241984472355) {
   out_4685965241984472355[0] = state[5];
   out_4685965241984472355[1] = state[6];
}
void H_28(double *state, double *unused, double *out_2686308391314673286) {
   out_2686308391314673286[0] = 0;
   out_2686308391314673286[1] = 0;
   out_2686308391314673286[2] = 0;
   out_2686308391314673286[3] = 0;
   out_2686308391314673286[4] = 0;
   out_2686308391314673286[5] = 1;
   out_2686308391314673286[6] = 0;
   out_2686308391314673286[7] = 0;
   out_2686308391314673286[8] = 0;
   out_2686308391314673286[9] = 0;
   out_2686308391314673286[10] = 0;
   out_2686308391314673286[11] = 0;
   out_2686308391314673286[12] = 0;
   out_2686308391314673286[13] = 0;
   out_2686308391314673286[14] = 1;
   out_2686308391314673286[15] = 0;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_3763381325127121296) {
  err_fun(nom_x, delta_x, out_3763381325127121296);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2374270180543176146) {
  inv_err_fun(nom_x, true_x, out_2374270180543176146);
}
void car_H_mod_fun(double *state, double *out_5587839573095323498) {
  H_mod_fun(state, out_5587839573095323498);
}
void car_f_fun(double *state, double dt, double *out_4671392523931821141) {
  f_fun(state,  dt, out_4671392523931821141);
}
void car_F_fun(double *state, double dt, double *out_1408835289975819796) {
  F_fun(state,  dt, out_1408835289975819796);
}
void car_h_25(double *state, double *unused, double *out_8167740060048081459) {
  h_25(state, unused, out_8167740060048081459);
}
void car_H_25(double *state, double *unused, double *out_136662046690260990) {
  H_25(state, unused, out_136662046690260990);
}
void car_h_24(double *state, double *unused, double *out_8028196242872299114) {
  h_24(state, unused, out_8028196242872299114);
}
void car_H_24(double *state, double *unused, double *out_4331524086921821580) {
  H_24(state, unused, out_4331524086921821580);
}
void car_h_30(double *state, double *unused, double *out_472137416670460425) {
  h_30(state, unused, out_472137416670460425);
}
void car_H_30(double *state, double *unused, double *out_8969507209450629326) {
  H_30(state, unused, out_8969507209450629326);
}
void car_h_26(double *state, double *unused, double *out_4661942715490933705) {
  h_26(state, unused, out_4661942715490933705);
}
void car_H_26(double *state, double *unused, double *out_5603131469989754690) {
  H_26(state, unused, out_5603131469989754690);
}
void car_h_27(double *state, double *unused, double *out_5277330746000624512) {
  h_27(state, unused, out_5277330746000624512);
}
void car_H_27(double *state, double *unused, double *out_7681925221614004014) {
  H_27(state, unused, out_7681925221614004014);
}
void car_h_29(double *state, double *unused, double *out_1415647269572122860) {
  h_29(state, unused, out_1415647269572122860);
}
void car_H_29(double *state, double *unused, double *out_5145708440321797373) {
  H_29(state, unused, out_5145708440321797373);
}
void car_h_28(double *state, double *unused, double *out_4685965241984472355) {
  h_28(state, unused, out_4685965241984472355);
}
void car_H_28(double *state, double *unused, double *out_2686308391314673286) {
  H_28(state, unused, out_2686308391314673286);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
