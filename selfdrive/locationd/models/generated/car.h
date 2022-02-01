#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_3763381325127121296);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2374270180543176146);
void car_H_mod_fun(double *state, double *out_5587839573095323498);
void car_f_fun(double *state, double dt, double *out_4671392523931821141);
void car_F_fun(double *state, double dt, double *out_1408835289975819796);
void car_h_25(double *state, double *unused, double *out_8167740060048081459);
void car_H_25(double *state, double *unused, double *out_136662046690260990);
void car_h_24(double *state, double *unused, double *out_8028196242872299114);
void car_H_24(double *state, double *unused, double *out_4331524086921821580);
void car_h_30(double *state, double *unused, double *out_472137416670460425);
void car_H_30(double *state, double *unused, double *out_8969507209450629326);
void car_h_26(double *state, double *unused, double *out_4661942715490933705);
void car_H_26(double *state, double *unused, double *out_5603131469989754690);
void car_h_27(double *state, double *unused, double *out_5277330746000624512);
void car_H_27(double *state, double *unused, double *out_7681925221614004014);
void car_h_29(double *state, double *unused, double *out_1415647269572122860);
void car_H_29(double *state, double *unused, double *out_5145708440321797373);
void car_h_28(double *state, double *unused, double *out_4685965241984472355);
void car_H_28(double *state, double *unused, double *out_2686308391314673286);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}