#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_3(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_19(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_8729358110465776939);
void live_err_fun(double *nom_x, double *delta_x, double *out_6985748779881356812);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_7782584527847184999);
void live_H_mod_fun(double *state, double *out_2807386512015406);
void live_f_fun(double *state, double dt, double *out_460231330011564285);
void live_F_fun(double *state, double dt, double *out_2331054952450132772);
void live_h_3(double *state, double *unused, double *out_8586752585035404355);
void live_H_3(double *state, double *unused, double *out_5698610554807408040);
void live_h_4(double *state, double *unused, double *out_8540102027610099);
void live_H_4(double *state, double *unused, double *out_1405536358766626176);
void live_h_9(double *state, double *unused, double *out_1977128371001968200);
void live_H_9(double *state, double *unused, double *out_3728041902637482146);
void live_h_10(double *state, double *unused, double *out_5757918426375340444);
void live_H_10(double *state, double *unused, double *out_2807340980368654320);
void live_h_12(double *state, double *unused, double *out_276606852640741133);
void live_H_12(double *state, double *unused, double *out_3512091870505835745);
void live_h_31(double *state, double *unused, double *out_552123329346621541);
void live_H_31(double *state, double *unused, double *out_837002095362415924);
void live_h_32(double *state, double *unused, double *out_3049938672845674483);
void live_H_32(double *state, double *unused, double *out_5094531309091300737);
void live_h_13(double *state, double *unused, double *out_4676866542797813623);
void live_H_13(double *state, double *unused, double *out_2474068813600300595);
void live_h_14(double *state, double *unused, double *out_1977128371001968200);
void live_H_14(double *state, double *unused, double *out_3728041902637482146);
void live_h_19(double *state, double *unused, double *out_5172964960218397550);
void live_H_19(double *state, double *unused, double *out_8160721753739913327);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}