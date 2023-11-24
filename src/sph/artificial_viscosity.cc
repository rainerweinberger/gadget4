/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  artificial_viscosity.cc
 *
 *  \brief Calculates time-dependent artificial viscosity parameter
 */

#include "gadgetconfig.h"

#include <gsl/gsl_linalg.h>

#include "../data/sph_particle_data.h"

/*! This file contains the function for the time-dependent artificial viscosity
 */

#ifdef IMPROVED_VELOCITY_GRADIENTS
void sph_particle_data::set_velocity_gradients(void)
{
#ifdef ONEDIMS
  if(fabs(dpos.dx_dx) > 0)
    {
      dvel[0][0] = dvel[0][0] / dpos.dx_dx;
      DivVel     = dvel[0][0];
    }
  else
    {
      DivVel = dvel[0][0] / Density;
    }
#elif defined(TWODIMS)
  double det = dpos.dx_dx * dpos.dy_dy - dpos.dx_dy * dpos.dx_dy;
  if(fabs(det) > 0)
    {
      double m11_inv = dpos.dy_dy / det;
      double m12_inv = -dpos.dx_dy / det;
      double m21_inv = -dpos.dx_dy / det;
      double m22_inv = dpos.dx_dx / det;

      double y11 = dvel[0][0];
      double y21 = dvel[0][1];
      double y12 = dvel[1][0];
      double y22 = dvel[1][1];

      dvel[0][0] = m11_inv * y11 + m12_inv * y21;
      dvel[0][1] = m21_inv * y11 + m22_inv * y21;
      dvel[1][0] = m11_inv * y12 + m12_inv * y22;
      dvel[1][1] = m21_inv * y12 + m22_inv * y22;
      DivVel     = dvel[0][0] + dvel[1][1];
      CurlVel    = fabs(dvel[0][1] - dvel[1][0]);
    }
  else
    {
      DivVel  = (dvel[0][0] + dvel[1][1]) / Density;
      CurlVel = fabs((dvel[1][0] - dvel[0][1]) / Density);
    }
#else
  gsl_matrix* distance_matrix = gsl_matrix_alloc(3, 3);
  gsl_matrix_set(distance_matrix, 0, 0, dpos.dx_dx);
  gsl_matrix_set(distance_matrix, 0, 1, dpos.dx_dy);
  gsl_matrix_set(distance_matrix, 0, 2, dpos.dx_dz);
  gsl_matrix_set(distance_matrix, 1, 0, dpos.dx_dy);
  gsl_matrix_set(distance_matrix, 1, 1, dpos.dy_dy);
  gsl_matrix_set(distance_matrix, 1, 2, dpos.dy_dz);
  gsl_matrix_set(distance_matrix, 2, 0, dpos.dx_dz);
  gsl_matrix_set(distance_matrix, 2, 1, dpos.dy_dz);
  gsl_matrix_set(distance_matrix, 2, 2, dpos.dz_dz);

  int sign                = 1;
  gsl_permutation* permut = gsl_permutation_alloc(3);
  gsl_linalg_LU_decomp(distance_matrix, permut, &sign);

  double det = gsl_linalg_LU_det(distance_matrix, sign);

  if(fabs(det) > 0)
    {
      gsl_matrix* inv_distance_matrix = gsl_matrix_alloc(3, 3);
      gsl_linalg_LU_invert(distance_matrix, permut, inv_distance_matrix);

      double m_inv[3][3];
      for(int i = 0; i < 3; i++)
        {
          for(int j = 0; j < 3; j++)
            {
              m_inv[i][j] = gsl_matrix_get(inv_distance_matrix, i, j);
            }
        }

      double y[3][3];
      for(int i = 0; i < 3; i++)
        {
          for(int j = 0; j < 3; j++)
            {
              y[i][j] = dvel[i][j];
            }
        }

      for(int i = 0; i < 3; i++)
        {
          for(int j = 0; j < 3; j++)
            {
              dvel[i][j] = 0;
              for(int k = 0; k < 3; k++)
                {
                  dvel[i][j] += m_inv[j][k] * y[k][i];
                }
            }
        }

      DivVel  = dvel[0][0] + dvel[1][1] + dvel[2][2];
      Rot[0]  = dvel[2][1] - dvel[1][2];
      Rot[1]  = dvel[0][2] - dvel[2][0];
      Rot[2]  = dvel[1][0] - dvel[0][1];
      CurlVel = sqrt(Rot[0] * Rot[0] + Rot[1] * Rot[1] + Rot[2] * Rot[2]);

      gsl_permutation_free(permut);
      gsl_matrix_free(distance_matrix);
      gsl_matrix_free(inv_distance_matrix);
    }
  else
    {
      DivVel  = (dvel[0][0] + dvel[1][1] + dvel[2][2]) / Density;
      Rot[0]  = (dvel[2][1] - dvel[1][2]) / Density;
      Rot[1]  = (dvel[0][2] - dvel[2][0]) / Density;
      Rot[2]  = (dvel[1][0] - dvel[0][1]) / Density;
      CurlVel = sqrt(Rot[0] * Rot[0] + Rot[1] * Rot[1] + Rot[2] * Rot[2]);
    }

#endif
}
#endif

#ifdef TIMEDEP_ART_VISC
/*UPS: This function is called in density.cc and loops over all (active) particles to forward integrate Alpha.
       I.e. Alpha below is really Alpha_p which is either set to alpha_tar, if the Alpha from the previous step
       is smaller than alpha_tar. Otherwise we want to let Alpha decay, but make sure it is never falling below zero. 
       G3 and Phantom actually trace DAlpha/dt and integrate it (G3 in kicks.c, Phantom with an 
       backward (implcit) Euler). Might be worth testing in a future iteration. */
void sph_particle_data::set_viscosity_coefficient(double dt)
{
  double dDivVel_dt = dt > 0 ? (DivVel - DivVelOld) / (dt) : 0;
  dDivVel_dt *= All.cf_a2inv; //now in physical coordinates
  double shockIndicator = -dDivVel_dt > 0 ? -dDivVel_dt : 0; //Cullen & Dehnen (2012), shock capturing
  double hsml_p           = Hsml * All.cf_atime;
  double hsml2_p          = hsml_p * hsml_p;
  double csnd_p           = Csnd * All.cf_afac3; 
  double decayVel2        = decayVel * decayVel * All.cf_afac3 * All.cf_afac3;
  double DivVel2       = (DivVel * All.cf_a2inv) * (DivVel * All.cf_a2inv);
  double CurlVel2      = (CurlVel * All.cf_a2inv) * (CurlVel * All.cf_a2inv);
  double CsndOverHsml2 = (csnd_p / hsml_p) * (csnd_p / hsml_p);
  double limiter       = DivVel2 / (DivVel2 + CurlVel2 + 0.00001 * CsndOverHsml2);
  /*UPS: Cullen and Dehnen (2010) take the limiter into the alpha_loc, which is alpha_tar here and take the the decayVel as 
         instead of the sound speed to estiamte alpha_loc. Thsi seems to reduce post hsock ringing quite a bit. */
  double alpha_tar        = (hsml2_p * limiter * shockIndicator) / (hsml2_p * limiter * shockIndicator + decayVel2) * All.ArtBulkViscConst;
#ifdef NO_SHEAR_VISCOSITY_LIMITER
  limiter = 1.;
#endif

  if(Alpha < alpha_tar)
    {
      Alpha = alpha_tar;
      return;
    }

  double devayVel_p = decayVel * All.cf_afac3;  //has the same a factor as sound speed

  //double DecayTime = 10. * hsml_p / devayVel_p; Price et al. (2018), eq. 48
  double DecayTime = 20. * hsml_p / devayVel_p; //Cullen & Dehnen for AlphaMin=0
  //Alpha            = limiter * (alpha_tar + (Alpha - alpha_tar) * exp(-dt / DecayTime));
  Alpha            = (alpha_tar + (Alpha - alpha_tar) * exp(-dt / DecayTime));
  if(Alpha < All.AlphaMin){
    Alpha = All.AlphaMin;
    printf("We should never arrive here");
  }
}

#endif
