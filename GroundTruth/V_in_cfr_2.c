/**********************************************************************************
UDF for specifying steady-state velocity profile boundary condition on a 2D inlet
***********************************************************************************/

#include "udf.h"

DEFINE_PROFILE(inlet_velocity, thread, position)
{
  real i[ND_ND]; /* this will hold the position vector */
  real x, y, r, R, Vmax;
  face_t f;

  Vmax = 0.161;  /* maximum velocity (m/s)*/
  R = 0.0023*0.0023; /* inlet radius (m) */

  begin_f_loop(f,thread)
  {
    F_CENTROID(i, f, thread);
    x = i[0];
    y = i[1];
    r = x*x + y*y;
    F_PROFILE(f, thread, position) = Vmax*(1.0 - (r/R) );
  }
  end_f_loop(f, thread)
}