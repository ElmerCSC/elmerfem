FUNCTION initbedrock(Model, nodenumber, inputarray) RESULT(elevation)
  USE ElementDescription
  USE DefUtils

  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL (KIND=dp) :: inputarray(*),x,y,elevation,elevation_x,&
       elevation_y,z0, slope, centrality, trunk_hwidth, &
       obs_length,obs_height,elevation_obstacle, bump_height, &
       bump_radius, bump_1_dist, bump_1_elev, bump_2_dist, bump_2_elev

  x = inputarray(1)
  y = inputarray(2)

  !Y varies from 5000.0 at input to 0 at calving front
  !X varies from 0 at right margin to 5000.0 at left

  !longitudinal component of bed function
  slope = 1.0_dp/40.0
  z0 = -550.0_dp
  elevation_y = z0 + slope*y

  trunk_hwidth = (3000.0 + 2000.0 * (y/5000.0)) * 0.5

  centrality = ABS(x - 2500.0) / trunk_hwidth !1 at margins, 0 at centerline
  elevation_x = centrality**2.0 * 200.0 !vary depth from -550.0 at centerline to -350.0 at edge

  obs_length = 1000.0_dp
  obs_height = 100.0_dp
  elevation_obstacle = obs_height * EXP(-(y/obs_length)**2.0)

  bump_height = 100.0_dp
  bump_radius = 500.0_dp

  bump_1_dist = ((x - 1800.0)**2.0 + (y - 0.0)**2.0)**0.5
  bump_1_elev = bump_height * EXP(-(bump_1_dist / bump_radius)**2.0)

  bump_2_dist = ((x - 3000.0)**2.0 + (y - 0.0)**2.0)**0.5
  bump_2_elev = bump_height * EXP(-(bump_2_dist / bump_radius)**2.0)
  
  elevation = elevation_y + elevation_x + elevation_obstacle + bump_1_elev + bump_2_elev

END FUNCTION initbedrock

FUNCTION initsurface(Model, nodenumber, inputarray) RESULT(elevation)
  USE ElementDescription
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL (KIND=dp) :: inputarray(*),x,y, elevation, elevation_y, slope, z0

  x = inputarray(1)
  y = inputarray(2)

  !Y varies from 5000.0 at input to 0 at calving front
  !X varies from 0 at right margin to 5000.0 at left

  !longitudinal component of bed function
  slope = 1.0_dp/40.0
  z0 = 50.0_dp
  elevation_y = z0 + slope*y
  elevation = elevation_y
END FUNCTION initsurface
