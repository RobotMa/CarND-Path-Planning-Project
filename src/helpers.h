#ifndef HELPERS_H
#define HELPERS_H

#include <math.h>
#include <string>
#include <vector>
#include <cmath>
#include "spline.h"
// for convenience
using std::string;
using std::vector;

const double MAX_SPEED_MPS = 22.0;

// States of ego vehicle
struct CarState
{
	CarState(double car_x_val, double car_y_val, double car_s_val, double car_d_val,
			double car_yaw_val, double car_speed_val):
		car_x(car_x_val), car_y(car_y_val), car_s(car_s_val), car_d(car_d_val),
		car_yaw(car_yaw_val), car_speed(car_speed_val) {}
	~CarState() {}

	double car_x;
	double car_y;
	double car_s;
	double car_d;
	double car_yaw;
	double car_speed;
};

// State of other vehicle on the street
struct StreetVehicleState
{
	StreetVehicleState(int id, double x, double y, double vx, double vy, double s, double d):
		id_n(id), x_m(x), y_m(y), vx_mps(vx), vy_mps(vy), s_m(s), d_m(d) {}
	~StreetVehicleState() {}

	int id_n;
	double x_m;
	double y_m;
	double vx_mps;
	double vy_mps;
	double s_m;
	double d_m;
};

// Lane ID
enum class Lane : int
{
	LEFTLANE = 0,
	MIDLANE,
	RIGHTLANE,
	INVALID
};

// Occupied status of the lanes
enum class LaneState : int
{
	MIDLANE_OPEN = 0,
	MIDLANE_CLOSED,
	LEFTLANE_OPEN,
	LEFTLANE_CLOSED,
	RIGHTLANE_OPEN,
	RIGHTLANE_CLOSED,
	INVALID
};

struct TrajectoryGoal
{
	TrajectoryGoal(double refSpeed_mps, Lane targetLane):
		ref_speed_mps(refSpeed_mps), target_lane(targetLane){}
	~TrajectoryGoal() {}

	double ref_speed_mps{};
	Lane   target_lane{};

};

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
//   else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

//
// Helper functions related to waypoints and converting from XY to Frenet
//   or vice versa
//

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Calculate distance between two points
double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

// Calculate closest waypoint to current x, y position
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, 
                    const vector<double> &maps_y) {
  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for (int i = 0; i < maps_x.size(); ++i) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x,y,map_x,map_y);
    if (dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }
  }

  return closestWaypoint;
}

// Returns next waypoint of the closest waypoint
int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, 
                 const vector<double> &maps_y) {
  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y-y),(map_x-x));

  double angle = fabs(theta-heading);
  angle = std::min(2*pi() - angle, angle);

  if (angle > pi()/2) {
    ++closestWaypoint;
    if (closestWaypoint == maps_x.size()) {
      closestWaypoint = 0;
    }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, 
                         const vector<double> &maps_x, 
                         const vector<double> &maps_y) {
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if (next_wp == 0) {
    prev_wp  = maps_x.size()-1;
  }

  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point
  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if (centerToPos <= centerToRef) {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for (int i = 0; i < prev_wp; ++i) {
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }

  frenet_s += distance(0,0,proj_x,proj_y);

  return {frenet_s,frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, 
                     const vector<double> &maps_x, 
                     const vector<double> &maps_y) {
  int prev_wp = -1;

  while (s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1))) {
    ++prev_wp;
  }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),
                         (maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};
}

/*
 * @brief Convert lane ID into lateral d distance
 * @input Lane ID
 * @output Distance along d
 */
double laneToD(const Lane lane)
{
	const double laneWidth = 4.0;
	return ((static_cast<double>(lane) + 0.5)*laneWidth);
}

/*
 * @brief Convert d distance into lane ID
 * @input d : lateral shift
 * @output Lane ID
 */
Lane dToLane(const double d)
{
	if(0 < d && d < 4.0)
	{
		return Lane::LEFTLANE;
	}
	else if (4.0 <= d && d < 8.0)
	{
		return Lane::MIDLANE;
	}
	else if(8.0 <= d && d < 12.0)
	{
		return Lane::RIGHTLANE;
	}
	else
	{
		return Lane::INVALID;
	}
}

/*
 * @brief Generate the final trajectory
 * @input
 * @output Trajectory to follow for next iteration
 */
// TODO: Passing too many arguments into one function is bad and this is due to further clean up
void generatePath(const CarState& carState, const TrajectoryGoal& goal,
				  const std::vector<double>& map_waypoints_s,
				  const std::vector<double>& map_waypoints_x,
				  const std::vector<double>& map_waypoints_y,
				  const std::vector<double>& previous_path_x,
				  const std::vector<double>& previous_path_y,
		          std::vector<double>& next_x_vals,
				  std::vector<double>& next_y_vals)
{
	static int counter = 0;
	double pos_x, pos_y, angle;
	double pos_x2, pos_y2; 			// Predicted history position of the ego vehicle
	std::vector<double> pos_set_x;	// x of the planned path in vehicle frame
	std::vector<double> pos_set_y;	// y of the planned path in vehicle frame
	next_x_vals.clear();
	next_y_vals.clear();

	// Set history points of planned path in vehicle frame
	int path_size = previous_path_x.size();
	if (path_size == 0)
	{
	  pos_x = carState.car_x;
	  pos_y = carState.car_y;
	  angle = deg2rad(carState.car_yaw);
	  pos_x2 = pos_x - cos(angle);
	  pos_y2 = pos_y - sin(angle);
	}
	else
	{
	  pos_x = previous_path_x.at(path_size - 1);
	  pos_y = previous_path_y.at(path_size - 1);

	  pos_x2 = previous_path_x.at(path_size - 2);
	  pos_y2 = previous_path_y.at(path_size - 2);
	  angle = atan2(pos_y-pos_y2,pos_x-pos_x2);
	}

	pos_set_x.push_back(pos_x2);
	pos_set_x.push_back(pos_x);
	pos_set_y.push_back(pos_y2);
	pos_set_y.push_back(pos_y);

	// Mock lane status
	counter++;
	Lane lane = Lane::LEFTLANE;
//	if (counter < 100)
//	{
//		lane = Lane::LEFTLANE;
//	}
//	else if (100 <= counter && counter < 200)
//	{
//		lane = Lane::MIDLANE;
//	}
//	else if (200 <= counter && counter < 300)
//	{
//		lane = Lane::RIGHTLANE;
//	}
//	else
//	{
//		counter = 0;
//	}

	// Create lookahead way points in global frame based on Frenet frame derivation
	std::vector<double> lookaheadDistance{30.0, 60.0, 90.0};
	for (const auto& e : lookaheadDistance)
	{
		std::vector<double> waypoint = getXY(carState.car_s + e, laneToD(lane),
											   map_waypoints_s,
											   map_waypoints_x,
											   map_waypoints_y);
		pos_set_x.push_back(waypoint.at(0));
		pos_set_y.push_back(waypoint.at(1));
	}

	// Transform lookahead points from global frame to vehicle frame for spline interpolation
	for(int i = 0; i < pos_set_x.size(); ++i)
	{
		double shift_x = pos_set_x.at(i) - pos_x;
		double shift_y = pos_set_y.at(i) - pos_y;
		pos_set_x.at(i) = shift_x*cos(0 - angle) - shift_y*sin(0 - angle);
		pos_set_y.at(i) = shift_x*sin(0 - angle) + shift_y*cos(0 - angle);
	}

	// Path interpolation using spline function
	tk::spline spline;
	spline.set_points(pos_set_x, pos_set_y);

	double target_x = 30.0;
	double target_y = spline(target_x);
	double targetDistance = std::hypot(target_x, target_y);

	for (int i = 0; i < path_size; ++i) {
	  next_x_vals.push_back(previous_path_x.at(i));
	  next_y_vals.push_back(previous_path_y.at(i));
	}

	double x_increment = 0.0;
	for (int i = 0; i < 50-path_size; ++i)
	{
		double n = targetDistance/(.02*goal.ref_speed_mps);
		double x_point = x_increment + target_x/n;
		double y_point = spline(x_point);

		x_increment = x_point;
		double tmp_x = x_point;
		double tmp_y = y_point;

		x_point = tmp_x*cos(angle) - tmp_y*sin(angle);
		y_point = tmp_x*sin(angle) + tmp_y*cos(angle);

		x_point += pos_x;
		y_point += pos_y;

		next_x_vals.push_back(x_point);
		next_y_vals.push_back(y_point);
	}
}

/*
 * @brief Predict the positions of other vehicles on the adjacent lanes
 * @input streetVehStateSet
 * @input timePredict
 * @output predicted positions of vehicles on the same side of the street
 */
void predictVehiclePosition(const double timePredict, std::vector<StreetVehicleState>& streetVehStateSet)
{


}

/*
 * @brief Calculate the goal position for the trajectory to be generated
 * @input predicated positions of street vehicles
 * @output goal position
 */
TrajectoryGoal planBehavior(const CarState& carState, const std::vector<StreetVehicleState>& streetVehStateSet)
{
	TrajectoryGoal goal(carState.car_speed, dToLane(carState.car_d));

	// Initialize all states of adjacent lanes to be open (valid to stay in or make lane change to)
	std::vector<LaneState> laneStateSet{LaneState::LEFTLANE_OPEN,
										LaneState::MIDLANE_OPEN,
										LaneState::RIGHTLANE_OPEN};

	Lane laneEgo = dToLane(carState.car_d);
	if (laneEgo == Lane::LEFTLANE)
	{
		laneStateSet.at(0) = LaneState::LEFTLANE_CLOSED;
	}
	else if (laneEgo == Lane::RIGHTLANE)
	{
		laneStateSet.at(2) = LaneState::RIGHTLANE_CLOSED;
	}

	for(const auto& e : streetVehStateSet)
	{
		Lane laneVeh = dToLane(e.d_m);
		double collisionDistance_m = e.s_m - carState.car_s;
		if(collisionDistance_m >= 0.0 && collisionDistance_m < 8.0 && laneVeh == laneEgo)
		{
			laneState.at(static_cast<int>(laneVeh))
		}

	}

	if(goal.ref_speed_mps < (MAX_SPEED_MPS - 0.1))
	{
		goal.ref_speed_mps += 0.01;
	}
	else
	{
		goal.ref_speed_mps = MAX_SPEED_MPS;
	}
	return goal;
}

#endif  // HELPERS_H
