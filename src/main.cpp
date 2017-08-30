#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"
#include "json.hpp"
#include "spline.h"
#include <algorithm>
using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
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

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
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

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}


vector<double> JMT(vector< double> start, vector <double> end, double T)
{
    /*
    Calculate the Jerk Minimizing Trajectory that connects the initial state
    to the final state in time T.

    INPUTS

    start - the vehicles start location given as a length three array
        corresponding to initial values of [s, s_dot, s_double_dot]

    end   - the desired end state for vehicle. Like "start" this is a
        length three array.

    T     - The duration, in seconds, over which this maneuver should occur.

    OUTPUT 
    an array of length 6, each value corresponding to a coefficent in the polynomial 
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

    */
    double si, si_dot,si_d_dot, sf, sf_dot, sf_d_dot;
    
    si = start[0];
    si_dot = start[1];
    si_d_dot = start[2];
    
    sf = end[0];
    sf_dot = end[1];
    sf_d_dot = end[2];
    
    double a0,a1, a2,a3,a4,a5;
    a0 = si;
    a1 = si_dot;
    a2 = si_d_dot/2;
    
    double C1, C2, C3;
    
    C1 = si + (si*T) + (si_d_dot*0.5*pow(T,2));
    C2 = si_dot + (si_d_dot*T);
    C3 = si_d_dot;
    
    MatrixXd T_mat(3,3);
    T_mat << pow(T,3), pow(T,4), pow(T,5),
             3*pow(T,2), 4*pow(T,3), 5*pow(T,4),
             6*T, 12*pow(T,2), 20*pow(T,3);
             
    VectorXd alpha(3);
    
    VectorXd Const(3);
    
    Const << sf - (si + (si_dot*T) + (si_d_dot*0.5*pow(T,2))),
             sf_dot - (si_dot + (si_d_dot*T)),
             sf_d_dot - si_d_dot;
    
    alpha = T_mat.inverse()*Const;
    
    a3 = alpha[0];
    a4 = alpha[1];
    a5 = alpha[2];
    //std::cout << "a3, a4, a5 = " << a3 << " " << a4 << " " << a5 << std::endl;
    return {a0,a1,a2,a3,a4,a5};
    
}

tk::spline Get_spline(vector<double> x_vec, vector<double> y_vec) {

    tk::spline s;
    
    s.set_points(x_vec, y_vec);

    return s;

}

int which_lane(double d, double lane_length) {

    int lane;
    lane = d/lane_length;

    return lane;

}

vector<double> realize_keep_lane(double curr_s, double curr_d, double car_speed, map< int, vector<vector<double>>> predictions) {

    std::map<int,vector < vector<double> > >::iterator it = predictions.begin();
    int curr_lane_num = which_lane(curr_d,4);
    double following_buffer = 15; // Lets give a buffer of 3 seconds behind the car in front  
    double curr_s_t_0_5 = curr_s + car_speed;
    double end_s = curr_s_t_0_5;
    vector<vector<vector<double>>> car_infront_sd;

    double accel = 10;
    double target_velo = 49/2.24;
    while (it != predictions.end()) {

        int id = it->first;
        vector<vector<double>> c_pred = it->second;

        int c_lane_num = which_lane(c_pred[0][1], 4);
        if ((c_lane_num == curr_lane_num) and (c_pred[0][0] > curr_s) and (fabs(c_pred[0][0]-curr_s) < 100)) {

            car_infront_sd.push_back(c_pred);

        }
        
    it++;

    }
    
    double min_s = 10000000;
    int car_in_front = -1;
    for (int i=0; i<car_infront_sd.size(); i++) {

        if(car_infront_sd[i][0][0] < min_s) {

            min_s = car_infront_sd[i][0][0];
            car_in_front = i;

        }

    }
    if (car_in_front != -1) {
        double car_in_front_t_0_5_s = car_infront_sd[car_in_front][50][0];
        accel = fmin(fmax((car_in_front_t_0_5_s - curr_s_t_0_5 - following_buffer), -5), 9);
        target_velo = fmin(car_infront_sd[car_in_front][1][2], 49/2.24); // Set same as this car's velo; 
        //std::cout << "car with id " << car_in_front << " , with s = " << min_s << " velo = " << car_infront_sd[car_in_front][1][2] << " in front, I am at " << curr_s << " accel = "  << accel <<  std::endl;
        //std::cout << "car_in_front will be at " << car_in_front_t_0_5_s << " I will be at << " << curr_s_t_0_5 << std::endl;
    }
    double end_d = curr_d;// 2 + (4*which_lane(curr_d, 4));
    //std::cout << " curr_d = " << curr_d <<  " translated to end_d = " << end_d << std::endl;




    vector<double> sda;
    sda.push_back(end_s);
    sda.push_back(end_d);
    sda.push_back(accel);
    sda.push_back(target_velo);

    
    return sda;


}


vector<double> realize_change_lane(double curr_s, double curr_d, double car_speed, map< int, vector<vector<double>>> predictions, string lane_direction)
 {
    int lane_num = which_lane(curr_d, 4);
    int delta = 1; // default to right
    if (lane_direction == "L") { 

        delta = -1;
    }
    
    vector<double> sda;
    double end_d = 2 + 4*(lane_num+delta);
    //std::cout << "chaging lane towards " << lane_direction << " curr d = " << curr_d << " end_d = " << end_d << std::endl;
    sda = realize_keep_lane(curr_s, end_d, car_speed, predictions);
    
    return sda;


}
// Lets define states as follows

/*
"KL" - Keep Lane
     - The vehicle will attempt to drive its target speed, unless there is 
       traffic in front of it, in which case it will slow down.

"LCL" or "LCR" - Lane Change Left / Right
     - The vehicle will IMMEDIATELY change lanes and then follow longitudinal
       behavior for the "KL" state in the new lane.

"PLCL" or "PLCR" - Prepare for Lane Change Left / Right
     - The vehicle will find the nearest vehicle in the adjacent lane which is
       BEHIND itself and will adjust speed to try to get behind that vehicle.

*/
vector<double> realize_state(double curr_s, double curr_d, double car_speed, string state, map< int, vector<vector<double>>> predictions) {

    vector<double> sda;
    if (state == "KL") {
        sda = realize_keep_lane(curr_s, curr_d, car_speed, predictions);
    }
    else if (state == "LCL") {
        //std::cout << "Lane change left, curr_d = "  << curr_d << std::endl;
        sda = realize_change_lane(curr_s, curr_d, car_speed, predictions, "L");
    }
    else if (state == "LCR") {
        //std::cout << "Lane change right, curr_d = "  << curr_d << std::endl;
        sda = realize_change_lane(curr_s, curr_d, car_speed, predictions, "R");
    }

    return sda;


}


map<int, vector< vector<double> > > generate_predictions(vector<vector<double>> sensor_fusion, double time_horizon,  vector<double> maps_x, vector<double> maps_y) {

    
    map< int, vector<vector<double>>> predictions;


    for (int i=0; i< sensor_fusion.size(); i++ ) {

        double c_id = sensor_fusion[i][0];
        double c_x = sensor_fusion[i][1];
        double c_y = sensor_fusion[i][2];
        double c_vx = sensor_fusion[i][3];
        double c_vy = sensor_fusion[i][4];
        double c_s = sensor_fusion[i][5];
        double c_d = sensor_fusion[i][6];
        vector<vector<double>> c_pred;
        vector<double> sdv;
        double car_velo = sqrt(pow(c_vx,2) + pow(c_vy,2));
        sdv.push_back(c_s);
        sdv.push_back(c_d);
        sdv.push_back(car_velo);
        c_pred.push_back(sdv);
        double c_x_curr = c_x;
        double c_y_curr = c_y;
        for (int t=1; t<time_horizon+1; t++) {

            //double c_x_t = c_x + (c_vx*t*0.02);
            //double c_y_t = c_y + (c_vy*t*0.02);
            //double angle = atan2(c_y_t - c_y_curr, c_x_t - c_x_curr);
            //c_x_curr = c_x_t;
            //c_y_curr = c_y_t;
            //vector<double> sd = getFrenet(c_x_t, c_y_t, angle, maps_x, maps_y);
            vector<double> sdv_t;
            sdv_t.push_back(c_s + (t*0.02*car_velo));
            sdv_t.push_back(c_d);
            sdv_t.push_back(car_velo);
            
            c_pred.push_back(sdv_t);
            

        }
        predictions[c_id] = c_pred;

    }
    
    return predictions;

}

double compute_speed_cost(double v, double speed_limit, double stop_cost) {
    
    double cost;
    double buffer;
    buffer = 0; // 1 m/s buffer
    double t_speed;
    t_speed = speed_limit - buffer;
    if (v > speed_limit) {
        cost = 1;
    } 
    else if (v < t_speed) {
        cost = stop_cost*(t_speed - v)/(double)t_speed;
    }
    else {
        cost = (v - t_speed)/(double)buffer;
    }
    //std::cout << " curr_speed : " << v << " speed limit = " << speed_limit << std::endl;
    return cost;
}

double compute_lane_cost(int target_lane, int lane_num) {
    // Cost is 1 - e^ (-dd/ds)
    double cost;
    int delta_lane;
    delta_lane = target_lane-lane_num;
    cost = 1 - exp(-1*abs(delta_lane));
    //std::cout << "delta lane , delta s : " << delta_lane << ", " <<  delta_s <<std::endl;
    //std::cout << " The cost calc is " << cost << std::endl;
    return cost;
}

bool check_collision(double s, double end_d, double start_d, double car_speed, int time_ahead, map< int, vector<vector<double>>> predictions) {

    bool collided = false;
    std::map<int,vector < vector<double> > >::iterator it;
    vector<vector<double>> loc; 
    for (int t=0; t<time_ahead;t++) {
        double s_t = s + (car_speed*t*0.02);
         // Check for collision
        for (it=predictions.begin();it!=predictions.end();it++) {
        
            loc = it->second;
             
            if (((which_lane(loc[t][1],4) == which_lane(end_d,4))) && (fabs(loc[t][0]-s_t) <= 20)) {
                collided = true;
                break;
            }
            
            
        }
    if (collided == true) { break; }
    }
    return collided;
} 
string update_state(double curr_s, double curr_d, double car_speed, string curr_state, map< int, vector<vector<double>>> predictions) {


    std::vector<string> possible_states;
    
    bool collided;
    string state;
    string curr_state_dummy = curr_state;

    if (curr_state == "KL") {
        possible_states = {"KL","LCL","LCR"};
        //possible_states = {"KL","LCL","LCR","PLCL","PLCR"};
    }
    else if (curr_state == "LCL") {
        possible_states = {"KL"};
    }
    else if (curr_state == "LCR") {
        possible_states = {"KL"};
    }
    else if (curr_state == "PLCL") {
        possible_states = {"KL", "LCL"};
    }
    else if (curr_state == "PLCR") {
        possible_states = {"KL","LCR"};
    }


    vector<double> cost_vec;
    double cost;
    double min_cost;
    int min_cost_index;
    for (string test_state: possible_states) {
        vector<double> sda_temp;
        cost = 0;
        // COMMENT: trial of the possible state
        sda_temp = realize_state(curr_s, curr_d, car_speed, test_state, predictions);
        //std::cout << " size of sda temp = " << sda_temp.size() << std::endl;
        //std::cout << " target speed = " << sda_temp[3] << " for poss state "<< test_state << std::endl;
        // COMMENT: Compute the speed-based cost based on the target velocity of the target lane
        cost = cost + 2*compute_speed_cost(sda_temp[3], 50/2.24, 1.0);  
        //std::cout << " speed cost = " << cost << std::endl;
        // COMMENT: Penalize lane-change more than keep-lane so that we change lanes only when required
        if (test_state != "KL") { 
            collided = check_collision(curr_s, sda_temp[1], curr_d, car_speed, 50, predictions);
            //std::cout << " collided = " << collided << std::endl;
            if (collided == true) { cost = cost + 100000; } 
            cost = 1.1*cost; 
        } // Penalize lane change
        //cost = cost + (compute_lane_cost(0,which_lane(sda_temp[1],4)));  
        // COMMENT: Cannot go out of the three lanes
        if (sda_temp[1]<0) { cost = cost + 100000; } // Cannot cross over the yellow line
        if (sda_temp[1]>11) { cost = cost + 100000; } // Cannot go outside of the lanes
        cost_vec.push_back(cost);

    }
    //std::cout << " size of cost vec = " << cost_vec.size() << std::endl;
    // COMMENT: Pick up the state with minimum cost
    min_cost = 99999;
    min_cost_index = 0;
    for (int i = 0; i < cost_vec.size(); i++) {
        std::cout << "cost[" << i << "]: " << cost_vec[i] << endl;
        if (cost_vec[i] < min_cost) {
            min_cost = cost_vec[i];
            min_cost_index = i;
            
        }
    }
    state = possible_states[min_cost_index];  
    cost_vec.clear();
    
    //sda = realize_state(curr_s, curr_d, car_speed, state, predictions);
    
    return state;


}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  string curr_state = "KL";
  double time_horizon = 0;
  bool changing_lane;
  double end_d;
  map< string, double > new_state_avg;
  new_state_avg["KL"] = 0;
  new_state_avg["LCL"] = 0;
  new_state_avg["LCR"] = 0;
  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  double ref_vel = 0; // In m/s
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &curr_state, &time_horizon, &ref_vel, &changing_lane, &end_d, &new_state_avg](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;
                vector<double> x_vec;
                vector<double> y_vec;
                double posX, posX2;
                double posY, posY2;
                double angle;
                int prev_path_size;
                double lookahead_dist = 30; // In meters
                // First figure out if lane needs to be changed based on a slow car in front. Keep moving fwd at target velo (Keep lane) till a slow car is encountered. Then generate target trajectories based on perturbed goals on the side lane (use Quintic polynomial to calculate coeff for each init-goal pair). T will be based on the target speed. On curves limit speed.  

                prev_path_size = previous_path_x.size();

                double end_s, prev_end_d;
                vector<double> sd_prev;
                    
                if (prev_path_size <2) {

                     posX = car_x;
                     posY = car_y;
                     angle = deg2rad(car_yaw);
                     posX2 = car_x - cos(angle);
                     posY2 = car_y - sin(angle);
                     x_vec.push_back(posX2);
                     x_vec.push_back(posX);
                     y_vec.push_back(posY2);
                     y_vec.push_back(posY);
                     sd_prev = getFrenet(posX, posY, angle, map_waypoints_x, map_waypoints_y);

                }

                else {
                    posX = previous_path_x[prev_path_size-1];
                    posY = previous_path_y[prev_path_size-1];
                    posX2 = previous_path_x[prev_path_size-2];
                    posY2 = previous_path_y[prev_path_size-2];
                    angle = atan2((posY-posY2),(posX-posX2));
                    x_vec.push_back(posX2);
                    x_vec.push_back(posX);
                    y_vec.push_back(posY2);
                    y_vec.push_back(posY);
                    sd_prev = getFrenet(previous_path_x[0], previous_path_y[0], angle, map_waypoints_x, map_waypoints_y);
                    car_s = end_path_s;
                    car_d = end_path_d;
  
                }
               vector<double> sd = getFrenet(posX, posY, angle, map_waypoints_x, map_waypoints_y);
               


               // Now that we have the initial reference position, we need to figure out our 'END' S and D values based on the sensor fusion data. Is there a car in our lane?
               //bool collision;
               //int car_lane = which_lane(sd[1], 4);

               map< int, vector<vector<double>> >predictions;
               bool wait;
               predictions = generate_predictions(sensor_fusion,50, map_waypoints_x, map_waypoints_y);
               if (time_horizon == 150) {
               time_horizon = 0;
               wait = false;
               new_state_avg["KL"] = 0;
               new_state_avg["LCL"] = 0;
               new_state_avg["LCR"] = 0;
               } else {
               time_horizon = time_horizon + 1;
               wait = true;
               }

               vector<double> sda;
               string new_state = update_state(car_s,car_d, car_speed/2.24, curr_state, predictions);
               //new_state_avg[new_state] += 1;
               //std::cout << "KL count = " << new_state_avg["KL"] << std::endl;
               //std::cout << "LCL count = " << new_state_avg["LCL"] << std::endl;
               //std::cout << "LCR count = " << new_state_avg["LCR"] << std::endl;
               std::cout << "curr_state = " << curr_state << " new state = " << new_state << " end d = " << end_d <<std::endl;


               // COMMENT: Realize a lane change only if there is no ongoing lane change. If initiated set a flag that lane change is going on. 
               if (((new_state == "LCL") || (new_state == "LCR")) & (changing_lane == false)) {
                   //std::cout << " I am in 1" << std::endl;
                   sda = realize_state(car_s, car_d, car_speed/2.24, new_state, predictions);
                   end_d = 2+ (4*which_lane(sda[1],4));
                   changing_lane = true;
                   curr_state = new_state;
               
               } else if (changing_lane == true) {
                   //std::cout << " I am in 2" << " end_d = " << end_d << " car_d = " << car_d << " changing lane = " << changing_lane << std::endl;

                     sda = realize_state(car_s, end_d, car_speed/2.24, new_state, predictions);
                     if (fabs(end_d-car_d) < 0.3) { changing_lane = false; } 
               } else {
                     //std::cout << " I am in 3" << std::endl;
                     sda = realize_state(car_s, car_d, car_speed/2.24, new_state, predictions);
                     end_d = 2+ (4*which_lane(sda[1],4));
                     curr_state = new_state;
               }
               ref_vel = fmin(ref_vel + (0.02*sda[2]), 49/2.24);

               
               //std::cout << "ref vel = " << ref_vel << " angle = " << angle << " (posX, posY) = " << posX << "," << posY << " (posX2,posY2) = " << posX2 << "," << posY2 << " --  "<< (posY-posY2)/(posX-posX2) << std::endl;
               //std::cout << " s, d = " << sd[0] << ", " << sd[1] << ", " " car speed = " << car_speed << std::endl;
               //if (sda[1] != 6)  { ref_vel = 0; }
               //end_s = sda[0];
               end_s = sd[0] + 30;
               vector<double> next_wp_0 = getXY(end_s, end_d, map_waypoints_s,map_waypoints_x, map_waypoints_y);
               vector<double> next_wp_1 = getXY(end_s+30, end_d, map_waypoints_s,map_waypoints_x, map_waypoints_y);
               vector<double> next_wp_2 = getXY(end_s+60, end_d, map_waypoints_s,map_waypoints_x, map_waypoints_y);

               x_vec.push_back(next_wp_0[0]);
               x_vec.push_back(next_wp_1[0]);
               x_vec.push_back(next_wp_2[0]);

               y_vec.push_back(next_wp_0[1]);
               y_vec.push_back(next_wp_1[1]);
               y_vec.push_back(next_wp_2[1]);


                     
                tk::spline poly;
                //std::cout << " ------- " << std::endl;
                //std::cout << "angle = " << angle << std::endl;

                // Transforming all the coordinates to the car domain. So that spline does not break
                for (int i=0; i < x_vec.size(); i++) {
                    double x_new, y_new;
                    x_new = ((x_vec[i] - posX)*cos(angle))  + ((y_vec[i] - posY)*sin(angle));
                    y_new = ((y_vec[i] - posY)*cos(angle))  - ((x_vec[i] - posX)*sin(angle));
                    //std::cout << "ground -> (" << x_vec[i] << "," << y_vec[i] << ")" << ",";
                    x_vec[i] = x_new;
                    y_vec[i] = y_new;
                    //std::cout << "(" << x_vec[i] << "," << y_vec[i] << ")" << ",";
                }
                //std::cout << std::endl;
                poly = Get_spline(x_vec, y_vec);
                //double target_x = 30;
                //double target_y = poly(target_x);
                //double target_dist = sqrt(pow(target_x,2) + pow(target_y,2));
                //double N = target_dist/(0.02*ref_vel/2.24);
                //std::cout << " N = " << N << " target_x = " << target_x << " posX = " << posX << std::endl;
                //double incr = (target_x)/N;
                //double next_x = posX;
                //double curr_x = posX;
                //double next_y = posY;
                //double curr_y = posY;
                //double direc = angle;
                double next_x;
                double curr_x = 0;
                double next_y;
                double curr_y = 0;
                double direc = 0;
                double next_x_ground;
                double next_y_ground;

                for (int i=0; i < prev_path_size; i++) {
                    next_x_vals.push_back(previous_path_x[i]);
                    next_y_vals.push_back(previous_path_y[i]);
                }
               

                for (int i=1; i<50 - prev_path_size;i++) {
                    
                    double incr = 0.02*ref_vel*cos(direc); // This ensures that we maintain the reference velocity
                    next_x = curr_x + incr;
                    next_y = poly(next_x);
                    
                    direc = atan2(next_y-curr_y,next_x-curr_x);
                    curr_x = next_x;
                    curr_y = next_y;



                    // Transform the relative coord to the global coord before pushing it in
                    next_x_ground = posX + (next_x*cos(angle)) - (next_y*sin(angle));
                    next_y_ground = posY + (next_x*sin(angle)) + (next_y*cos(angle));
                    
                    //std::cout << " next_x = " << next_x << " next_y = " << next_y << " direc = " << direc << std::endl;
                    //std::cout << " next_x_g = " << next_x_ground << " next_y_g = " << next_y_ground << " direc = " <<  std::endl;
                    //vector<double> xy =  getXY(s, d, map_waypoints_s,map_waypoints_x, map_waypoints_y);
                    next_x_vals.push_back(next_x_ground);
                    next_y_vals.push_back(next_y_ground);
                }

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































