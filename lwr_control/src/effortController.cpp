#include "ros/ros.h"
#include "boost/thread.hpp"
#include "sensor_msgs/JointState.h"
#include "geometry_msgs/PoseStamped.h"
#include <std_msgs/Float64.h>

#include <kdl_parser/kdl_parser.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolverpos_nr.hpp>
#include <kdl/chaindynparam.hpp>

#include "../include/lwr_control/planner.h"

using namespace std;


class KUKA_INVDYN {
	public:
		KUKA_INVDYN();
		void run();
		bool init_robot_model();
		void get_dirkin();

		void joint_states_cb( sensor_msgs::JointState );
		void ctrl_loop();


	private:
		ros::NodeHandle _nh;
		KDL::Tree iiwa_tree;

		KDL::ChainFkSolverPos_recursive *_fksolver; //Forward position solver
		KDL::ChainIkSolverVel_pinv *_ik_solver_vel;   	//Inverse velocity solver
		KDL::ChainIkSolverPos_NR *_ik_solver_pos;

		KDL::Chain _k_chain;

		ros::Subscriber _js_sub;
		ros::Publisher _cartpose_pub;
		KDL::JntArray *_initial_q;
		KDL::JntArray *_q_in;
		KDL::JntArray *_dq_in;
		bool _first_js;
		bool _first_fk;
		bool _sync;
		ros::Publisher _cmd_pub[7];
		KDL::	Frame _p_out;
		KDL::ChainDynParam *_dyn_param;
		geometry_msgs::PoseStamped _pose;
};


bool KUKA_INVDYN::init_robot_model() {
	std::string robot_desc_string;
	_nh.param("robot_description", robot_desc_string, std::string());
	if (!kdl_parser::treeFromString(robot_desc_string, iiwa_tree)){
		ROS_ERROR("Failed to construct kdl tree");
		return false;
	}
	else {
		ROS_INFO("Robot tree found!");
	}

	std::string base_link = "lwr_base_link";
	std::string tip_link  = "tool_link_ee";
	if ( !iiwa_tree.getChain(base_link, tip_link, _k_chain) ) return false;

	_fksolver = new KDL::ChainFkSolverPos_recursive( _k_chain );
	_ik_solver_vel = new KDL::ChainIkSolverVel_pinv( _k_chain );
	_ik_solver_pos = new KDL::ChainIkSolverPos_NR( _k_chain, *_fksolver, *_ik_solver_vel, 100, 1e-6 );

	_q_in = new KDL::JntArray( _k_chain.getNrOfJoints() );
	_dq_in = new KDL::JntArray( _k_chain.getNrOfJoints() );
	_initial_q = new KDL::JntArray( _k_chain.getNrOfJoints() );
	_dyn_param = new KDL::ChainDynParam(_k_chain,KDL::Vector(0,0,-9.81));

	return true;
}


KUKA_INVDYN::KUKA_INVDYN() {

	if (!init_robot_model()) exit(1);
	ROS_INFO("Robot tree correctly loaded from parameter server!");

	cout << "Joints and segments: " << iiwa_tree.getNrOfJoints() << " - " << iiwa_tree.getNrOfSegments() << endl;

	_js_sub = _nh.subscribe("/lwr/joint_states", 0, &KUKA_INVDYN::joint_states_cb, this);

	_cartpose_pub = _nh.advertise<geometry_msgs::PoseStamped>("/lwr/eef_pose", 0);

	_cmd_pub[0] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J1_controller/command", 0);
	_cmd_pub[1] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J2_controller/command", 0);
	_cmd_pub[2] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J3_controller/command", 0);
	_cmd_pub[3] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J4_controller/command", 0);
	_cmd_pub[4] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J5_controller/command", 0);
	_cmd_pub[5] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J6_controller/command", 0);
	_cmd_pub[6] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J7_controller/command", 0);

	_first_js = false;
	_first_fk = false;
}


void KUKA_INVDYN::joint_states_cb( sensor_msgs::JointState js ) {

	for(int i=0; i<2; i++ ) {
		_q_in->data[i] = js.position[i];
		_dq_in->data[i] = js.velocity[i];
		if( !_first_js ) {
			_initial_q->data[i] = js.position[i];
		}
	}

	_q_in->data[2] = js.position[6];
	_dq_in->data[2] = js.velocity[6];
	if( !_first_js ) _initial_q->data[2] = js.position[6];

	for(int i=3; i<7; i++ ) {
		_q_in->data[i] = js.position[i-1];
		_dq_in->data[i] = js.velocity[i-1];
		if( !_first_js ) _initial_q->data[i] = js.position[i-1];
	}

	get_dirkin();

	_first_js = true;
	_sync = true;
}

void KUKA_INVDYN::ctrl_loop() {

	std_msgs::Float64 cmd[7];
  KDL::JntArray coriol_(7);
  KDL::JntArray grav_(7);
	KDL::JntArray q_out(_k_chain.getNrOfJoints());
	KDL::JntArray qd_out(_k_chain.getNrOfJoints());

  KDL::JntSpaceInertiaMatrix jsim_;
  jsim_.resize(_k_chain.getNrOfJoints());

	ros::Rate r(500);
	double Kp = 50;
	double Kd = 20;

	while( !_first_js ) usleep(0.1);

	std::vector<geometry_msgs::PoseStamped> waypoints;
	geometry_msgs::PoseStamped p=_pose;
	waypoints.push_back(_pose);
	p.pose.position.x -= 0.2;
	p.pose.position.y -= 0.2;
	p.pose.position.z -= 0.2;
	waypoints.push_back(p);
	waypoints.push_back(_pose);

	std::vector<double> times;
	times.push_back(0);
	times.push_back(5);
	times.push_back(10);

	CARTESIAN_PLANNER	cplanner(500);
	cplanner.set_waypoints(waypoints,times);
	cplanner.compute();

//	cout<<_pose.pose.orientation.x<<" "<<_pose.pose.orientation.y<<" "<<_pose.pose.orientation.z<<" "<<_pose.pose.orientation.w<<endl;
//	cout<<cplanner._x[0].pose.orientation.x<<" "<<cplanner._x[0].pose.orientation.y<<" "<<cplanner._x[0].pose.orientation.z<<" "<<cplanner._x[0].pose.orientation.w<<endl<<endl;
	int i = 0;
	while( ros::ok() ) {

    while( !_sync ) usleep(0.1);
	//	cout<<cplanner._x[i].pose.orientation.x<<" "<<cplanner._x[i].pose.orientation.y<<" "<<cplanner._x[i].pose.orientation.z<<" "<<cplanner._x[i].pose.orientation.w<<endl<<endl;

		KDL::Rotation Rot;
		Rot.Quaternion(cplanner._x[i].pose.orientation.x,cplanner._x[i].pose.orientation.y,cplanner._x[i].pose.orientation.z,cplanner._x[i].pose.orientation.w);
		KDL::Vector pos(cplanner._x[i].pose.position.x,cplanner._x[i].pose.position.y,cplanner._x[i].pose.position.z);
		KDL::Frame F_dest(Rot,pos);
		//cout<<pos.data<<endl<<endl;

		KDL::Vector pvel(cplanner._xd[i].twist.linear.x,cplanner._xd[i].twist.linear.y,cplanner._xd[i].twist.linear.z);
		KDL::Vector avel(cplanner._xd[i].twist.angular.x,cplanner._xd[i].twist.angular.y,cplanner._xd[i].twist.angular.z);
		KDL::Twist des_vel(pvel,avel);

		if( _ik_solver_pos->CartToJnt(*_q_in, F_dest, q_out) != KDL::SolverI::E_NOERROR )
			cout << "failing in ik!" << endl;

		if( _ik_solver_vel->CartToJnt(*_q_in, des_vel, qd_out) != KDL::SolverI::E_NOERROR )
			cout << "failing in ik vel!" << endl;

  	Eigen::VectorXd e = q_out.data - _q_in->data; //Position Error
    Eigen::VectorXd de = qd_out.data - _dq_in->data; //Desired velocity
		_dyn_param->JntToMass(*_q_in, jsim_);
		_dyn_param->JntToCoriolis(*_q_in, *_dq_in, coriol_);
    _dyn_param->JntToGravity(*_q_in, grav_);

		Eigen::VectorXd q_out = jsim_.data * ( Kd*de + Kp*e ) + coriol_.data + grav_.data;
		//Eigen::VectorXd q_out = grav_.data + coriol_.data;

		for(int i=0; i<7; i++ ) {
			cmd[i].data = q_out(i);
      //cout<<"Joint "<<i<<" :"<<q_out(i)<<endl;
		}
		for(int i=0; i<7; i++ ) {
			_cmd_pub[i].publish( cmd[i] );
		}

		_sync = false;

		if(i<(cplanner._x.size()-1)) i++;

		r.sleep();
	}

}

void KUKA_INVDYN::get_dirkin() {

	_fksolver->JntToCart(*_q_in, _p_out);
	_pose.pose.position.x = _p_out.p.x();
	_pose.pose.position.y = _p_out.p.y();
	_pose.pose.position.z = _p_out.p.z();

	double qx, qy, qz, qw;
	_p_out.M.GetQuaternion( qx, qy, qz, qw);
	_pose.pose.orientation.w = qw;
	_pose.pose.orientation.x = qx;
	_pose.pose.orientation.y = qy;
	_pose.pose.orientation.x = qz;

	_cartpose_pub.publish( _pose );
	_first_fk = true;

}


void KUKA_INVDYN::run() {


	boost::thread ctrl_loop_t ( &KUKA_INVDYN::ctrl_loop, this);
	ros::spin();

}


int main(int argc, char** argv) {
	ros::init(argc, argv, "iiwa_kdl");
	KUKA_INVDYN ik;
	ik.run();

	return 0;
}
