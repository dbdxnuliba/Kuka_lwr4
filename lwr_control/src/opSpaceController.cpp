#include "ros/ros.h"
#include "boost/thread.hpp"
#include "sensor_msgs/JointState.h"
#include "geometry_msgs/PoseStamped.h"
#include "geometry_msgs/TwistStamped.h"
#include "geometry_msgs/AccelStamped.h"
#include "gazebo_msgs/ContactsState.h"
#include <std_msgs/Float64.h>

#include <kdl_parser/kdl_parser.hpp>
#include <kdl/chainfksolvervel_recursive.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolverpos_nr.hpp>
#include <kdl/chainjnttojacsolver.hpp>
#include <kdl/chainjnttojacdotsolver.hpp>
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
		void interaction_wrench_cb(const gazebo_msgs::ContactsStateConstPtr&);
		void ctrl_loop();
		void compute_errors(const geometry_msgs::PoseStamped& p_des, const geometry_msgs::TwistStamped& v_des, const geometry_msgs::AccelStamped& a_des);

	private:
		ros::NodeHandle _nh;
		KDL::Tree iiwa_tree;

		KDL::ChainFkSolverPos_recursive *_fksolver; //Forward position solver
		KDL::ChainFkSolverVel_recursive *_fk_solver_pos_vel; //Forward position and velocity solver
		KDL::ChainIkSolverVel_pinv *_ik_solver_vel;   	//Inverse velocity solver
		KDL::ChainIkSolverPos_NR *_ik_solver_pos;
		KDL::ChainJntToJacSolver *_J_solver;
		KDL::ChainJntToJacDotSolver *_Jdot_solver;

		KDL::Chain _k_chain;

		ros::Subscriber _js_sub;
		ros::Subscriber _wrench_sub;
		ros::Publisher _cartpose_pub;
		ros::Publisher _plannedpose_pub;
		KDL::JntArray *_initial_q;
		KDL::JntArray *_q_in;
		KDL::JntArray *_q_in_old;
		KDL::JntArray *_dq_in;
		bool _first_js;
		bool _first_fk;
		bool _sync;
		ros::Publisher _cmd_pub[7];
		KDL::FrameVel _dirkin_out;
		KDL::Frame _p_out;
		KDL::Twist _v_out;
		KDL::ChainDynParam *_dyn_param;
		geometry_msgs::PoseStamped _pose;
		geometry_msgs::TwistStamped _vel;
		Eigen::VectorXd x_t;
		Eigen::VectorXd xDot_t;
		Eigen::VectorXd xDotDot;
		Eigen::MatrixXd _J;
		Eigen::MatrixXd _Jold;
		Eigen::MatrixXd _JDot;
		Eigen::VectorXd _gradManMeas;
		Eigen::VectorXd _extWrench;
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
	_fk_solver_pos_vel = new KDL::ChainFkSolverVel_recursive( _k_chain );
	_ik_solver_vel = new KDL::ChainIkSolverVel_pinv( _k_chain );
	_ik_solver_pos = new KDL::ChainIkSolverPos_NR( _k_chain, *_fksolver, *_ik_solver_vel, 100, 1e-6 );
	_J_solver = new KDL::ChainJntToJacSolver( _k_chain );
	_Jdot_solver = new KDL::ChainJntToJacDotSolver( _k_chain );
	//_Jdot_solver->setRepresentation(2);//INERTIAL
	//_Jdot_solver->setRepresentation(0);//HYBRID

	_q_in = new KDL::JntArray( _k_chain.getNrOfJoints() );
	_q_in_old = new KDL::JntArray( _k_chain.getNrOfJoints() );
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
	_wrench_sub = _nh.subscribe("/tool_contact_sensor_state", 0, &KUKA_INVDYN::interaction_wrench_cb, this);

	_cartpose_pub = _nh.advertise<geometry_msgs::PoseStamped>("/lwr/eef_pose", 0);
	_plannedpose_pub = _nh.advertise<geometry_msgs::PoseStamped>("/lwr/planned_pose", 0);

	_cmd_pub[0] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J1_controller/command", 0);
	_cmd_pub[1] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J2_controller/command", 0);
	_cmd_pub[2] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J3_controller/command", 0);
	_cmd_pub[3] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J4_controller/command", 0);
	_cmd_pub[4] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J5_controller/command", 0);
	_cmd_pub[5] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J6_controller/command", 0);
	_cmd_pub[6] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J7_controller/command", 0);

	x_t.resize(6);x_t.resize(6);
	xDot_t.resize(6);
	xDotDot.resize(6);
	_extWrench.resize(6);
	_extWrench = Eigen::VectorXd::Zero(6);
	_J.resize(6,_k_chain.getNrOfJoints());
	_Jold.resize(6,_k_chain.getNrOfJoints());
	_JDot.resize(6,_k_chain.getNrOfJoints());
	_J = MatrixXd::Zero(6,_k_chain.getNrOfJoints());
	_Jold = MatrixXd::Zero(6,_k_chain.getNrOfJoints());
	_gradManMeas.resize(7);
	_gradManMeas = Eigen::VectorXd::Zero(7);

	_first_js = false;
	_first_fk = false;
}

void KUKA_INVDYN::interaction_wrench_cb(const gazebo_msgs::ContactsStateConstPtr& message) {

	int nContacts=message->states.size();

	if(nContacts==0) {
		_extWrench = Eigen::VectorXd::Zero(6);
	}
	else {
		for (int i=0; i<nContacts; i++) {
			_extWrench(0)=message->states[i].total_wrench.force.x;
			_extWrench(1)=message->states[i].total_wrench.force.y;
			_extWrench(2)=message->states[i].total_wrench.force.z;
			_extWrench(3)=message->states[i].total_wrench.torque.x;
			_extWrench(4)=message->states[i].total_wrench.torque.y;
			_extWrench(5)=message->states[i].total_wrench.torque.z;
		}

	}
//cout<<_extWrench<<endl<<endl;
}

void KUKA_INVDYN::joint_states_cb( sensor_msgs::JointState js ) {

	_q_in_old->data=_q_in->data;

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
		if( !_first_js) _initial_q->data[i] = js.position[i-1];
	}

	get_dirkin();

	if(_first_js) {
		Eigen::MatrixXd man = _J*_J.transpose();
		double manMeas = sqrt(man.determinant());
		//cout<<manMeas<<endl<<endl;
		for(int i=0; i<7; i++)
			_gradManMeas(i) = manMeas/(_q_in->data[i] - _q_in_old->data[i]);
	}

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
	double Kp = 500;
	double Kd = 100;

	while( !_first_js ) usleep(0.1);

	std::vector<geometry_msgs::PoseStamped> waypoints;
	geometry_msgs::PoseStamped p=_pose;
	waypoints.push_back(_pose);

	tf::Quaternion qorient(p.pose.orientation.x,p.pose.orientation.y,p.pose.orientation.z,p.pose.orientation.w);
	tf::Matrix3x3 orient(qorient);
	double roll,pitch,yaw;
	orient.getRPY(roll,pitch,yaw);
	qorient.setRPY(roll,pitch-90,yaw);
	//p.pose.position.x -= 0.4;
	p.pose.position.y -= 0.4;
	p.pose.position.z -= 0.5;
	p.pose.orientation.z = qorient.z();
	p.pose.orientation.w = qorient.w();
	p.pose.orientation.x = qorient.x();
	p.pose.orientation.y = qorient.y();
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
		_plannedpose_pub.publish(cplanner._x[i]);
		//cout<<cplanner._x[i].pose.orientation.x<<" "<<cplanner._x[i].pose.orientation.y<<" "<<cplanner._x[i].pose.orientation.z<<" "<<cplanner._x[i].pose.orientation.w<<endl<<endl;
		compute_errors(cplanner._x[i],cplanner._xd[i],cplanner._xdd[i]); //Calcolo errori spazio operativo

		Eigen::MatrixXd J_pinv_partial = _J*_J.transpose();
		double manMeas = sqrt(J_pinv_partial.determinant());
		while(manMeas<0.01) {
			J_pinv_partial += 0.01*Eigen::MatrixXd::Identity(6,6);
			manMeas = sqrt(J_pinv_partial.determinant());
		}
		Eigen::MatrixXd J_pinv = _J.transpose()*J_pinv_partial.inverse();
		Eigen::MatrixXd JDot_pinv_partial = _JDot*_JDot.transpose();
		Eigen::MatrixXd JDot_pinv = _JDot.transpose()*JDot_pinv_partial.inverse();
		Eigen::MatrixXd Id(7,7);
		Id.setIdentity();
		Eigen::MatrixXd nullP = Id - (J_pinv*_J);
		double KN = 1;
		Eigen::VectorXd eDot_N = nullP*(-10*_dq_in->data);
		Eigen::VectorXd yNull(7);
		yNull = KN*eDot_N;// + (J_pinv*_JDot*J_pinv + JDot_pinv)*_J*(_dq_in->data);

		//cout<< _J*yNull<<endl<<endl;

  	Eigen::VectorXd y = J_pinv * ( xDotDot + Kd*xDot_t + Kp*x_t - _JDot*(_dq_in->data)) + yNull;
		_dyn_param->JntToMass(*_q_in, jsim_);
		_dyn_param->JntToCoriolis(*_q_in, *_dq_in, coriol_);
    _dyn_param->JntToGravity(*_q_in, grav_);

		Eigen::VectorXd q_out = jsim_.data * y + coriol_.data + grav_.data;
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
	KDL::JntArrayVel q_qdot(*_q_in,*_dq_in);
	_fk_solver_pos_vel->JntToCart(q_qdot, _dirkin_out);
	_p_out = _dirkin_out.GetFrame();
	_v_out = _dirkin_out.GetTwist();
	_pose.pose.position.x = _p_out.p.x();
	_pose.pose.position.y = _p_out.p.y();
	_pose.pose.position.z = _p_out.p.z();

	double qx, qy, qz, qw;
	_p_out.M.GetQuaternion( qx, qy, qz, qw);
	//tf::Quaternion quat(qx, qy, qz, qw);
	//quat.normalize()
	_pose.pose.orientation.w = qw;
	_pose.pose.orientation.x = qx;
	_pose.pose.orientation.y = qy;
	_pose.pose.orientation.z = qz;

	_vel.twist.linear.x = _v_out.vel.x();
	_vel.twist.linear.y = _v_out.vel.y();
	_vel.twist.linear.z = _v_out.vel.z();
	_vel.twist.angular.x = _v_out.rot.x();
	_vel.twist.angular.y = _v_out.rot.y();
	_vel.twist.angular.z = _v_out.rot.z();

	KDL::Jacobian Jac(_k_chain.getNrOfJoints());
	KDL::Jacobian JacDot(_k_chain.getNrOfJoints());
	if( _J_solver->JntToJac(*_q_in, Jac) != KDL::ChainJntToJacSolver::E_NOERROR )
		cout << "failing in Jacobian computation!" << endl;

	_Jold = _J;
	_J = Jac.data;
	if( _Jdot_solver->JntToJacDot(q_qdot, JacDot) != KDL::ChainJntToJacDotSolver::E_NOERROR )
		cout << "failing in JacobianDot computation!" << endl;

	_JDot = JacDot.data;
/*	for(int i=0; i<6;i++)
		for(int j=0; j<7; j++) {
			_JDot(i,j) = (_J(i,j)-_Jold(i,j))*500;
		} */

	Eigen::VectorXd vel(6);
	vel = _J*(_dq_in->data);

	_vel.twist.linear.x = vel(0);
	_vel.twist.linear.y = vel(1);
	_vel.twist.linear.z = vel(2);
	_vel.twist.angular.x = vel(3);
	_vel.twist.angular.y = vel(4);
	_vel.twist.angular.z = vel(5);

	_cartpose_pub.publish( _pose );
	_first_fk = true;
}

void KUKA_INVDYN::compute_errors(const geometry_msgs::PoseStamped& p_des, const geometry_msgs::TwistStamped& v_des, const geometry_msgs::AccelStamped& a_des) {
	x_t(0) = _pose.pose.position.x - p_des.pose.position.x;
	x_t(1) = _pose.pose.position.y - p_des.pose.position.y;
	x_t(2) = _pose.pose.position.z - p_des.pose.position.z;

	x_t = -1*x_t; //inverti segno

	tf::Quaternion qe(_pose.pose.orientation.x,_pose.pose.orientation.y,_pose.pose.orientation.z,_pose.pose.orientation.w);
	tf::Quaternion qd(p_des.pose.orientation.x,p_des.pose.orientation.y,p_des.pose.orientation.z,p_des.pose.orientation.w);
	tf::Matrix3x3 Re_tf, Rd_tf;
	Eigen::Matrix3d Re,Rd;
	Re_tf.setRotation(qe);
	Rd_tf.setRotation(qd);
	tf::matrixTFToEigen(Re_tf,Re);
	tf::matrixTFToEigen(Rd_tf,Rd);

	Eigen::Matrix3d Rerr = Re.transpose()*Rd;
	tf::Matrix3x3 Rerr_tf;
	tf::matrixEigenToTF(Rerr,Rerr_tf);
	tf::Quaternion qerr;
	Rerr_tf.getRotation(qerr);
	double angle = qerr.getAngle();
	tf::Vector3 axis = qerr.getAxis();
	Eigen::Vector3d eps;
	tf::vectorTFToEigen(axis,eps);
	eps = Re*(sin(angle/2.0)*eps);
	x_t(3) = eps(0);x_t(4) = eps(1);x_t(5) = eps(2);

	//cout<<x_t<<endl<<endl;

	xDot_t(0) = _vel.twist.linear.x - v_des.twist.linear.x;
	xDot_t(1) = _vel.twist.linear.y - v_des.twist.linear.y;
	xDot_t(2) = _vel.twist.linear.z - v_des.twist.linear.z;
	xDot_t(3) = _vel.twist.angular.x - v_des.twist.angular.x;
	xDot_t(4) = _vel.twist.angular.y - v_des.twist.angular.y;
	xDot_t(5) = _vel.twist.angular.z - v_des.twist.angular.z;
/*
	xDot_t(0) += v_des.twist.linear.x;
	xDot_t(1) +=v_des.twist.linear.y;
	xDot_t(2) += v_des.twist.linear.z;
	xDot_t(3) += v_des.twist.angular.x;
	xDot_t(4) +=  v_des.twist.angular.y;
	xDot_t(5) +=v_des.twist.angular.z;
*/
	xDot_t = -1*xDot_t; //inverti segno

	xDotDot(0) = a_des.accel.linear.x;
	xDotDot(1) = a_des.accel.linear.y;
	xDotDot(2) = a_des.accel.linear.z;
	xDotDot(3) = a_des.accel.angular.x;
	xDotDot(4) = a_des.accel.angular.y;
	xDotDot(5) = a_des.accel.angular.z;

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
