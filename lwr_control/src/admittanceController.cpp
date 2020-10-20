#include "ros/ros.h"
#include "boost/thread.hpp"
#include "sensor_msgs/JointState.h"
#include "geometry_msgs/PoseStamped.h"
#include "geometry_msgs/TwistStamped.h"
#include "geometry_msgs/AccelStamped.h"
#include "gazebo_msgs/ContactsState.h"
#include <std_msgs/Float64.h>
#include <std_msgs/Float64MultiArray.h>

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

class ETank {
	public:
		ETank(double Einit, double Emin, double Emax, double dt) {_Et=Einit;_Emin=Emin; _Emax=Emax; _xt=sqrt(2*_Et);_dt=dt;};
		void update(const Eigen::VectorXd input, const Eigen::MatrixXd Kd, const Eigen::VectorXd x_dot, const Eigen::VectorXd Kpxtilde);
		double getEt() {return _Et;};
		double _alpha;
	private:
		double _Et, _Emin, _Emax, _xt;
		double _beta;
		double _dt;
};

void ETank::update(const Eigen::VectorXd input, const Eigen::MatrixXd Kd, const Eigen::VectorXd x_dot, const Eigen::VectorXd Kpxtilde) {
	if(_Et<=_Emax) _beta=1;
	else _beta=0;

	double w = input.dot(Kpxtilde);

	double f_energy = 0.5*( 1 - cos(M_PI*(_Et-_Emin)/(_Emax-_Emin)) );
	double g_input;
	if(w<=0) g_input=0;
	else	g_input=1;

	_alpha = f_energy*g_input + (1-g_input);

	double gamma;
	if( (_Et>=_Emin) && w>=0 ) gamma=_alpha;
	else gamma=0;

	w=gamma*w;

	double eta=0.5;
	double Diss = x_dot.transpose()*Kd*x_dot;

	double ut = -(1.0/_xt)*w;
	cout<<"Ut: "<<ut<<endl;
	double xt_dot = (_beta*eta/_xt)*Diss + ut;
	_xt += xt_dot*_dt;
	_Et = 0.5*_xt*_xt;

}

class ETankGen {
	public:
		ETankGen(double Einit, double Emin, double Emax, double dt, int inputSize) {_Et=Einit;_Emin=Emin; _Emax=Emax; _xt=sqrt(2*_Et);_dt=dt;_alpha.resize(inputSize);};
		void update(const std::vector<Eigen::VectorXd> inputs, const std::vector<double> dissInputs, const std::vector<Eigen::VectorXd> products);
		double getEt() {return _Et;};
		std::vector<double> _alpha;
	private:
		double _Et, _Emin, _Emax, _xt;
		double _beta;
		double _dt;
};

void ETankGen::update(const std::vector<Eigen::VectorXd> inputs, const std::vector<double> dissInputs, const std::vector<Eigen::VectorXd> products) {
	if(_Et<=_Emax) _beta=1;
	else _beta=0;

	double f_energy = 0.5*( 1 - cos(M_PI*(_Et-_Emin)/(_Emax-_Emin)) );
	double eta=0.5;
	double Diss=0;
	double wtot = 0;

	for (int i=0; i<_alpha.size(); i++)
		Diss+=dissInputs[i];

	for (int i=0; i<_alpha.size(); i++) {
		double w = inputs[i].dot(products[i]);
		double g_input;
		if(w<=0) g_input=0;
		else	g_input=1;

		_alpha[i] = f_energy*g_input + (1-g_input);
		double gamma;
		if( (_Et>=_Emin) && w>=0 ) gamma=_alpha[i];
		else gamma=0;

		wtot+=gamma*w;
	}

	double ut = -(1.0/_xt)*wtot;
	cout<<"Ut: "<<ut<<endl;
	double xt_dot = (_beta*eta/_xt)*Diss + ut;
	_xt += xt_dot*_dt;
	_Et = 0.5*_xt*_xt;
	if(ut>0)
		ROS_ERROR("Diss: %f\n",Diss);
}



class DERIV {
	public:
		DERIV(double freq=500,double gain=100) {_f=freq;_dt=1.0/_f;_integral=Eigen::VectorXd::Zero(6);_gain=gain;};
		void update(Eigen::VectorXd x) {_xd=_gain*(x-_integral); _integral+=_gain*_dt*(x-_integral);};

		Eigen::VectorXd _xd;
	private:
		double _f,_dt;
		Eigen::VectorXd _integral;
		double _gain;

};

class KUKA_INVDYN {
	public:
		KUKA_INVDYN();
		void run();
		bool init_robot_model();
		void get_dirkin();

		void joint_states_cb( sensor_msgs::JointState );
		void interaction_wrench_cb(const gazebo_msgs::ContactsStateConstPtr&);
		void ctrl_loop();
		void compute_force_errors(const Eigen::VectorXd h, const Eigen::VectorXd hdot, const Eigen::VectorXd mask);
		void compute_errors(const geometry_msgs::PoseStamped& p_des, const geometry_msgs::TwistStamped& v_des, const geometry_msgs::AccelStamped& a_des);
		void compute_compliantFrame(const geometry_msgs::PoseStamped& p_des, const geometry_msgs::TwistStamped& v_des, const geometry_msgs::AccelStamped& a_des);
		void compute_compliantFrame(const geometry_msgs::PoseStamped& p_des, const geometry_msgs::TwistStamped& v_des, const geometry_msgs::AccelStamped& a_des, const std::vector<double> alpha);
		bool newTrajectory(const std::vector<geometry_msgs::PoseStamped> waypoints, const std::vector<double> times);
		bool newTrajectory(const std::vector<geometry_msgs::PoseStamped> waypoints, const std::vector<double> times, const Eigen::VectorXd xdi, const Eigen::VectorXd xdf, const Eigen::VectorXd xddi, const Eigen::VectorXd xddf);
		bool newForceTrajectory(const std::vector<Eigen::VectorXd> waypoints, const std::vector<double> times, const Eigen::VectorXd mask);
		bool getPose(geometry_msgs::PoseStamped& p_des);
		bool getDesPose(geometry_msgs::PoseStamped& p_des);
		bool getWrench(Eigen::VectorXd& _wrench);
		bool robotReady() {return _first_fk;};
		void exitForceControl() {_fControl=false;};
	private:
		void updatePose();
		void updateForce();
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
		ros::Publisher _plannedpose_pub,_plannedtwist_pub,_plannedacc_pub,_plannedwrench_pub;
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
		Eigen::VectorXd _acc;
		Eigen::VectorXd x_t;
		Eigen::VectorXd xDot_t;
		Eigen::VectorXd xDotDot;
		Eigen::MatrixXd _J;
		Eigen::MatrixXd _Jold;
		Eigen::MatrixXd _JDot;
		Eigen::VectorXd _gradManMeas;
		Eigen::VectorXd _extWrench;
		Eigen::VectorXd z_t,zDot_t,zDotDot_t;
		geometry_msgs::PoseStamped _complPose;
		geometry_msgs::TwistStamped _complVel;
		geometry_msgs::AccelStamped _complAcc;
		geometry_msgs::PoseStamped _desPose;
		geometry_msgs::TwistStamped _desVel;
		geometry_msgs::AccelStamped _desAcc;
		bool _fControl;
		bool _trajEnd;
		bool _newPosReady;
		geometry_msgs::PoseStamped _nextdesPose;
		geometry_msgs::TwistStamped _nextdesVel;
		geometry_msgs::AccelStamped _nextdesAcc;
		Eigen::MatrixXd _Mt;
		Eigen::MatrixXd _Kdt;
		Eigen::MatrixXd _Kpt;
		Eigen::VectorXd xf,xf_dot,xf_dotdot;
		Eigen::VectorXd _h_des,_hdot_des, _nexth_des,_nexthdot_des, _forceMask;
		DERIV numericAcc;
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
	_plannedtwist_pub = _nh.advertise<geometry_msgs::TwistStamped>("/lwr/planned_twist", 0);
	_plannedacc_pub = _nh.advertise<geometry_msgs::AccelStamped>("/lwr/planned_acc", 0);
	_plannedwrench_pub = _nh.advertise<std_msgs::Float64>("/lwr/planned_wrench", 0);

	_cmd_pub[0] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J1_controller/command", 0);
	_cmd_pub[1] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J2_controller/command", 0);
	_cmd_pub[2] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J3_controller/command", 0);
	_cmd_pub[3] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J4_controller/command", 0);
	_cmd_pub[4] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J5_controller/command", 0);
	_cmd_pub[5] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J6_controller/command", 0);
	_cmd_pub[6] = _nh.advertise< std_msgs::Float64 > ("/lwr/JointEffortController_J7_controller/command", 0);

	x_t.resize(6);
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

	z_t.resize(6); z_t=Eigen::VectorXd::Zero(6);
	zDot_t.resize(6); zDot_t=Eigen::VectorXd::Zero(6);
	zDotDot_t.resize(6); zDotDot_t=Eigen::VectorXd::Zero(6);

	xf.resize(7); xf=Eigen::VectorXd::Zero(7);
	xf_dot.resize(6); xf_dot=Eigen::VectorXd::Zero(6);
	xf_dotdot.resize(6); xf_dotdot=Eigen::VectorXd::Zero(6);

	_Mt = 7*Eigen::MatrixXd::Identity(6,6);
	_Mt.bottomRightCorner(3,3) = 70*Eigen::MatrixXd::Identity(3,3);
	_Kdt = 70*Eigen::MatrixXd::Identity(6,6);
	_Kpt = 100*Eigen::MatrixXd::Identity(6,6);
	_Kpt.bottomRightCorner(3,3) = 1000*Eigen::MatrixXd::Identity(3,3);

	_h_des.resize(6); _h_des=Eigen::VectorXd::Zero(6);
	_hdot_des.resize(6); _hdot_des=Eigen::VectorXd::Zero(6);
	_acc.resize(6);_acc=Eigen::VectorXd::Zero(6);

	_first_js = false;
	_first_fk = false;
	_fControl = false;
	_trajEnd = true;
	_newPosReady = false;
}

bool KUKA_INVDYN::getPose(geometry_msgs::PoseStamped& p_des) {
	if(!_first_fk) return false;

	p_des = _pose;
	return true;
}

bool KUKA_INVDYN::getDesPose(geometry_msgs::PoseStamped& p_des) {
	if(!_first_fk) return false;

	p_des = _desPose;
	return true;
}

bool KUKA_INVDYN::getWrench(Eigen::VectorXd& wrench) {
	if(!_first_fk) return false;

	wrench = _extWrench;
	return true;
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
		tf::Quaternion qe(_pose.pose.orientation.x,_pose.pose.orientation.y,_pose.pose.orientation.z,_pose.pose.orientation.w);
		tf::Matrix3x3 Re_tf;
		Eigen::Matrix3d Re;
		Re_tf.setRotation(qe);
		tf::matrixTFToEigen(Re_tf,Re);
		_extWrench.head(3) = Re*_extWrench.head(3);
		_extWrench.tail(3) = Re*_extWrench.tail(3);
		//cout<<_extWrench<<endl<<endl;
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
	double Kp = 700;
	double Kd = 100;

	//ETank tank(1.0,0.01,1.0,0.002);
	ETankGen tankGen(2.0,0.01,2.0,0.002,2);

	while( !_first_js ) usleep(0.1);

	while( ros::ok() ) {

    while( !_sync ) usleep(0.1);

		if(_fControl) {
			updateForce();
			compute_force_errors(_h_des, _hdot_des,_forceMask);
			_newPosReady=true;
		}
		updatePose();

		/*
		Eigen::VectorXd desVelEigen, desAccEigen, complVelEigen;
		twist2Vector(_desVel,desVelEigen);
		accel2Vector(_desAcc,desAccEigen);
		twist2Vector(_complVel,complVelEigen);
		std::vector<Eigen::VectorXd> tankInputs, tankProds;
		std::vector<double> tankDiss;
		tankDiss.push_back(complVelEigen.transpose()*_Kdt*complVelEigen);
		tankInputs.push_back(-desVelEigen);
		tankProds.push_back(_Kpt*z_t);
		tankInputs.push_back(_Mt*desAccEigen + _Kdt*desVelEigen);
		tankProds.push_back(complVelEigen);
		tankGen.update(tankInputs,tankDiss,tankProds);
		cout<<tankGen.getEt()<<endl;
		*/
		//compute_compliantFrame(_desPose,_desVel,_desAcc,tankGen._alpha);
		compute_compliantFrame(_desPose,_desVel,_desAcc);
		compute_errors(_complPose,_complVel,_complAcc); //Calcolo errori spazio operativo

		_plannedpose_pub.publish(_complPose);
		_plannedtwist_pub.publish(_complVel);
		_plannedacc_pub.publish(_complAcc);

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
		Eigen::VectorXd eDot_N = nullP*(-100*_dq_in->data);
		Eigen::VectorXd yNull(7);
		yNull = KN*eDot_N;// + (J_pinv*_JDot*J_pinv + JDot_pinv)*_J*(_dq_in->data);

		//cout<< _J*yNull<<endl<<endl;

  	Eigen::VectorXd y = J_pinv * ( xDotDot + Kd*xDot_t + Kp*x_t - _JDot*(_dq_in->data)) + yNull;
		_dyn_param->JntToMass(*_q_in, jsim_);
		_dyn_param->JntToCoriolis(*_q_in, *_dq_in, coriol_);
    _dyn_param->JntToGravity(*_q_in, grav_);

		Eigen::VectorXd q_out = jsim_.data * y + coriol_.data + grav_.data;

		for(int i=0; i<7; i++ ) {
			cmd[i].data = q_out(i);
		}
		for(int i=0; i<7; i++ ) {
			_cmd_pub[i].publish( cmd[i] );
		}

		_sync = false;

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

	if(!_first_fk) _desPose = _pose;

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
	numericAcc.update(vel);
	_acc = numericAcc._xd;

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

void KUKA_INVDYN::compute_compliantFrame(const geometry_msgs::PoseStamped& p_des, const geometry_msgs::TwistStamped& v_des, const geometry_msgs::AccelStamped& a_des, const std::vector<double> alpha) {
	geometry_msgs::TwistStamped vmod_des;
	geometry_msgs::AccelStamped amod_des;

	if (alpha.size()<2) {
		vmod_des.twist.linear.x = alpha[0]*v_des.twist.linear.x;
		vmod_des.twist.linear.y = alpha[0]*v_des.twist.linear.y;
		vmod_des.twist.linear.z = alpha[0]*v_des.twist.linear.z;
		vmod_des.twist.angular.x = alpha[0]*v_des.twist.angular.x;
		vmod_des.twist.angular.y = alpha[0]*v_des.twist.angular.y;
		vmod_des.twist.angular.z = alpha[0]*v_des.twist.angular.z;
	}
	else {
		double alphaMin = min(alpha[0],alpha[1]);
		vmod_des.twist.linear.x = alphaMin*v_des.twist.linear.x;
		vmod_des.twist.linear.y = alphaMin*v_des.twist.linear.y;
		vmod_des.twist.linear.z = alphaMin*v_des.twist.linear.z;
		vmod_des.twist.angular.x = alphaMin*v_des.twist.angular.x;
		vmod_des.twist.angular.y = alphaMin*v_des.twist.angular.y;
		vmod_des.twist.angular.z = alphaMin*v_des.twist.angular.z;

		amod_des.accel.linear.x = alpha[1]*a_des.accel.linear.x;
		amod_des.accel.linear.y = alpha[1]*a_des.accel.linear.y;
		amod_des.accel.linear.z = alpha[1]*a_des.accel.linear.z;
		amod_des.accel.angular.x = alpha[1]*a_des.accel.angular.x;
		amod_des.accel.angular.y = alpha[1]*a_des.accel.angular.y;
		amod_des.accel.angular.z = alpha[1]*a_des.accel.angular.z;
		cout<<"alpha2: "<<alpha[1]<<endl;
	}

	compute_compliantFrame(p_des,vmod_des,amod_des);
}

void KUKA_INVDYN::compute_compliantFrame(const geometry_msgs::PoseStamped& p_des, const geometry_msgs::TwistStamped& v_des, const geometry_msgs::AccelStamped& a_des) {

	zDotDot_t = _Mt.inverse() * ( _extWrench - _Kdt*zDot_t - _Kpt*z_t);
	zDot_t += zDotDot_t*0.002;
	z_t += zDot_t*0.002;

	_complAcc.accel.linear.x = a_des.accel.linear.x + zDotDot_t(0);
	_complAcc.accel.linear.y = a_des.accel.linear.y + zDotDot_t(1);
	_complAcc.accel.linear.z = a_des.accel.linear.z + zDotDot_t(2);
	_complAcc.accel.angular.x = a_des.accel.angular.x + zDotDot_t(3);
	_complAcc.accel.angular.y = a_des.accel.angular.y + zDotDot_t(4);
	_complAcc.accel.angular.z = a_des.accel.angular.z + zDotDot_t(5);

	_complVel.twist.linear.x = v_des.twist.linear.x + zDot_t(0);
	_complVel.twist.linear.y = v_des.twist.linear.y + zDot_t(1);
	_complVel.twist.linear.z = v_des.twist.linear.z + zDot_t(2);
	_complVel.twist.angular.x = v_des.twist.angular.x + zDot_t(3);
	_complVel.twist.angular.y = v_des.twist.angular.y + zDot_t(4);
	_complVel.twist.angular.z = v_des.twist.angular.z + zDot_t(5);

	_complPose.pose.position.x = p_des.pose.position.x + z_t(0);
	_complPose.pose.position.y = p_des.pose.position.y + z_t(1);
	_complPose.pose.position.z = p_des.pose.position.z + z_t(2);

	tf::Quaternion qe(p_des.pose.orientation.x,p_des.pose.orientation.y,p_des.pose.orientation.z,p_des.pose.orientation.w);
	tf::Quaternion qd(p_des.pose.orientation.x,p_des.pose.orientation.y,p_des.pose.orientation.z,p_des.pose.orientation.w);
	tf::Matrix3x3 Re_tf, Rd_tf;
	Eigen::Matrix3d Re,Rd;
	Re_tf.setRotation(qe);
	//Rd_tf.setRotation(qd);
	tf::matrixTFToEigen(Re_tf,Re);
	//	tf::matrixTFToEigen(Rd_tf,Rd);
	Eigen::Vector3d eps;
	eps << z_t(3),z_t(4),z_t(5);
	eps = Re.transpose()*eps;
	double eta = sqrt(1-eps(0)*eps(0)-eps(1)*eps(1)-eps(2)*eps(2));
	if(eta>1) eta=1;
	else if (eta<-1) eta=-1;
	double theta = 2*acos(eta);
	if(theta!=0) { //qd actually different from qe
		Eigen::Vector3d axis = (1.0/sin(theta*0.5))*eps;
		tf::Vector3 axis_tf;
		tf::vectorEigenToTF(axis,axis_tf);
		tf::Quaternion qerr(axis_tf,theta);
		tf::Matrix3x3 Rerr_tf(qerr);
		Eigen::Matrix3d Rerr;
		tf::matrixTFToEigen(Rerr_tf,Rerr);
		Rd = Re*Rerr;
		tf::matrixEigenToTF(Rd,Rd_tf);
		Rd_tf.getRotation(qd);
	}

	_complPose.pose.orientation.x = qd.x();
	_complPose.pose.orientation.y = qd.y();
	_complPose.pose.orientation.z = qd.z();
	_complPose.pose.orientation.w = qd.w();

}

void KUKA_INVDYN::updatePose() {

	if(!_newPosReady)
		return;

	_desPose=_nextdesPose;
	_desVel = _nextdesVel;
	_desAcc = _nextdesAcc;

	_newPosReady = false;
}

void KUKA_INVDYN::updateForce() {
	if(!_newPosReady) return;

	_h_des=_nexth_des;
	_hdot_des = _nexthdot_des;
	_newPosReady = false;
}

bool KUKA_INVDYN::newTrajectory(const std::vector<geometry_msgs::PoseStamped> waypoints, const std::vector<double> times, const Eigen::VectorXd xdi, const Eigen::VectorXd xdf, const Eigen::VectorXd xddi, const Eigen::VectorXd xddf) {
	if(!_trajEnd) return false;

	_trajEnd=false;
	CARTESIAN_PLANNER	cplanner(500);
	cplanner.set_waypoints(waypoints,times,xdi,xdf,xddi,xddf);
	cplanner.compute();

	_fControl = false;

	while(cplanner.isReady() && ros::ok()) {
		while(_newPosReady && ros::ok()) usleep(1);
		cplanner.getNext(_nextdesPose,_nextdesVel,_nextdesAcc);
		_newPosReady=true;
	}

	_trajEnd=true;
	return true;
}

bool KUKA_INVDYN::newTrajectory(const std::vector<geometry_msgs::PoseStamped> waypoints, const std::vector<double> times) {
	Eigen::VectorXd dummy(6);
	dummy = Eigen::VectorXd::Zero(6);
	newTrajectory(waypoints,times,dummy,dummy,dummy,dummy);
}

bool KUKA_INVDYN::newForceTrajectory(const std::vector<Eigen::VectorXd> waypoints, const std::vector<double> times, const Eigen::VectorXd mask) {

	if(!_trajEnd) return false;
	_trajEnd=false;

	_forceMask = mask;
	SPLINE_PLANNER* w[6];

	for(int i=0; i<6; i++) {
		w[i] = new SPLINE_PLANNER(500);
		std::vector<double> componentPoints;
		for(int j=0; j<waypoints.size(); j++) {
			componentPoints.push_back(waypoints[j](i));
		}
		w[i]->set_waypoints(componentPoints,times);
		w[i]->compute_traj();
	}

	xf(0) = _desPose.pose.position.x;
	xf(1) = _desPose.pose.position.y;
	xf(2) = _desPose.pose.position.z;
	xf(3) = _desPose.pose.orientation.x;
	xf(4) = _desPose.pose.orientation.y;
	xf(5) = _desPose.pose.orientation.z;
	xf(6) = _desPose.pose.orientation.w;
	xf_dot(0) = _desVel.twist.linear.x;
	xf_dot(1) = _desVel.twist.linear.y;
	xf_dot(2) = _desVel.twist.linear.z;
	xf_dot(3) = _desVel.twist.angular.x;
	xf_dot(4) = _desVel.twist.angular.y;
	xf_dot(5) = _desVel.twist.angular.z;

	_fControl=true;

	while(w[0]->isReady() && ros::ok()) {
		Eigen::VectorXd h(6), hdot(6);
		for(int i=0; i<6; i++) {
			double f, fdot, fdotdot;
			w[i]->getNext(f, fdot, fdotdot);
			h(i) = f;
			hdot(i) = fdot;
		}
		while(_newPosReady && ros::ok()) usleep(1);
		_nexth_des=h;
		_nexthdot_des=hdot;
		_newPosReady=true;
	}

	for(int i=0; i<6; i++)
		delete w[i];

	_trajEnd=true;
	return true;
}

void KUKA_INVDYN::compute_force_errors(const Eigen::VectorXd h, const Eigen::VectorXd hdot, const Eigen::VectorXd mask) {

	Eigen::VectorXd vel(6),ht(6);
	double Kh = 10.0;

	ht = _extWrench-h;
	std_msgs::Float64 data;
	data.data=ht(1);
	_plannedwrench_pub.publish(data);
	//cout<<"Error: "<<ht(1)<<" / "<<_extWrench(1)<<endl<<endl;
	//cout<<ht(1)<<endl<<endl;
	vel(0) = _desVel.twist.linear.x;
	vel(1) = _desVel.twist.linear.y;
	vel(2) = _desVel.twist.linear.z;
	vel(3) = _desVel.twist.angular.x;
	vel(4) = _desVel.twist.angular.y;
	vel(5) = _desVel.twist.angular.z;

	xf_dotdot = -_Kdt.inverse()*( -_Kpt*(vel-xf_dot) - Kh*ht + hdot);
	//xf_dotdot <<0,xf_dotdot(1),0,0,0,0;
	for (int i=0; i<6;i++) {
		if(mask(i)!=0)
			xf_dot(i) += xf_dotdot(i)*0.002;
		else
			xf_dotdot(i) = 0.0;
	}
	for (int i=0; i<3;i++) {
		if(mask(i)!=0)
			xf(i) += xf_dot(i)*0.002;
	}

	if((mask(3)!=0)||(mask(4)!=0)||(mask(5)!=0)) {
		Eigen::VectorXd wdes(3);
		wdes = xf_dot.tail(3);
		Eigen::Quaterniond quat(xf(6),xf(3),xf(4),xf(5));
		Eigen::AngleAxisd qax(quat);
		double theta = qax.angle();
		Eigen::Vector3d r = qax.axis();
		Eigen::Vector3d eps = sin(theta/2.0)*r;
		double eta = cos(theta/2.0);

		double eta_dot = -0.5*eps.transpose()*wdes;
		Eigen::Vector3d eps_dot = 0.5*(eta*Eigen::Matrix3d::Identity() - Skew(eps))*wdes;

		eta += eta_dot*0.002;
		eps += eps_dot*0.002;

		if(eta>1) eta=1;
		else if (eta<-1) eta=-1;
		theta = 2*acos(eta);
		if(theta!=0) { //qd actually different from qe
			r = (1.0/sin(theta*0.5))*eps;
			Eigen::AngleAxisd qax_new(theta,r);
			Eigen::Quaterniond q_new(qax_new);
			xf(3)=q_new.x();
			xf(4)=q_new.y();
			xf(5)=q_new.z();
			xf(6)=q_new.w();
		}
	}


	_nextdesPose.pose.position.x = xf(0);
	_nextdesPose.pose.position.y = xf(1);
	_nextdesPose.pose.position.z = xf(2);
	_nextdesPose.pose.orientation.x = xf(3);
	_nextdesPose.pose.orientation.y = xf(4);
	_nextdesPose.pose.orientation.z = xf(5);
	_nextdesPose.pose.orientation.w = xf(6);

	//xf_dot == Eigen::VectorXd::Zero(6);
	_nextdesVel.twist.linear.x = xf_dot(0);
	_nextdesVel.twist.linear.y = xf_dot(1);
	_nextdesVel.twist.linear.z = xf_dot(2);
	_nextdesVel.twist.angular.x = xf_dot(3);
	_nextdesVel.twist.angular.y = xf_dot(4);
	_nextdesVel.twist.angular.z = xf_dot(5);

	_nextdesAcc.accel.linear.x = xf_dotdot(0);
	_nextdesAcc.accel.linear.y = xf_dotdot(1);
	_nextdesAcc.accel.linear.z = xf_dotdot(2);
	_nextdesAcc.accel.angular.x = xf_dotdot(3);
	_nextdesAcc.accel.angular.y = xf_dotdot(4);
	_nextdesAcc.accel.angular.z = xf_dotdot(5);

}

void KUKA_INVDYN::run() {
	boost::thread ctrl_loop_t( &KUKA_INVDYN::ctrl_loop, this);
	//ros::spin();
}


int main(int argc, char** argv) {
	ros::init(argc, argv, "iiwa_kdl");

	ros::AsyncSpinner spinner(1); // Use 1 thread
	spinner.start();

	KUKA_INVDYN lwr;
	lwr.run();

	geometry_msgs::PoseStamped pose;
	while(!lwr.getPose(pose) && ros::ok()) sleep(2);

	std::vector<geometry_msgs::PoseStamped> waypoints;
	geometry_msgs::PoseStamped p=pose;
	waypoints.push_back(pose);

	tf::Quaternion qorient(p.pose.orientation.x,p.pose.orientation.y,p.pose.orientation.z,p.pose.orientation.w);
	tf::Matrix3x3 orient(qorient);
	double roll,pitch,yaw;
	orient.getRPY(roll,pitch,yaw);
	//cout<<roll<<" "<<pitch<<" "<<yaw<<" "<<endl;
	//qorient.setRPY(roll-90,pitch,yaw);
	p.pose.position.x = 0;
	p.pose.position.y -= 0.47;
	p.pose.position.z -= 0.3;
	/*p.pose.orientation.z = qorient.z();
	p.pose.orientation.w = qorient.w();
	p.pose.orientation.x = qorient.x();
	p.pose.orientation.y = qorient.y();*/
	p.pose.orientation.z = -0.5;
	p.pose.orientation.w = 0.5;
	p.pose.orientation.x = -0.5;
	p.pose.orientation.y = 0.5;
	waypoints.push_back(p);
	std::vector<double> times;
	times.push_back(0);
	times.push_back(5);

	lwr.newTrajectory(waypoints,times); //Compute new trajectory
	cout<<"Trajectory complete!"<<endl<<endl;
	sleep(5);

	cout<<"Force control!"<<endl<<endl;
	Eigen::VectorXd he(6),mask(6);
	lwr.getWrench(he);
	std::vector<Eigen::VectorXd> force_wp;
	force_wp.push_back(he);
	mask=Eigen::VectorXd::Zero(6);
	mask(1)=1;
	he(1)+=5;
	force_wp.push_back(he);
	times[1]=10;

	lwr.newForceTrajectory(force_wp,times,mask);

	cout<<"Force control complete!"<<endl<<endl;
	lwr.exitForceControl();

	waypoints.clear(); times.clear();
	lwr.getDesPose(p);
	waypoints.push_back(p); times.push_back(0);
	waypoints.push_back(pose); times.push_back(5);

	lwr.newTrajectory(waypoints,times); //Compute new trajectory
	cout<<"Trajectory complete!"<<endl<<endl;

	ros::waitForShutdown();

	return 0;
}
