#pragma once
#include"Eigen\Dense"
#include"Eigen/QR"
#include<list>
#include<vector>
#include<string>
#include<map>
#include<ctime>
#include<fstream>
#include<sstream>
#include<iostream>
#include"WinSock2.h"
#include <stdio.h>   
#include< Ws2tcpip.h>  
#include<thread>
#pragma comment(lib,"WS2_32.lib")

namespace LEEBEN {
	class UWB
	{
	private:
		int _AnchorNum;
		Eigen::Vector3d* _Anchors;
		double* _AnchorsBlunders;//基站的系统误差
		std::map<std::string,Eigen::Vector3d>  _TagPos;//不仅包括标签的二维坐标，第三个表示中心基站到标签的距离r1，单位cm
//		std::vector<std::string> _TagIDs;//
		std::map<std::string, std::list<double>* > _deltaRs;//存放各标签得到的deltaR
		std::ofstream _LOC_record;//用于记录经典chan解算结果的文件流
		std::ofstream _LOC_record_LS_TLS;//用于记录混合最小二乘chan的结果文件流
		std::ofstream _LOC_record_pro;//用于记录改进Chan算法结果的文件流
		std::ofstream _Filtered_Tdoa;//用于记录粗差剔除后的TDOA
		std::string _log;//这里是日志文件的文件名，存放系统误差，在程序结束后将新系统误差保存
		class Soket_Handle {
		public:
			SOCKET sclient;
			sockaddr_in recv_add;
			sockaddr_in send_add;
			Soket_Handle();
			void start_socket();
			void end_socket();
			~Soket_Handle() {}
		};
		Soket_Handle _sochl;
		char _buffer[1024];
		bool _receiveOn;//接收的开关

		//开始解算TDOA
		friend void solve_start(UWB& uwb);

		//总的计算位置，既有控制台输出，也有文件输出,既有chan直接得出的估值结果，也有
		friend bool GetLoc(UWB& u, const std::string&id,std::vector<double>& deltaR,double timestamp,std::vector<double>& Q);

		//经典Chan算法
		friend Eigen::Vector3d chan(UWB& u, std::vector<double>& deltaR,std::vector<double>& Q);

		//改进chan算法,需要给入初值
		friend void chan_pro(UWB& u, Eigen::Vector3d& ,std::vector<double>& deltaR,std::vector<double>& Q);

		//最小二乘迭代
		//friend bool MS(UWB& u, std::vector<double>*deltaR,std::vector<double>& Q);
		
	public:
		UWB();
		//初始化基站，给入数量和坐标、基站系统误差；
		void Initialize_Anchors(int num, const Eigen::Vector3d* points);

		////初始化标签，给入标签id用于筛选
		//void Initialize_Tag(const char* id);

		//将粗差剔除结果写到文件中去
		void Filter2file(const std::string& id,std::vector<double>& R, std::vector<double>& Q);

		//把数据从txt中读入，并进行输出和解算
		void fromfile(const char* filename,const char* tagid);
		

		//显示控制菜单
		void show_menu();


		//实时定位
		void realtime();

		//向55510发送TDOA广播请求
		bool broadcast_begin();

		//停止55510的TDOA广播
		int broadcast_stop();

		//改变接收TDOA的开关
		bool switch_receive();

		//输入校准标签和位置，对系统误差进行校正
		void calibrate(const std::string& tagid, const Eigen::Vector3d& pos);

		//如果选择了实时解算，那么控制菜单应该为
		void show_menu2();

		~UWB();
	private:
		//给入一行文本，push进入deltaR
		void push(std::istringstream& iss,const std::string&);

		//获得时间\标签ID和TDOA个数，三者作为筛选条件,返回的是TDOA个数
		int get_sift(std::istringstream& iss, std::string& id, double& time);

		//剔除粗差，返回num-1组数据，每一组代表去掉粗差之后的剩余距离差数据，给入的是未剔除粗差的距离差数据
		void Gross_Free(std::string&);

		//给入剔除粗差的每组数据，返回AnchorNum-1个均值数字，以及均值的方差
		void Filter_Result(const std::string&, std::vector<double>& R, std::vector<double>& Q);

		//时间的转换，得到一个三位小数
		double StringToDatetime(std::string str);

		//double 转日期
		std::string LEEBEN::UWB::Double2Datetime(time_t timestamp);

		//释放map<string,list*>的内存
		void release_deltaRs();

		//从文件中读取系统误差
		void blunderfromfile();

		
	};
	
	//普通最小二乘法，线性
	Eigen::VectorXd Normal_LS(const Eigen::MatrixXd& A, const Eigen::MatrixXd& b,const Eigen::MatrixXd& P,Eigen::MatrixXd& Q);

	//整体最小二乘法
	Eigen::VectorXd Total_LS(const Eigen::MatrixXd& A, const Eigen::MatrixXd& b);

	//混合最小二乘法
	Eigen::VectorXd LS_TLS(const Eigen::MatrixXd& A1,const Eigen::MatrixXd& A2, const Eigen::MatrixXd& b);

	//QR分解
	void Maqr(const Eigen::MatrixXd& A, Eigen::MatrixXd& Q, Eigen::MatrixXd& R);
}

//this is the trial of github usage.
//let me check what changes have happened.

