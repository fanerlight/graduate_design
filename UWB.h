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
		double* _AnchorsBlunders;//��վ��ϵͳ���
		std::map<std::string,Eigen::Vector3d>  _TagPos;//����������ǩ�Ķ�ά���꣬��������ʾ���Ļ�վ����ǩ�ľ���r1����λcm
//		std::vector<std::string> _TagIDs;//
		std::map<std::string, std::list<double>* > _deltaRs;//��Ÿ���ǩ�õ���deltaR
		std::ofstream _LOC_record;//���ڼ�¼����chan���������ļ���
		std::ofstream _LOC_record_LS_TLS;//���ڼ�¼�����С����chan�Ľ���ļ���
		std::ofstream _LOC_record_pro;//���ڼ�¼�Ľ�Chan�㷨������ļ���
		std::ofstream _Filtered_Tdoa;//���ڼ�¼�ֲ��޳����TDOA
		std::string _log;//��������־�ļ����ļ��������ϵͳ���ڳ����������ϵͳ����
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
		bool _receiveOn;//���յĿ���

		//��ʼ����TDOA
		friend void solve_start(UWB& uwb);

		//�ܵļ���λ�ã����п���̨�����Ҳ���ļ����,����chanֱ�ӵó��Ĺ�ֵ�����Ҳ��
		friend bool GetLoc(UWB& u, const std::string&id,std::vector<double>& deltaR,double timestamp,std::vector<double>& Q);

		//����Chan�㷨
		friend Eigen::Vector3d chan(UWB& u, std::vector<double>& deltaR,std::vector<double>& Q);

		//�Ľ�chan�㷨,��Ҫ�����ֵ
		friend void chan_pro(UWB& u, Eigen::Vector3d& ,std::vector<double>& deltaR,std::vector<double>& Q);

		//��С���˵���
		//friend bool MS(UWB& u, std::vector<double>*deltaR,std::vector<double>& Q);
		
	public:
		UWB();
		//��ʼ����վ���������������ꡢ��վϵͳ��
		void Initialize_Anchors(int num, const Eigen::Vector3d* points);

		////��ʼ����ǩ�������ǩid����ɸѡ
		//void Initialize_Tag(const char* id);

		//���ֲ��޳����д���ļ���ȥ
		void Filter2file(const std::string& id,std::vector<double>& R, std::vector<double>& Q);

		//�����ݴ�txt�ж��룬����������ͽ���
		void fromfile(const char* filename,const char* tagid);
		

		//��ʾ���Ʋ˵�
		void show_menu();


		//ʵʱ��λ
		void realtime();

		//��55510����TDOA�㲥����
		bool broadcast_begin();

		//ֹͣ55510��TDOA�㲥
		int broadcast_stop();

		//�ı����TDOA�Ŀ���
		bool switch_receive();

		//����У׼��ǩ��λ�ã���ϵͳ������У��
		void calibrate(const std::string& tagid, const Eigen::Vector3d& pos);

		//���ѡ����ʵʱ���㣬��ô���Ʋ˵�Ӧ��Ϊ
		void show_menu2();

		~UWB();
	private:
		//����һ���ı���push����deltaR
		void push(std::istringstream& iss,const std::string&);

		//���ʱ��\��ǩID��TDOA������������Ϊɸѡ����,���ص���TDOA����
		int get_sift(std::istringstream& iss, std::string& id, double& time);

		//�޳��ֲ����num-1�����ݣ�ÿһ�����ȥ���ֲ�֮���ʣ���������ݣ��������δ�޳��ֲ�ľ��������
		void Gross_Free(std::string&);

		//�����޳��ֲ��ÿ�����ݣ�����AnchorNum-1����ֵ���֣��Լ���ֵ�ķ���
		void Filter_Result(const std::string&, std::vector<double>& R, std::vector<double>& Q);

		//ʱ���ת�����õ�һ����λС��
		double StringToDatetime(std::string str);

		//double ת����
		std::string LEEBEN::UWB::Double2Datetime(time_t timestamp);

		//�ͷ�map<string,list*>���ڴ�
		void release_deltaRs();

		//���ļ��ж�ȡϵͳ���
		void blunderfromfile();

		
	};
	
	//��ͨ��С���˷�������
	Eigen::VectorXd Normal_LS(const Eigen::MatrixXd& A, const Eigen::MatrixXd& b,const Eigen::MatrixXd& P,Eigen::MatrixXd& Q);

	//������С���˷�
	Eigen::VectorXd Total_LS(const Eigen::MatrixXd& A, const Eigen::MatrixXd& b);

	//�����С���˷�
	Eigen::VectorXd LS_TLS(const Eigen::MatrixXd& A1,const Eigen::MatrixXd& A2, const Eigen::MatrixXd& b);

	//QR�ֽ�
	void Maqr(const Eigen::MatrixXd& A, Eigen::MatrixXd& Q, Eigen::MatrixXd& R);
}


