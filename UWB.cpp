#include "stdafx.h"
#include "UWB.h"
#include<algorithm>

using namespace LEEBEN;
using namespace Eigen;

UWB::UWB()
{
	_log = "ϵͳ�����־�ļ�";
}

void LEEBEN::UWB::Initialize_Anchors(int num, const Eigen::Vector3d * points)
{
	_AnchorNum = num;
	_Anchors = new Eigen::Vector3d[_AnchorNum];
	_AnchorsBlunders = new double[num-1];
	for (int i = 0; i < _AnchorNum; i++) {
		_Anchors[i] = points[i];
	}
}

//void LEEBEN::UWB::Initialize_Tag(const char * id)
//{
//	_TagID = id;
//}

void LEEBEN::UWB::Gross_Free(std::string& id)
{
	for (int i = 0; i < _AnchorNum - 1; i++) {
		int n = _deltaRs[id][i].size();
		double ave = 0., stdv = 0.;
		for(auto iter= _deltaRs[id][i].begin();iter!= _deltaRs[id][i].end();iter++)
			ave += *iter;
		ave /= n;
		//1sigma remove
		for (auto iter = _deltaRs[id][i].begin(); iter != _deltaRs[id][i].end(); iter++)
			stdv += (*iter - ave)*(*iter - ave);
		stdv = sqrt(stdv / (n - 1));
		auto iter = _deltaRs[id][i].begin();
		while(iter != _deltaRs[id][i].end()) {
			if (abs(*iter - ave) >  2*stdv) {
				iter = _deltaRs[id][i].erase(iter);
				continue;
			}
			iter++;
		}
	}
}

void LEEBEN::UWB::Filter_Result(const std::string& id, std::vector<double>& R, std::vector<double>& Q)
{//R Q��ʼʱ�ǿյ�
	for (int i = 0; i < _AnchorNum - 1; i++) {
		R.push_back(0.); Q.push_back(0.);
		for (auto iter = _deltaRs[id][i].begin(); iter != _deltaRs[id][i].end(); iter++)
			R[i] += *iter;
		R[i]/= _deltaRs[id][i].size();
		for (auto iter = _deltaRs[id][i].begin(); iter != _deltaRs[id][i].end(); iter++)
			Q[i] += (*iter - R[i])*(*iter - R[i]);
		if (_deltaRs[id][i].size() > 1) Q[i] /= _deltaRs[id][i].size()-1;
	}
}

void LEEBEN::UWB::Filter2file(const std::string& id,std::vector<double>& R, std::vector<double>& Q)
{
	for (int i = 0; i < _AnchorNum - 1; i++) {
		_Filtered_Tdoa <<id<< ":R" << i + 2 << " average:" << R[i] << '\t' << "variance:" << Q[i] << '\n';
		for (auto iter = _deltaRs[id][i].begin(); iter != _deltaRs[id][i].end(); iter++) {
			for (int j = 0; j < i + 3; j++)
				_Filtered_Tdoa << '\t';
			_Filtered_Tdoa << *iter << '\n';
		}
	}
}

void LEEBEN::UWB::fromfile(const char * filename,const char* tagid)
{
	blunderfromfile();
	using std::cout;
	cout << "��ȡ�ļ��У�\n";
	//������ֻ��ȡ1s���ڵ�TDOA����
	std::ifstream fin(filename);
	if (!fin.is_open()) {
		std::cout << "�ļ�������\n";
		return;
	}
	std::string str,id,log;
	std::istringstream iss;
	int tdoas;
	//std::map<std::string, bool> start;
	std::map<std::string, double> marktable1,marktable2;
	double mark;
	//mark1��ʾ��һ�е�ʱ�� mark2��ʾ��һ�е�ʱ��
	log = "F:\\������\\����׼��\\UWB\\����\\UWB_LOC�������λ��.txt";
	//log.append(to_string(clock()));
	_LOC_record.open(log.c_str(), std::ios::out);

	log = "F:\\������\\����׼��\\UWB\\����\\UWB_LOC�������λ��(LS_TLS).txt";
	_LOC_record_LS_TLS.open(log.c_str(), std::ios::out);

	log = "F:\\������\\����׼��\\UWB\\����\\UWB_LOC�������λ��(Chan_pro).txt";
	_LOC_record_pro.open(log,std::ios::out);

	log = "F:\\������\\����׼��\\UWB\\����\\ȥ���ֲ���.txt";
	_Filtered_Tdoa.open(log.c_str(), std::ios::out);
	std::vector<double> R, Q;
	while (getline(fin,str)){
		iss.str(str);
		tdoas = get_sift(iss, id, mark);
		if (tdoas != _AnchorNum)
			continue;
		if (_deltaRs.find(id) == _deltaRs.end())//��������µı�ǩ���򴴽�deltaR
		{
			_deltaRs[id] = new std::list<double>[_AnchorNum - 1];
			marktable1[id] = mark;
			marktable2[id] = mark;
		}
		else
			marktable2[id] = mark;
		if (marktable2[id]- marktable1[id]>1. 
			//&& start
			) {
			//��������ͽ���
			Gross_Free(id);
			Gross_Free(id);
			Filter_Result(id,R, Q);
			Filter2file(id,R, Q);
			for (int i = 0; i < _AnchorNum - 1; i++) {
				R[i] -= _AnchorsBlunders[i];
				_deltaRs[id][i].clear();
			}
			GetLoc(*this,id,R, mark,Q);
			marktable1[id] = marktable2[id];
		}/*
		if (!start) {
			mark1 = mark2;
			start = true;
		}*/
		push(iss,id);
	}
	fin.close();
	_LOC_record.close();
	_LOC_record_LS_TLS.close();
	_LOC_record_pro.close();
	_Filtered_Tdoa.close();
	release_deltaRs();
	std::cout << "�������\n";
}

void LEEBEN::UWB::push(std::istringstream& iss,const std::string& id)
{
	std::string str;
	double dt;
	for (int i = 0; i < _AnchorNum; i++)
		getline(iss, str, ',');//����num��
	getline(iss, str, ',');//����0
	for (int i = 0; i < _AnchorNum - 1; i++) {
		iss >> dt;
		_deltaRs[id][i].push_back(dt);//�ڵõ�����ֵ��ͬʱ��ȥϵͳ���
		getline(iss, str, ',');
	}
}

int LEEBEN::UWB::get_sift(std::istringstream& iss, std::string& id, double& t)
{
	std::string str;
	getline(iss, str, ',');//TDOA
	getline(iss, str, ',');//unit
	getline(iss, str, ',');//datetime
	str = str.substr(1);
	t=StringToDatetime(str);
	int num;
	getline(iss, str, ',');//Tagid
	id = str;
	iss >> num;
	getline(iss, str, ',');
	//std::cout << iss. << std::endl;
	return num;
}


//�ַ���תʱ��,���������Ҫ�Ѽ��ʱ����̣����Կ��Ƿ���һ��double
double LEEBEN::UWB::StringToDatetime(std::string str)
{
	char *cha = (char*)str.data();             // ��stringת����char*��
	tm tm_;                                    // ����tm�ṹ�塣
	int year, month, day, hour, minute;// ����ʱ��ĸ���int��ʱ������
	double second;
	sscanf_s(cha, "%d-%d-%d %d:%d:%lf", &year, &month, &day, &hour, &minute, &second);// ��string�洢������ʱ�䣬ת��Ϊint��ʱ������
	tm_.tm_year = year - 1900;                 // �꣬����tm�ṹ��洢���Ǵ�1900�꿪ʼ��ʱ�䣬����tm_yearΪint��ʱ������ȥ1900��
	tm_.tm_mon = month - 1;                    // �£�����tm�ṹ����·ݴ洢��ΧΪ0-11������tm_monΪint��ʱ������ȥ1��
	tm_.tm_mday = day;                         // �ա�
	tm_.tm_hour = hour;                        // ʱ��
	tm_.tm_min = minute;                       // �֡�
	tm_.tm_sec =int(second);                       // �롣
	tm_.tm_isdst = 0;                          // ������ʱ��
	time_t t_ = mktime(&tm_);                  // ��tm�ṹ��ת����time_t��ʽ��
	t_ -= int(second);
	return t_+second;                                 // ����ֵ�� 
}

void LEEBEN::UWB::release_deltaRs()
{//�ͷ�mapֵ���ڴ�
	for (auto iter = _deltaRs.begin(); iter != _deltaRs.end(); iter++) {
		for (int i = 0; i < _AnchorNum - 1; i++) {
			if((*iter).second!=nullptr)
				(*iter).second[i].clear();
		}
		delete[](*iter).second;
		(*iter).second = nullptr;
	}
}

void LEEBEN::UWB::blunderfromfile()
{
	std::ifstream logfile;
	logfile.open(_log.c_str());
	std::cout << "Blunders\n";
	for (int i = 0; i < _AnchorNum - 1; i++) {
		logfile >> _AnchorsBlunders[i];
		std::cout << "have read:" << _AnchorsBlunders[i] << std::endl;
	}
	logfile.close();
}

void LEEBEN::UWB::show_menu2()
{
	std::cout << "you choose the real-time mode, please enter your choice:\n";
	std::cout << "a)start solving\n";
	std::cout << "b)stop solving\n";
	std::cout << "c)exit the application\n";
}

//double ת����
std::string LEEBEN::UWB::Double2Datetime(time_t timestamp) {
	tm tp;
	char c[30];
	localtime_s(&tp, &timestamp);
	asctime_s(c, &tp);
	return std::string(c);
}



void LEEBEN::UWB::show_menu()
{
	std::cout << "this is the menu of UWB_LOC application, please enter your choice:\n";
	std::cout << "a) real-time positioning\tb)solve position from previous file\n";
}

void LEEBEN::UWB::realtime()//ʵʱ����ĵ�һ��
{
	std::cout << "Do you want to calibrate the system error?[y/n]\n";
	char choice; std::cin >> choice;
	while (choice != 'y' && choice != 'n') {
		while (std::cin.get() != '\n')continue;
		std::cout << "[y/n]\n";
		std::cin.clear();
		std::cin >> choice;
	}
	if (choice == 'y') {
		std::string id;
		Eigen::Vector3d pos;
		std::cout << "input the positioning tag id:\n";
		std::cin >> id;
		std::cout << "input the GCP(3d):\n";
		for (int i = 0; i < 3; i++)
			std::cin >> pos(i);
		calibrate(id, pos);
	}
	else {
		blunderfromfile();
	}
}

bool LEEBEN::UWB::broadcast_begin()
{
	if (_receiveOn){
		std::cout << "receiving\n";
		return false;
	}
	else {
		_sochl.start_socket();
		_receiveOn = true;
		if (sendto(_sochl.sclient, "#TDOATRAN,START,__127.0.0.1_5510 ", 33, 0, (sockaddr*)&_sochl.send_add, sizeof(_sochl.send_add)) < 0) {
			std::cout <<"UDP�����룺"<< WSAGetLastError() << std::endl;
		}
		return true;
	}
}

int LEEBEN::UWB::broadcast_stop()
{
	//_sochl.sclient;
	_sochl.send_add.sin_family = AF_INET;
	inet_pton(AF_INET, "127.0.0.1", &_sochl.send_add.sin_addr);
	_sochl.send_add.sin_port = htons(55510);
	std::cout<<sendto(_sochl.sclient, "#TDOATRAN,STOP,__127.0.0.1_5510 ", 32, 0, (sockaddr*)&_sochl.send_add, sizeof(_sochl.send_add));
	//Sleep(1000);
	std::cout << "��ֹͣ�㲥\n";
	return 1;
}

void LEEBEN::solve_start(LEEBEN::UWB& u)
{
	std::string str,id,log;
	std::istringstream iss;
	double mark;
	std::map<std::string, double> marktable1, marktable2;
	int tdoas;
	std::vector<double> R, Q;//R��ÿ����������֮��ľ�ֵ��Q�ǹ��˺�ķ�����ߴ�С��Ϊanchor.num-1
	log = "F:\\������\\����׼��\\UWB\\����\\UWB_LOC�������λ��.txt";
	//log.append(to_string(clock()));
	u._LOC_record.open(log.c_str(), std::ios::out);

	log = "F:\\������\\����׼��\\UWB\\����\\UWB_LOC�������λ��(LS_TLS).txt";
	u._LOC_record_LS_TLS.open(log.c_str(), std::ios::out);

	log = "F:\\������\\����׼��\\UWB\\����\\UWB_LOC�������λ��(Chan_pro).txt";
	u._LOC_record_pro.open(log, std::ios::out);

	log = "F:\\������\\����׼��\\UWB\\����\\ȥ���ֲ���.txt";
	u._Filtered_Tdoa.open(log.c_str(), std::ios::out);
	//log.append
	int send_add_size;
	send_add_size = sizeof(u._sochl.send_add);
	while (true) {
		memset(u._buffer, 0, sizeof(u._buffer));
		recvfrom(u._sochl.sclient, u._buffer, sizeof(u._buffer), 0, (SOCKADDR*)&u._sochl.send_add, &send_add_size);
		str = u._buffer;
		if (str == "leave") {
			std::cout << "��ֹͣ����\n";
			u._receiveOn ^= u._receiveOn;
			u._sochl.end_socket();//�ͷ���Դ
			std::cout << "�ѹر��׽���\n";
			break;
		}
		iss.str(str);
		tdoas = u.get_sift(iss, id, mark);
		if (tdoas != u._AnchorNum )
			continue;
		if (u._deltaRs.find(id) == u._deltaRs.end()) {
			marktable1[id] = mark;
			marktable2[id] = mark;
			u._deltaRs[id] = new std::list<double>[u._AnchorNum - 1];
		}
		else
			marktable2[id] = mark;
		if (marktable2[id] - marktable1[id]>1	//���ƽ���һ�ε�ʱ���
			) {
			//��������ͽ���
			u.Gross_Free(id);
			u.Gross_Free(id);
			u.Filter_Result(id,R, Q);
			for (int i = 0; i <u. _AnchorNum - 1; i++) {
				R[i] -= u._AnchorsBlunders[i];
				u._deltaRs[id][i].clear();
			}
			std::thread t(GetLoc, std::ref(u), id, R, mark, Q);//��ʼ�����߳�
			t.join();
			marktable1[id] = marktable2[id];
		}
		u.push(iss,id);
	}
	u._LOC_record.close();
	u._LOC_record_LS_TLS.close();
	u._LOC_record_pro.close();
	u._Filtered_Tdoa.close();
	std::cout << "�ѹر�����ļ���\n";
	u.release_deltaRs();
	u._deltaRs.clear();
}

bool LEEBEN::GetLoc(UWB& u,const std::string&id, std::vector<double>& deltaR,double timestamp,std::vector<double>& Q)
{
	Eigen::Vector3d ans;
	std::string t = u.Double2Datetime((time_t)timestamp);
		//chan�㷨
	ans = chan(u, deltaR, Q);
	u._LOC_record <<id<<','<< t ;
	u._LOC_record_LS_TLS << id << ',' << t ;
	chan_pro(u,ans, deltaR, Q);
	u._LOC_record_pro << id << ',' << t ;
	u._TagPos[id] = ans;
	//��Ļ��ʾ
	std::cout << id << ',';
	for (int i = 0; i < ans.size() - 1; i++)
		std::cout << ans(i) << ',';
	std::cout << t;
		//std::cout << "chan's result: \n" << ans << std::endl << std::endl;
		/*if (MS(u, deltaR)) {
			std::cout << "��С���˵���������\n";
			u._Tag = ans;
		}
		delete[] deltaR;*/
		/*std::cout << "time: " << timestamp << ", LOC: ";
		for (int i = 0; i < u._Tag.size() - 1; i++)
			std::cout << u._Tag(i) << ',';*/
		//std::cout << u._Tag(3);
		/*u._LOC_record << "time: " << timestamp << ", LOC: ";
		for (int i = 0; i < u._Tag.size() - 1; i++)
			u._LOC_record << u._Tag(i) << ',';
		u._LOC_record << u._Tag(3)<<std::endl;*/
	deltaR.clear(); Q.clear();
	return true;
}

Vector3d LEEBEN::chan(UWB& u, std::vector<double>& deltaR,std::vector<double>& Q)
{
	MatrixXd A(u._AnchorNum-1, 3), b(u._AnchorNum-1, 1),D=MatrixXd::Zero(u._AnchorNum-1,u._AnchorNum-1);
	Vector3d ans;
	double k1;
	k1 = u._Anchors[0].squaredNorm();
	for (int i = 0; i < u._AnchorNum - 1; i++) {
		double ki, xi, yi;
		ki= u._Anchors[i+1].squaredNorm();
		xi = u._Anchors[i + 1](0) - u._Anchors[0](0);
		yi = u._Anchors[i + 1](1) - u._Anchors[0](1);
		A(i, 0) = xi; A(i, 1) = yi; A(i, 2) = deltaR[i];
		b(i, 0) = ki - k1 - deltaR[i] * deltaR[i];
		D(i, i) = 2 * deltaR[i] * Q[i];
	}
	A *= 2;
	//std::cout << A << std::endl << b << std::endl;
	//MatrixXd ans(3, 1);
	//normal���
	ans = (A.transpose()*D*A).inverse()*A.transpose()*D*b; 
	/*double z=u._Anchors[0](2)
		- sqrt(
			ans(2, 0) * ans(2, 0)
			- (u._Anchors[0](0) - ans(0))*(u._Anchors[0](0) - ans(0))
			- (u._Anchors[0](1) - ans(1))*(u._Anchors[0](1) - ans(1))
		);
	MatrixXd G(4, 3), za(3, 1), ha(4, 1),B(4,4),fai(4,4);
	G << 1, 0, 0,
		0, 1, 0,
		0, 0, 1,
		1, 1, 1;
	ha << (ans(0) - u._Anchors[0](0))*(ans(0) - u._Anchors[0](0)),
		(ans(1) - u._Anchors[0](1))*(ans(1) - u._Anchors[0](1)),
		(z - u._Anchors[0](2))*(z - u._Anchors[0](2)),
		ans(2)*ans(2);
	B << ans(0) - u._Anchors[0](0), 0, 0, 0,
		0, ans(1) - u._Anchors[0](1), 0, 0,
		0, 0, z - u._Anchors[0](2), 0,
		0, 0, 0, ans(2);
	fai=4*B*(A.transpose()*D*A).inverse()
	za=(G.transpose()**///����z���겻̫����ʵ���������Ϊz�������Է���ֱ�ӽ�õģ��ٽ������ƽ���н����������ã�����ƽ������֮ǰ��Ȳ��ᷢ���ı�
	//u.ans(0, 0) = ans(0,0); u.ans(1,0) = ans(1, 0); u.ans(2,0) = ans(2,0);
	for (int i = 0; i < ans.size(); i++)
		u._LOC_record << ans(i) << ',';
	//LS_TLS���
	ans = LS_TLS(A.block(0, 0, A.rows(), 2), A.block(0, 2, A.rows(), 1), b);
	/*u.ans(2, 0) = u._Anchors[0](2, 0)
		- sqrt(
			ans(2, 0) * ans(2, 0)
			- (u._Anchors[0](0, 0) - ans(0, 0))*(u._Anchors[0](0, 0) - ans(0, 0))
			- (u._Anchors[0](1, 0) - ans(1, 0))*(u._Anchors[0](1, 0) - ans(1, 0))
		);*/
	//u.ans(0, 0) = ans(0, 0); u.ans(1, 0) = ans(1, 0); u.ans(2, 0) = ans(2, 0);
	for (int i = 0; i < ans.size(); i++)
		u._LOC_record_LS_TLS << ans (i)<< ',';
	return ans;
}

void LEEBEN::chan_pro(UWB & u,Eigen::Vector3d& ans, std::vector<double>& deltaR, std::vector<double>& Q)
{
	MatrixXd A=MatrixXd::Zero(u._AnchorNum - 1, u._AnchorNum-1), v(u._AnchorNum - 1, 1),B(u._AnchorNum-1,3),x(3,1),W(u._AnchorNum-1,1),
		D = MatrixXd::Zero(u._AnchorNum - 1, u._AnchorNum - 1);
	double k1, k[4]{0,0,0,0};
	k1 = u._Anchors[0].squaredNorm();
	for (int i = 0; i < u._AnchorNum - 1; i++) {
		double xi, yi;
		k[i] = u._Anchors[i + 1].squaredNorm();
		xi = u._Anchors[i + 1](0) - u._Anchors[0](0);
		yi = u._Anchors[i + 1](1) - u._Anchors[0](1);
		A(i, i) = deltaR[i] + ans(2); 
		B(i, 0) = xi; B(i, 1) = yi; B(i, 2) = deltaR[i];
		W(i, 0) = k1 - k[i] + deltaR[i] * deltaR[i] + 2 * (xi*ans(0) + yi*ans(1) + deltaR[i] * ans(2));
		D(i, i) = Q[i];
	}
	A *= 2; B *= 2; W *= -1;
	//std::cout << A << std::endl << b << std::endl;
	MatrixXd N(u._AnchorNum - 1, u._AnchorNum - 1), Nb(3, 3), K(u._AnchorNum - 1, 1);
	N = A*D*A.transpose(); Nb = B.transpose()*N.inverse()*B;
	x = Nb.inverse()*B.transpose()*N.inverse()*W;
	K = -N.inverse()*(B*x - W);
	v = D*A.transpose()*K;
	int count = 0;
	while (v.squaredNorm() + x.squaredNorm() > 0.1 && count<20) {
		ans(0) += x(0); ans(1) += x(1); ans(2) += x(2);
		for (int i = 0; i < u._AnchorNum - 1; i++)
			deltaR[i] += v(i);
		for (int i = 0; i < u._AnchorNum - 1; i++) {
			A(i, i) = deltaR[i] + ans(2);
			B(i, 2) = 2 * deltaR[i];
			W(i, 0) = k1 - k[i] 
				+ deltaR[i] * deltaR[i] +(B(i,0)*ans(0) + B(i,1)*ans(1) + B(i,2) * ans(2));
		}
		A *= 2; W *= -1;
		N = A*D*A.transpose(); Nb = B.transpose()*N.inverse()*B;
		x = Nb.inverse()*B.transpose()*N.inverse()*W;
		K = -N.inverse()*(B*x - W);
		v = D*A.transpose()*K;
		count++;
	}
	if (count == 20)
		u._LOC_record_pro << "������,";
	else {
		for (int i = 0; i < ans.size(); i++)
			u._LOC_record_pro << ans(i) << ',';
	}
}
//
//bool LEEBEN::MS(UWB& u, std::vector<double>* deltaR,std::vector<double>&Q)
//{
//	int count = 0,min=deltaR[0].size();
//	MatrixXd deltax(4, 1), A(1 + min*(u._AnchorNum - 1), 4), b(A.rows(), 1);
//	Eigen::MatrixXd Qx(A.cols(), A.cols());
//	A(0, 0) = u._Anchors[0](0) - u._Tag(0); A(0, 1) = u._Anchors[0](1,0) - u._Tag(1,0); A(0, 2) = u._Anchors[0](2,0) - u._Tag(2,0);
//	A(0, 3) = u._Tag(2); b(0, 0) = (u._Anchors[0] - u._Tag.block(0, 0, 3, 1)).squaredNorm() - u._Tag(3,0) * u._Tag(3,0);
//	for (int i = 0; i < u._AnchorNum - 1; i++) {
//		double xi, yi, zi, ki;
//		xi = u._Anchors[i + 1](0) - u._Tag(0);
//		yi = u._Anchors[i + 1](1) - u._Tag(1);
//		zi = u._Anchors[i + 1](2) - u._Tag(2);
//		ki = xi*xi + yi*yi + zi*zi;
//		for (int j = 0; j < min; j++) {
//			A(j + i*min + 1, 0) = xi;
//			A(j + i*min + 1, 1) = yi;
//			A(j + i*min + 1, 2) = zi;
//			A(j + i*min + 1, 3) = u._Tag(2) + deltaR[i][j];
//			b(j + i*min + 1, 0) = ki - (deltaR[i][j] - u._Tag(3))*(deltaR[i][j] - u._Tag(2));
//		}
//	}
//	A *= 2;
//	//std::cout << A << std::endl << b << std::endl<<std::endl;
//	deltax = Normal_LS(A, b, Eigen::MatrixXd::Identity(b.rows(), b.rows()), Qx);
//	//std::cout << deltax << std::endl;
//	while (count<15 && deltax.squaredNorm()>10e-10) {
//		u._Tag += deltax;
//		A(0, 0) = u._Anchors[0](0) - u._Tag(0); A(0, 1) = u._Anchors[0](1) - u._Tag(1); A(0, 2) = u._Anchors[0](2) - u._Tag(2);
//		A(0, 3) = u._Tag(3); b(0, 0) = (u._Anchors[0] - u._Tag.block(0, 0, 3, 1)).squaredNorm() - u._Tag(3) * u._Tag(3);
//		for (int i = 0; i < u._AnchorNum - 1; i++) {
//			double xi, yi, zi, ki;
//			xi = u._Anchors[i + 1](0) - u._Tag(0);
//			yi = u._Anchors[i + 1](1) - u._Tag(1);
//			zi = u._Anchors[i + 1](2) - u._Tag(2);
//			ki = xi*xi + yi*yi + zi*zi;
//			for (int j = 0; j < min; j++) {
//				A(j + i*min + 1, 0) = xi;
//				A(j + i*min + 1, 1) = yi;
//				A(j + i*min + 1, 2) = zi;
//				A(j + i*min + 1, 3) = u._Tag(3) + deltaR[i][j];
//				b(j + i*min + 1, 0) = ki - (deltaR[i][j] - u._Tag(3))*(deltaR[i][j] - u._Tag(3));
//			}
//		}
//		A *= 2;
//		deltax = (A.transpose()*A).inverse()*A.transpose()*b;
//		count++;
//	}
//	return (count == 15);
//}

Eigen::VectorXd LEEBEN::Normal_LS(const Eigen::MatrixXd & A, const Eigen::MatrixXd & b, const Eigen::MatrixXd & P,Eigen::MatrixXd& Q)
{
	//std::cout << A << std::endl;
	Eigen::VectorXd x(A.cols(), 1);
	Q.resize(x.rows(), x.rows());
	Q = (A.transpose()*P*A).inverse();
	x = Q*A.transpose()*P*b;
	return x;
}

Eigen::VectorXd LEEBEN::Total_LS(const Eigen::MatrixXd & A, const Eigen::MatrixXd & b)
{
	Eigen::MatrixXd C(b.rows(),A.cols() + 1);
	C.block(0, 0, A.rows(), A.cols()) = A;
	C.block(0, A.cols(), A.rows(),1) = b;
	Eigen::JacobiSVD<MatrixXd> svd(C, Eigen::ComputeFullU | ComputeFullV);
	Eigen::MatrixXd U(C.rows(), C.rows()), V(C.cols(), C.cols());
	U = svd.matrixU(); V = svd.matrixV();
	//S = U.inverse()*C*V.transpose().inverse();
	//std::cout << C << std::endl<< U << std::endl << V << std::endl;
	/*for (int i = 0; i < C.rows(); i++)
		std::cout << U(0, i) << "  ";*/
	Eigen::VectorXd x(A.cols());
	x = -1. / V(A.cols(), A.cols())*V.block(A.cols(),0,1,A.cols()).transpose();
	return x;
}

Eigen::VectorXd LEEBEN::LS_TLS(const Eigen::MatrixXd & A1, const Eigen::MatrixXd & A2, const Eigen::MatrixXd & b)
{
	Eigen::MatrixXd C(b.rows(), A1.cols() + A2.cols() + 1);
	C.block(0, 0, A1.rows(), A1.cols()) = A1;
	C.block(0, A1.cols(),A2.rows(),A2.cols()) = A2;
	C.col(A1.cols() + A2.cols()) = b;
	Eigen::MatrixXd Q(C.rows(),C.rows()), R(C.rows(),C.cols());
	Maqr(C,Q,R);
	//std::cout << C << std::endl << Q << std::endl << R << std::endl;
	Eigen::VectorXd x1,x2;
	x2 = Total_LS(R.block(A1.cols(), A1.cols(), A2.cols() + 1, A2.cols()), R.block(A1.cols(), R.cols()-1, A2.cols() + 1, 1));
	MatrixXd Qx;
	//std::cout << A1 << std::endl << A2 << std::endl << x2 << std::endl;
	//std::cout << x2 << std::endl;
	x1 = Normal_LS(R.block(0, 0, A1.cols(), A1.cols()), R.block(0, R.cols() - 1, A1.cols(), 1)- R.block(0, A1.cols(), A1.cols(), A2.cols())*x2
		, Eigen::MatrixXd::Identity(A1.cols(), A1.cols()), Qx);
	Eigen::VectorXd x(x1.rows() + x2.rows());
	x.block(0, 0, x1.rows(), 1) = x1;
	x.block(x1.rows(), 0, x2.rows(), 1) = x2;
	return x;
}

void LEEBEN::Maqr(const Eigen::MatrixXd & CA, Eigen::MatrixXd & Q, Eigen::MatrixXd & A)
{//����һ��ʵ����QR�ֽ�ĺ���
	using std::cout;
	using std::endl;
	int i, j, k, nn, jj,m,n;
	double u, alpha, w, t;
	m = CA.rows(); n = CA.cols();
	A = CA;
	if (m < n) {
		cout << "\nQR�ֽ�ʧ�ܣ�" << endl; exit(1);
	} //��֤����>��������ʵ��QR�ֽ⣬������ 
	Q = Eigen::MatrixXd::Identity(m, m);
	nn = n;
	if (m == n) nn = m - 1;
	for (k = 0; k <= nn - 1; k++)//�ڴ�ѭ��k��0~m���У�����H�������⣬���Q���Լ����A 
	{
		u = 0.0;
		for (i = k; i <= m - 1; i++) {
			w = fabs(A(i,k));
			if (w > u) u = w;
		}
		alpha = 0.0;
		for (i = k; i <= m - 1; i++) {
			t = A(i,k) / u; alpha = alpha + t * t;
		}
		if (A(k,k) > 0.0) u = -u;
		alpha = u * sqrt(alpha);
		if (fabs(alpha) + 1.0 == 1.0) {
			cout << "\nQR�ֽ�ʧ�ܣ�" << endl;
			exit(1);
		} u = sqrt(2.0*alpha*(alpha - A(k,k)));
		if ((u + 1.0) != 1.0) {
			A(k,k) = (A(k,k) - alpha) / u;
			for (i = k + 1; i <= m - 1; i++)
				A(i,k) = A(i,k) / u;
			//���Ͼ���H�������ã�ʵ���ϳ���û�������κ����ݽṹ���洢H�� 
			//�󣬶���ֱ�ӽ�u������Ԫ�ظ�ֵ��ԭA�����ԭ��������Ӧ��λ�� 
			//��������Ϊ�˼�����˾���Q��A 
			for (j = 0; j <= m - 1; j++) {
				t = 0.0;
				for (jj = k; jj <= m - 1; jj++)
					t = t + A(jj,k) * Q(jj,j);
				for (i = k; i <= m - 1; i++)
					Q(i,j) = Q(i,j) - 2.0*t*A(i,k);
			} //��˾���Q��ѭ��������õ�һ�������ٽ��������ת��һ�¾͵õ�QR�ֽ��е�Q���� 
				//Ҳ������������ 
			for (j = k + 1; j <= n - 1; j++) {
				t = 0.0;
				for (jj = k; jj <= m - 1; jj++)
					t = t + A(jj,k) * A(jj,j);
				for (i = k; i <= m - 1; i++)
					A(i,j) = A(i,j) - 2.0*t*A(i,k);
			} //H�������A����ѭ�����֮���������ǲ��ֵ����ݾ��������Ǿ���R
			A(k,k) = alpha;
			for (i = k + 1; i <= m - 1; i++)
				A(i,k)= 0.0;
		}
	}
	Q.transposeInPlace();//QR�ֽ����
}

bool LEEBEN::UWB::switch_receive()
{
	if (_receiveOn) {
		broadcast_stop();//���ù㲥ͣ����
		//_receiveOn = false;
		SOCKET stoprecv = socket(AF_INET, SOCK_DGRAM, 0);
		if(sendto(stoprecv, "leave", 6, 0, (sockaddr*)&_sochl.recv_add, sizeof(_sochl.recv_add))<0)
			std::cout<<"����㲥ֹͣ���ִ��󣬴����룺"<<WSAGetLastError()<<std::endl;
		closesocket(stoprecv);
		return true;
	}
	else {
		std::cout << "�׽����ѹر�\n";
		return false;
	}
}

void LEEBEN::UWB::calibrate(const std::string & tagid, const Eigen::Vector3d & pos)
{
	//_sochl.start_socket();
	_deltaRs[tagid]=new std::list<double>[_AnchorNum - 1];
	int count = 0, tdoas, sendsize;
	double timestamp;
	std::string str,id;
	std::istringstream iss;
	sendsize = sizeof(_sochl.send_add);
	std::cout << "����ɣ�\n";
	while (count < 100) {
		memset(_buffer, 0, sizeof(_buffer));
		recvfrom(_sochl.sclient, _buffer, sizeof(_buffer), 0, (sockaddr*)&_sochl.send_add, &sendsize);
		str = _buffer;
		iss.str(str);
		tdoas = get_sift(iss, id, timestamp);
		if (tdoas == _AnchorNum && id==tagid) {
			count++;
			push(iss,id);
			std::cout << count << "\b\b";
		}
	}//���չ�100���������ݺ󣬿�ʼ����ϵͳ���
	double* real = new double[_AnchorNum ];
	for (int i = 0; i < _AnchorNum; i++) {
		real[i] = sqrt(
			(_Anchors[i](0) - pos(0))*(_Anchors[i](0) - pos(0))
			+ (_Anchors[i](1) - pos(1))*(_Anchors[i](1) - pos(1))
			+ (_Anchors[i](2) - pos(2))*(_Anchors[i](2) - pos(2))
		);
	}
	Gross_Free(id);
	Gross_Free(id);
	std::vector<double>R, Q;
	Filter_Result(id,R, Q);
	std::cout << "\n������ɣ�ϵͳ���ȡֵ��";
	for (int i = 0; i < _AnchorNum - 1; i++) {
		_AnchorsBlunders[i] = R[i] - (real[i + 1] - real[0]);
		std::cout << _AnchorsBlunders[i] << '\t';
		_deltaRs[id][i].clear();
	}
	delete[]real;
	delete[] _deltaRs[id];
	_deltaRs.erase(id);
	//solve_start(*this);
	//_sochl.end_socket();
}

UWB::~UWB()
{
	//������ʹ�õ�ϵͳ��������־�ļ�
	std::ofstream logfile;
	logfile.open(_log.c_str(), std::ios::out);
	for (int i = 0; i < _AnchorNum - 1; i++)
		logfile << _AnchorsBlunders[i] << std::endl;
	logfile.close();
	if (_Anchors != nullptr)
		delete[] _Anchors;
	if(_AnchorsBlunders!=nullptr)
		delete[] _AnchorsBlunders;
	/*for(auto iter=_deltaRs.begin();iter!=_deltaRs.end();iter++)
		if((*iter).second!=nullptr)*/
	release_deltaRs();
}

LEEBEN::UWB::Soket_Handle::Soket_Handle()
{
	//���ñ���IP�Ͷ˿�
	recv_add.sin_family = AF_INET;
	//recv_add.sin_addr.S_un.S_addr = INADDR_ANY;
	inet_pton(AF_INET, "127.0.0.1", &recv_add.sin_addr);
	recv_add.sin_port = htons(5510);
}

void LEEBEN::UWB::Soket_Handle::start_socket()
{
	//�󶨱��ص�ַ��ӿ�
	WSADATA wsadata;
	WSAStartup(MAKEWORD(2, 2), &wsadata);
	sclient = socket(AF_INET, SOCK_DGRAM, 0);
	//����Ŀ��IP��˿�
	send_add.sin_family = AF_INET;
	inet_pton(AF_INET, "127.0.0.1", &send_add.sin_addr);
	send_add.sin_port = htons(55510);
	if (bind(sclient, (SOCKADDR*)&recv_add, sizeof(SOCKADDR)) <0)
		std::cout << WSAGetLastError() << std::endl;
}

void LEEBEN::UWB::Soket_Handle::end_socket()
{
	closesocket(sclient);
	WSACleanup();
}
