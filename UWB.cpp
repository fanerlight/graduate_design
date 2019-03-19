#include "stdafx.h"
#include "UWB.h"
#include<algorithm>

using namespace LEEBEN;
using namespace Eigen;

UWB::UWB()
{
	_log = "系统误差日志文件";
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
{//R Q开始时是空的
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
	cout << "读取文件中：\n";
	//以下是只获取1s以内的TDOA数据
	std::ifstream fin(filename);
	if (!fin.is_open()) {
		std::cout << "文件不存在\n";
		return;
	}
	std::string str,id,log;
	std::istringstream iss;
	int tdoas;
	//std::map<std::string, bool> start;
	std::map<std::string, double> marktable1,marktable2;
	double mark;
	//mark1表示上一行的时间 mark2表示这一行的时间
	log = "F:\\电子书\\毕设准备\\UWB\\数据\\UWB_LOC程序解算位置.txt";
	//log.append(to_string(clock()));
	_LOC_record.open(log.c_str(), std::ios::out);

	log = "F:\\电子书\\毕设准备\\UWB\\数据\\UWB_LOC程序解算位置(LS_TLS).txt";
	_LOC_record_LS_TLS.open(log.c_str(), std::ios::out);

	log = "F:\\电子书\\毕设准备\\UWB\\数据\\UWB_LOC程序解算位置(Chan_pro).txt";
	_LOC_record_pro.open(log,std::ios::out);

	log = "F:\\电子书\\毕设准备\\UWB\\数据\\去除粗差结果.txt";
	_Filtered_Tdoa.open(log.c_str(), std::ios::out);
	std::vector<double> R, Q;
	while (getline(fin,str)){
		iss.str(str);
		tdoas = get_sift(iss, id, mark);
		if (tdoas != _AnchorNum)
			continue;
		if (_deltaRs.find(id) == _deltaRs.end())//如果遇到新的标签，则创建deltaR
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
			//考虑输出和解算
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
	std::cout << "解算完毕\n";
}

void LEEBEN::UWB::push(std::istringstream& iss,const std::string& id)
{
	std::string str;
	double dt;
	for (int i = 0; i < _AnchorNum; i++)
		getline(iss, str, ',');//跳过num次
	getline(iss, str, ',');//跳过0
	for (int i = 0; i < _AnchorNum - 1; i++) {
		iss >> dt;
		_deltaRs[id][i].push_back(dt);//在得到测量值的同时减去系统误差
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


//字符串转时间,这里如果想要把间隔时间调短，可以考虑返回一个double
double LEEBEN::UWB::StringToDatetime(std::string str)
{
	char *cha = (char*)str.data();             // 将string转换成char*。
	tm tm_;                                    // 定义tm结构体。
	int year, month, day, hour, minute;// 定义时间的各个int临时变量。
	double second;
	sscanf_s(cha, "%d-%d-%d %d:%d:%lf", &year, &month, &day, &hour, &minute, &second);// 将string存储的日期时间，转换为int临时变量。
	tm_.tm_year = year - 1900;                 // 年，由于tm结构体存储的是从1900年开始的时间，所以tm_year为int临时变量减去1900。
	tm_.tm_mon = month - 1;                    // 月，由于tm结构体的月份存储范围为0-11，所以tm_mon为int临时变量减去1。
	tm_.tm_mday = day;                         // 日。
	tm_.tm_hour = hour;                        // 时。
	tm_.tm_min = minute;                       // 分。
	tm_.tm_sec =int(second);                       // 秒。
	tm_.tm_isdst = 0;                          // 非夏令时。
	time_t t_ = mktime(&tm_);                  // 将tm结构体转换成time_t格式。
	t_ -= int(second);
	return t_+second;                                 // 返回值。 
}

void LEEBEN::UWB::release_deltaRs()
{//释放map值的内存
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

//double 转日期
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

void LEEBEN::UWB::realtime()//实时解算的第一步
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
			std::cout <<"UDP错误码："<< WSAGetLastError() << std::endl;
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
	std::cout << "已停止广播\n";
	return 1;
}

void LEEBEN::solve_start(LEEBEN::UWB& u)
{
	std::string str,id,log;
	std::istringstream iss;
	double mark;
	std::map<std::string, double> marktable1, marktable2;
	int tdoas;
	std::vector<double> R, Q;//R是每个距离差过滤之后的均值，Q是过滤后的方差，两者大小均为anchor.num-1
	log = "F:\\电子书\\毕设准备\\UWB\\数据\\UWB_LOC程序解算位置.txt";
	//log.append(to_string(clock()));
	u._LOC_record.open(log.c_str(), std::ios::out);

	log = "F:\\电子书\\毕设准备\\UWB\\数据\\UWB_LOC程序解算位置(LS_TLS).txt";
	u._LOC_record_LS_TLS.open(log.c_str(), std::ios::out);

	log = "F:\\电子书\\毕设准备\\UWB\\数据\\UWB_LOC程序解算位置(Chan_pro).txt";
	u._LOC_record_pro.open(log, std::ios::out);

	log = "F:\\电子书\\毕设准备\\UWB\\数据\\去除粗差结果.txt";
	u._Filtered_Tdoa.open(log.c_str(), std::ios::out);
	//log.append
	int send_add_size;
	send_add_size = sizeof(u._sochl.send_add);
	while (true) {
		memset(u._buffer, 0, sizeof(u._buffer));
		recvfrom(u._sochl.sclient, u._buffer, sizeof(u._buffer), 0, (SOCKADDR*)&u._sochl.send_add, &send_add_size);
		str = u._buffer;
		if (str == "leave") {
			std::cout << "已停止解算\n";
			u._receiveOn ^= u._receiveOn;
			u._sochl.end_socket();//释放资源
			std::cout << "已关闭套接字\n";
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
		if (marktable2[id] - marktable1[id]>1	//控制解算一次的时间差
			) {
			//考虑输出和解算
			u.Gross_Free(id);
			u.Gross_Free(id);
			u.Filter_Result(id,R, Q);
			for (int i = 0; i <u. _AnchorNum - 1; i++) {
				R[i] -= u._AnchorsBlunders[i];
				u._deltaRs[id][i].clear();
			}
			std::thread t(GetLoc, std::ref(u), id, R, mark, Q);//开始解算线程
			t.join();
			marktable1[id] = marktable2[id];
		}
		u.push(iss,id);
	}
	u._LOC_record.close();
	u._LOC_record_LS_TLS.close();
	u._LOC_record_pro.close();
	u._Filtered_Tdoa.close();
	std::cout << "已关闭输出文件流\n";
	u.release_deltaRs();
	u._deltaRs.clear();
}

bool LEEBEN::GetLoc(UWB& u,const std::string&id, std::vector<double>& deltaR,double timestamp,std::vector<double>& Q)
{
	Eigen::Vector3d ans;
	std::string t = u.Double2Datetime((time_t)timestamp);
		//chan算法
	ans = chan(u, deltaR, Q);
	u._LOC_record <<id<<','<< t ;
	u._LOC_record_LS_TLS << id << ',' << t ;
	chan_pro(u,ans, deltaR, Q);
	u._LOC_record_pro << id << ',' << t ;
	u._TagPos[id] = ans;
	//屏幕显示
	std::cout << id << ',';
	for (int i = 0; i < ans.size() - 1; i++)
		std::cout << ans(i) << ',';
	std::cout << t;
		//std::cout << "chan's result: \n" << ans << std::endl << std::endl;
		/*if (MS(u, deltaR)) {
			std::cout << "最小二乘迭代不收敛\n";
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
	//normal结果
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
	za=(G.transpose()**///解算z坐标不太符合实际情况，因为z是由线性方程直接解得的，再将其放入平差中将不再起作用，所以平差结果与之前相比不会发生改变
	//u.ans(0, 0) = ans(0,0); u.ans(1,0) = ans(1, 0); u.ans(2,0) = ans(2,0);
	for (int i = 0; i < ans.size(); i++)
		u._LOC_record << ans(i) << ',';
	//LS_TLS结果
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
		u._LOC_record_pro << "不收敛,";
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
{//进行一般实矩阵QR分解的函数
	using std::cout;
	using std::endl;
	int i, j, k, nn, jj,m,n;
	double u, alpha, w, t;
	m = CA.rows(); n = CA.cols();
	A = CA;
	if (m < n) {
		cout << "\nQR分解失败！" << endl; exit(1);
	} //保证列数>行数，才实现QR分解，所以在 
	Q = Eigen::MatrixXd::Identity(m, m);
	nn = n;
	if (m == n) nn = m - 1;
	for (k = 0; k <= nn - 1; k++)//在大循环k：0~m当中，进行H矩阵的求解，左乘Q，以及左乘A 
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
			cout << "\nQR分解失败！" << endl;
			exit(1);
		} u = sqrt(2.0*alpha*(alpha - A(k,k)));
		if ((u + 1.0) != 1.0) {
			A(k,k) = (A(k,k) - alpha) / u;
			for (i = k + 1; i <= m - 1; i++)
				A(i,k) = A(i,k) / u;
			//以上就是H矩阵的求得，实际上程序并没有设置任何数据结构来存储H矩 
			//阵，而是直接将u向量的元素赋值给原A矩阵的原列向量相应的位置 
			//这样做是为了计算左乘矩阵Q和A 
			for (j = 0; j <= m - 1; j++) {
				t = 0.0;
				for (jj = k; jj <= m - 1; jj++)
					t = t + A(jj,k) * Q(jj,j);
				for (i = k; i <= m - 1; i++)
					Q(i,j) = Q(i,j) - 2.0*t*A(i,k);
			} //左乘矩阵Q，循环结束后得到一个矩阵，再将这个矩阵转置一下就得到QR分解中的Q矩阵 
				//也就是正交矩阵 
			for (j = k + 1; j <= n - 1; j++) {
				t = 0.0;
				for (jj = k; jj <= m - 1; jj++)
					t = t + A(jj,k) * A(jj,j);
				for (i = k; i <= m - 1; i++)
					A(i,j) = A(i,j) - 2.0*t*A(i,k);
			} //H矩阵左乘A矩阵，循环完成之后，其上三角部分的数据就是上三角矩阵R
			A(k,k) = alpha;
			for (i = k + 1; i <= m - 1; i++)
				A(i,k)= 0.0;
		}
	}
	Q.transposeInPlace();//QR分解完毕
}

bool LEEBEN::UWB::switch_receive()
{
	if (_receiveOn) {
		broadcast_stop();//先让广播停下来
		//_receiveOn = false;
		SOCKET stoprecv = socket(AF_INET, SOCK_DGRAM, 0);
		if(sendto(stoprecv, "leave", 6, 0, (sockaddr*)&_sochl.recv_add, sizeof(_sochl.recv_add))<0)
			std::cout<<"请求广播停止出现错误，错误码："<<WSAGetLastError()<<std::endl;
		closesocket(stoprecv);
		return true;
	}
	else {
		std::cout << "套接字已关闭\n";
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
	std::cout << "已完成：\n";
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
	}//接收过100个以上数据后，开始解算系统误差
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
	std::cout << "\n解算完成，系统误差取值：";
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
	//将本轮使用的系统误差存入日志文件
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
	//设置本地IP和端口
	recv_add.sin_family = AF_INET;
	//recv_add.sin_addr.S_un.S_addr = INADDR_ANY;
	inet_pton(AF_INET, "127.0.0.1", &recv_add.sin_addr);
	recv_add.sin_port = htons(5510);
}

void LEEBEN::UWB::Soket_Handle::start_socket()
{
	//绑定本地地址与接口
	WSADATA wsadata;
	WSAStartup(MAKEWORD(2, 2), &wsadata);
	sclient = socket(AF_INET, SOCK_DGRAM, 0);
	//设置目标IP与端口
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
