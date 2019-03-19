// UWB_LOC.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include"Eigen\Dense"
#include<string>
#include<thread>
#include"WinSock2.h"
#include <stdio.h>   
#include< Ws2tcpip.h> 
#pragma comment(lib,"WS2_32.lib")
#include"UWB.h"

using namespace LEEBEN;
using namespace Eigen;


int main()
{
	Vector3d p[5];
	p[0] << -73., 2., 290.;
	p[1] << 110., -181., 290.;
	p[2] << -258., 292., 290.;
	p[3] << -239.7, -228.7, 290.;
	p[4] << 148.6, 296.7, 290.;
	UWB mytest1;
	mytest1.Initialize_Anchors(5, p);
	//以下是多线程动态解算
	mytest1.show_menu();
	char x; std::cin >> x;
	if (x == 'a') {
		mytest1.broadcast_begin();
		std::cout << "broadcast starts.\n";
		mytest1.realtime();
		mytest1.show_menu2();
		std::thread t(solve_start, std::ref(mytest1));
		t.detach();
		while (std::cin >> x) {
			if (x == 'a') {
				if (mytest1.broadcast_begin()) {
					std::thread t(solve_start, std::ref(mytest1));
					t.detach();
				}
			}
			else if (x == 'b') {
				mytest1.switch_receive();
			}
			else if (x == 'c') {
				break;
			}
			else
				mytest1.show_menu2();
		}
	}
	else {
		std::string filename, id;
		std::cout << "enter the input file:\n";
		std::cin >> filename;
		mytest1.fromfile(filename.c_str(), id.c_str());
	}
	system("pause");
    return 0;
}

