#include <iostream>
#include <string>
#include <Windows.h>
#include <GL/glut.h>
#include <yzLib/yz_lib.h>

//namespace yz {	namespace windows {
//
///**
//	Windows Serial Communication Abstraction
//
//	Communication through serial port, commonly used for communicate with micro controllers
//
//	how to use:	\n
//	1, call Open(port_name, baud_rate), such as ("COM1", 9600)	\n
//	2, call Read(), Write()	\n
//	3, call Purge() to clear buffer	\n
//*/
//class WinSerialPort {
//public:
//	WinSerialPort() {
//		h_comm = NULL;
//	}
//
//	~WinSerialPort() {
//		if (h_comm)
//			CloseHandle(h_comm);
//	}
//
//	/**
//		Open the serial port
//
//		\param	com_port_name	name of the port, typical format: COM1
//		\param	baud_rate		baud_rate of this port
//		\return					whether open port succeed
//	*/
//	int Open(const char* com_port_name, int baud_rate) {
//		//	open the port
//		h_comm = CreateFileA(
//			com_port_name, 
//			GENERIC_READ | GENERIC_WRITE, 
//			0, 
//			NULL, 
//			OPEN_EXISTING,
//			0, 
//			NULL);
//
//		if (h_comm == INVALID_HANDLE_VALUE) {
//			std::cout << "error: WinSerialComPort::Open, cannot open port: " << com_port_name << std::endl;
//			h_comm = NULL;
//			return 0;
//		}
//
//		// set timeouts
//		COMMTIMEOUTS cto = { MAXDWORD, 0, 0, 0, 0 };
//		if (!SetCommTimeouts(h_comm, &cto)) {
//			std::cout << "error: WinSerialComPort::Open, cannot set timeout" << std::endl;
//			CloseHandle(h_comm);
//			h_comm = NULL;
//			return 0;
//		}
//
//		// set baud rate by DCB
//		memset(&dcb, 0, sizeof(dcb));
//		dcb.DCBlength = sizeof(dcb);
//		dcb.BaudRate = baud_rate;
//		dcb.fBinary = 1;
//		dcb.fDtrControl = DTR_CONTROL_ENABLE;
//		dcb.fRtsControl = RTS_CONTROL_ENABLE;
//		dcb.Parity = NOPARITY;
//		dcb.StopBits = ONESTOPBIT;
//		dcb.ByteSize = 8;
//
//		if (!SetCommState(h_comm, &dcb)) {
//			std::cout << "error: WinSerialComPort::Open, cannot set DCB" << std::endl;
//			CloseHandle(h_comm);
//			h_comm = NULL;
//			return 0;
//		}
//
//		//	record port name
//		port_name = com_port_name;
//
//		return 1;
//	}
//
//	/**
//		Write data to the port
//
//		\param	data		the data buffer
//		\param	data_len	length of the data, if length is negative, data will be treated as a string and calculate length
//		\return				whether write succeed
//	*/
//	int Write(const char* data, int data_len = -1) {
//		if (!h_comm)
//			return 0;
//
//		//	calculate data length (if not input)
//		if (data_len < 0)
//			data_len = strlen(data);
//		if (data_len == 0)
//			return 1;
//
//		DWORD dwNumBytesWritten;
//		BOOL ret = WriteFile(
//			h_comm, 
//			data, 
//			data_len, 
//			&dwNumBytesWritten, 
//			NULL);
//
//		if (!ret) {
//			std::cout << "error: WinSerialPort::Write, write failed" << std::endl;
//			return 0;
//		}
//
//		return 1;
//	}
//
//	/**
//		Read data from the port
//
//		\param	buffer			the buffer to hold data
//		\param	buffer_len		size of the buffer
//		\param	zero_pending	whether append '\0' to the end of buffer, so it can be used as a string directly
//		\return					whether read succeed
//	*/
//	int Read(char *buffer, int buffer_len, int zero_pending = 1) {
//		if (!h_comm)
//			return 0;
//
//		if (zero_pending)
//			buffer_len--;
//
//		DWORD numRead;
//		BOOL ret = ReadFile(
//			h_comm, 
//			buffer, 
//			buffer_len, 
//			&numRead, 
//			NULL);
//
//		if (!ret)
//			return 0;
//
//		if (zero_pending)
//			buffer[numRead] = '\0';
//
//		return numRead;
//	}
//
//	/**
//		Clear buffer
//	*/
//	void Purge() {
//		char buf[8];
//		while (Read(buf, 8));
//	}
//
//public:
//	std::string port_name;
//	HANDLE		h_comm;
//
//protected:
//	DCB			dcb;
//};
//
//
//}}

yz::windows::WinSerialPort		com_port;

yz::opengl::DemoWindowManager	manager;
yz::opengl::DemoWindow2D		win3d;

std::string						message;

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		exit(0);
	default:
		com_port.Write((char*)&key, 1);
	}
}

void print() {
	glColor3f(0, 0, 0);
	yz::opengl::printInfo(0, 0, "echo: %s", message.c_str());
}

void idle() {
	char buf[128];
	int len = com_port.Read(buf, 128);
	if (len) {
		message = buf;
	}
}

int main() {
	com_port.Open("COM3", 19200);

	win3d.keyboardFunc = keyboard;
	win3d.SetDrawAppend(print);
	win3d.CreateGLUTWindow();

	manager.AddIdleFunc(win3d.idleFunc);
	manager.AddIdleFunc(idle);
	manager.EnterMainLoop();
}