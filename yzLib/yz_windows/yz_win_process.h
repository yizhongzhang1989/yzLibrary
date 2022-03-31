/***********************************************************/
/**	\file
	\brief		Process in Windows
	\details	All string are treated as ASCII, we don't use unicode.
				This file include :
				windows process abstraction,
				shared memory
				pipe
	\author		Yizhong Zhang
	\date		5/28/2012
*/
/***********************************************************/
#ifndef __YZ_WIN_PROCESS_H__
#define __YZ_WIN_PROCESS_H__

#include "yzLib/yz_setting.h"

#ifndef YZ_windows_h
#	error yz_win_process.h must be included after windows.h
#endif

#include <iostream>
#include "yzLib/yz_windows/yz_win_utils.h"

namespace yz{

/**
	namespace contain windows related classes and functions
*/
namespace windows{
/**
	Windows Process Abstraction

	This class will make it easy for us to create another process
	in our program. To do this, you should

	1,	Create an instance of WinProcess class \n
	2,	Call Init to pass the process name to WinProcess \n
	3,	Call Start() to start the process

	The above three steps can also be achieved by the constructor

	If we want to wait the process to terminate, call Wait() and 
	the function will not return until this process terminates.
*/
class WinProcess{
public:
	char*				name;

public:
	WinProcess(){
		name = NULL;
		Reset();
	}
	/**
		Constructor, can start the process directly

		\param	process_name	[directory] + process_name + [command line]
		\param	start_flag		whether start the process, true: start now; false: start later
	*/
	WinProcess(const char* process_name, bool start_flag=true){
		name = NULL;
		Reset();
		Init(process_name);
		if( start_flag )
			Start();
	}
	/**
		Initialize WinProcess by giving process name and directory

		\param	process_name	[directory] + process_name + [command line]
		\return					1: succeed;		2: failed	
	*/
	inline int Init(const char* process_name){
		if( status & WINPROCESS_STATUS_BIT_INITIALIZED ){
			#ifndef BE_QUIET
				std::cout << "win32process : " << name << " already initialized" << std::endl;
			#endif
			return 0;
		}

		Reset();
		int len = strlen(process_name);
		if( len <= 0 ){
			#ifndef BE_QUIET
				std::cout << "invalid process name" << std::endl;
			#endif
			return 0;
		}
		name = new char[len+1];
		strcpy(name, process_name);

		status |= WINPROCESS_STATUS_BIT_INITIALIZED;
		return 1;
	}

	/**
		Start the process immediantly

		\return			1: succeed;		2: failed	
	*/
	inline int	Start(){
		if( ! (status & WINPROCESS_STATUS_BIT_INITIALIZED) ){
			#ifndef BE_QUIET
				std::cout << "win32process : " << name << " not initialized" << std::endl;
			#endif
			return 0;
		}
		if( status & WINPROCESS_STATUS_BIT_STARTED ){
			#ifndef BE_QUIET
				std::cout << "win32process : " << name << " has started, cannot start again" << std::endl;
			#endif
			return 0;
		}

		if( !CreateProcessA( NULL,	// No module name (use command line)
			name,			// Command line
			NULL,           // Process handle not inheritable
			NULL,           // Thread handle not inheritable
			FALSE,          // Set handle inheritance to FALSE
			0,              // No creation flags
			NULL,           // Use parent's environment block
			NULL,           // Use parent's starting directory 
			&si,            // Pointer to STARTUPINFO structure
			&pi )           // Pointer to PROCESS_INFORMATION structure
			) {
				#ifndef BE_QUIET
					std::cout << "start process " << name << " failed" << std::endl;
				#endif
				PrintLastError();
				return 0;
		}

		status |= WINPROCESS_STATUS_BIT_STARTED;

		return 1;
	}

	/**
		Wait until the process terminate

		\return			1: succeed;		2: failed
	*/
	inline int Wait(){
		if( ! (status & WINPROCESS_STATUS_BIT_INITIALIZED) ){
			#ifndef BE_QUIET
				std::cout << "win32process : " << name << " not initialized" << std::endl;
			#endif
			return 0;
		}
		if( ! (status & WINPROCESS_STATUS_BIT_STARTED) ){
			#ifndef BE_QUIET
				std::cout << "win32process : " << name << " has started, cannot start again" << std::endl;
			#endif
			return 0;
		}

		// Wait until child process exits.
		if( WaitForSingleObject( pi.hProcess, INFINITE ) == WAIT_FAILED )
			PrintLastError();

		status &= ~WINPROCESS_STATUS_BIT_STARTED;
		return 1;
	} 

	/**
		Close the process, clear everything

		If the process is running, then terminate the process.
	*/
	inline void Close(){
		if( status & WINPROCESS_STATUS_BIT_STARTED ){
			//	the process has started, so we check whether it has terminated
			DWORD exit_code;

			if(!GetExitCodeProcess(pi.hProcess, &exit_code))	PrintLastError();

			if( exit_code == STILL_ACTIVE )
				if(!TerminateProcess(pi.hProcess, 0))	PrintLastError();

			if(!CloseHandle( pi.hProcess ))	PrintLastError();
			if(!CloseHandle( pi.hThread ))	PrintLastError();
		}

		Reset();
	}

	/**
		Print help information
	*/
	inline void Help(){
		#ifndef BE_QUIET
			std::cout
				<< "WinProcess Help:\n"
				<< "create a new process, just construct a WinProcess:\n"
				<< "\tWinProcess pros(\"reciever.exe\", 1);    // 1 indicate start directly\n"
				<< std::endl;
		#endif
	}

protected:
	STARTUPINFOA		si;
	PROCESS_INFORMATION	pi;
	int					status;
	/**
		Status Bit of WinProcess, default 0
	*/
	enum{
		WINPROCESS_STATUS_BIT_INITIALIZED	= 0x01,
		WINPROCESS_STATUS_BIT_STARTED		= 0x02,
	};

protected:
	/**
		Reset the class, all values are clear
	*/
	inline void Reset(){
		if( name )	delete[] name;
		name = NULL;

		ZeroMemory( &si, sizeof(si) );
		si.cb = sizeof(si);
		ZeroMemory( &pi, sizeof(pi) );

		status = 0;
	}

};


/**
	Windows Shared Memory Abstraction

	This class can create a named shared memory, to be used by 
	different processes to transfer data.

	1,	Process 1 call CreateSharedMemory() to create a shared 
		memory space with name \n
	2,	Process 2 call OpenSharedMemory() to get the shared memory, 
		their name and size must match 1 \n
	3,	The two processes can transfer data now by calling CopyToSharedMemory()
*/
class WinSharedMemory{
public:
	char*			name;	///<name of the shared mamory
	int				size;	///<size of buffer in bytes
	unsigned char*	buffer;	///<buffer

public:
	WinSharedMemory(){
		name = NULL;
		Reset();
	}
	/**
		Create a Shared Memory, called by Process 1

		\param	shared_memory_name	name of the shared memory
		\param	buffer_size			size of the shared memory in bytes
		\return						1: succeed;		2: failed
	*/
	inline int CreateSharedMemory(const char* shared_memory_name, int buffer_size){
		if( status & SHARED_MEMORY_STATUS_BIT_CREATED ){
			#ifndef BE_QUIET
				std::cout << "this class has been used to create, cannot create again" << std::endl;
			#endif
			return 0;
		}
		if( status & SHARED_MEMORY_STATUS_BIT_OPENED ){
			#ifndef BE_QUIET
				std::cout << "this class has been used to open, cannot create" << std::endl;
			#endif
			return 0;
		}

		int len = strlen(shared_memory_name);
		if( len <= 0 ){
			#ifndef BE_QUIET
				std::cout << "invalid shared memory name" << std::endl;
			#endif
			return 0;
		}
		name = new char[len];
		strcpy(name, shared_memory_name);

		size = buffer_size;

		//	create mapping and check
		hMapFile = CreateFileMappingA(
			INVALID_HANDLE_VALUE,	// use paging file
			NULL,					// default security
			PAGE_READWRITE,			// read/write access
			0,						// maximum object size (high-order DWORD)
			size,					// maximum object size (low-order DWORD)
			(LPSTR)name);			// name of mapping object
		if(hMapFile == NULL){
			#ifndef BE_QUIET
				std::cout << "Could not create file mapping object creating shared memory: " << name << ", error: " << GetLastError() << std::endl;
			#endif
			PrintLastError();
			Reset();
			return 0;
		}

		//	map buffer
		buffer = (unsigned char*) MapViewOfFile(hMapFile,   // handle to map object
			FILE_MAP_ALL_ACCESS, // read/write permission
			0,
			0,
			size);
		if(buffer == NULL){
			#ifndef BE_QUIET
				std::cout << "Could not map view of file creating shared memory: " << name << ", error: " << GetLastError() << std::endl;
			#endif
			PrintLastError();
			if(!CloseHandle(hMapFile))	PrintLastError();
			Reset();
			return 0;
		}

		status |= SHARED_MEMORY_STATUS_BIT_CREATED;
		return 1;
	}

	/**
		Open an existing Shared Memory, called by Process 2

		\param	shared_memory_name	name of the shared memory
		\param	buffer_size			size of the shared memory in bytes
		\return						1: succeed;		2: failed
	*/
	inline int OpenSharedMemory(const char* shared_memory_name, int buffer_size){
		if( status & SHARED_MEMORY_STATUS_BIT_OPENED ){
			#ifndef BE_QUIET
				std::cout << "this class has been used to open, cannot open again" << std::endl;
			#endif
			return 0;
		}
		if( status & SHARED_MEMORY_STATUS_BIT_CREATED ){
			#ifndef BE_QUIET
				std::cout << "this class has been used to create, cannot open" << std::endl;
			#endif
			return 0;
		}

		int len = strlen(shared_memory_name);
		if( len <= 0 ){
			#ifndef BE_QUIET
				std::cout << "invalid shared memory name" << std::endl;
			#endif
			return 0;
		}
		name = new char[len];
		strcpy(name, shared_memory_name);
		size = buffer_size;

		//	open and check
		hMapFile = OpenFileMappingA(
			FILE_MAP_ALL_ACCESS,	// read/write access
			FALSE,					// do not inherit the name
			(LPSTR)name);			// name of mapping object
		if(hMapFile == NULL){
			#ifndef BE_QUIET
				std::cout << "Could not open file mapping object open shared memory: " << name << ", error: " << GetLastError() << std::endl;
			#endif
			PrintLastError();
			Reset();
			return 0;
		}

		//	map buffer
		buffer = (unsigned char*) MapViewOfFile(hMapFile, // handle to map object
			FILE_MAP_ALL_ACCESS,  // read/write permission
			0,
			0,
			size);
		if(buffer == NULL){
			#ifndef BE_QUIET
				std::cout << "Could not map view of file open shared memory: " << name << ", error: " << GetLastError() << std::endl;
			#endif
			PrintLastError();
			if(!CloseHandle(hMapFile))	PrintLastError();
			Reset();
			return 0;
		}

		status |= SHARED_MEMORY_STATUS_BIT_OPENED;
		return 1;

	}

	/**
		Copy Data to this Shared Memory. Size of this memory cannot
		exceed size of buffer, or print error information and return failed

		\param	ptr			data buffer pointer
		\param	mem_size	size of data in bytes
		\return				1: succeed;		2: failed
	*/
	inline int CopyToSharedMemory(const void* ptr, int mem_size){
		if( mem_size > size ){
			#ifndef BE_QUIET
				std::cout << "size: " << mem_size << " bigger than shared memory size: " << size << std::endl;
			#endif
			return 0;
		}

		CopyMemory((PVOID)buffer, ptr, size);
		return 1;
	}

	/**
		Remove this Shared Memory
	*/
	inline void RemoveSharedMemory(){
		if( buffer != NULL )
			if(!UnmapViewOfFile(buffer))	PrintLastError();
		if( hMapFile != INVALID_HANDLE_VALUE )
			if(!CloseHandle(hMapFile))	PrintLastError();
		Reset();
	}

	/**
		Print Help Information
	*/
	inline void Help(){
		#ifndef BE_QUIET
			std::cout
				<< "SharedMemory Help:\n"
				<< "In client process code, create shared memory:\n"
				<< "\tSharedMemory sm;\n"
				<< "\tsm.CreateSharedMemory(\"shared_memory_name\", 512);\n"
				<< "\tsm.CopyToSharedMemory(buf, strlen(buf)+1);\n"
				<< std::endl;

			std::cout
				<< "In server process code, open shared memory:\n"
				<< "\tSharedMemory sm;\n"
				<< "\tsm.OpenSharedMemory(\"shared_memory_name\", 512);\n"
				<< std::endl;
		#endif
	}

protected:
	HANDLE			hMapFile;
	int				status;
	/**
		Status Bit of SharedMemory, default 0
	*/	
	enum{
		SHARED_MEMORY_STATUS_BIT_CREATED	= 0x01,
		SHARED_MEMORY_STATUS_BIT_OPENED		= 0x02
	};

protected:
	/**
		Reset the class, all values are clear
	*/
	inline void Reset(){
		if(name) delete[] name;
		name		= NULL;
		size		= 0;
		buffer		= NULL;
		hMapFile	= INVALID_HANDLE_VALUE;
		status		= 0;
	}
};


/**
	Windows Named Pipe Abstraction

	Named pipe is used to communicate between processes. Pipe
	use Client-Server communication mode. This class set the 
	pipe to be bi-directional, which means that data can transfer
	both direction. So Client, Server is just logical.

	Follow the following steps to set a typical C-S pipe model

	1,	Process 1 call CreatePipeServer() to set process 1 as 
		pipe server \n
	2,	Process 2 call CreatePipeClient() to set process 2 as
		pipe client \n
	3,	Process 1 start a loop, wait to Read() then Write() 
		according to read \n
	4,	Process 2 call Write() when necessary then wait to Read()

	Read() Write() must match on a single pipe, because they will
	wait until success. So if two processes call Read() at the same
	time, dead lock will appear because pipe has been taken by Read(),
	but data will never come.
*/
class WinNamedPipe{
public:
	char*	name;	///<The unique pipe name. This string must have the following form: "\\.\pipe\pipename"

public:
	WinNamedPipe(){
		name = NULL;
		Reset();
	}
	/**
		Create a named pipe on the server side and wait for client to connect

		"\\.\pipe\" will be added to pipe_name automatically, so no need to include
		this prefix in pipe_name

		\param	pipe_name		the name of the pipe
		\return					whether the pipe has been created successfully
								1: success; 0: failed
	*/
	inline int CreatePipeServer(const char* pipe_name){
		if( status & NAMED_PIPE_STATUS_BIT_SET_AS_SERVER ){
			#ifndef BE_QUIET
				std::cout << "This pipe class has been set as Pipe Server" << std::endl;
			#endif
			return 0;
		}
		if( status & NAMED_PIPE_STATUS_BIT_SET_AS_CLIENT ){
			#ifndef BE_QUIET
				std::cout << "This pipe class has been set as Pipe Client" << std::endl;
			#endif
			return 0;
		}

		int len = strlen(pipe_name);
		if( len <= 0 ){
			#ifndef BE_QUIET
				std::cout << "invalid pipe_name" << std::endl;
			#endif
			return 0;
		}
		name = new char[len+9];
		strcpy(name, "\\\\.\\pipe\\");	//add prefix \\.\pipe\ to the name
		strcat(name, pipe_name);

		//	create named pipe
		pipe = CreateNamedPipeA(
			name,						//	name of pipe, \\.\pipe\pipe_name
			PIPE_ACCESS_DUPLEX,			//	open mode, bi-directional
			PIPE_TYPE_MESSAGE | PIPE_READMODE_MESSAGE | PIPE_WAIT,	//	pipe mode, message
			PIPE_UNLIMITED_INSTANCES,	//	maximal instances
			4096,						//	out buffer size, the sized used in MSDN demo
			4096,						//	in buffer size, the sized used in MSDN demo
			0,							//	time out, default (zero value): 50ms
			NULL );						//	security attribute

		if( pipe == INVALID_HANDLE_VALUE ){
			#ifndef BE_QUIET
				std::cout << "create named pipe: " << name << " failed" << std::endl;
			#endif
			PrintLastError();
			Reset();
			return 0;
		}

		//	wait for client to connect
		BOOL result = ConnectNamedPipe(pipe, NULL);
		if( !result && GetLastError()!=ERROR_PIPE_CONNECTED ){	//	it is possible that the client has connected already
			#ifndef BE_QUIET
				std::cout << "connect named pipe: " << name << " failed" << std::endl;
			#endif
			PrintLastError();
			if(!CloseHandle(pipe))	PrintLastError();
			Reset();
			return 0;
		}

		status |= NAMED_PIPE_STATUS_BIT_SET_AS_SERVER;
		return 1;
	}

	/**
		Create a named pipe client and connect to the named pipe.

		It is possible that the pipe doesn't exist, so try several times 

		\param	pipe_name		the name of the pipe, added prefix automatically
		\param	max_try_count	how many times trying to connect to the pipe
		\param	try_delay_ms	delay time between each try
		\return					how many times tried to connect before succeed, 0: failed
	*/
	inline int CreatePipeClient(const char* pipe_name, int max_try_count=50, int try_delay_ms=100){
		if( status & NAMED_PIPE_STATUS_BIT_SET_AS_SERVER ){
			#ifndef BE_QUIET
				std::cout << "This pipe class has been set as Pipe Server" << std::endl;
			#endif
			return 0;
		}
		if( status & NAMED_PIPE_STATUS_BIT_SET_AS_CLIENT ){
			#ifndef BE_QUIET
				std::cout << "This pipe class has been set as Pipe Client" << std::endl;
			#endif
			return 0;
		}

		int len = strlen(pipe_name);
		if( len <= 0 ){
			#ifndef BE_QUIET
				std::cout << "invalid pipe_name" << std::endl;
			#endif
			return 0;
		}
		name = new char[len+9];
		strcpy(name, "\\\\.\\pipe\\");	//add prefix \\.\pipe\ to the name
		strcat(name, pipe_name);

		//	try to connect to an existing pipe
		int try_count = 0;
		while(try_count < max_try_count){	//	we try max_try_count times at most
			pipe = CreateFileA(
				name,			//	pipe name
				GENERIC_READ | GENERIC_WRITE,	//	R/W access
				0,				//	no sharing
				NULL,			//	default security attributes
				OPEN_EXISTING,	//	open existing pipe
				0,				//	default attributes
				NULL );			//	no template file

			try_count ++;

			if( pipe==INVALID_HANDLE_VALUE ){	//	if didn't succeed this time, delay then retry
				if(try_count < max_try_count){
					if( GetLastError() == ERROR_PIPE_BUSY )
						WaitNamedPipeA(name, try_delay_ms);
					else
						Sleep(try_delay_ms);
				}
			}
			else
				break;
		}

		if( pipe == INVALID_HANDLE_VALUE ){
			#ifndef BE_QUIET
				std::cout << "open named pipe: " << name << " failed" << std::endl;
			#endif
			PrintLastError();
			Reset();
			return 0;
		}

		//	set pipe status
		DWORD mode = PIPE_READMODE_MESSAGE;
		BOOL success = SetNamedPipeHandleState(
			pipe,	//	pipe handle
			&mode,	//	new pipe mode
			NULL,	//	don't set maximum bytes
			NULL );	//	don't set maximum times
		if( !success ){
			#ifndef BE_QUIET
				std::cout << "set named pipe: " << name << " handle state failed" << std::endl;
			#endif
			PrintLastError();
			if( !CloseHandle(pipe) )	PrintLastError();
			Reset();
			return 0;
		}

		status |= NAMED_PIPE_STATUS_BIT_SET_AS_CLIENT;
		return try_count;	//	at least 1
	}

	/**
		Read data from pipe to buffer

		If the pipe has been connected, this function will wait until
		the other process calls Write()

		Since one read cannot exceed the limit of buffer_size, more_data_flag
		is set to 1 if there are more data to be read. So we don't know the size 
		of data to read, we should call Read in a loop, and break when more_data_falg is zero.

		\param	buffer			the buffer to hold data
		\param	buffer_size		size of the buffer
		\param	more_data_flag	if more data are remain to read, set this flag to 1
		\return					how many bytes read
	*/
	inline int Read(char* buffer, int buffer_size, int& more_data_flag){
		more_data_flag = 0;

		if( !(status & NAMED_PIPE_STATUS_BIT_SET_AS_SERVER) && 
			!(status & NAMED_PIPE_STATUS_BIT_SET_AS_CLIENT) ){
				#ifndef BE_QUIET
					std::cout << "This pipe class has not setup yet, cannot read" << std::endl;
				#endif
				return 0;
		}

		DWORD read_bytes = 0;
		BOOL success = ReadFile(
			pipe,
			buffer,
			buffer_size,
			&read_bytes,
			NULL );

		if( !success ){
			if( GetLastError()==ERROR_MORE_DATA )
				more_data_flag = 1;
			else{
				#ifndef BE_QUIET
					std::cout << "failed to read from pipe" << std::endl;
				#endif
				//PrintLastError();
				return 0;
			}
		}

		return read_bytes;
	}

	/**
		Write data buffer to pipe

		If the pipe has been connected, this function will wait until 
		the other process calls Read()

		\param	buffer		data buffer
		\param	size		size of data in bytes
		\return				how many bytes really written
	*/
	inline int Write(char* buffer, int size){
		if( size<=0 )	return 0;

		if( !(status & NAMED_PIPE_STATUS_BIT_SET_AS_SERVER) && 
			!(status & NAMED_PIPE_STATUS_BIT_SET_AS_CLIENT) ){
				#ifndef BE_QUIET
					std::cout << "This pipe class has not setup yet, cannot write" << std::endl;
				#endif
				return 0;
		}

		DWORD written_bytes = 0;
		BOOL success = WriteFile(
			pipe,
			buffer,
			size,
			&written_bytes,
			NULL );

		if( !success ){
			#ifndef BE_QUIET
				std::cout << "failed to write to pipe" << std::endl;
			#endif
			PrintLastError();
			return 0;
		}

		return written_bytes;
	}

	/**
		close the pipe

		call this function when the pipe is no longer needed
	*/
	inline void ClosePipe(){
		if( !(status & NAMED_PIPE_STATUS_BIT_SET_AS_SERVER) && 
			!(status & NAMED_PIPE_STATUS_BIT_SET_AS_CLIENT) ){
				#ifndef BE_QUIET
					std::cout << "This pipe class has not setup yet, cannot close" << std::endl;
				#endif
				return;
		}

		if( pipe != INVALID_HANDLE_VALUE ){
			if( !FlushFileBuffers(pipe) )		PrintLastError();
			if( !DisconnectNamedPipe(pipe) )	PrintLastError();
			if( !CloseHandle(pipe) )			PrintLastError();
		}
		Reset();
	}

	/**
		Print help information of WinNamedPipe class
	*/
	inline void Help(){
		#ifndef BE_QUIET
			std::cout
				<< "Pipe is used to transfer command between processes.\n"
				<< "In process 1, call CreatePipeServer() to set the process as pipe server.\n"
				<< "In process 2, call CreatePipeClient() to set the process as pipe client.\n"
				<< "After that, the processes can transfer data by Read() and Write()\n"
				<< "Read() Write() must match, or the program will reach a dead lock.\n"
				<< "Call ClosePipe() when the pipe is no longer needed"
				<< std::endl;
		#endif
	}

protected:
	HANDLE	pipe;	///<handle to the pipe
	int status;		///<Status bit of the class
	/**
		Status Bit of Named Pipe, default 0
	*/	
	enum{
		NAMED_PIPE_STATUS_BIT_SET_AS_SERVER	= 0x01,
		NAMED_PIPE_STATUS_BIT_SET_AS_CLIENT	= 0x02
	};

protected:
	/**
		Reset all member variables
	*/
	inline void Reset(){
		if(name) delete[] name;
		name		= NULL;
		pipe		= INVALID_HANDLE_VALUE;
		status		= 0;
	}

};

/**
	Windows Serial Communication Abstraction

	Communication through serial port, commonly used for communicate with micro controllers

	how to use:	\n
	1, call Open(port_name, baud_rate), such as ("COM1", 9600)	\n
	2, call Read(), Write()	\n
	3, call Purge() to clear buffer	\n
*/
class WinSerialPort {
public:
	WinSerialPort() {
		h_comm = NULL;
	}

	~WinSerialPort() {
		if (h_comm)
			CloseHandle(h_comm);
	}

	/**
	Open the serial port

	\param	com_port_name	name of the port, typical format: COM1
	\param	baud_rate		baud_rate of this port
	\return					whether open port succeed
	*/
	int Open(const char* com_port_name, int baud_rate) {
		//	open the port
		h_comm = CreateFileA(
			com_port_name,
			GENERIC_READ | GENERIC_WRITE,
			0,
			NULL,
			OPEN_EXISTING,
			0,
			NULL);

		if (h_comm == INVALID_HANDLE_VALUE) {
			std::cout << "error: WinSerialComPort::Open, cannot open port: " << com_port_name << std::endl;
			h_comm = NULL;
			return 0;
		}

		// set timeouts
		COMMTIMEOUTS cto = { MAXDWORD, 0, 0, 0, 0 };
		if (!SetCommTimeouts(h_comm, &cto)) {
			std::cout << "error: WinSerialComPort::Open, cannot set timeout" << std::endl;
			CloseHandle(h_comm);
			h_comm = NULL;
			return 0;
		}

		// set baud rate by DCB
		memset(&dcb, 0, sizeof(dcb));
		dcb.DCBlength = sizeof(dcb);
		dcb.BaudRate = baud_rate;
		dcb.fBinary = 1;
		dcb.fDtrControl = DTR_CONTROL_ENABLE;
		dcb.fRtsControl = RTS_CONTROL_ENABLE;
		dcb.Parity = NOPARITY;
		dcb.StopBits = ONESTOPBIT;
		dcb.ByteSize = 8;

		if (!SetCommState(h_comm, &dcb)) {
			std::cout << "error: WinSerialComPort::Open, cannot set DCB" << std::endl;
			CloseHandle(h_comm);
			h_comm = NULL;
			return 0;
		}

		//	record port name
		port_name = com_port_name;

		return 1;
	}

	/**
	Write data to the port

	\param	data		the data buffer
	\param	data_len	length of the data, if length is negative, data will be treated as a string and calculate length
	\return				whether write succeed
	*/
	int Write(const char* data, int data_len = -1) {
		if (!h_comm)
			return 0;

		//	calculate data length (if not input)
		if (data_len < 0)
			data_len = strlen(data);
		if (data_len == 0)
			return 1;

		DWORD dwNumBytesWritten;
		BOOL ret = WriteFile(
			h_comm,
			data,
			data_len,
			&dwNumBytesWritten,
			NULL);

		if (!ret) {
			std::cout << "error: WinSerialPort::Write, write failed" << std::endl;
			return 0;
		}

		return 1;
	}

	/**
	Read data from the port

	\param	buffer			the buffer to hold data
	\param	buffer_len		size of the buffer
	\param	zero_pending	whether append '\0' to the end of buffer, so it can be used as a string directly
	\return					whether read succeed
	*/
	int Read(char *buffer, int buffer_len, int zero_pending = 1) {
		if (!h_comm)
			return 0;

		if (zero_pending)
			buffer_len--;

		DWORD numRead;
		BOOL ret = ReadFile(
			h_comm,
			buffer,
			buffer_len,
			&numRead,
			NULL);

		if (!ret)
			return 0;

		if (zero_pending)
			buffer[numRead] = '\0';

		return numRead;
	}

	/**
	Clear buffer
	*/
	void Purge() {
		char buf[8];
		while (Read(buf, 8));
	}

	/**
		Get current baud rate

		\return		current baud rate. -1 if the port is not set
	*/
	int GetBaudRate() {
		if (h_comm)
			return dcb.BaudRate;
		return -1;
	}

public:
	std::string port_name;
	HANDLE		h_comm;

protected:
	DCB			dcb;
};


}	//	namespace windows
}	//	namespace yz

#endif	//	__YZ_WIN_PROCESS_H__