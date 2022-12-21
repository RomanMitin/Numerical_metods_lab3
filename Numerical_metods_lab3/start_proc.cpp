#include <Windows.h>
#include <iostream>


void start_proc()
{
    system("start excel out.csv");
    //LPCWSTR app_name = L"notepad";
    //auto tmp_arg = L"Excel2013.lnk";
    //auto arg = new wchar_t[128];


    //memcpy(arg, tmp_arg, 14 * sizeof(wchar_t));

    ////printf("arg: %s", arg);

    ////system(arg);

    ////LPSECURITY_ATTRIBUTES
    //STARTUPINFO si;
    //PROCESS_INFORMATION pi;

    ////memset(&si, 0, sizeof(si));
    ////memset(&pi, 0, sizeof(pi));
    //ZeroMemory(&si, sizeof(si));
    //si.cb = sizeof(si);
    //ZeroMemory(&pi, sizeof(pi));

    //if (!CreateProcessW(NULL, arg, NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi))
    //    std::cerr << "Error on excel start\n";

    //auto err = GetLastError();
    //std::cout << err;

    //CloseHandle(pi.hProcess);
    //CloseHandle(pi.hThread);

    ////if (!CreateProcess(NULL,   // No module name (use command line)
    ////    argv[1],        // Command line
    ////    NULL,           // Process handle not inheritable
    ////    NULL,           // Thread handle not inheritable
    ////    FALSE,          // Set handle inheritance to FALSE
    ////    0,              // No creation flags
    ////    NULL,           // Use parent's environment block
    ////    NULL,           // Use parent's starting directory 
    ////    &si,            // Pointer to STARTUPINFO structure
    ////    &pi)           // Pointer to PROCESS_INFORMATION structure
    ////    )
}