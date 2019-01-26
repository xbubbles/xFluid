
#include "stopwatch.h"

StopWatch::StopWatch()
{
}

void StopWatch::start() {
    
    #if defined(__linux__) || defined(__APPLE__) || defined(__MACOSX)
        struct timeval tp;
        gettimeofday(&tp, nullptr);
        _tbegin = (double)tp.tv_sec + (double)tp.tv_usec / 1000000.0;
    #elif defined(_WIN32)
        _tbegin = (double)GetTickCount() / 1000.0;
    #else
    #endif
    
    _isStarted = true;
}


void StopWatch::stop() {
    if (!_isStarted) {
        return;
    }

    #if defined(__linux__) || defined(__APPLE__) || defined(__MACOSX)
        struct timeval tp;
        gettimeofday(&tp, nullptr);
        _tend = (double)tp.tv_sec + (double)tp.tv_usec / 1000000.0;
    #elif defined(_WIN32)
        _tend = (double)GetTickCount() / 1000.0;
    #else
    #endif
    
    double time = _tend - _tbegin;
    _timeRunning += time;
}

void StopWatch::reset() {
    _isStarted = false;
    _timeRunning = 0.0;
}

double StopWatch::getTime() {
    return _timeRunning >= 0.0 ? _timeRunning : 0.0;
}