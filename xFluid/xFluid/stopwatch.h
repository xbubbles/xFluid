
#ifndef STOPWATCH_H
#define STOPWATCH_H

#if defined(__linux__) || defined(__APPLE__) || defined(__MACOSX)
    #include <sys/time.h>
#elif defined(_WIN32)
    #include <Windows.h>
    #include <Winbase.h>
#else
#endif

class StopWatch
{
public:
    StopWatch();
    void start();
    void stop();
    void reset();
    double getTime();    // in seconds

private:
    bool _isStarted = false;
    double _tbegin, _tend;
    double _timeRunning = 0.0;
};

#endif