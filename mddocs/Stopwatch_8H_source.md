
# File Stopwatch.H

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**Stopwatch.H**](Stopwatch_8H.md)

[Go to the documentation of this file.](Stopwatch_8H.md) 


````cpp
#ifndef _STOPWATCH__
#define _STOPWATCH__

#define USE_STOPWATCH
//#define DEBUG_STOPWATCH_STOPLAP
#include "AMReX_ParallelDescriptor.H"
#include <fstream>

class Stopwatch{
public:
    static void init(int plevel=0,int wlevel=0,std::string logname="stopwatch_log.dat");
    static void start(int level=0);
    static void stop(std::string msg,int level=0);
    static void startlap(std::string msg="", int level=0, bool do_print=false);
    static void stoplap(std::string msg="");
    static void off() {all_off=true;}
    static void on() {all_off=false;}
    static void starts_on() {print_starts=true;}
    
private:
    static double starttime;
    static std::vector<double> startlaptime;
    static std::vector<int> levels;
    static std::vector<std::string> msgs;
    static void do_write(std::ostream &s, double time, double delta, std::string traceback, std::string message);
    
    static int printlevel;
    static int writelevel;
    static bool all_off;
    static bool print_starts;
    static std::ofstream of;
    static void write(double time, double delta, std::string message, int wlevel);
    static void print(double time, double delta, std::string message, int plevel);
    
};


#endif
````

