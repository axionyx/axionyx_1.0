#ifndef _STOPWATCH__
#define _STOPWATCH__

//! undef to disable logging.
#define USE_STOPWATCH

#include "AMReX_ParallelDescriptor.H"
#include <fstream>
//! this is a little helper class that implements some profiling capabilities. The main useful routines are start/stoplap. They can be used to track the amount of time spent in a subroutine/code block, and the output also comes with a traceback of the laps already running.

//! For example, let us say you have a task A that consists of substeps B and C. If you put startlaps at the beginning of each (sub)task and stoplaps at the end of them, it will give a output the timing of A,B, and C, with a backtrace of the relation between the two (i.e. B and C being substeps of A)
class Stopwatch{
public:
    //! initialize object, store start time and verbosity level. plevel is the level of a lap to be mentioned when printing on the screen. wlevel is the same, but for writing the logfile. logname is the name of file in which to store the output.
    static void init(int plevel=0,int wlevel=0,std::string logname="log.dat");
    //! start the clock
    static void start(int level=0);
    //! stop the clock
    static void stop(std::string msg,int level=0);
    //! start a new lap with name msg; level sets the level of importance of the message.
    static void startlap(std::string msg="", int level=0, bool do_print=false);
    //! stop a running lap
    static void stoplap();
    //! turn off all messages
    static void off() {all_off=true;}
    //! turn on messages
    static void on() {all_off=false;}
    //! print information if any lap is started 
    static void starts_on() {print_starts=true;}
    
private:
    //! actual start time of the clock
    static double starttime;
    //! holds the start times of the active laps (works like a stack)
    static std::vector<double> startlaptime;
    //! holds the levels of the active laps (works like a stack)
    static std::vector<int> levels;
    //! holds the messages belonging to the active laps (works like a stack)
    static std::vector<std::string> msgs;
    //! do the actual write of a message into a stream buffer
    static void do_write(std::ostream &s, double time, double delta, std::string traceback, std::string message);
    
    //! the level above which we will print things on the screen
    static int printlevel;
    //! the level above which we will print things into a file
    static int writelevel;
    //! if true, nothing will be printed
    static bool all_off;
    //! if true, all starts will be printed
    static bool print_starts;
    //! the output log file
    static std::ofstream of;
    //! short hand to write to logfile
    static void write(double time, double delta, std::string message, int wlevel);
    //! short hand to write to screen
    static void print(double time, double delta, std::string message, int plevel);
    
};


#endif
