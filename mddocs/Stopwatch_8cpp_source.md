
# File Stopwatch.cpp

[**File List**](files.md) **>** [**Source**](dir_74389ed8173ad57b461b9d623a1f3867.md) **>** [**Stopwatch.cpp**](Stopwatch_8cpp.md)

[Go to the documentation of this file.](Stopwatch_8cpp.md) 


````cpp
#include "Stopwatch.H"
#include <ios>
#include <iostream>
#include <iomanip>
using namespace amrex;

double Stopwatch::starttime=0.0;
std::vector<double> Stopwatch::startlaptime;
std::vector<int> Stopwatch::levels;
std::vector<std::string> Stopwatch::msgs;
bool Stopwatch::all_off = false;
bool Stopwatch::print_starts = false;

int Stopwatch::printlevel=0;
int Stopwatch::writelevel=0;

std::ofstream Stopwatch::of;

#ifdef USE_STOPWATCH

void Stopwatch::init(int plevel,int wlevel,std::string logname)
{
    printlevel = plevel;
    writelevel = wlevel;
    
    starttime = 0;
    of.open(logname);
    
}
void Stopwatch::start(int level)
{
    starttime = ParallelDescriptor::second();
    write(starttime,-1,"starting clock",level);
    print(starttime,-1,"starting clock",level);
}
void Stopwatch::stop(std::string msg,int level)
{
    double stoptime = ParallelDescriptor::second();
    write(stoptime,stoptime-starttime,"stopping clock;",level);
    print(stoptime,stoptime-starttime,"stopping clock;",level);
}
void Stopwatch::startlap(std::string msg, int level, bool do_print)
{
    double s = ParallelDescriptor::second();
    startlaptime.push_back(s);
    levels.push_back(level);
    msgs.push_back(msg);
    if(do_print || print_starts)
    {
        write(startlaptime.back(),-1,"starting lap",level);
        print(startlaptime.back(),-1,"starting lap",level);
    }
}
void Stopwatch::stoplap(std::string msg)
{
    double stoplaptime = ParallelDescriptor::second();
    
    write(stoplaptime,stoplaptime-startlaptime.back(),"finished in",levels.back());
    print(stoplaptime,stoplaptime-startlaptime.back(),"finished in",levels.back());
#ifdef DEBUG_STOPWATCH_STOPLAP
    if(msg!="" && msg!=msgs.back() && ParallelDescriptor::IOProcessor())
        std::cout << "Stopwatch::stoplap(): WARNING; should be stopping " << msg << " but found " << msgs.back() << std::endl;
#endif        
    startlaptime.pop_back();
    levels.pop_back();
    msgs.pop_back();
}

void Stopwatch::do_write(std::ostream &s, double time, double delta, std::string traceback, std::string message)
{
    std::ios_base::fmtflags f( std::cout.flags() );
    if(delta<0.0)
        s << "[" << std::setw(5) << ((int)(time)) << "s]: " << traceback << message << std::setprecision(12) << std::endl;
    else
        s << "[" << std::setw(5) << ((int)(time)) << "s]: " << traceback << message << " " << std::setprecision(5)  << delta << "s" << std::setprecision(12) << std::endl; 
    std::cout.flags( f );
}


void Stopwatch::write(double time, double delta,std::string message, int wlevel)
{
    if(writelevel<wlevel)
        return;
    
    std::string traceback;
    for(auto it=msgs.begin();it!=msgs.end();++it)
    {
        traceback += *it+" >> ";
        
    }
    traceback += "; ";
    if(ParallelDescriptor::IOProcessor())
    {
        do_write(of, time-starttime, delta, traceback, message);
    }
}

void Stopwatch::print(double time, double delta, std::string message, int plevel)
{
    if(printlevel<plevel)
        return;
    if(all_off)
        return;
    std::string traceback;
    for(unsigned int i=0;i<msgs.size();i++)
    {
        traceback += msgs[i];
        if(i!=msgs.size()-1)
            traceback += " >> ";
        
    }
    traceback += "; ";
    if(ParallelDescriptor::IOProcessor())
    {
        do_write(std::cout, time-starttime, delta, traceback, message);
    }
}
#else
void Stopwatch::init(int plevel,int wlevel,std::string logname)
{
    printlevel = plevel;
    writelevel = wlevel;
    
    starttime = 0;
}
void Stopwatch::start(int level)
{
}
void Stopwatch::stop(std::string msg,int level)
{
}
void Stopwatch::startlap(std::string msg,bool do_print, int level)
{
}
void Stopwatch::stoplap()
{
}
void Stopwatch::write(double time, double delta,std::string message, int wlevel)
{
}
void Stopwatch::print(double time, double delta, std::string message, int plevel)
{
}
#endif

````

