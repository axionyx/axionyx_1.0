
# Class Stopwatch


[**Class List**](annotated.md) **>** [**Stopwatch**](classStopwatch.md)



_def for debugging missing stoplap calls._ [More...](#detailed-description)

* `#include <Stopwatch.H>`
















## Public Static Functions

| Type | Name |
| ---: | :--- |
|  void | [**init**](classStopwatch.md#function-init) (int plevel=0, int wlevel=0, std::string logname="stopwatch\_log.dat") <br>_initialize object, store start time and verbosity level. plevel is the level of a lap to be mentioned when printing on the screen. wlevel is the same, but for writing the logfile. logname is the name of file in which to store the output._  |
|  void | [**off**](classStopwatch.md#function-off) () <br>_turn off all messages_  |
|  void | [**on**](classStopwatch.md#function-on) () <br>_turn on messages_  |
|  void | [**start**](classStopwatch.md#function-start) (int level=0) <br>_start the clock_  |
|  void | [**startlap**](classStopwatch.md#function-startlap) (std::string msg="", int level=0, bool do\_print=false) <br>_start a new lap with name msg; level sets the level of importance of the message._  |
|  void | [**starts\_on**](classStopwatch.md#function-starts-on) () <br>_print information if any lap is started_  |
|  void | [**stop**](classStopwatch.md#function-stop) (std::string msg, int level=0) <br>_stop the clock_  |
|  void | [**stoplap**](classStopwatch.md#function-stoplap) (std::string msg="") <br>_stop a running lap - msg is not used unless you define DEBUG\_STOPWATCH\_STOPLAP above for debugging purposes._  |







## Public Static Functions Documentation


### <a href="#function-init" id="function-init">function init </a>


```cpp
static void Stopwatch::init (
    int plevel=0,
    int wlevel=0,
    std::string logname="stopwatch_log.dat"
) 
```



### <a href="#function-off" id="function-off">function off </a>


```cpp
static inline void Stopwatch::off () 
```



### <a href="#function-on" id="function-on">function on </a>


```cpp
static inline void Stopwatch::on () 
```



### <a href="#function-start" id="function-start">function start </a>


```cpp
static void Stopwatch::start (
    int level=0
) 
```



### <a href="#function-startlap" id="function-startlap">function startlap </a>


```cpp
static void Stopwatch::startlap (
    std::string msg="",
    int level=0,
    bool do_print=false
) 
```



### <a href="#function-starts-on" id="function-starts-on">function starts\_on </a>


```cpp
static inline void Stopwatch::starts_on () 
```



### <a href="#function-stop" id="function-stop">function stop </a>


```cpp
static void Stopwatch::stop (
    std::string msg,
    int level=0
) 
```



### <a href="#function-stoplap" id="function-stoplap">function stoplap </a>


```cpp
static void Stopwatch::stoplap (
    std::string msg=""
) 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/Stopwatch.H`