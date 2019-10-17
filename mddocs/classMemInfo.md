
# Class MemInfo


[**Class List**](annotated.md) **>** [**MemInfo**](classMemInfo.md)




















## Public Functions

| Type | Name |
| ---: | :--- |
|  int | [**BelowThreshold**](classMemInfo.md#function-belowthreshold) (float threshold) <br> |
|  bool | [**GetMemInfo**](classMemInfo.md#function-getmeminfo) (float \* availMem, float \* totalMem) <br> |
|  void | [**Init**](classMemInfo.md#function-init) (const std::string & fname) <br> |
|  void | [**LogSummary**](classMemInfo.md#function-logsummary) (const char \* info) <br> |
|  void | [**PrintAll**](classMemInfo.md#function-printall) (FILE \* fname) <br> |

## Public Static Functions

| Type | Name |
| ---: | :--- |
|  [**MemInfo**](classMemInfo.md) \* | [**GetInstance**](classMemInfo.md#function-getinstance) () <br> |







## Public Functions Documentation


### <a href="#function-belowthreshold" id="function-belowthreshold">function BelowThreshold </a>


```cpp
int MemInfo::BelowThreshold (
    float threshold
) 
```



### <a href="#function-getmeminfo" id="function-getmeminfo">function GetMemInfo </a>


```cpp
bool MemInfo::GetMemInfo (
    float * availMem,
    float * totalMem
) 
```



### <a href="#function-init" id="function-init">function Init </a>


```cpp
void MemInfo::Init (
    const std::string & fname
) 
```



### <a href="#function-logsummary" id="function-logsummary">function LogSummary </a>


```cpp
void MemInfo::LogSummary (
    const char * info
) 
```



### <a href="#function-printall" id="function-printall">function PrintAll </a>


```cpp
void MemInfo::PrintAll (
    FILE * fname
) 
```


## Public Static Functions Documentation


### <a href="#function-getinstance" id="function-getinstance">function GetInstance </a>


```cpp
static inline MemInfo * MemInfo::GetInstance () 
```



------------------------------
The documentation for this class was generated from the following file `/home/uni06/cosmo/cbehren2/AXIONYX/axionyx/Source/Monitors/MemInfo.H`