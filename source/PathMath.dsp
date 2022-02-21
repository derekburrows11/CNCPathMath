# Microsoft Developer Studio Project File - Name="PathMath" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=PathMath - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "PathMath.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "PathMath.mak" CFG="PathMath - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "PathMath - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "PathMath - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "PathMath - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD BASE RSC /l 0xc09
# ADD RSC /l 0xc09
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "PathMath - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /MDd /W4 /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /FR /YX"stdafx.h" /FD /c
# ADD BASE RSC /l 0xc09
# ADD RSC /l 0xc09
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "PathMath - Win32 Release"
# Name "PathMath - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\CNCControl\PathSpeed.cpp
# End Source File
# Begin Source File

SOURCE=.\PointSumInfo.cpp
# End Source File
# Begin Source File

SOURCE=..\Common\PolyFunc.cpp
# End Source File
# Begin Source File

SOURCE=.\PolySegDblFit.cpp
# End Source File
# Begin Source File

SOURCE=.\PolySegFit.cpp
# End Source File
# Begin Source File

SOURCE=.\PolySegFitPoints.cpp
# End Source File
# Begin Source File

SOURCE=.\PolySegInfo.cpp
# End Source File
# Begin Source File

SOURCE=.\StdAfx.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h"
# Begin Source File

SOURCE=.\PointSumInfo.h
# End Source File
# Begin Source File

SOURCE=..\Common\PolyFunc.h
# End Source File
# Begin Source File

SOURCE=.\PolySegDblFit.h
# End Source File
# Begin Source File

SOURCE=.\PolySegFit.h
# End Source File
# Begin Source File

SOURCE=.\PolySegFitPoints.h
# End Source File
# Begin Source File

SOURCE=.\PolySegInfo.h
# End Source File
# Begin Source File

SOURCE=.\StdAfx.h
# End Source File
# End Group
# Begin Source File

SOURCE=".\Formula Path Math.txt"
# End Source File
# End Target
# End Project
