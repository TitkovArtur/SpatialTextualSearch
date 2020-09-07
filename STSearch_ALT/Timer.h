//
// Project:    Spatially Combined Text Searches
// Filename:   Timer.cpp
// Created on: 11.06.20.
//
// Developer(s):
//    I) Artur Titkov
//       ArturTitkov@icloud.com
//       https://github.com/TitkovArtur
//    II) unknown
//
// Description:
//
//
// Copyright © 2020 Artur Titkov. All rights reserved.
//

#ifndef Timer_h
#define Timer_h

#include <stdio.h>
#include <chrono>

class Timer{
private:
    using Clock = std::chrono::high_resolution_clock;
    Clock::time_point start_time, stop_time;
    
public:
    Timer();
    void start();
    double getElapsedTimeInSeconds();
    double stop();
};







#endif /* Timer_hpp */
