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

#include "Timer.h"
#include <chrono>

using namespace std;


Timer::Timer(){
    start();
}

void Timer::start(){
    start_time = Clock::now();
}


double Timer::getElapsedTimeInSeconds(){
    return std::chrono::duration<double>(stop_time - start_time).count();
}


double Timer::stop(){
    stop_time = Clock::now();
    double tmp = getElapsedTimeInSeconds();
    start();
    return tmp;
}
