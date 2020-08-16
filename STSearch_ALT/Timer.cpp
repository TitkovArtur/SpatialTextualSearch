//
//  Timer.cpp
//  STSearch_ALT
//
//  Created by Artur Titkov on 11.06.20.
//  Copyright © 2020 Artur Titkov. All rights reserved.
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
