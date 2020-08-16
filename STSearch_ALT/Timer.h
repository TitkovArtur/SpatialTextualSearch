//
//  Timer.hpp
//  STSearch_ALT
//
//  Created by Artur Titkov on 11.06.20.
//  Copyright © 2020 Artur Titkov. All rights reserved.
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
