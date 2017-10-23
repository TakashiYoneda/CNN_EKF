# CNN_EKF


here attached src files for EKF project.  
In src,
- main.cpp
- FusionEKF.h
- FusionEKF.cpp
- kalman_filter.h
- kalman_filter.cpp
- tools.h
- tools.cpp
- measurement_package.h
- json.hpp
   
  
1. Compile
I succeeded to compile those files through process Udacity provided.  

2. Accuracy
By setting P with specific set of numbers, the accuracy was less than [px:0.11, py:0.11, vx:0.52, vy:0.52] for both Dataset1 and Dataset2.  

3. Follows the Correct Algorithm  
I followed instructions in the lessons for KF/EKF algorithm.  
I set processing of initialization and predict&update, and both can handle laser/radar data.  
  
4. Code Efficiency
I tried to avoid redundancy before accomplishing coding.  

