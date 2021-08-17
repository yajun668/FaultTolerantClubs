# Code for Finding Fault-tolerant s-Clubs


This code accompanies the paper "On fault-tolerant low-diameter clusters in graphs" and is written in C++. If you wish to use or cite this code, please cite the paper: 

    @article{LSBB2021fault-tolerant, 
      title={On fault-tolerant low-diameter clusters in graphs}, 
      author={Lu, Yajun and Salemi, Hosseinali and Balasundaram, Balabhaskar and Buchanan, Austin}, 
      year={2021}}

This repository includes two folders:
1. r_robust_s_club: used to find r-robust s-clubs
2. h_hereditary_s_club: used to find h-hereditary s-clubs

## Compiling the code
The following steps show how to compile and run the code to find r-robust s-clubs in a Linux environment using a makefile (you can also run the code in Mac or Windows environment by configuring your IDE appropriately). In the folder of r_robust_s_club, parameter.txt is used to configure parameters r and s, and the "data" folder includes InputFile.txt (used to determine which instances you want to test) and 10th DIMACS graph instances (you can also downlaod these instances from the website: https://www.cc.gatech.edu/dimacs10/archive/clustering.shtml). It is similar to run the code to find h-hereditary s-clubs.


### Steps to run the code to find r-robust s-clubs in Linux environment:
1. Download or clone the repository to your machine.
2. Go to the folder r_robust_s_club.
3. Open the "Makefile" and set GUROBI_HOME to the directory of your Gurobi installation, e.g.: /opt/gurobi/9.0.1/linux64.
4. From the terminal, go to the folder of r_robust_s_club.
5. Type "make" and hit enter to compile. After successful complilation, type "./main" to run the code.


## Terms and Use:

MIT License

Copyright (c) 2021 Yajun Lu, Hosseinali Salemi, Balabhaskar Balasundaram, and Austin Buchanan.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
