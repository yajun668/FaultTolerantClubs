#include <iostream>
#include "model.h"
int main(int argc, const char * argv[]) {
    ifstream parameter_file("parameter.txt");// Parameters r and s are set up in this txt file
    string line;
    vector<string> parameters;// Instruction for parameters: 1st line: R; 2nd: S; 3rd: graph format-- 1 for DIMACS 10 and 2 for Edge list format; 4th: Instance folder name
    while (getline(parameter_file, line)){
        if (line[0]=='%')
            continue;
        parameters.push_back(line);
    }
    Model new_model;

    //change working directory
    chdir("..");; //go back to parent folder
    string wd = "./Dataset/"+parameters[3] + "/";
    chdir(wd.c_str());
    new_model.read_masterfile("InputFile.txt");
    cout<<"-------------Solving the maximum "<<parameters[0] + "-Hereditary " + parameters[1] + "-Club------------------"<<endl;
    string outputfile;
    outputfile=parameters[0] + "_Hereditary_" + parameters[1] + "_Club"+".csv"; // Result file
    ofstream fout;
    fout.open(outputfile.c_str(), ios::out);

    //FILE OUTPUT: columns headers
    fout<<"Graph Name,#vertex,#edge,h,s,PreproTime,BuildTime,HeuTime,grbSolveTime,WallTime,";
    fout<<"Heuristic, Best Obj, Best UB, MIPgap, #BC nodes, #Upfront Constr, #Lazy Cut, Status";
    fout<<"\n";
    fout.close();
    //Run the model for each instance
    for(long i=0;i<new_model.num_instances;i++)
    {
        fout.open(outputfile.c_str(), ios::app);
        Graph inputG;
        inputG.graphname = new_model.filenames[i]; //i-th instance name
        inputG.param_r = atol(parameters[0].c_str());
        inputG.param_s = atol(parameters[1].c_str());

        if (atol(parameters[2].c_str()) == 1) {
            cout << "Attempting to read DIMACS10 clustering format instance ...\n";
            inputG.ReadDIMACS10cluster();//read graph file
        }
        else if(atol(parameters[2].c_str()) == 2) {
            cout << "Attempting to read edge list format instance ...\n";
            inputG.ReadVBgraph();//read graph file
        }
        else{
            cout<<"Caution!!! Please enter 1 or 2 for the third line about graph format"<<endl;
            return 0;
        }

        fout<<inputG.graphname<<","<< inputG.nverts<<","<<inputG.nedges<<",";
        
        //Wal-clock time starts
        chrono_time_point Walltime = chrono_clock::now();
        //All results are saved in Sol; Initialization
        SolutionStats Sol;
        Sol.WallTime = 0;
        Sol.preProcessTime=0;
        Sol.BuildTime = 0;
        Sol.grbSolveTime = 0;
        Sol.SolHeuTime = 0;
        Sol.UpperBound = -1;
        Sol.nodes_exp=0;
        Sol.lazy_cut = vector<long>(2,0); //index 0 for #upfront constraints; index 1 for #lazy cuts
        Sol.currentBestSol.clear();
        Sol.BestObj = 0;
        Sol.heuristic_size=0;
        Sol.termistatus = 1; //1 optimal; 2 suboptimal;
        double MIPgap=0.0;

        inputG.KCore(inputG.param_r); //first preprocess graph using k-core
        //obtain wall-time before calling RecursiveBlockDecom
        chrono::duration<double>WallTime_span1 = chrono::duration_cast<std::chrono::duration<double> >(chrono_clock::now() - Walltime);
        Sol.WallTime += WallTime_span1.count();
        new_model.RecursiveBlockDecom(inputG,Sol); // run Recursive block decomposition
        Sol.BestObj = Sol.currentBestSol.size();

        //print Pre-processing time and Wall-clock time
        cout<<"The preprocessing time and  Wall-clock time is: "<<Sol.preProcessTime<<"->"<<Sol.WallTime<<endl;

        // print the best solution:
        cout<<"The best " <<inputG.param_r<<"_hereditary "<<inputG.param_s<<"_club objective and solution in graph "<<inputG.graphname<< " is: "<<Sol.BestObj<<endl;
        for (long i3=0; i3<Sol.BestObj; i3++)
            cout<<Sol.currentBestSol[i3]+1<<"->";
        cout<<endl;

        //save results in the csv file
        if( (Sol.UpperBound > -1) && (Sol.BestObj >= 1) )
            MIPgap = 100*(Sol.UpperBound-Sol.BestObj)/Sol.BestObj;
        //save results in the csv file
        fout<<inputG.param_r<<","<< inputG.param_s<<","<<Sol.preProcessTime<<","<<Sol.BuildTime<< ",";
        fout<<Sol.SolHeuTime<< ","<<Sol.grbSolveTime<< ","<<Sol.WallTime<< ",";
        fout<<Sol.heuristic_size<<","<< Sol.BestObj<<","<<Sol.UpperBound<<","<<MIPgap<<",";
        fout<<Sol.nodes_exp<<","<<Sol.lazy_cut[0]<<","<<Sol.lazy_cut[1]<<","<<Sol.termistatus;
        fout<<"\n";
        fout.close();
    }
    return 0;
}
