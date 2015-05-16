/// The definition of the class Task

#ifndef TASK_H
#define TASK_H
#include "stdafx.h"
#include "f3cmax.h"
    
using namespace std;
class Task{
public:
    int p1;
    int p2;
    int p3;

public:
    Task(int p1, int p2, int p3){
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
    }

    Task(){
        this->p1 = this->p2 = this->p3 = 0;
    }

    void toString(){
        cout<<"\n Task pij: "<<p1<<" "<<p2<<" "<<p3<<" ==== "<<endl;
    }

    // Transforms Task data to a long array in order to be used by "branch and bound"
    static long* tasksToArray( vector<Task> & tasks){
        long * arrTask = new long[tasks.size()*3];
        int i=0;
        for_each(tasks.begin(), tasks.end(), [&arrTask,&i](Task & t){
            arrTask[i++]=t.p1;
            arrTask[i++]=t.p2;
            arrTask[i++]=t.p3;
        });
        return arrTask;
    }
};

#endif
