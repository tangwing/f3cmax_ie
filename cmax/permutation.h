#ifndef PERMUTATION_H
#define PERMUTATION_H
#include "stdafx.h"
#include "f3cmax.h"
#include "task.h"
using namespace::std;
class Permutation{
public:
    vector<int> permutation;
    vector< vector<int> > cmaxMatrix;
    int cmax;

    Permutation(){cmax = -1;}
    Permutation( const vector<int>& perm){
        permutation = perm;
    }

    void setPermutation(const vector<int> & perm){
        permutation = perm;
        cmax = -1;
        cmaxMatrix.clear();
    }

    int calculateCmax(const vector<Task> & tasks){
        assert(tasks.size() > *std::max_element(permutation.begin(), permutation.end()))  ;
        cmaxMatrix.clear();
        for(auto it = permutation.begin(); it!= permutation.end(); it++)
        {
            vector<int> cmaxi(3);
            if(it == permutation.begin()){
                cmaxi[0] = tasks[*it].p1;
                cmaxi[1] = cmaxi[0] + tasks[*it].p2;
                cmaxi[2] = cmaxi[1]+tasks[*it].p3;
            }else{
                vector<int> & last = cmaxMatrix.back();
                cmaxi[0] = tasks[*it].p1 + last[0];
                cmaxi[1] = std::max(cmaxi[0] , last[1])+ tasks[*it].p2;
                cmaxi[2] = std::max(cmaxi[1] , last[2])+ tasks[*it].p3;
            }
            //cout<< cmaxi[0] << " " <<cmaxi[1] <<" " <<cmaxi[2] <<endl;
            cmaxMatrix.push_back(cmaxi);
        }
        cmax = cmaxMatrix.back()[2];
        //print cmax matrix
        //for (int i=0; i<3; i++)
        //{
        //    for (auto it=cmaxMatrix.cbegin(); it!=cmaxMatrix.cend(); it++)
        //    {
        //        cout<< (*it)[i] << " ";
        //    }
        //    cout<<endl;
        //}
        return cmax;
    }
};
#endif
