#ifndef F3CMAX_H
#define F3CMAX_H
#include "stdafx.h"
#include "permutation.h"
#include "task.h"
#include "branch_and_bound.h"
#define N 15 //Max number of tasks
using namespace std;
typedef boost::multiprecision::cpp_int bigint;

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

class F3Cmax{
public:
    vector<Task> tasks;

    map<string, bigint> nsDB;

    bigint tmpCountDicUsage;

    F3Cmax( const vector<Task> & tasks){
        tmpCountDicUsage=0;
        this->tasks = tasks;
    }

    /// Return an upper bound of the permutation
    int upperBound( const vector<int> & indTasks){
        assert( !indTasks.empty());
        if(indTasks.size()==1)
            return tasks[indTasks.front()].p1 + tasks[indTasks.front()].p2 + tasks[indTasks.front()].p3;
        // Calculate upper bound! To search
        // Currently we choose a random permu
        vector<int> indPermu = indTasks;
        std::random_shuffle(indPermu.begin(), indPermu.end());
        Permutation perm(indPermu);
        int cmax1 = perm.calculateCmax(tasks);
        std::random_shuffle(indPermu.begin(), indPermu.end());
        perm.setPermutation(indPermu);
        int cmax2 = perm.calculateCmax(tasks);
        return std::min(cmax1, cmax2);
    }

    /// Create a unique string key for the ns function call
    string makeKey(const bitset<N> & S, int t1, int t2, int t3){
        stringstream key;
        key<<S.to_string()<<","<<t1<<","<<t2<<","<<t3;
        //cout<<key.str();
        return key.str();
    }

    /// The function to calculate N(S): number of affectation of intervals excluding jobs from S, with cmax < t3, c2<t2, c1<t1
    bigint ns(const bitset<N> & S, int t1, int t2, int t3){
        // Ensure that ti<t(i+1)
        t2=min(t2, t3-1);
        t1=min(t1, t2-1);

        // Halting conditions
        if(t1<0) return 0;
        else if (t1==0) return 1;
        // Excluding all jobs: one way no schedule
        if(S.all())
            return 1;

        // Try to find if it's already precomputed
        string key=makeKey(S,t1,t2,t3);
        map<string,bigint>::iterator itVal = nsDB.find(key);
        if(itVal != nsDB.end()){
            tmpCountDicUsage++;
            //cout<<"yes!";
            return itVal->second;
        }
        // Recurrence conditions
        bigint result = ns(S, t1,t2, t3-1);
        for (unsigned int i=0; i<tasks.size();  i++)
        {
            if(S.test(i) == false){           // job i is in S_bar
                int arg3=t3-tasks[i].p3;
                int arg2=min(t2-tasks[i].p2, arg3-tasks[i].p2);
                int arg1=min(t1-tasks[i].p1, arg2-tasks[i].p1);
                result += ns(S, arg1, arg2, arg3);
            }
        }

        nsDB[key] = result;
        return result;
    }

    /// 
    bigint calc_InclusionExclusion( const int k){
        // (-1)^|W|N(W)
        bigint count=0;
        int n = tasks.size();   // Problem size
        bigint nW = int(std::pow(double(2), int(n)));    // Number of subsets
        // Consider n<64
        for(unsigned __int64 i=0; i<nW; i++){
            // We flip the selection by "nW - 1 -i"
            // So selection is all properties that we could have, not essentially all of them.
            std::bitset<N> selection(i);
            bigint partialcount=ns(selection, k,k,k);
            cout<<"========"<< i<<endl;
            cout<<selection<<endl;
            cout <<partialcount<<endl;
            count += selection.count()%2==0? partialcount: -partialcount;
        }

        return count;
    }


    /// Solve instance by force
    pair<int,string> solve_enu(){
        int cmax=numeric_limits<int>::max();
        stringstream ssPerm;

        vector<int> vecPerm;
        for(unsigned int i=0; i<tasks.size(); i++)
            vecPerm.push_back(i);
        Permutation perm;
        do{
            perm.setPermutation(vecPerm);
            //cout << "\nPerm = ";
            //std::copy(vecPerm.begin(), vecPerm.end(), std::ostream_iterator<int>(std::cout, ","));
            //cout<<"Cmax Matrix:===="<<endl;
            perm.calculateCmax(tasks);
            //cout << "C2 = " <<perm.cmaxMatrix[1][2];
            //cout << "; Cmax = " <<perm.cmax<< endl;
            if(perm.cmax < cmax){
                cmax = perm.cmax;
                ssPerm.str("");
                ssPerm<<"(";
                for(auto it = perm.permutation.cbegin(); it!=perm.permutation.cend(); it++)
                    ssPerm<<*it<<",";
            }

        }while(std::next_permutation(vecPerm.begin(), vecPerm.end()));

        return make_pair(cmax, ssPerm.str());
    }

    /// Solve by calling branch and bound
    unsigned int solve_bb(){
        long* arrtasks = Task::tasksToArray(tasks);
        unsigned int r= branchAndBound(3,tasks.size(),arrtasks);
        return r;
    }



    /// (Old version)Inclusion-Exclusion return the number of permutation of indTasks with cmax < k
    int calc_InclusionExclusion( const vector<int> & indTasks, const int k){
        assert(!indTasks.empty());
        cout<<">>>> Step into I-E, k = "<< k <<endl;
        if(indTasks.size()==1)    //Halting condition
            return tasks[indTasks.front()].p1 + tasks[indTasks.front()].p2 + tasks[indTasks.front()].p3 > k ? 0: 1;
        int counter=0;
        // (-1)^|W|N(W)
        int n = indTasks.size();   // Problem size
        int nW = int(std::pow(double(2), int(n)));    // Number of subsets
        for(int i=1; i<nW-1; i++){
            // We flip the selection by "nW - 1 -i"
            // So selection is all properties that we could have, not essentially all of them.
            std::bitset<N> * selection = new std::bitset<N>(nW-1 - i);
            vector<int> indTasksSubset;
            for(int j=0; j<n; j++){
                if(selection->test(j))
                    indTasksSubset.push_back(indTasks[j]);
            }
            //indTasksSubset is the tasks we COULD have in the permutation
            int ub = upperBound(indTasksSubset);
            cout<< "=====\nindTasksSubset=";
            std::copy(indTasksSubset.begin(), indTasksSubset.end(), std::ostream_iterator<int>(cout, " "));
            cout<<endl <<"Upper bound = "<<ub<<endl;
            int addition = 1;
            if( (n-indTasksSubset.size())%2 == 1) addition = -1;
            if(ub<=k) {
                // Add all possible permutations that formed by subsets of indTasksSubset
                // (Could be precalculated)
                addition = addition * int(factorial(indTasksSubset.size() * E -1 ));
            }
            else if(indTasksSubset.size() <= 1){
                addition = 0;
                cerr<<"k is too small even for only one task. Try a bigger cmax."<<endl;
                exit(1);
            }
            // If ub > k,,,,,,then we fail......................................no way to continue

            addition = addition * calc_InclusionExclusion(indTasksSubset, k);
            counter+=addition;
            cout<<"factor = "<<addition<<", counter = "<<counter<<endl;
            delete selection;
        }
        cout<<"<<<< Step out I-E"<<endl;
        return counter + factorial(n);
    }


}; //Class

#endif
