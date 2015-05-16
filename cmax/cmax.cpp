// cmax.cpp? définit le point d'entrée pour l'application console.
//

#include "stdafx.h"
#include "f3cmax.h"
#include "task.h"
#include "permutation.h"
#include "branch_and_bound.h"
using namespace std;

void genRandomTask(vector<Task>& tasks,int min, int max, int nbTask);
void mainIE();

int _tmain(int argc, _TCHAR* argv[])
{
    
    mainIE();
    return 1;
}

void mainIE(){

//while(true)
{
     
    vector<Task> tasks;
    genRandomTask(tasks, 1,5, N);

    F3Cmax f(tasks);

    auto resPair = f.solve_bb();    cout<<resPair<<endl;

    int count = f.calc_InclusionExclusion(resPair-1);    cout<<endl<<count<<"!"<<endl;    cout<<"Db usage:"<<f.tmpCountDicUsage<<"!"<<endl;
    
    /*char c;
    if(count<=0) 
        cin.get();

    count = f.calc_InclusionExclusion(resPair.first-1);
    if(count!=0) 
        cin.get();*/
}
}


void genRandomTask(vector<Task>& tasks,int min, int max, int nbTask){
    default_random_engine generator(time(NULL));
    uniform_int_distribution<int> dist(min,max);
    auto getRandomP = [&dist, &generator](){ return dist(generator);};
    
    tasks.resize(nbTask);
    for(int i=0; i<nbTask; i++){
        tasks[i].p1=getRandomP();
        tasks[i].p2=getRandomP();
        tasks[i].p3=getRandomP();
        tasks[i].toString();
    }

    /*for_each(tasks.cbegin(), tasks.cend(), [](Task& t){
        t.toString();
    });
    Task t1(43,  18,  25);
    Task t2(40,  60,  68);
    Task t3(9,  37,  50);
    Task t4(81,  95,  96);
    Task t5(26,  70,  85);*/
}


void genRandomTask()
{
    default_random_engine generator;
    uniform_int_distribution<int> dist(1,5);
    auto getRandomP = [&dist, &generator](){ return dist(generator);};
    vector<int> Pij(9);
    generate(Pij.begin(), Pij.end(), getRandomP);

    Task &t1 = *(new Task(Pij[0],  Pij[1],  Pij[2]));
    Task &t2 = *(new Task(Pij[3],  Pij[4],  Pij[5]));
    Task &t3 = *(new Task(Pij[6],  Pij[7],  Pij[8]));
    t1.toString();
    t2.toString();
    t3.toString();


    vector<Task> tasks;
    tasks.push_back(t1);
    tasks.push_back(t2);
    tasks.push_back(t3);
    //tasks.push_back(t4);

    // Initial order
    vector<int> vecPerm;
    vecPerm.push_back(0);
    vecPerm.push_back(1);
    vecPerm.push_back(2);
    //vecPerm.push_back(3);

    Permutation perm;
    do{
        perm.setPermutation(vecPerm);
        cout << "\nPerm = ";
        std::copy(vecPerm.begin(), vecPerm.end(), std::ostream_iterator<int>(std::cout, ","));
        cout<<"Cmax Matrix:===="<<endl;
        perm.calculateCmax(tasks);
        cout << "C2 = " <<perm.cmaxMatrix[1][2];
        cout << "; Cmax = " <<perm.cmax<< endl;
    }while(std::next_permutation(vecPerm.begin(), vecPerm.end()));

    //_______________________________

    //F3Cmax f(tasks);
    //vector<int> indTasks(tasks.size());
    //for(int i=0; i<tasks.size(); i++)
    //    indTasks[i]=i;

    //int count = f.calc_InclusionExclusion(indTasks, 12);
    //cout<<endl<<count<<"!"<<endl;

    delete &t1;
    delete &t2;
    delete &t3;
}