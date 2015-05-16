#include "f3cmax.h"
#include "task.h"
#include "permutation.h"

using namespace std;

int main()
{
    Task t1(1,2,3);
    Task t2(1,2,2);
    Task t3(4,2,1);

    vector<Task> tasks;
    tasks.push_back(t1);
    tasks.push_back(t2);
    tasks.push_back(t3);
//
//    vector<int> vecPerm({0,1,2});
//    Permutation perm;
//    do{
//        perm.setPermutation(vecPerm);
//        perm.calculateCmax(tasks);
//        cout << "Perm = ";
//        std::copy(vecPerm.begin(), vecPerm.end(), std::ostream_iterator<int>(std::cout, ","));
//        cout << " Cmax = " <<perm.cmax<< endl;
//    }while(std::next_permutation(vecPerm.begin(), vecPerm.end()));
///_______________________________

    F3Cmax f(tasks);
    vector<int> indTasks(tasks.size());
    for(int i=0; i<tasks.size(); i++)
        indTasks[i]=i;

    int count = f.calc_InclusionExclusion(indTasks, 12);
    cout<<endl<<count<<"!"<<endl;
    return 1;
}
