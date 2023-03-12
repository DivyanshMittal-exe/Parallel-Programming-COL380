#include <iostream>

using namespace std;

int main()
{
    int i=0;
    int j;
    do
    {
        i++;
        j= 2*i;
    }
    while(j<10);
    cout<<i<<endl;
}