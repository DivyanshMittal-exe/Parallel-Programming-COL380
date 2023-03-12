#include<bits/stdc++.h>
// #include<mpi.h>

using namespace std;

int main(int argc, char** argv)
{
    fstream outdata("test1.dat", ios::out | ios::binary);
    fstream outdata_gra("test1.gra", ios::out | ios::binary);
    int n = 7;
    int m = 12;

    // outdata_gra.write(reinterpret_cast<char *>(&n), 4);
    // outdata_gra.write(reinterpret_cast<char *>(&m), 4);

    int arr[] = {7,12,    0,5,2,3,4,1,5    ,1,5,6,4,3,2,0     ,2,3,0,1,3      ,3,4,2,0,1,4     ,4,3,3,0,1        ,5,2,0,6,  6,2,1,5};
    int arr_head[] = {8,36,64,84,108,128,144};
    for(int i=0; i<40; i++)
    {
        outdata_gra.write(reinterpret_cast<char *>(&(arr[i])), 4);
    }
    for(int i=0; i<7; i++)
    {
        outdata.write(reinterpret_cast<char *>(&(arr_head[i])), 4);
    }

    return(0);
}