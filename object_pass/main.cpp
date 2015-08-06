/// object passer
// a lesson in passing objects to processors via mpi

#include <iostream>
#include <stdio.h>
#include <mpi.h>
#include <vector>

using namespace std;

struct A
{
  int a;
  int b;
  int c;
};

struct Point2
{
  int x,y;
};

class Point3
{
public:
  struct Point2 xy;
  int z;
};


int main(int argc, char ** argv)
{
  MPI_Init(NULL,NULL);
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  //    // pass via send receive
  //    struct A obj;
  //    if(rank == 0)
  //    {
  //
  //      obj.a = 1;
  //      obj.b = 2;
  //      obj.c = 3;
  //      for(int i = 0; i < nproc; i++)
  //        if(i != rank)
  //            MPI_Send(&obj,sizeof(struct A),MPI_BYTE,i,0,MPI_COMM_WORLD);
  //
  //    }else{
  //        MPI_Recv(&obj,sizeof(struct A),MPI_BYTE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  //    }
  //
  //    MPI_Barrier(MPI_COMM_WORLD);
  //    printf("Proc %d: a = %d, b = %d, c = %d\n",rank,obj.a, obj.b, obj.c);
  
  Point3 p;
  vector<Point3> pts;
  if(rank == 0)
  {
    for(int i = 0; i < nproc; i++)
    {
      Point3 tmp{{i,2*i},3*i};
      pts.push_back(tmp);
    }
    
    for (Point3 p : pts) printf("%d, %d, %d\n", p.xy.x,p.xy.y,p.z);
  }
//  
  MPI_Scatter(pts.data(),sizeof(Point3),MPI_BYTE,
              &p,sizeof(Point3), MPI_BYTE,
              0,MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  printf("PROC %d: {%d, %d, %d}\n", rank,p.xy.x,p.xy.y,p.z);
//  delete obj2;
  MPI_Finalize();
}
