#include<math.h>
#include<stdio.h>
#include <sys/time.h>
#include <string.h>
#include <stdlib.h>


void stream(double s[101][101]);
void vorticity(double w[101][101]);
void velocity(double m[101][101]);
void prs(double pr[101][101]);

int main()
{
  int Nx=101,Ny=101,it=0,xx; 
  double dx,dy,dt = 0.00025,l=1,h=1,h1,h2,r,nu = 0.000001,bta,max1,max2,lt,pr[Ny][Nx],p[Ny][Nx],m[Ny][Nx],d1[Ny][Nx],y[Ny],d2[Ny][Nx],s[Ny][Nx],s1[Ny][Nx],u[Ny][Nx],u1[Ny][Nx],v[Ny][Nx],v1[Ny][Nx],w[Ny][Nx],w1[Ny][Nx];
  struct timeval stop, start;


  
/*  //printf("\nEnter number of grid points in x-direction:");
  //scanf("%d",&Nx);

  //printf("\nEnter number of grid points in y-direction:");
 // scanf("%d",&Ny);
  
  //printf("\nEnter the total length:");
  //scanf("%lf",&l);

  //printf("\nEnter the total height:");
  //scanf("%lf",&h);
 */ 
  
  //delx calculation
  dx=(l/(Nx-1)); 
  dy=(h/(Ny-1));
  bta=(dx/dy);
  h1=h;


  printf("\nThe Difference b/w two grid points in x-direction is: %lf" ,dx);
  printf("\nThe Difference b/w two grid points in y-direction is: %lf" ,dy);
  printf("\nThe Value of beta is: %lf" ,bta);
 
//initialization 

for(int i = 1; i <= Ny; i++)
  {
      for(int j = 1; j <= Nx; j++)
      {
        s[i][j]=0; 
        s1[i][j]=0;               
      } 
  }

for(int i = 1; i <= Ny; i++)
  {
      for(int j = 1; j <= Nx; j++)
      {
        w[i][j]=0;
        w1[i][j]=0;
      }
  }

  for(int i = 1; i <= Ny; i++)
  {
      for(int j = 1; j <= Nx; j++)
      {
        u[i][j]=0; 
        u1[i][j]=0;               
      } 
  }

  for(int i = 1; i <= Ny; i++)
  {
      for(int j = 1; j <= Nx; j++)
      {
        v[i][j]=0; 
        v1[i][j]=0;               
      } 
  }

//bc

  printf("\nEnter the value of Reynolds Number :\n"); //get inlet Velocity from this
  scanf("%lf",&r);

//left

  for(int j = 1; j <= Ny; j++)
    {
        u[j][1]=0;
        u1[j][1]=0;
    }

    for(int j = 1; j <= Ny; j++)
    {
        v[j][1]=0;
        v1[j][1]=0;
    }

  for(int j = 1; j <= Ny; j++)
  {
        s[j][1] = 0;
        s1[j][1] = 0;

  }

for(int j = 1; j <= Ny; j++)
  {
        w[j][1] = (s[j][1]-s[j][2])*(2/(dx*dx));
        w1[j][1] = (s1[j][1]-s1[j][2])*(2/(dx*dx));

  }

//right

for(int j = 1; j <= Ny; j++)
   {
       u[j][Nx]=0;
       u1[j][Nx]=0;
 }

   for(int j = 1; j <= Ny; j++)
   {
    v[j][Nx]=0;
    v1[j][Nx]=0;
   }

  for(int j = 1; j <= Ny; j++)
  {
       s[j][Nx]=0;
       s1[j][Nx]=0;
  } 
for(int j = 1; j <= Ny; j++)
  {
        w[j][Nx] = (s[j][Nx]-s[j][Nx-1])*(2/(dx*dx));
        w1[j][Nx] = (s1[j][Nx]-s1[j][Nx-1])*(2/(dx*dx));

  }


//top

  for(int j = 1; j <= Nx; j++)
    {
        u[1][j]=1;
        u1[1][j]=1;
    }

    for(int j = 1; j <= Nx; j++)
    {
        v[1][j]=0;
        v1[1][j]=0;
    }


  for(int i = 1; i <= Nx; i++)
  {
        s[1][i]= 0;
        s1[1][i]= 0;
  }
for(int i = 1; i <= Nx; i++)
  {
        w[1][i] = (s[1][i]-s[2][i])*(2/(dy*dy))+(2*1/dy);
        w1[1][i] = (s1[1][i]-s1[2][i])*(2/(dy*dy))+(2*1/dy);

  }

//bottom

  for(int j = 1; j <= Nx; j++)
    {
        u[Ny][j]=0;
        u1[Ny][j]=0;
    }

    for(int j = 1; j <= Nx; j++)
    {
        v[Ny][j]=0;
        v1[Ny][j]=0;
    }
  for(int i = 1; i <= Nx; i++)
  {
        s[Ny][i]= 0;
        s1[Ny][i]= 0;
  }

  for(int i = 1; i <= Nx; i++)
  {
        w[Ny][i] = (s[Ny][i]-s[Ny-1][i])*(2/(dy*dy));
        w1[Ny][i] = (s1[Ny][i]-s1[Ny-1][i])*(2/(dy*dy));
 }


for ( int i = 2 ; i < Ny ; i++ )
{
    for ( int j = 2 ; j < Nx ; j++ )

    {
         w[i][j] = 0;
         w1[i][j] = 0; //interior

    }
} 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  printf(" \nsolve\n\n"); 

  gettimeofday(&start, NULL);
    do
   {
     
    it=(it+1);

//stream fn  

      for(int i = 2; i < Ny; i++)
        {
            for(int j = 2; j < Nx; j++)
            {

                s1[i][j] = (s[i+1][j]+s[i-1][j]+s[i][j-1]+s[i][j+1]+(w[i][j]*(dx*dx)))*(0.25); //interior

            }
        }
  //swap

      for(int i = 1; i <= Ny; i++)
        {
            for(int j = 1; j <= Nx; j++)
            {

                s[i][j] = s1[i][j];

            }
        }   
 

//new velocities

      for(int j = 2; j < Nx; j++)
        {
            for(int i = 2; i < Ny; i++)
            {

                v1[i][j]= (-s[i][j+1]+s[i][j-1])*(2*dx);
                u1[i][j]=(s[i+1][j]-s[i-1][j])*(2*dy); 

            }
        }
        
      //swap

      for(int i = 1; i <= Ny; i++)
        {
            for(int j = 1; j <= Nx; j++)
            {

                u[i][j] = u1[i][j];

            }
        }  
//swap

      for(int i = 1; i <= Ny; i++)
        {
            for(int j = 1; j <= Nx; j++)
            {

                v[i][j] = v1[i][j];

            }
        } 
//BC 

//left

  for(int j = 1; j <= Ny; j++)
    {
        u[j][1]=0;
        u1[j][1]=0;
    }

    for(int j = 1; j <= Ny; j++)
    {
        v[j][1]=0;
        v1[j][1]=0;
    }

  for(int j = 1; j <= Ny; j++)
  {
        s[j][1] = 0;
        s1[j][1] = 0;

  }

for(int j = 1; j <= Ny; j++)
  {
        w[j][1] = (s[j][1]-s[j][2])*(2/(dx*dx));
        w1[j][1] = (s1[j][1]-s1[j][2])*(2/(dx*dx));

  }

//right

for(int j = 1; j <= Ny; j++)
   {
       u[j][Nx]=0;
       u1[j][Nx]=0;
 }

   for(int j = 1; j <= Ny; j++)
   {
    v[j][Nx]=0;
    v1[j][Nx]=0;
   }

  for(int j = 1; j <= Ny; j++)
  {
       s[j][Nx]=0;
       s1[j][Nx]=0;
  } 
for(int j = 1; j <= Ny; j++)
  {
        w[j][Nx] = (s[j][Nx]-s[j][Nx-1])*(2/(dx*dx));
        w1[j][Nx] = (s1[j][Nx]-s1[j][Nx-1])*(2/(dx*dx));

  }


//top

  for(int j = 1; j <= Nx; j++)
    {
        u[1][j]=1;
        u1[1][j]=1;
    }

    for(int j = 1; j <= Nx; j++)
    {
        v[1][j]=0;
        v1[1][j]=0;
    }


  for(int i = 1; i <= Nx; i++)
  {
        s[1][i]= 0;
        s1[1][i]= 0;
  }
for(int i = 1; i <= Nx; i++)
  {
        w[1][i] = (s[1][i]-s[2][i])*(2/(dy*dy))+(2*1/dy);
        w1[1][i] = (s1[1][i]-s1[2][i])*(2/(dy*dy))+(2*1/dy);

  }

//bottom

  for(int j = 1; j <= Nx; j++)
    {
        u[Ny][j]=0;
        u1[Ny][j]=0;
    }

    for(int j = 1; j <= Nx; j++)
    {
        v[Ny][j]=0;
        v1[Ny][j]=0;
    }
  for(int i = 1; i <= Nx; i++)
  {
        s[Ny][i]= 0;
        s1[Ny][i]= 0;
  }

  for(int i = 1; i <= Nx; i++)
  {
        w[Ny][i] = (s[Ny][i]-s[Ny-1][i])*(2/(dy*dy));
        w1[Ny][i] = (s1[Ny][i]-s1[Ny-1][i])*(2/(dy*dy));
 }

//vorticity solver

     for(int i = 2; i < Ny; i++)
        {
            for(int j = 2; j < Nx; j++)
            {
              w1[i][j] = (((0.01)*(w[i+1][j]-4*w[i][j]+w[i-1][j]+w[i][j+1]+w[i][j-1])*(1/(dx*dx)))-(v[i][j]*(w[i+1][j]-w[i-1][j])-(u[i][j]*(w[i][j+1]-w[i][j-1])))*(1/(2*dx)))*dt + w[i][j];
            
            }
        } 

        //creating the difference matrix 1 for w

        for(int i = 2; i < Ny; i++)
        {
           for(int j = 2; j < Nx; j++) 
            {
               d1[i][j] = fabs(w1[i][j]-w[i][j]); 
                
            }
        }

//swap

      for(int i = 1; i <= Ny; i++)
        {
            for(int j = 1; j <= Nx; j++)
            {

                w[i][j] = w1[i][j];

            }
        }  

  //finding the minimum value w

       max1 = d1[1][1];
        
             
       for (int i = 2; i < Ny; i++) 
        {
            for (int j = 2; j < Nx; j++)
        
               {
                   if(d1[i][j] > max1)  
                      {
                      max1 = d1[i][j];  
                      }
               }

       } 
       
       printf("\nTotal Variation  : 1st = %lf, 2nd = %lf \n",max1, max2);     

    }while(max1>0.001);

  gettimeofday(&stop, NULL);

    printf("\n");

    printf("\n\nValues of stream function after iteration Number: %d\n\n", it);

///resultant velocity

    for(int i = 1; i <= Ny; i++)
        {
            for(int j = 1; j <= Nx; j++)
            {

                p[i][j] = (u[i][j]*u[i][j])+(v[i][j]*v[i][j]);

            }
        }  

        for(int i = 1; i <= Ny; i++)
        {
            for(int j = 1; j <= Nx; j++)
            {

                m[i][j] = pow((p[i][j]), 0.5);

            }
        } 


  for(int i = 1; i <= Ny; i++)
 {
     for(int j = 1; j <= Nx; j++)
    {
      
      printf(" %lf ",s[i][j]);
    } 

    printf("\n");
      

  } 

   printf("\n\nValues of vorticity after iteration Number: %d\n\n", it);


  for(int i = 1; i <= Ny; i++)
 {
     for(int j = 1; j <= Nx; j++)
    {
      
      printf(" %lf ",w[i][j]);
    } 

    printf("\n");
 }
printf("\n\nValues of u velocity after iteration Number: %d\n\n", it);


  for(int i = 1; i <= Ny; i++)
 {
     for(int j = 1; j <= Nx; j++)
    {
      
      printf(" %lf ",u[i][j]);
    } 

    printf("\n");
 }   

 printf("\n\nValues of v velocity after iteration Number: %d\n\n", it);


  for(int i = 1; i <= Ny; i++)
 {
     for(int j = 1; j <= Nx; j++)
    {
      
      printf(" %lf ",v[i][j]);
    } 

    printf("\n");   
}

    printf("\nComputational Time  : %lu microsecond\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);

    printf("\nNumber of Iteration : %d\n\n", it);


//excel

stream(s);
vorticity(w);
velocity(m);
prs(pr);

 return 0;

}

void stream(double s[101][101])
{
  FILE *fp;
    fp = fopen("stream_function.csv","w");

    for(int i = 1; i <= 101; i++)
 {
     for(int j = 1; j <= 101; j++)
    {
    fprintf( fp, "%lf,",s[i][j]);
    }
    fprintf( fp,"\n");
 }

 fclose(fp);
 fp =0;
 system("pause");
  
}

void vorticity(double w[101][101])
{
  FILE *fp;
    fp = fopen("vorticity.csv","w");

    for(int i = 1; i <= 101; i++)
 {
     for(int j = 1; j <= 101; j++)
    {
    fprintf( fp, "%lf,",w[i][j]);
    }
    fprintf( fp,"\n");
 }

 fclose(fp);
 fp =0;
 system("pause");
  
}

void velocity(double m[101][101])
{
  FILE *fp;
    fp = fopen("velocity.csv","w");

    for(int i = 1; i <= 101; i++)
 {
     for(int j = 1; j <= 101; j++)
    {
    fprintf( fp, "%lf,",m[i][j]);
    }
    fprintf( fp,"\n");
 }

 fclose(fp);
 fp =0;
 system("pause");
  
}

void prs(double s[101][101])
{
  FILE *fp;
    //fp = fopen("pressure.csv","w");

    for(int i = 1; i <= 101; i++)
 {
     for(int j = 1; j <= 101; j++)
    {
    //fprintf( fp, "%lf,",pr[i][j]);
    }
   // fprintf( fp,"\n");
 }

 fclose(fp);
 fp =0;
 system("pause");
  
}