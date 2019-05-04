#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/timeb.h>
#include <string.h>
#include <string> 
#include "knapsack_template.h" 
#include <bits/stdc++.h>
#define TIMEB struct timeb
#define FTIME ftime
#define UDT_TIME long
#define MAX_SIZE_TO_PRINT 10

UDT_TIME gRefTime = 0;
int Maxs=0;
UDT_TIME GetMilliSecondTime(TIMEB timeBuf);
void SetTime(void);
UDT_TIME GetTime(void);

bool byratio(std::pair<int,int> a,
              std::pair<int,int> b)
{
  double r1=(double) a.second/a.first;
  double r2=(double) b.second/b.first;
    return (r1 > r2);
}
int main(int argc, char* argv[])
{
	UDT_TIME time, BFTime,B3;
    float speedup;
	int itemCnt;
	KnapsackInstance* inst; //a Knapsack instance object
	KnapsackDPSolver DPSolver; //dynamic programming solver
	KnapsackBFSolver BFSolver; //brute-force solver
	KnapsackBTSolver BTSolver; //backtracking solver
	KnapsackBBSolver BBSolver1(UB1); //branch-and-bound solver with UB1
	KnapsackBBSolver BBSolver2(UB2); //branch-and-bound solver with UB2
	KnapsackBBSolver BBSolver3(UB3); //branch-and-bound solver with UB3
        KnapsackBBSolver BBSolver4(UB4); //branch-and-bound solver with UB3
	KnapsackSolution *DPSoln, *BFSoln, *BTSoln, *BBSoln1, *BBSoln2, *BBSoln3,*BBSoln4;

    if(argc != 2) {
        printf("Invalid Number of command-line arguments\n");
	    exit(1);
    }
    itemCnt = atoi(argv[1]);
    if(itemCnt < 1) {
        printf("Invalid number of items\n");
        exit(1);
    }

	inst   =   new KnapsackInstance(itemCnt);
	DPSoln = new KnapsackSolution(inst);
	BFSoln = new KnapsackSolution(inst);
	BTSoln = new KnapsackSolution(inst);
	BBSoln1= new KnapsackSolution(inst);
	BBSoln2= new KnapsackSolution(inst);
	BBSoln3= new KnapsackSolution(inst);
	BBSoln4= new KnapsackSolution(inst);

	inst->Generate();
	inst->Print();

	SetTime();
	DPSolver.Solve(inst,DPSoln);
	time = GetTime();
	printf("\n\nSolved using dynamic programming (DP) in %ld ms. Optimal value = %d",time, DPSoln->GetValue());
	if(itemCnt <= MAX_SIZE_TO_PRINT)
		DPSoln->Print("Dynamic Programming Solution");


	SetTime();
	BFSolver.Solve(inst,BFSoln);
	BFTime = time = GetTime();
	printf("\n\nSolved using brute-force enumeration (BF) in %ld ms. Optimal value = %d",time, BFSoln->GetValue());
	if(itemCnt <= MAX_SIZE_TO_PRINT)
		BFSoln->Print("Brute-Force Solution");
	if(*DPSoln == *BFSoln)
		printf("\nSUCCESS: DP and BF solutions match");
	else
		printf("\nERROR: DP and BF solutions mismatch");


	SetTime();
	BTSolver.Solve(inst,BTSoln);
	time = GetTime();
	printf("\n\nSolved using backtracking (BT) in %ld ms. Optimal value = %d",time, BTSoln->GetValue());
	if(itemCnt <= MAX_SIZE_TO_PRINT)
		BTSoln->Print("Backtracking Solution");
	if(*BFSoln == *BTSoln)
		printf("\nSUCCESS: BF and BT solutions match");
	else
		printf("\nERROR: BF and BT solutions mismatch");
    speedup = time==0? 0 : 100.0*(BFTime-time)/(float)BFTime;
	printf("\nSpeedup of BT relative to BF is %.2f%c",speedup,'%');


	SetTime();
	BBSolver1.Solve(inst,BBSoln1);
	time = GetTime();
	printf("\n\nSolved using branch-and-bound (BB) with UB1 in %ld ms. Optimal value = %d",time, BBSoln1->GetValue());
	if(itemCnt <= MAX_SIZE_TO_PRINT)
		BBSoln1->Print("BB-UB1 Solution");
	if(*BFSoln == *BBSoln1)
		printf("\nSUCCESS: BF and BB-UB1 solutions match");
	else
		printf("\nERROR: BF and BB-UB1 solutions mismatch");
    speedup = time==0? 0 : 100.0*(BFTime-time)/(float)BFTime;
	printf("\nSpeedup of BB-UB1 relative to BF is %.2f%c",speedup,'%');


	SetTime();
	BBSolver2.Solve(inst,BBSoln2);
	time = GetTime();
	printf("\n\nSolved using branch-and-bound (BB) with UB2 in %ld ms. Optimal value = %d",time, BBSoln2->GetValue());
	if(itemCnt <= MAX_SIZE_TO_PRINT)
		BBSoln2->Print("BB-UB2 Solution");
	if(*BFSoln == *BBSoln2)
		printf("\nSUCCESS: BF and BB-UB2 solutions match");
	else
		printf("\nERROR: BF and BB-UB2 solutions mismatch");
    speedup = time==0? 0 : 100.0*(BFTime-time)/(float)BFTime;
	printf("\nSpeedup of BB-UB2 relative to BF is %.2f%c",speedup,'%');
       


int n=inst->GetItemCnt();
std::pair <int,int> ps[n+1];
ps[0].first=0;
ps[0].second=0;

for(int i=1;i<=n;i++)
{
ps[i].first=inst->GetItemWeight(i);
ps[i].second=inst->GetItemValue(i);
}

sort(ps+1,ps+n+1,byratio);

for(int i=1;i<=n;i++)
{
inst->SetItemWeight(i,ps[i].first);
inst->SetItemValue(i,ps[i].second);
}      

	SetTime();
	BBSolver3.Solve(inst,BBSoln3);
	B3=time = GetTime();
	printf("\n\nSolved using branch-and-bound (BB) with UB3 in %ld ms. Optimal value = %d",time, BBSoln3->GetValue());
	if(itemCnt <= MAX_SIZE_TO_PRINT)
		BBSoln3->Print("BB-UB3 Solution");
	if(*BFSoln == *BBSoln3)
		printf("\nSUCCESS: BF and BB-UB3 solutions match");
	else
		printf("\nERROR: BF and BB-UB3 solutions mismatch");
    speedup = time==0? 0 : 100.0*(BFTime-time)/(float)BFTime;
	printf("\nSpeedup of BB-UB3 relative to BF is %.2f%c",speedup,'%');

SetTime();
        BBSolver4.Solve(inst,BBSoln4);
        time = GetTime();
        printf("\n\nSolved using branch-and-bound (BB) with UB4 in %ld ms. Optimal value = %d",time, BBSoln4->GetValue());
        if(itemCnt <= MAX_SIZE_TO_PRINT)
                BBSoln4->Print("BB-UB4 Solution");
        if(*BFSoln == *BBSoln4)
                printf("\nSUCCESS: BF and BB-UB4 solutions match");
        else
            	printf("\nERROR: BF and BB-UB4 solutions mismatch");
    speedup = time==0? 0 : 100.0*(BFTime-time)/(float)BFTime;
        printf("\nSpeedup of BB-UB4 relative to BF is %.2f%c",speedup,'%');
    speedup = time==0? 0 : 100.0*(B3-time)/(float)B3;
      printf("\nSpeedup of BB-UB4 relative to BB-UB3 is %.2f%c",speedup,'%');


	
	delete inst;
	delete DPSoln;
	delete BFSoln;
	delete BTSoln;
	delete BBSoln1;
	delete BBSoln2;
	delete BBSoln3;
        delete BBSoln4; 
	printf("\n\nProgram Completed Successfully\n");

	return 0;
}
/********************************************************************/

KnapsackInstance::KnapsackInstance(int itemCnt_)
{
	itemCnt = itemCnt_;

	weights = new int[itemCnt+1];
	values = new int[itemCnt+1];
	cap = 0;
}
/********************************************************************/
	
KnapsackInstance::~KnapsackInstance()
{
	delete [] weights;
	delete [] values;
}
/********************************************************************/

void KnapsackInstance::Generate()
{
    int i, wghtSum;  
        
    weights[0] = 0;
    values[0] = 0;
        
    wghtSum = 0;
    for(i=1; i<= itemCnt; i++)
    {
        weights[i] = rand()%100 + 1;
        values[i] = weights[i] + 10;
        wghtSum += weights[i]; 
    }
    cap = wghtSum/2;
}
/********************************************************************/

int KnapsackInstance::GetItemCnt()
{
	return itemCnt;
}
/********************************************************************/

int KnapsackInstance::GetItemWeight(int itemNum)
{
	return weights[itemNum];
}
/********************************************************************/

int KnapsackInstance::GetItemValue(int itemNum)
{
	return values[itemNum];
}
void KnapsackInstance::SetItemWeight(int itemNum,int wt)
{
        weights[itemNum]=wt;
}
/********************************************************************/

void KnapsackInstance::SetItemValue(int itemNum,int val)
{
        values[itemNum]=val;
}
/********************************************************************/

int KnapsackInstance::GetCapacity()
{
	return cap;
}
/********************************************************************/

void KnapsackInstance::Print()
{
	int i;

	printf("Number of items = %d, Capacity = %d\n",itemCnt, cap);
	printf("Weights: ");
	for(i=1; i<=itemCnt; i++)
	{
		printf("%d ",weights[i]);
	}
	printf("\nValues: ");
	for(i=1; i<=itemCnt; i++)
	{
		printf("%d ",values[i]);
	}
	printf("\n");
}
/*****************************************************************************/

KnapsackSolution::KnapsackSolution(KnapsackInstance* inst_)
{
	int i, itemCnt = inst_->GetItemCnt();

	inst = inst_;
	isTaken = new bool[itemCnt+1];
	value = INVALID_VALUE;

	for(i=1; i<=itemCnt; i++)
	{
		isTaken[i] = false;
	}
}
/********************************************************************/

KnapsackSolution::~KnapsackSolution()
{
	delete [] isTaken;
}
/********************************************************************/

bool KnapsackSolution::operator == (KnapsackSolution& otherSoln)
{
	return value == otherSoln.value;
}
/********************************************************************/

void KnapsackSolution::TakeItem(int itemNum)
{
	isTaken[itemNum] = true;
}
/********************************************************************/
	
void KnapsackSolution::DontTakeItem(int itemNum)
{
	isTaken[itemNum] = false;	
}
/********************************************************************/

int KnapsackSolution::ComputeValue()
{
	int i, itemCnt = inst->GetItemCnt(), weight = 0;

	value = 0;
	for(i=1; i<=itemCnt; i++)
	{
		if(isTaken[i] == true)
		{
			weight += inst->GetItemWeight(i);
			if(weight > inst->GetCapacity())
			{
				value = INVALID_VALUE;
				break;
			}
			value += inst->GetItemValue(i);
		}
	}
	return value;
}
/********************************************************************/

int KnapsackSolution::GetValue()
{
	return value;
}

/********************************************************************/

void KnapsackSolution::Copy(KnapsackSolution* otherSoln)
{
	int i, itemCnt = inst->GetItemCnt();

	for(i=1; i<=itemCnt; i++)
	{
		isTaken[i] = otherSoln->isTaken[i];
	}
	value = otherSoln->value;
}
/********************************************************************/

void KnapsackSolution::Print(std::string title)
{
	int i, itemCnt = inst->GetItemCnt();

	printf("\n%s: ",title.c_str());
	for(i=1; i<=itemCnt; i++)
	{
		if(isTaken[i] == true)
			printf("%d ",i);
	}
	printf("\nValue = %d\n",value);
	
}
/*****************************************************************************/

KnapsackBFSolver::KnapsackBFSolver()
{
	crntSoln = NULL;
}
/********************************************************************/

KnapsackBFSolver::~KnapsackBFSolver()
{
	if(crntSoln != NULL)
		delete crntSoln;
}
/********************************************************************/

void KnapsackBFSolver::Solve(KnapsackInstance* inst_, KnapsackSolution* soln_)
{
	inst = inst_;	
	bestSoln = soln_;
	crntSoln = new KnapsackSolution(inst);
	FindSolns(1);
}
/********************************************************************/

void KnapsackBFSolver::FindSolns(int itemNum)
{
	int itemCnt = inst->GetItemCnt();

	if(itemNum == itemCnt + 1)
	{
		CheckCrntSoln();
		return;
	}
	crntSoln->DontTakeItem(itemNum);
	FindSolns(itemNum+1);
	crntSoln->TakeItem(itemNum);
	FindSolns(itemNum+1);
}
/********************************************************************/

void KnapsackBFSolver::CheckCrntSoln()
{
	int crntVal = crntSoln->ComputeValue();

#ifdef KNAPSACK_DEBUG
	printf("\nChecking solution ");
	crntSoln->Print(" ");
#endif

	if(crntVal == INVALID_VALUE)
		return;

	if(bestSoln->GetValue() == INVALID_VALUE) //The first solution is initially the best solution
		bestSoln->Copy(crntSoln);
	else
	{
		if(crntVal > bestSoln->GetValue())
			bestSoln->Copy(crntSoln);
	}
}
/********************************************************************/

// Write code below to implement the DP solver, the backtracking (BT) solver
// and the Branch-and-Bound (BB) solver. 
// Note that the BT and BB solvers inherit from the brute-force solver.
// You may add any private data members or functions that you need.
// Your solve() function takes an object of class KnapsackInstance as input
// and produces an object of class KnapsackSolution as output. 
// See how the given KnapsackBFSolver::Solve() writes its result into the 
// KnapsackSolution object and make the solvers that you write do the same.

KnapsackDPSolver::KnapsackDPSolver()
{

}

KnapsackDPSolver::~KnapsackDPSolver()
{

}

void KnapsackDPSolver::Solve(KnapsackInstance* inst_, KnapsackSolution* soln_)
{
	inst = inst_;
	soln = soln_;
    int n=inst->GetItemCnt()+1;
    int c=inst->GetCapacity()+1;
    long temp[n][c];
    for(int i=0;i<n;i++)
    {
        temp[i][0]=0;
    }
    for(int i=0;i<c;i++)
    {
        temp[0][i]=0;
    }
    for(int i=1;i<n;i++)
    {
        for(int j=1;j<c;j++)
        {
            if((inst->GetItemWeight(i))>j)
                temp[i][j]=temp[i-1][j];
            else
            {
            if((inst->GetItemValue(i))+temp[i-1][j-(inst->GetItemWeight(i))]>temp[i-1][j])
            {
                temp[i][j]=(inst->GetItemValue(i))+temp[i-1][j-(inst->GetItemWeight(i))];
            }
            else
            {
                temp[i][j]=temp[i-1][j];
            }
          }
       }
    }
    int i=n-1,j=c-1;
    for(;;)
    {
        if(i==0 || j==0) break;
        else if(temp[i][j]>temp[i-1][j])
        {
            soln->TakeItem(i);
            j=j-(inst->GetItemWeight(i));
            i=i-1;
            
        }
        else
            i=i-1;
        
    }
    
    
    soln->ComputeValue();    
    
    
    
}
/*****************************************************************************/

KnapsackBTSolver::KnapsackBTSolver()
{
crntSoln=NULL;
}
/********************************************************************/

KnapsackBTSolver::~KnapsackBTSolver()
{
if(crntSoln!=NULL) delete crntSoln;
}
/********************************************************************/

void KnapsackBTSolver::Solve(KnapsackInstance* inst_, KnapsackSolution* soln_)
{
inst=inst_;
bestSoln=soln_;
crntSoln=new KnapsackSolution(inst);

FindSolns(1,0);	
}
void KnapsackBTSolver::CheckCrntSoln()
{
        int crntVal = crntSoln->ComputeValue();

#ifdef KNAPSACK_DEBUG
        printf("\nChecking solution ");
        crntSoln->Print(" ");
#endif

      	if(crntVal == INVALID_VALUE)
                return;

        if(bestSoln->GetValue() == INVALID_VALUE) //The first solution is initially the best solution
                bestSoln->Copy(crntSoln);
        else
	{
                if(crntVal > bestSoln->GetValue())
                        bestSoln->Copy(crntSoln);
        }
}

void KnapsackBTSolver::FindSolns(int itemNum,int wt)
{        if(wt>inst->GetCapacity())
            return;       
          
        int itemCnt = inst->GetItemCnt();
        
        if(itemNum == itemCnt + 1)
        {    CheckCrntSoln();
               return;
        }
	crntSoln->DontTakeItem(itemNum);
        FindSolns(itemNum+1,wt);
        crntSoln->TakeItem(itemNum);
        FindSolns(itemNum+1,wt+(inst->GetItemWeight(itemNum)));
}

/*****************************************************************************/

KnapsackBBSolver::KnapsackBBSolver(enum UPPER_BOUND ub_)
{
	ub = ub_;
       crntSoln=NULL;
}
/********************************************************************/

KnapsackBBSolver::~KnapsackBBSolver()
{
if(crntSoln==NULL) delete crntSoln;

}
/********************************************************************/
void KnapsackBBSolver::Solve(KnapsackInstance* inst_, KnapsackSolution* soln_)
{

inst=inst_;
bestSoln=soln_;
crntSoln=new KnapsackSolution(inst);

int sumv=0;
for(int i=1;i<=inst->GetItemCnt();i++)
 sumv+=inst->GetItemValue(i);



if(ub==UB1)
FindSolns1(1,0,0,sumv,0);
else if(ub==UB2)
{
FindSolns2(1,0,0,sumv,0);
}
else if(ub==UB3)
 {
FindSolns3(1,0,0);
}	
else
{
FindSolns4(1,0,0,sumv,0);
}

}
void KnapsackBBSolver::FindSolns1(int itemNum,int valueTaken,int notvalueTaken,int sum,int wt)
{   
     
       int itemCnt = inst->GetItemCnt();
     
if((sum-notvalueTaken)<=bestSoln->GetValue())
             return;


   if(itemNum == itemCnt + 1)
        {
              if(wt<=inst->GetCapacity())
                {CheckCrntSoln();
                }
                

       return;
        }


	crntSoln->DontTakeItem(itemNum);
        FindSolns1(itemNum+1,valueTaken,notvalueTaken+inst->GetItemValue(itemNum),sum,wt);
        crntSoln->TakeItem(itemNum);
        FindSolns1(itemNum+1,valueTaken+inst->GetItemValue(itemNum),notvalueTaken,sum,wt+inst->GetItemWeight(itemNum));
}

/************/
int KnapsackBBSolver::checkweight(int itemNum,int valueTaken,int wt)
{
    int n =inst->GetItemCnt();
     int capc=inst->GetCapacity();
     
int newcap=capc-wt;
    for(int i=itemNum;i<=n;i++)
      { if(inst->GetItemWeight(i)<=newcap)
           {
             crntSoln->TakeItem(i);
            valueTaken+=inst->GetItemValue(i);
            }
         
       }

return valueTaken;
}
void KnapsackBBSolver::FindSolns2(int itemNum,int valueTaken,int notvaluetaken,int sum,int wt)
{
   if(checkweight(itemNum,valueTaken,wt)<=bestSoln->GetValue())
     return;
   
     int itemCnt = inst->GetItemCnt();
        if(itemNum == itemCnt + 1)
        {    if(wt<=inst->GetCapacity())
               	CheckCrntSoln();
              
       return;
        }     


	crntSoln->DontTakeItem(itemNum);
        FindSolns2(itemNum+1,valueTaken,notvaluetaken+inst->GetItemValue(itemNum),sum,wt);
        crntSoln->TakeItem(itemNum);
        FindSolns2(itemNum+1,valueTaken+inst->GetItemValue(itemNum),notvaluetaken,sum,wt+inst->GetItemWeight(itemNum));
}

/************/

/******/

/************/
int KnapsackBBSolver::checkweightfractional(int itemNum,int valueTaken,int wt)
{
    int n =inst->GetItemCnt();
     int capc=inst->GetCapacity();

int newcap=capc-wt;
    for(int i=itemNum;i<=n;i++)
      { if(inst->GetItemWeight(i)<=newcap)
           {
             newcap-=inst->GetItemWeight(i);
            valueTaken+=inst->GetItemValue(i);
            }
           else
               	{valueTaken+=(int )((newcap*inst->GetItemValue(i))/(inst->GetItemWeight(i)));

              newcap=0;
                 break;
                }

       }

return valueTaken;
}
void KnapsackBBSolver::FindSolns3(int itemNum,int valueTaken,int wt)
{ 
 
    if(checkweightfractional(itemNum,valueTaken,wt)<=bestSoln->GetValue())
     return;

   int itemCnt = inst->GetItemCnt();
        if(itemNum == itemCnt + 1)
        {  if(wt<=inst->GetCapacity()) 
             CheckCrntSoln();

             return;
        }

	crntSoln->DontTakeItem(itemNum);
        FindSolns3(itemNum+1,valueTaken,wt);
        crntSoln->TakeItem(itemNum);
        FindSolns3(itemNum+1,valueTaken+inst->GetItemValue(itemNum),wt+inst->GetItemWeight(itemNum));




}
/*******/
void KnapsackBBSolver::FindSolns4(int itemNum,int valueTaken,int notvaluetaken,int sum,int wt)
{


int remainingcap=inst->GetCapacity()-wt;
int res=remainingcap*inst->GetItemValue(itemNum);
  if(remainingcap>0 && inst->GetItemWeight(itemNum)>0)
  {   res=(int )(res/inst->GetItemWeight(itemNum));
  
}

if(valueTaken+res<=bestSoln->GetValue())
     return;

 

     int itemCnt = inst->GetItemCnt();
        if(itemNum == itemCnt + 1)
        {   if(wt<=inst->GetCapacity())
                CheckCrntSoln();

         return;
        }

      crntSoln->DontTakeItem(itemNum);
     FindSolns4(itemNum+1,valueTaken,notvaluetaken+inst->GetItemValue(itemNum),sum,wt);
        crntSoln->TakeItem(itemNum);
      FindSolns4(itemNum+1,valueTaken+inst->GetItemValue(itemNum),notvaluetaken,sum,wt+inst->GetItemWeight(itemNum));
}
/*****************************************************************************/

UDT_TIME GetCurrentTime(void)
{
	UDT_TIME crntTime=0;

	TIMEB timeBuf;
	FTIME(&timeBuf);
    crntTime = GetMilliSecondTime(timeBuf);

	return crntTime;
}
/********************************************************************/

void SetTime(void)
{
	gRefTime = GetCurrentTime();
}
/********************************************************************/

UDT_TIME GetTime(void)
{
	UDT_TIME crntTime = GetCurrentTime();

	return (crntTime - gRefTime);
}
/********************************************************************/

UDT_TIME GetMilliSecondTime(TIMEB timeBuf)
{
	UDT_TIME mliScndTime;

	mliScndTime = timeBuf.time;
	mliScndTime *= 1000;
	mliScndTime += timeBuf.millitm;
	return mliScndTime;
}
/*****************************************************************************/
