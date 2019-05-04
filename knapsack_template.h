
//#define KNAPSACK_DEBUG

#define INVALID_VALUE -1

enum UPPER_BOUND { UB1, UB2, UB3,UB4};

class KnapsackInstance
{
private:
	int itemCnt; //Number of items
	int cap; //The capacity
	int* weights;//An array of weights
	int* values;//An array of values

public:
	KnapsackInstance(int itemCnt_);
	~KnapsackInstance();

	void Generate();

	int GetItemCnt();
	int GetItemWeight(int itemNum);
	int GetItemValue(int itemNum);
        void SetItemValue(int itemNum,int value);
        void SetItemWeight(int itemNum,int weight); 
	int GetCapacity();
	void Print();
};

class KnapsackSolution
{
private:
	bool* isTaken;
        
	int value;
	int wght;
	KnapsackInstance* inst;

public:
	KnapsackSolution(KnapsackInstance* inst);
	~KnapsackSolution();

	bool operator == (KnapsackSolution& otherSoln);
	void TakeItem(int itemNum);
	void DontTakeItem(int itemNum);
	int ComputeValue();
	int GetValue();
	int GetWeight();
	void Print(std::string str);
	void Copy(KnapsackSolution* otherSoln);
        
       
           
};

// Dynamic programming solver
class KnapsackDPSolver
{
private:
	KnapsackInstance* inst;
	KnapsackSolution* soln;

public:
	KnapsackDPSolver();
	~KnapsackDPSolver();
	void Solve(KnapsackInstance* inst, KnapsackSolution* soln);
};


// Brute-force solver
class KnapsackBFSolver
{
protected:
	KnapsackInstance* inst;
	KnapsackSolution* crntSoln;
	KnapsackSolution* bestSoln;

	virtual void FindSolns(int itemNum);
	virtual void CheckCrntSoln();

public:
	KnapsackBFSolver();
	~KnapsackBFSolver();
	virtual void Solve(KnapsackInstance* inst, KnapsackSolution* soln);
};


// Backtracking solver
class KnapsackBTSolver: public KnapsackBFSolver
{

public:
	KnapsackBTSolver();
	~KnapsackBTSolver();
        void FindSolns(int itemNum,int wt);
        void CheckCrntSoln();
	void Solve(KnapsackInstance* inst, KnapsackSolution* soln);
};

// Branch-and-Bound solver
class KnapsackBBSolver: public KnapsackBFSolver
{
protected:
enum UPPER_BOUND ub;

public:
      
	KnapsackBBSolver(enum UPPER_BOUND ub_);
	~KnapsackBBSolver();
        void FindSolns1(int a,int b,int c,int d,int e);
        void FindSolns2(int a,int b,int	c,int d,int e);
        void FindSolns3(int a,int b,int c);
        void FindSolns4(int a,int b,int c,int d,int e);
         int checkweight(int a,int b,int c);
         int checkweightfractional(int a,int b,int c);
       
	void Solve(KnapsackInstance* inst, KnapsackSolution* soln);
};

