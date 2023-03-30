#include<stdio.h>
#include <math.h>
#include <cstdlib>
#include <ctime>

#define numdigit 9
#define square 3
#define Ne 120
#define Nc 200
#define Ng 100000

int indexOf(int x,int y){
	return (numdigit*x+y);
}

using namespace std;

typedef struct{
	int A[10];
	int size;
}List;

void initList(List *L){
	L->size = 0;
}

void initListWithZero(List *L, int size){
	L->size = size;
	int i;
	for(i = 0; i < size; i++)
		L->A[i] = 0;
}

void appendList(List *L, int x){
	L->A[L->size++] = x;
}

void erase(int x, List *L){
	int pos;
	for(pos = 0; pos < L->size; pos++)
		if(L->A[pos] == x) break;
	int i;
	for(i = pos; i < L->size-1; i++)
		L->A[i] = L->A[i+1];
	L->size--;
}

typedef struct{
	int values[numdigit][numdigit];
	float fitness;
}Candidate;

void initCandidate(Candidate* candidate){
	int i, j;
	for(i = 0; i < numdigit; i++)
		for(j = 0; j < numdigit; j++)
			candidate->values[i][j] = 0;
	candidate->fitness = 0.0;
}

int isRowDuplicate(Candidate candidate, int row, int value){
	int col;
	for(col = 0; col < numdigit; col++)
		if(candidate.values[row][col] == value)	return 1;
	return 0;
}

int isColDuplicate(Candidate candidate, int col, int value){
	int row;
	for(row = 0; row < numdigit; row++)
		if(candidate.values[row][col] == value)	return 1;
	return 0;
}

int isBlockDuplicate(Candidate candidate, int row, int col, int value){
	int i, j;
	int areaX = (row/square)*square;
	int areaY = (col/square)*square;
	for(i = 0; i < square; i++)
		for(j = 0; j < square; j++)
			if(candidate.values[areaX+i][areaY+j] == value)
				return 1;
	return 0;
}

void updateFitness(Candidate *candidate){
	List columnCount;
	initListWithZero(&columnCount, numdigit);
	List blockCount;
	initListWithZero(&blockCount, numdigit);
	float columnSum = 0;
	float blockSum = 0;
	int i, j;
	for(i = 0; i < numdigit; i++)
	{
		float nonzero = 0;
		for(j = 0; j < numdigit; j++)
			columnCount.A[candidate->values[j][i] - 1] += 1;
		int k;
		for(k = 0; k < numdigit; k++)
			if(columnCount.A[k] != 0)
				nonzero += 1;
		nonzero = nonzero/numdigit;
		columnSum = columnSum + nonzero;
		initListWithZero(&columnCount, numdigit);
	}
	columnSum = columnSum/numdigit;
	for(i = 0; i < numdigit; i += square)
	{
		for(j = 0; j < numdigit; j += square)
		{
	        blockCount.A[candidate->values[i][j]-1] += 1;
    	    blockCount.A[candidate->values[i][j+1]-1] += 1;
        	blockCount.A[candidate->values[i][j+2]-1] += 1;
                
        	blockCount.A[candidate->values[i+1][j]-1] += 1;
        	blockCount.A[candidate->values[i+1][j+1]-1] += 1;
        	blockCount.A[candidate->values[i+1][j+2]-1] += 1;
                
        	blockCount.A[candidate->values[i+2][j]-1] += 1;
        	blockCount.A[candidate->values[i+2][j+1]-1] += 1;
        	blockCount.A[candidate->values[i+2][j+2]-1] += 1;		
			float nonzero = 0;
			int k;
			for(k = 0; k < numdigit; k++)
				if(blockCount.A[k] != 0)
					nonzero += 1;
			nonzero = nonzero/numdigit;
			blockSum = blockSum + nonzero;
			initListWithZero(&blockCount, numdigit);
		}
	}
	blockSum = blockSum/numdigit;
	float fitness;
	if((int)columnSum == 1 && (int)blockSum == 1)
		fitness = 1.0;
	else	fitness = columnSum*blockSum;
	candidate->fitness = fitness;
}

float randFloat(float min, float max){
    float scale = rand()/(float)RAND_MAX;
    return min + scale * (max - min); 
}

void printMap(Candidate sudoku){
	int i,j;
	printf("Sudoku:\n");
	for(i = 0; i < numdigit; i++)
	{
		if(i%square == 0) printf("\n");
		for(j = 0; j < numdigit; j++)
		{
			if(j%square == 0) printf("| ");
			printf("%d ", sudoku.values[i][j]);
		}
		printf("|\n");
	}
	printf("-------------------------\n");
}

int mutate(Candidate* candidate, float mutationRate, Candidate given){
	int r = randFloat(0, 1.1);
	while(r > 1)	r = randFloat(0, 1.1);
	int success = 0;
	if(r < mutationRate)
	{
		while(!success)
		{
			int row = rand()%numdigit;
			int fromCol = rand()%numdigit;
			int toCol = rand()%numdigit;
			while(fromCol == toCol)
			{
				fromCol = rand()%numdigit;
				toCol = rand()%numdigit;
			}
			if(given.values[row][fromCol] == 0 && given.values[row][toCol] == 0)
			{
				if(!isColDuplicate(given, toCol, candidate->values[row][fromCol])
				&& !isColDuplicate(given, fromCol, candidate->values[row][toCol])
				&& !isBlockDuplicate(given, row, toCol, candidate->values[row][fromCol])
				&& !isBlockDuplicate(given, row, fromCol, candidate->values[row][toCol])
				)
				{
					int temp = candidate->values[row][toCol];
					candidate->values[row][toCol] = candidate->values[row][fromCol];
					candidate->values[row][fromCol] = temp;
					success = 1;
				}
			}
		}
	}
	return success;
}

typedef struct{
	Candidate candidates[Nc];
	int numOfGen;
}Population;

void initPopulation(Population* population){
	population->numOfGen = 0;
}

int set(List row){
	int availables[numdigit+1];
	int i;
	for(i = 1; i < numdigit + 1; i++) availables[i] = 1;
	for(i = 0; i < row.size; i++)
		availables[row.A[i]] = 0;
	for(i = 1; i < numdigit+1; i++)
		if(availables[i]) return 0;
	return 1;
}

void rowAppend(Candidate *candidate, List x, int row){
	int col;
	for(col = 0; col < numdigit; col++)
		candidate->values[row][col] = x.A[col];
}

void candidateAppend(Population *population, Candidate given){
	population->candidates[population->numOfGen++] = given;
}

void updatePopulationFitness(Population *population){
	int i;
	for(i = 0; i < population->numOfGen; i++)
		updateFitness(&population->candidates[i]);
}

void seed(Population *population, Candidate given){
	initPopulation(population);
	int i;
	List available[numdigit*numdigit];
	for(i = 0; i < numdigit*numdigit; i++)
		initList(&available[i]);
	int row, col, value;		
	for(row = 0; row < numdigit; row++)
	{
		for(col = 0; col < numdigit; col++)
		{
			for(value = 1; value <= numdigit; value++)
			{
				if((given.values[row][col] == 0) &&
				!(isColDuplicate(given, col, value) || isRowDuplicate(given, row, value) || isBlockDuplicate(given, row, col, value)))
					appendList(&available[numdigit*row + col], value);
				else if(given.values[row][col] != 0)
				{
					appendList(&available[numdigit*row + col], given.values[row][col]);
					break;
				}
			}
		}
	}
	int p;
	for(p = 0; p < Nc; p++)
	{
		Candidate g;
		initCandidate(&g);
		int i;
		for(i = 0; i < numdigit; i++)
		{
			List row;
			initListWithZero(&row, numdigit);
			int j;
			for(j = 0; j < numdigit; j++)
			{
				if(given.values[i][j] != 0)
					row.A[j] = given.values[i][j];
				else if(given.values[i][j] == 0)
					row.A[j] = available[i*numdigit + j].A[rand()%available[i*numdigit + j].size];
			}
			while(!set(row))
			{
				int j;
				for(j = 0; j < numdigit; j++)
					row.A[j] = available[i*numdigit + j].A[rand()%available[i*numdigit + j].size];				
			}
			rowAppend(&g, row, i);
		}
		candidateAppend(population, g);
	}
	updatePopulationFitness(population);
	printf("Seeding complete.");
}

void sort(Population *population){
	int i, j;
	for(i = 0; i < population->numOfGen-1; i++)
	{
		int max = i;
		for(j = i+1; j < population->numOfGen; j++)
		{
			if(population->candidates[max].fitness < population->candidates[j].fitness)
			{
				max = j;
				Candidate temp = population->candidates[i];
				population->candidates[i] = population->candidates[max];
				population->candidates[max] = temp;				
			}

		}
	}
}

Candidate compete(Population population){
	Candidate c1 = population.candidates[rand()%population.numOfGen];
	Candidate c2 = population.candidates[rand()%population.numOfGen];
	float f1 = c1.fitness;
	float f2 = c2.fitness;
	Candidate fittest;
	Candidate weakest;
	if(f1 > f2)
	{
		fittest = c1;
		weakest = c2;
	}
	else
	{
		fittest = c2;
		weakest = c1;		
	}
	float selection_rate = 0.85;
	float r = randFloat(0, 1.1);
	while(r > 1)	r = randFloat(0, 1.1);
	if(r < selection_rate)	return fittest;
	return weakest;
}

void copy(Candidate parent, Candidate *child){
	int i, j;
	for(i = 0; i < numdigit; i++)
		for(j = 0; j < numdigit; j++)
			child->values[i][j] = parent.values[i][j];
	child->fitness = parent.fitness;
}

int findValue(List parent_row, int value){
	int i;
	for(i = 0; i < parent_row.size; i++)
		if(parent_row.A[i] == value)
			return i + 1;
	return 0;
}

int findUnused(List parent_row, List remaining){
	int i;
	for(i = 0; i < parent_row.size; i++)
		if(findValue(remaining, parent_row.A[i]))
			return i + 1;
	return 0;
}

void crossoverRows(List row1, List row2, Candidate *child1, Candidate *child2, int p){
	List childRow1;
	List childRow2;
	initListWithZero(&childRow1, numdigit);
	initListWithZero(&childRow2, numdigit);
	List remaining;	
	initList(&remaining);
	int i;
	for(i = 1; i < numdigit+1; i++)
		appendList(&remaining, i);
	int cycle = 0;
	while(findValue(childRow1, 0) && findValue(childRow2, 0))
	{
		if(cycle%2 == 0)
		{
			int index = findUnused(row1, remaining) - 1;
			int start = row1.A[index];
			erase(row1.A[index], &remaining);
			childRow1.A[index] = row1.A[index];
			childRow2.A[index] = row2.A[index];
			int next = row2.A[index];
			while(next != start)
			{
				index = findValue(row1, next) - 1;
				childRow1.A[index] = row1.A[index];
				erase(row1.A[index], &remaining);
				childRow2.A[index] = row2.A[index];
				next = row2.A[index];
			}
			cycle++;
		}
		else
		{
			int index = findUnused(row1, remaining) - 1;
			int start = row1.A[index];
			erase(row1.A[index], &remaining);
			childRow1.A[index] = row2.A[index];
			childRow2.A[index] = row1.A[index];
			int next = row2.A[index];
			while(next != start)
			{
				index = findValue(row1, next) - 1;
				childRow1.A[index] = row2.A[index];
				erase(row1.A[index], &remaining);
				childRow2.A[index] = row1.A[index];
				next = row2.A[index];
			}
			cycle++;
		}
	}
	rowAppend(child1, childRow1, p);
	rowAppend(child2, childRow2, p);
}

List rowGetter(Candidate candidate, int row){
	int col;
	List a;
	initList(&a);
	for(col = 0; col < numdigit; col++)
		appendList(&a, candidate.values[row][col]);
	return a;
}

void crossover(Candidate parent1, Candidate parent2, Candidate *child1, Candidate *child2, float crossoverRate){
	copy(parent1, child1);
	copy(parent2, child2);
	float r = randFloat(0, 1.1);
	while(r > 1)	r = randFloat(0, 1.1);
	if(r < crossoverRate)
	{
		int crossoverPoint1 = rand()%9;
		int crossoverPoint2 = 1+rand()%9;
		while(crossoverPoint1 == crossoverPoint2)
		{
			crossoverPoint1 = rand()%9;
			crossoverPoint2 = 1+rand()%9;			
		}
		if(crossoverPoint1 > crossoverPoint2)
		{
			int temp = crossoverPoint1;
			crossoverPoint1 = crossoverPoint2;
			crossoverPoint2 = temp;
		}
		int i;
		for(i = crossoverPoint1; i < crossoverPoint2; i++)
		{
			List row1 = rowGetter(*child1, i);
			List row2 = rowGetter(*child2, i);
			crossoverRows(row1, row2, child1, child2, i);
		}
	}
}

Candidate solve(Candidate given){
	int Nm = 0;
	//dem so lan the he dang theo doi
	int staleCount = 0;
	float prevFitness = 0;
	//dem so lan con tot hon cha, me
	//duoc su dung de cap nhat ty le dot bien
	float mutationRate = 0.5;
	//tao quan the ban dau HOAC seed
	Population population;
	initPopulation(&population);
	seed(&population, given);
	Candidate bestSolution;
	int g;//generation
	for(g = 0; g < Ng; ++g){
		printf("\nGeneration %d", g);

		float bestFitness = 0.0;
		bestSolution = given;
		// doi voi moi the he, duyet qua tat ca
		//cac ung vien hoac nhiem sac the de kiem tra giai phap
		int c;
		for(c=0; c<Nc; ++c){
			float fitness=population.candidates[c].fitness;
			if(int(fitness) == 1){
				printf("\nSolution found at generation %d!", g);
				printMap(population.candidates[c]);
				return population.candidates[c];
			}

			if(fitness > bestFitness){
				bestFitness = fitness;
				bestSolution = population.candidates[c];
			}
		}
		printf("\nBest fitness: %f", bestFitness);

		Population nextPopulation;
		initPopulation(&nextPopulation);
		sort(&population);

		Population elites;
		initPopulation(&elites);
		int e;
		for(e = 0; e < Ne; ++e){
			Candidate elite;
			initCandidate(&elite);
			elite = population.candidates[e];
			candidateAppend(&elites, elite);
		}
		
		int count;
		for(count = Ne; count < Nc; count += 2){
			Candidate parent2 = compete(population);
						
			Candidate child1;
			initCandidate(&child1);
			Candidate child2;
			initCandidate(&child2);
			crossover( parent1, parent2, &child1, &child2, 1.0);

			
			updateFitness(&child1);
			float oldFitness = child1.fitness;
			//Dot bien
			int success = mutate(&child1, mutationRate, given);
			float newFitness = child1.fitness;
			updateFitness(&child1);
			if(success){
				Nm += 1;
			}					
			//# Mutate child2.
			updateFitness(&child2);
			oldFitness = child2.fitness;		
			//Dot bien
			success = mutate(&child2, mutationRate, given);
			newFitness = child2.fitness;
			updateFitness(&child2);			
			if(success){
				Nm += 1;
			}		
            candidateAppend(&nextPopulation, child1);
        	candidateAppend(&nextPopulation, child2);
		}

	// Ket nap elites vao phan con lai cua nextPopulation
	//chung se khong bi anh huong boi  crossover hoac dot bien
	//hoac dot bien
	
	for(e=0; e<Ne; ++e)
        candidateAppend(&nextPopulation, elites.candidates[e]);
    initPopulation(&population);
	population = nextPopulation;
    updatePopulationFitness(&population);
	mutationRate = randFloat(0, 1.1);
	while(mutationRate > 1)
		mutationRate = randFloat(0, 1.1);	
		if(prevFitness == bestFitness) staleCount += 1;else
			if(prevFitness!=bestFitness){
				staleCount = 0;
                prevFitness = bestFitness;
				}
            if(staleCount >= 100){
            	printf("\nThe population has gone stale. Re-seeding...");
                seed(&population,given);
                staleCount = 0;
                mutationRate = 0.5;
			}		
	}
	printf("\nNo solution found.\n");
	printMap(bestSolution);
}

void readMap(Candidate* sudoku){
//	freopen("easy1.txt","r",stdin);
	freopen("medium1.txt","r",stdin);
//	freopen("hard1.txt","r",stdin);
//	freopen("AIEscarcot.txt", "r", stdin);
	initCandidate(sudoku);
	int i,j;
	for(i = 0; i < numdigit; i++)
		for(j = 0; j < numdigit; j++)
			scanf("%d", &sudoku->values[i][j]);
}

int main(){
	Candidate sudoku;
	readMap(&sudoku);
	printMap(sudoku);
	srand(time(0));	
	Candidate solution = solve(sudoku);
	return 0;
}
