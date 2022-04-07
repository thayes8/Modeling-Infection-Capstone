/**
Looks at all of a current cell's neighbors and gets a random num for each
infected neighbor and generates a random num from 0 to 1, then compares that
num to the infection rate to determine if it gets infected.
Variable Definitions:
i = current x pos
j = current y pos
nRows = numRows
nCols = numColumns
tau = Transmission Rate
**/
bool infect(int Pop[][], int i, int j, float tau, int nRows){
    //Tracks whether current cell has been infected
    bool t = 0;
    //if not the leftmost wall
    if(i > 1) {
        //if left neighbor is sick
        if(Pop[i-1][j] > 0){
            t = (rand < tau);
        }
    }
    //if i is not the rightmost wall
    if(i < n) {
        //if left neighbor is sick
        if(Pop[i+1][j] > 0){
            t = t + (rand < tau);
        }
    }
    if(j > 1) {
        //if left neighbor is sick
        if(Pop[i][j-1] > 0){
            t = t + (rand < tau);
        }
    }
    if(j < n) {
        //if left neighbor is sick
        if(Pop[i+1][j] > 0){
            t = t + (rand < tau);
        }
    }
    int p = 0;
    if(t > 0){
        p = 1;
    }
}