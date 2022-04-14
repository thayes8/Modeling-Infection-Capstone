void newPosition(float delta, int n, float rand, int** pop, int** npop){
    if(delta > 0){
        for(i = 0; i<n; i++){
            for(j=0; j<n; j++){
                if(rand<delta){
                    int inew = floor(rand*n+1);
                    int jnew = floor(rand*n+1);
                    int tt = npop[i][j];
                    npop[i][j] = npop[inew][jnew];
                    npop[inew][jnew] = tt;
                }
            }
        }
    }
    pop = npop;
}